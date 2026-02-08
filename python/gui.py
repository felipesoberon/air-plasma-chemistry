"""Tkinter GUI that runs the C++ airGM solver and plots output.csv live."""

from __future__ import annotations

import math
import shlex
import subprocess
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

from matplotlib import colormaps
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from .constants import NO_SPECIES, SPECIES_FORMULAS

DEFAULT_PLASMA_TIME_S = 1e-9
SUM_POSITIVE_IONS_KEY = 1001
SUM_NEGATIVE_IONS_KEY = 1002
POSITIVE_ION_INDICES = tuple(range(1, 17))
NEGATIVE_ION_INDICES = tuple(range(18, 28))


class AirGMGui:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("airGM GUI (C++ backend)")
        self.root.geometry("1400x900")

        self.repo_root = Path(__file__).resolve().parent.parent
        self.output_csv_path = self.repo_root / "output.csv"

        self.process: subprocess.Popen | None = None
        self.running = False
        self.last_csv_mtime: float | None = None
        self.last_csv_size: int | None = None
        self.latest_species_values: list[float] = []

        self.plot_times: list[float] = []
        self.plot_data: dict[int, list[float]] = {i: [] for i in range(1, NO_SPECIES + 1)}
        self.plot_data[SUM_POSITIVE_IONS_KEY] = []
        self.plot_data[SUM_NEGATIVE_IONS_KEY] = []

        self.lines = {}
        self.species_vars: dict[int, tk.BooleanVar] = {}
        self.endpoint_label_artists = []

        self._build_layout()
        self._init_plot()
        self._set_defaults()
        self._schedule_poll()

    def _build_layout(self) -> None:
        paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        self.left = ttk.Frame(paned, padding=8)
        self.right = ttk.Frame(paned, padding=8)

        paned.add(self.left, weight=1)
        paned.add(self.right, weight=3)

        self._build_inputs(self.left)
        self._build_species_controls(self.left)
        self._build_plot_panel(self.right)

    def _build_inputs(self, parent: ttk.Frame) -> None:
        inputs = ttk.LabelFrame(parent, text="Model Inputs", padding=8)
        inputs.pack(fill=tk.X, pady=(0, 8))

        self.te_var = tk.StringVar()
        self.total_time_var = tk.StringVar()
        self.rh_var = tk.StringVar()
        self.metric_min_var = tk.StringVar()
        self.metric_max_var = tk.StringVar()

        self.show_all_var = tk.BooleanVar(value=False)
        self.show_latest_nonzero_var = tk.BooleanVar(value=False)
        self.show_latest_zero_var = tk.BooleanVar(value=False)
        self.current_time_var = tk.StringVar(value="t = 0 s")
        self.current_metric_var = tk.StringVar(value="max rel dC/C = 0")
        self.current_dt_var = tk.StringVar(value="dt = n/a")

        self._add_labeled_entry(inputs, 0, "Te [eV]", self.te_var)
        self._add_labeled_entry(inputs, 1, "Total time [s]", self.total_time_var)
        self._add_labeled_entry(inputs, 2, "RH [%]", self.rh_var)
        self._add_labeled_entry(inputs, 3, "Metric min", self.metric_min_var)
        self._add_labeled_entry(inputs, 4, "Metric max", self.metric_max_var)

        ttk.Separator(inputs, orient=tk.HORIZONTAL).grid(row=5, column=0, columnspan=2, sticky="ew", pady=6)

        controls = ttk.Frame(inputs)
        controls.grid(row=6, column=0, columnspan=2, sticky="ew", pady=(8, 2))
        ttk.Button(controls, text="Start", command=self._start).pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(controls, text="Stop", command=self._stop).pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(controls, text="Reset", command=self._reset).pack(side=tk.LEFT)

        ttk.Label(inputs, textvariable=self.current_time_var).grid(row=7, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(inputs, textvariable=self.current_metric_var).grid(row=8, column=0, columnspan=2, sticky="w", pady=(2, 0))
        ttk.Label(inputs, textvariable=self.current_dt_var).grid(row=9, column=0, columnspan=2, sticky="w", pady=(2, 0))

        self.status_var = tk.StringVar(value="Idle")
        ttk.Label(inputs, textvariable=self.status_var).grid(row=10, column=0, columnspan=2, sticky="w", pady=(6, 0))

        for col in (0, 1):
            inputs.columnconfigure(col, weight=1)

    @staticmethod
    def _add_labeled_entry(parent: ttk.Frame, row: int, label: str, var: tk.StringVar) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(parent, textvariable=var, width=18).grid(row=row, column=1, sticky="ew", pady=2)

    def _build_species_controls(self, parent: ttk.Frame) -> None:
        species_frame = ttk.LabelFrame(parent, text="Species Visibility", padding=8)
        species_frame.pack(fill=tk.BOTH, expand=True)

        ttk.Checkbutton(
            species_frame,
            text="Show all species",
            variable=self.show_all_var,
            command=self._toggle_show_all,
        ).pack(anchor="w", pady=(0, 4))
        ttk.Checkbutton(
            species_frame,
            text="Show non-zero at latest time",
            variable=self.show_latest_nonzero_var,
            command=self._toggle_show_latest_nonzero,
        ).pack(anchor="w", pady=(0, 2))
        ttk.Checkbutton(
            species_frame,
            text="Show zero at latest time",
            variable=self.show_latest_zero_var,
            command=self._toggle_show_latest_zero,
        ).pack(anchor="w", pady=(0, 4))

        grid_body = ttk.Frame(species_frame)
        grid_body.pack(fill=tk.BOTH, expand=True)

        columns = 4
        for col in range(columns):
            grid_body.columnconfigure(col, weight=1)

        for i in range(1, NO_SPECIES + 1):
            var = tk.BooleanVar(value=False)
            self.species_vars[i] = var
            ttk.Checkbutton(
                grid_body,
                text=f"{i:2d}: {SPECIES_FORMULAS[i]}",
                variable=var,
                command=self._update_line_visibility,
            ).grid(row=(i - 1) // columns, column=(i - 1) % columns, sticky="w", padx=(0, 8), pady=1)

        next_row = (NO_SPECIES - 1) // columns + 1
        self.species_vars[SUM_POSITIVE_IONS_KEY] = tk.BooleanVar(value=False)
        self.species_vars[SUM_NEGATIVE_IONS_KEY] = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            grid_body,
            text="SUM + ions",
            variable=self.species_vars[SUM_POSITIVE_IONS_KEY],
            command=self._update_line_visibility,
        ).grid(row=next_row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Checkbutton(
            grid_body,
            text="SUM - ions",
            variable=self.species_vars[SUM_NEGATIVE_IONS_KEY],
            command=self._update_line_visibility,
        ).grid(row=next_row, column=1, sticky="w", padx=(0, 8), pady=3)

    def _build_plot_panel(self, parent: ttk.Frame) -> None:
        self.figure = Figure(figsize=(10, 7), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel("Time [s]")
        self.ax.set_ylabel("Density [m^-3]")
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.grid(True, which="both", alpha=0.3)

        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

    def _init_plot(self) -> None:
        colors = colormaps["tab20"]
        for i in range(1, NO_SPECIES + 1):
            color = colors((i - 1) % 20)
            line, = self.ax.plot([], [], linewidth=1.0, color=color, label=SPECIES_FORMULAS[i])
            self.lines[i] = line
        self.lines[SUM_POSITIVE_IONS_KEY], = self.ax.plot(
            [], [], linewidth=1.8, linestyle="--", color="black", label="SUM + ions"
        )
        self.lines[SUM_NEGATIVE_IONS_KEY], = self.ax.plot(
            [], [], linewidth=1.8, linestyle=":", color="dimgray", label="SUM - ions"
        )
        self._update_line_visibility()
        self.canvas.draw_idle()

    def _set_defaults(self) -> None:
        self.te_var.set("5.5")
        self.total_time_var.set("1")
        self.rh_var.set("50")
        self.metric_min_var.set("0.01")
        self.metric_max_var.set("0.05")

    def _toggle_show_all(self) -> None:
        value = self.show_all_var.get()
        self.show_latest_nonzero_var.set(False)
        self.show_latest_zero_var.set(False)
        for var in self.species_vars.values():
            var.set(value)
        self._update_line_visibility()

    def _toggle_show_latest_nonzero(self) -> None:
        if self.show_latest_nonzero_var.get():
            self.show_latest_zero_var.set(False)
            self._apply_latest_presence_filter(show_nonzero=True)
            return
        self._update_line_visibility()

    def _toggle_show_latest_zero(self) -> None:
        if self.show_latest_zero_var.get():
            self.show_latest_nonzero_var.set(False)
            self._apply_latest_presence_filter(show_nonzero=False)
            return
        self._update_line_visibility()

    def _apply_latest_presence_filter(self, show_nonzero: bool) -> None:
        if not self.latest_species_values:
            self.status_var.set("No data yet for latest-time species filter")
            if show_nonzero:
                self.show_latest_nonzero_var.set(False)
            else:
                self.show_latest_zero_var.set(False)
            self._update_line_visibility()
            return

        for i in range(1, NO_SPECIES + 1):
            value = self.latest_species_values[i - 1]
            is_nonzero = value > 0.0
            self.species_vars[i].set(is_nonzero if show_nonzero else not is_nonzero)
        self._update_line_visibility()

    def _update_line_visibility(self) -> None:
        visible_count = 0
        for i, line in self.lines.items():
            visible = self.species_vars[i].get()
            line.set_visible(visible)
            if visible:
                visible_count += 1

        self.show_all_var.set(visible_count == len(self.lines))
        self._autoscale_axes()
        self._update_endpoint_labels()
        self.canvas.draw_idle()

    def _clear_endpoint_labels(self) -> None:
        for artist in self.endpoint_label_artists:
            artist.remove()
        self.endpoint_label_artists.clear()

    def _update_endpoint_labels(self) -> None:
        self._clear_endpoint_labels()
        if not self.plot_times:
            return

        candidates = []
        for key, line in self.lines.items():
            if not line.get_visible():
                continue
            series = self.plot_data.get(key)
            if not series:
                continue

            for idx in range(len(series) - 1, -1, -1):
                y_val = series[idx]
                if not math.isnan(y_val) and y_val > 0.0:
                    x_val = self.plot_times[idx]
                    label = SPECIES_FORMULAS[key] if 1 <= key <= NO_SPECIES else line.get_label()
                    candidates.append((key, x_val, y_val, label, line.get_color()))
                    break

        if not candidates:
            return

        candidates.sort(key=lambda item: item[2])
        y_min_px = self.ax.bbox.ymin + 4.0
        y_max_px = self.ax.bbox.ymax - 4.0
        min_gap_px = 11.0
        dpi = self.figure.dpi

        adjusted_display_y: list[float] = []
        for _, x_val, y_val, _, _ in candidates:
            display_y = self.ax.transData.transform((x_val, y_val))[1]
            display_y = max(y_min_px, min(y_max_px, display_y))
            if adjusted_display_y:
                display_y = max(display_y, adjusted_display_y[-1] + min_gap_px)
                display_y = min(display_y, y_max_px)
            adjusted_display_y.append(display_y)

        for idx in range(len(adjusted_display_y) - 2, -1, -1):
            if adjusted_display_y[idx + 1] - adjusted_display_y[idx] < min_gap_px:
                adjusted_display_y[idx] = max(y_min_px, adjusted_display_y[idx + 1] - min_gap_px)

        for (key, x_val, y_val, label, color), final_display_y in zip(candidates, adjusted_display_y):
            raw_display_y = self.ax.transData.transform((x_val, y_val))[1]
            y_offset_points = (final_display_y - raw_display_y) * 72.0 / dpi
            text = self.ax.annotate(
                label,
                xy=(x_val, y_val),
                xytext=(6.0, y_offset_points),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=8,
                color=color,
                clip_on=True,
            )
            self.endpoint_label_artists.append(text)

    @staticmethod
    def _windows_path_to_wsl(path: Path) -> str:
        drive = path.drive.rstrip(":").lower()
        tail = path.as_posix().split(":", 1)[1]
        return f"/mnt/{drive}{tail}"

    def _collect_inputs(self) -> dict[str, float] | None:
        try:
            te = float(self.te_var.get())
            total_time = float(self.total_time_var.get())
            rh = float(self.rh_var.get())
            metric_min = float(self.metric_min_var.get())
            metric_max = float(self.metric_max_var.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Please provide numeric values for all model inputs.")
            return None

        if metric_min <= 0.0 or metric_max <= 0.0 or metric_min >= metric_max:
            messagebox.showerror("Invalid limits", "Require 0 < metric min < metric max.")
            return None

        return {
            "te": te,
            "total_time": total_time,
            "rh": max(0.0, min(100.0, rh)),
            "metric_min": metric_min,
            "metric_max": metric_max,
        }

    def _start(self) -> None:
        if self.running:
            return

        params = self._collect_inputs()
        if params is None:
            return

        repo_wsl = self._windows_path_to_wsl(self.repo_root)
        args = [
            "./src/airGM2.1",
            "-Te",
            f"{params['te']}",
            "-totaltime",
            f"{params['total_time']}",
            "-plasmatime",
            f"{DEFAULT_PLASMA_TIME_S}",
            "-RH",
            f"{params['rh']}",
            "-metricmin",
            f"{params['metric_min']}",
            "-metricmax",
            f"{params['metric_max']}",
        ]
        quoted_args = " ".join(shlex.quote(item) for item in args)
        bash_cmd = (
            f"cd {shlex.quote(repo_wsl)} && "
            f"if [ ! -x ./src/airGM2.1 ]; then (cd src && make); fi && "
            f"{quoted_args}"
        )

        try:
            self.process = subprocess.Popen(
                ["wsl.exe", "--", "bash", "-lc", bash_cmd],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except FileNotFoundError:
            messagebox.showerror("WSL not found", "wsl.exe is required to run the C++ backend from this GUI.")
            return

        self.running = True
        self.status_var.set("Running (C++ backend)")

    def _stop(self) -> None:
        if self.process is None:
            return

        self.status_var.set("Stopping...")
        self.process.terminate()
        try:
            self.process.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            self.process.kill()
        self.process = None
        self.running = False
        self.status_var.set("Stopped")

    def _reset(self) -> None:
        self._stop()
        if self.output_csv_path.exists():
            self.output_csv_path.unlink()
        self.last_csv_mtime = None
        self.last_csv_size = None
        self._clear_plot_data()
        self.status_var.set("Reset complete")

    def _clear_plot_data(self) -> None:
        self.plot_times.clear()
        self.latest_species_values = []
        for key in self.lines:
            self.plot_data[key].clear()
            self.lines[key].set_data([], [])
        self._clear_endpoint_labels()

        self.current_time_var.set("t = 0 s")
        self.current_metric_var.set("max rel dC/C = 0")
        self.current_dt_var.set("dt = n/a")
        self._autoscale_axes()
        self.canvas.draw_idle()

    def _load_csv_rows(self) -> list[list[float]]:
        rows: list[list[float]] = []
        if not self.output_csv_path.exists():
            return rows

        with self.output_csv_path.open("r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rows.append([float(x) for x in line.split(",") if x != ""])

        return rows

    def _refresh_from_csv(self) -> None:
        if not self.output_csv_path.exists():
            return

        stat = self.output_csv_path.stat()
        if self.last_csv_mtime == stat.st_mtime and self.last_csv_size == stat.st_size:
            return
        self.last_csv_mtime = stat.st_mtime
        self.last_csv_size = stat.st_size

        rows = self._load_csv_rows()
        if not rows:
            return

        self.plot_times.clear()
        for key in self.lines:
            self.plot_data[key].clear()

        for row in rows:
            time_s = row[53]
            if time_s <= 0:
                continue
            self.plot_times.append(time_s)
            for i in range(1, NO_SPECIES + 1):
                value = row[i - 1]
                self.plot_data[i].append(value if value > 0 else math.nan)
            pos_sum = sum(row[i - 1] for i in POSITIVE_ION_INDICES)
            neg_sum = sum(row[i - 1] for i in NEGATIVE_ION_INDICES)
            self.plot_data[SUM_POSITIVE_IONS_KEY].append(pos_sum if pos_sum > 0 else math.nan)
            self.plot_data[SUM_NEGATIVE_IONS_KEY].append(neg_sum if neg_sum > 0 else math.nan)

        for key in self.lines:
            self.lines[key].set_data(self.plot_times, self.plot_data[key])

        last = rows[-1]
        self.latest_species_values = [last[i - 1] for i in range(1, NO_SPECIES + 1)]
        current_time = last[53]
        self.current_time_var.set(f"t = {current_time:.6g} s")

        if len(rows) >= 2:
            prev = rows[-2]
            dt_saved = last[53] - prev[53]
            if dt_saved > 0:
                self.current_dt_var.set(f"dt = {dt_saved:.3e} s (saved)")

            max_rel = 0.0
            for idx in range(53):
                prev_value = prev[idx]
                if prev_value > 0.0:
                    rel = abs((last[idx] - prev_value) / prev_value)
                    if rel > max_rel:
                        max_rel = rel
            self.current_metric_var.set(f"max rel dC/C = {max_rel:.3e} (saved)")
        else:
            self.current_metric_var.set("max rel dC/C = 0")
            self.current_dt_var.set("dt = n/a")

        if self.show_latest_nonzero_var.get():
            self._apply_latest_presence_filter(show_nonzero=True)
        elif self.show_latest_zero_var.get():
            self._apply_latest_presence_filter(show_nonzero=False)

        self._autoscale_axes()
        self._update_endpoint_labels()
        self.canvas.draw_idle()

    def _autoscale_axes(self) -> None:
        if not self.plot_times:
            self.ax.set_xlim(1e-12, 1e-3)
            self.ax.set_ylim(1e-3, 1e26)
            return

        x_min = min(self.plot_times)
        x_max = max(self.plot_times)
        if x_min == x_max:
            x_min *= 0.9
            x_max *= 1.1

        y_min = float("inf")
        y_max = 0.0
        for key in self.lines:
            if not self.lines[key].get_visible():
                continue
            for y in self.plot_data[key]:
                if not math.isnan(y) and y > 0:
                    y_min = min(y_min, y)
                    y_max = max(y_max, y)

        if y_min == float("inf"):
            y_min, y_max = 1e-3, 1e26
        if y_min == y_max:
            y_min *= 0.9
            y_max *= 1.1

        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

    def _poll(self) -> None:
        self._refresh_from_csv()

        if self.process is not None and self.running:
            code = self.process.poll()
            if code is not None:
                self.running = False
                self.process = None
                if code == 0:
                    self.status_var.set("Completed")
                else:
                    self.status_var.set(f"C++ process failed (exit {code})")

    def _schedule_poll(self) -> None:
        self._poll()
        self.root.after(200, self._schedule_poll)


def main() -> None:
    root = tk.Tk()
    AirGMGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()
