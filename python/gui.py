"""Tkinter GUI for the airGM Python model."""

from __future__ import annotations

import math
import queue
import threading
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import cm
from matplotlib.figure import Figure

from .constants import EV_TO_KELVIN, NO_SPECIES, SPECIES_FORMULAS
from .globalmodel import GlobalModel

BOLTZMANN_CONSTANT = 1.380649e-23
DEFAULT_GAS_TEMPERATURE_K = 298.0
DEFAULT_PRESSURE_PA = 101325.0
DEFAULT_PLASMA_TIME_S = 1e-9


class AirGMGui:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("airGM GUI")
        self.root.geometry("1400x900")

        self.model_thread: threading.Thread | None = None
        self.stop_event = threading.Event()
        self.data_queue: queue.Queue = queue.Queue()
        self.running = False

        self.times: list[float] = []
        self.plot_times: list[float] = []
        self.plot_data: dict[int, list[float]] = {i: [] for i in range(1, NO_SPECIES + 1)}

        self.lines = {}
        self.species_vars: dict[int, tk.BooleanVar] = {}
        self.pending_points = 0

        self._build_layout()
        self._init_plot()
        self._set_defaults()
        self._schedule_queue_poll()

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
        self.current_time_var = tk.StringVar(value="t = 0 s")
        self.current_metric_var = tk.StringVar(value="max rel dC/C = 0")
        self.current_dt_var = tk.StringVar(value="dt = 5e-11 s")

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

        ttk.Label(inputs, textvariable=self.current_time_var).grid(
            row=7, column=0, columnspan=2, sticky="w", pady=(6, 0)
        )
        ttk.Label(inputs, textvariable=self.current_metric_var).grid(
            row=8, column=0, columnspan=2, sticky="w", pady=(2, 0)
        )
        ttk.Label(inputs, textvariable=self.current_dt_var).grid(
            row=9, column=0, columnspan=2, sticky="w", pady=(2, 0)
        )

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
        colors = cm.get_cmap("tab20")
        for i in range(1, NO_SPECIES + 1):
            color = colors((i - 1) % 20)
            line, = self.ax.plot([], [], linewidth=1.0, color=color, label=SPECIES_FORMULAS[i])
            self.lines[i] = line
        self._update_line_visibility()
        self.canvas.draw_idle()

    def _set_defaults(self) -> None:
        self.te_var.set("2.6")
        self.total_time_var.set("1e-6")
        self.rh_var.set("50")
        self.metric_min_var.set("0.05")
        self.metric_max_var.set("0.5")

    def _toggle_show_all(self) -> None:
        value = self.show_all_var.get()
        for var in self.species_vars.values():
            var.set(value)
        self._update_line_visibility()

    def _update_line_visibility(self) -> None:
        visible_count = 0
        for i, line in self.lines.items():
            visible = self.species_vars[i].get()
            line.set_visible(visible)
            if visible:
                visible_count += 1

        if visible_count != NO_SPECIES:
            self.show_all_var.set(False)
        elif visible_count == NO_SPECIES:
            self.show_all_var.set(True)

        self._autoscale_axes()
        self.canvas.draw_idle()

    def _compute_saturation_pressure_pa(self, temperature_k: float) -> float:
        # Tetens approximation (water over liquid), valid near room conditions.
        temp_c = temperature_k - 273.15
        return 610.94 * math.exp((17.625 * temp_c) / (temp_c + 243.04))

    def _rh_to_density(self, rh: float, temperature_k: float) -> float:
        rh_fraction = max(0.0, min(100.0, rh)) / 100.0
        p_sat = self._compute_saturation_pressure_pa(temperature_k)
        p_h2o = rh_fraction * p_sat
        return p_h2o / (BOLTZMANN_CONSTANT * temperature_k)

    def _clear_plot_data(self) -> None:
        self.times.clear()
        self.plot_times.clear()
        for i in range(1, NO_SPECIES + 1):
            self.plot_data[i].clear()
            self.lines[i].set_data([], [])
        self.pending_points = 0
        self.current_time_var.set("t = 0 s")
        self.current_metric_var.set("max rel dC/C = 0")
        self.current_dt_var.set("dt = 5e-11 s")
        self._autoscale_axes()
        self.canvas.draw_idle()

    def _autoscale_axes(self) -> None:
        has_points = len(self.plot_times) > 0
        if not has_points:
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
        for i in range(1, NO_SPECIES + 1):
            if not self.lines[i].get_visible():
                continue
            for y in self.plot_data[i]:
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

    def _collect_inputs(self) -> dict[str, float] | None:
        try:
            return {
                "te_ev": float(self.te_var.get()),
                "total_time": float(self.total_time_var.get()),
                "rh_percent": float(self.rh_var.get()),
                "metric_min": float(self.metric_min_var.get()),
                "metric_max": float(self.metric_max_var.get()),
            }
        except ValueError:
            messagebox.showerror("Invalid input", "Please provide numeric values for all model inputs.")
            return None

    def _start(self) -> None:
        if self.running:
            return

        params = self._collect_inputs()
        if params is None:
            return
        if params["metric_min"] <= 0.0 or params["metric_max"] <= 0.0 or params["metric_min"] >= params["metric_max"]:
            messagebox.showerror("Invalid limits", "Require 0 < metric min < metric max.")
            return

        self.stop_event.clear()
        self.running = True
        self.status_var.set("Running")

        model = GlobalModel()
        model.set_species_formula()
        model.set_default_species_densities()

        model.peak_electron_temperature_kelvin = params["te_ev"] * EV_TO_KELVIN
        model.total_time = params["total_time"]
        model.plasma_time = DEFAULT_PLASMA_TIME_S
        model.set_gas_temperature_kelvin(DEFAULT_GAS_TEMPERATURE_K)
        model.set_h2o_density(self._rh_to_density(params["rh_percent"], DEFAULT_GAS_TEMPERATURE_K))
        model.relative_change_min_limit = params["metric_min"]
        model.relative_change_max_limit = params["metric_max"]

        model.read_species_density_data_file()
        model.set_reaction_rates()
        model.set_reaction_reactant_and_product_species()
        model.set_balance_equations()

        self.model_thread = threading.Thread(target=self._run_model, args=(model,), daemon=True)
        self.model_thread.start()

    def _run_model(self, model: GlobalModel) -> None:
        try:
            model.process_main_loop(
                stop_requested=self.stop_event.is_set,
                on_saved_row=lambda t, step, dens: self.data_queue.put(
                    ("row", t, step, dens, model.last_max_relative_change, model.dt)
                ),
            )
            self.data_queue.put(("done",))
        except Exception as exc:  # pragma: no cover - UI runtime guard
            self.data_queue.put(("error", str(exc)))

    def _stop(self) -> None:
        if self.running:
            self.stop_event.set()
            self.status_var.set("Stopping...")

    def _reset(self) -> None:
        self._stop()
        if self.model_thread is not None and self.model_thread.is_alive():
            self.model_thread.join(timeout=2.0)

        output_path = Path("output.csv")
        if output_path.exists():
            output_path.unlink()

        self._clear_plot_data()
        self.status_var.set("Reset complete")

    def _append_point(
        self, time_s: float, densities: list[float], max_relative_change: float, dt: float
    ) -> None:
        self.times.append(time_s)
        self.current_time_var.set(f"t = {time_s:.6g} s")
        self.current_metric_var.set(f"max rel dC/C = {max_relative_change:.3e}")
        self.current_dt_var.set(f"dt = {dt:.3e} s")
        if time_s <= 0:
            return

        self.plot_times.append(time_s)
        for i in range(1, NO_SPECIES + 1):
            value = densities[i - 1]
            self.plot_data[i].append(value if value > 0 else math.nan)

        for i in range(1, NO_SPECIES + 1):
            self.lines[i].set_data(self.plot_times, self.plot_data[i])

        self.pending_points += 1

    def _drain_queue(self) -> None:
        refreshed = False
        while True:
            try:
                item = self.data_queue.get_nowait()
            except queue.Empty:
                break

            kind = item[0]
            if kind == "row":
                _, time_s, _step_no, densities, max_relative_change, dt = item
                self._append_point(time_s, densities, max_relative_change, dt)
                refreshed = True
            elif kind == "done":
                self.running = False
                self.status_var.set("Completed")
                refreshed = True
            elif kind == "error":
                self.running = False
                self.status_var.set("Error")
                messagebox.showerror("Simulation error", item[1])
                refreshed = True

        if refreshed and self.pending_points > 0:
            self._autoscale_axes()
            self.canvas.draw_idle()
            self.pending_points = 0

    def _schedule_queue_poll(self) -> None:
        self._drain_queue()
        self.root.after(100, self._schedule_queue_poll)


def main() -> None:
    root = tk.Tk()
    AirGMGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()
