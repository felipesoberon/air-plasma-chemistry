"""Global model port of src/globalmodel.cpp."""

from __future__ import annotations

import math
from pathlib import Path

from .commandline import CommandLine
from .constants import EV_TO_KELVIN, NO_REACTIONS, NO_SPECIES, SPECIES_CONTAINING_H
from .reaction import Reaction
from .species import Species


class GlobalModel:
    def __init__(self) -> None:
        self.species = [Species(i) for i in range(NO_SPECIES + 1)]
        self.reactions = [Reaction() for _ in range(NO_REACTIONS + 1)]

        self.peak_electron_temperature_kelvin = 298.0
        self.electron_temperature_kelvin = 298.0
        self.gas_temperature_kelvin = 298.0

        self.step_count = 0
        self.dt = 50e-12
        self.simulation_time = 0.0
        self.last_saved_simulation_time = -1.0
        self.plasma_time = 1e-9
        self.total_time = 1e-6

        self.command_line = CommandLine()

    def set_parameters_from_command_line_input(self, argv: list[str]) -> None:
        print("\nCOMMAND_LINE_INPUT_PARAMETERS\n")

        self.command_line.set_argument_list(len(argv), argv)
        self.command_line.print_argument_list()

        self.command_line.set_flag_name("-Te", "Electron temperature in eV")
        self.command_line.set_flag_name("-[H2O]", "Density of water in m-3")
        self.command_line.set_flag_name("-totaltime", "Total simulation time is s")
        self.command_line.set_flag_name("-plasmatime", "Plasma pulse time in s")
        self.command_line.set_flag_name("-dt", "Simulation time step in s")
        self.command_line.print_flag_name_list()

        self.command_line.set_flag_values()
        self.command_line.print_flag_values()

        if self.command_line.flag_value_is_number(0):
            self.peak_electron_temperature_kelvin = EV_TO_KELVIN * self.command_line.return_float_flag_value(0)
        if self.command_line.flag_value_is_number(1):
            self.set_h2o_density(self.command_line.return_float_flag_value(1))
        if self.command_line.flag_value_is_number(2):
            self.total_time = self.command_line.return_float_flag_value(2)
        if self.command_line.flag_value_is_number(3):
            self.plasma_time = self.command_line.return_float_flag_value(3)
        if self.command_line.flag_value_is_number(4):
            self.dt = self.command_line.return_float_flag_value(4)

        print()

    def set_species_formula(self) -> None:
        for i in range(0, NO_SPECIES + 1):
            self.species[i].set_formula(i)

    def set_default_species_densities(self) -> None:
        for i in range(0, NO_SPECIES + 1):
            self.species[i].set_density(0.0)

        self.species[0].set_density(2.40e25)
        self.species[17].set_density(1.00e3)
        self.species[36].set_density(0.00e0)
        self.species[39].set_density(0.00e0)
        self.species[51].set_density(1.92e25)
        self.species[52].set_density(4.80e24)
        self.species[53].set_density(1.20e24)

    def set_h2o_density(self, density_value: float) -> None:
        self.species[53].set_density(density_value)

    def set_no2_density(self, density_value: float) -> None:
        self.species[39].set_density(density_value)

    def set_o3_density(self, density_value: float) -> None:
        self.species[36].set_density(density_value)

    def set_electron_temperature_ev(self, temperature_value_ev: float) -> None:
        self.electron_temperature_kelvin = temperature_value_ev * EV_TO_KELVIN

    def set_electron_temperature_kelvin(self, temperature_value: float) -> None:
        self.electron_temperature_kelvin = temperature_value

    def set_gas_temperature_kelvin(self, temperature_value: float) -> None:
        self.gas_temperature_kelvin = temperature_value

    def set_reaction_rates(self) -> None:
        for j in range(1, NO_REACTIONS + 1):
            self.reactions[j].set_reaction_rate(j, self.gas_temperature_kelvin, self.electron_temperature_kelvin)

    def set_reaction_reactant_and_product_species(self) -> None:
        for j in range(1, NO_REACTIONS + 1):
            self.reactions[j].set_reactant_and_product_species(j)

    def set_balance_equations(self) -> None:
        for i in range(1, NO_SPECIES + 1):
            for j in range(1, NO_REACTIONS + 1):
                no_reactants = self.reactions[j].return_number_of_reactants()
                repeat_loss = 0
                for k in range(1, no_reactants + 1):
                    if self.reactions[j].return_reactant(k) == i:
                        repeat_loss += 1
                if repeat_loss > 0:
                    self.species[i].increment_number_of_reaction_no_losses()
                    item_index = self.species[i].return_number_of_reaction_no_losses()
                    self.species[i].set_reaction_no_losses_item(item_index, j)
                    self.species[i].set_reaction_no_losses_multiplier_item(item_index, repeat_loss)

                no_products = self.reactions[j].return_number_of_products()
                repeat_source = 0
                for k in range(1, no_products + 1):
                    if self.reactions[j].return_product(k) == i:
                        repeat_source += 1
                if repeat_source > 0:
                    self.species[i].increment_number_of_reaction_no_sources()
                    item_index = self.species[i].return_number_of_reaction_no_sources()
                    self.species[i].set_reaction_no_sources_item(item_index, j)
                    self.species[i].set_reaction_no_sources_multiplier_item(item_index, repeat_source)

            self.species[i].process_reduce_lists_of_reaction_no()

    def process_balance_equations(self) -> None:
        for i in range(1, 51):
            self.species[i].set_loss(0.0)
            self.species[i].set_source(0.0)

            if self.species[53].return_density() == 0.0 and i in SPECIES_CONTAINING_H:
                continue

            if self.species[i].return_density() > 0.0:
                no_reactions = self.species[i].return_number_of_reaction_no_losses()
                for h in range(1, no_reactions + 1):
                    j = self.species[i].return_reaction_no_losses(h)
                    no_reactants = self.reactions[j].return_number_of_reactants()
                    aux_density = 1.0
                    for ri in range(1, no_reactants + 1):
                        reactant_index = self.reactions[j].return_reactant(ri)
                        aux_density *= self.species[reactant_index].return_density()

                    multiplier = self.species[i].return_reaction_no_losses_multiplier(h)
                    reaction_rate = self.reactions[j].return_reaction_rate()
                    aux_loss = self.species[i].return_loss()
                    aux_loss += reaction_rate * aux_density * multiplier
                    self.species[i].set_loss(aux_loss)

            no_reactions = self.species[i].return_number_of_reaction_no_sources()
            for h in range(1, no_reactions + 1):
                j = self.species[i].return_reaction_no_sources(h)
                no_reactants = self.reactions[j].return_number_of_reactants()
                aux_density = 1.0
                for ri in range(1, no_reactants + 1):
                    reactant_index = self.reactions[j].return_reactant(ri)
                    aux_density *= self.species[reactant_index].return_density()

                multiplier = self.species[i].return_reaction_no_sources_multiplier(h)
                reaction_rate = self.reactions[j].return_reaction_rate()
                aux_source = self.species[i].return_source()
                aux_source += reaction_rate * aux_density * multiplier
                self.species[i].set_source(aux_source)

    def process_time_step_species_densities(self) -> None:
        for i in range(1, 51):
            if self.species[i].return_loss() > 0.0 or self.species[i].return_source() > 0.0:
                sources_minus_losses = self.species[i].return_source() - self.species[i].return_loss()
                aux_density = sources_minus_losses * self.dt

                if (
                    sources_minus_losses > 0.0
                    and aux_density > self.species[i].return_density()
                    and self.species[i].return_density() > 1e5
                ):
                    print(
                        f"WARNING: species [{self.species[i].return_formula()}] > x2 density at time {self.simulation_time}"
                    )

                aux_density = self.species[i].return_density() + aux_density
                if sources_minus_losses != 0.0:
                    self.species[i].set_density(aux_density)

    @staticmethod
    def _fmt_float(value: float) -> str:
        return f"{value:.6g}"

    def process_main_loop(self) -> None:
        output_path = Path("output.csv")
        with output_path.open("a", encoding="utf-8") as dumpfile:
            if self.simulation_time == 0.0:
                header = ["#" + self.species[i].return_formula() for i in range(1, NO_SPECIES + 1)]
                dumpfile.write(",".join(header) + ",Time(s),StepNo\n")

            print(f"\nPLASMA PULSE: (duration = {self.plasma_time})")
            while True:
                if self.simulation_time >= self.total_time or self.simulation_time >= 10.0 * self.plasma_time:
                    break

                self.electron_temperature_kelvin = self.return_electron_temperature_kelvin_at_time()
                self.set_reaction_rates()

                self.process_balance_equations()
                self.process_time_step_species_densities()

                if (
                    self.step_count % self.return_save_interval_step() == 0
                    and self.simulation_time > self.last_saved_simulation_time
                ):
                    row = [self._fmt_float(self.species[i].return_density()) for i in range(1, NO_SPECIES + 1)]
                    row.append(self._fmt_float(self.simulation_time))
                    row.append(str(self.step_count))
                    dumpfile.write(",".join(row) + "\n")
                    print(f"{self.simulation_time}\t{self.step_count}\t{self.electron_temperature_kelvin}")

                self.simulation_time += self.dt
                self.step_count += 1

            self.electron_temperature_kelvin = self.return_electron_temperature_kelvin_at_time()
            self.set_reaction_rates()

            print("\nAFTERGLOW:\n")
            while True:
                if self.simulation_time >= self.total_time:
                    break

                self.process_balance_equations()
                self.process_time_step_species_densities()

                if (
                    self.step_count % self.return_save_interval_step() == 0
                    and self.simulation_time > self.last_saved_simulation_time
                ):
                    row = [self._fmt_float(self.species[i].return_density()) for i in range(1, NO_SPECIES + 1)]
                    row.append(self._fmt_float(self.simulation_time))
                    row.append(str(self.step_count))
                    dumpfile.write(",".join(row) + "\n")
                    print(f"{self.simulation_time}\t{self.step_count}\t{self.electron_temperature_kelvin}")

                self.simulation_time += self.dt
                self.step_count += 1

    def print_species_formula_and_density(self) -> None:
        print("\nSPECIES\n")
        for i in range(0, NO_SPECIES + 1):
            print(f"{i},{self.species[i].return_formula()},{self.species[i].return_density()}")

    def print_list_of_reactions(self) -> None:
        print("\nREACTION_LIST\n")
        print("No.,Rate,r1,r2,r3,r4,--->,p1,p2,p3,p4")
        for i in range(1, NO_REACTIONS + 1):
            reaction_rate_value = self.reactions[i].return_reaction_rate()
            no_reactants = self.reactions[i].return_number_of_reactants()
            no_products = self.reactions[i].return_number_of_products()

            reactants = []
            for j in range(1, 5):
                if j <= no_reactants:
                    reactant_index = self.reactions[i].return_reactant(j)
                    reactants.append(self.species[reactant_index].return_formula())
                else:
                    reactants.append("")

            products = []
            for j in range(1, 5):
                if j <= no_products:
                    product_index = self.reactions[i].return_product(j)
                    products.append(self.species[product_index].return_formula())
                else:
                    products.append("")

            print(f"{i},{reaction_rate_value},{','.join(reactants)},--->,{','.join(products)},")

    def read_species_density_data_file(self) -> None:
        path = Path("output.csv")
        if not path.exists():
            print("WARNING: Problem opening <output.csv> file.")
            print("NOTICE: Program will default to initial species density values.")
            return

        data_line = ""
        with path.open("r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if len(line) > 50:
                    data_line = line

        if not data_line:
            return

        values = []
        for token in data_line.split(","):
            token = token.strip()
            if token:
                values.append(float(token))
            if len(values) >= NO_SPECIES + 2:
                break

        if len(values) < NO_SPECIES + 2:
            return

        print(f"\nNOTICE: <output.csv> file read... {len(values)}")

        density_of_species_m = 0.0
        for i in range(1, NO_SPECIES + 1):
            density_of_species_m += values[i - 1]
            self.species[i].set_density(values[i - 1])

        self.species[0].set_density(density_of_species_m)
        self.simulation_time = values[53]
        self.last_saved_simulation_time = self.simulation_time
        self.step_count = int(values[54])

    def return_electron_temperature_kelvin_at_time(self) -> float:
        temperature_value = self.gas_temperature_kelvin
        if self.simulation_time < 10.0 * self.plasma_time:
            temperature_value = self.gas_temperature_kelvin + (
                self.peak_electron_temperature_kelvin - self.gas_temperature_kelvin
            ) * math.exp(-0.5 * ((self.simulation_time - 5.0 * self.plasma_time) / self.plasma_time) ** 2)
        return temperature_value

    def return_save_interval_step(self) -> int:
        result = 1
        ratio_factor = self.dt / 50e-12

        if self.simulation_time <= 1e3:
            result = 1000000000000
        if self.simulation_time <= 1e2:
            result = 100000000000
        if self.simulation_time <= 1e1:
            result = 10000000000
        if self.simulation_time <= 1e0:
            result = 1000000000
        if self.simulation_time <= 1e-1:
            result = 100000000
        if self.simulation_time <= 1e-2:
            result = 10000000
        if self.simulation_time <= 1e-3:
            result = 1000000
        if self.simulation_time <= 1e-4:
            result = 100000
        if self.simulation_time <= 1e-5:
            result = 10000
        if self.simulation_time <= 1e-6:
            result = 1000
        if self.simulation_time <= 1e-7:
            result = 100
        if self.simulation_time <= 1e-8:
            result = 10
        if self.simulation_time <= 1e-9:
            result = 1

        ratio_factor = result / ratio_factor
        result = int(ratio_factor)
        if result < 1:
            result = 1

        return result
