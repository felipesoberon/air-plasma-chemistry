"""Python entry point for airGM port."""

from __future__ import annotations

import sys

from .globalmodel import GlobalModel


def print_program_description() -> None:
    print("**********************************************")
    print("*                                            *")
    print("* Felipe Soberon (felipe.soberon@gmail.com)  *")
    print("* 2023                                       *")
    print("*                                            *")
    print("* Global_Model_2.1 of atmospheric pressure   *")
    print("* plasma discharge in (humid) air;           *")
    print("* using Sakiyama et al., 2012 reaction data. *")
    print("*                                            *")
    print("**********************************************")
    print("Usage:")
    print("$ airGM2.1 <-flag> <flag value>")


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv

    print_program_description()

    microdischarge_model = GlobalModel()

    microdischarge_model.set_species_formula()
    microdischarge_model.set_default_species_densities()
    microdischarge_model.set_parameters_from_command_line_input(argv)

    microdischarge_model.read_species_density_data_file()

    microdischarge_model.set_reaction_rates()
    microdischarge_model.set_reaction_reactant_and_product_species()
    microdischarge_model.set_balance_equations()

    microdischarge_model.process_main_loop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
