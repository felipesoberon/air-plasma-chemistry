"""Minimal command-line parser matching src/commandline.cpp behavior."""

from __future__ import annotations

import struct


class CommandLine:
    def __init__(self) -> None:
        self.argument: list[str] = []
        self.flag_name: list[str] = []
        self.flag_value: list[str] = []
        self.flag_description: list[str] = []

    def set_argument_list(self, input_number_of_arguments: int, value_of_argument: list[str]) -> None:
        self.argument = [value_of_argument[i] for i in range(input_number_of_arguments)]

    def print_argument_list(self) -> None:
        print("\nECHO COMMAND:", " ".join(self.argument))

    def set_flag_name(self, flag_tag: str, flag_description_value: str) -> None:
        self.flag_name.append(flag_tag)
        self.flag_description.append(flag_description_value)
        self.flag_value.append("EMPTY")

    def print_flag_name_list(self) -> None:
        print()
        for i, name in enumerate(self.flag_name):
            print(f"{i}   {name:>12}\t{self.flag_description[i]}")

    def set_flag_values(self) -> None:
        print()
        for i, name in enumerate(self.flag_name):
            for j, arg in enumerate(self.argument):
                if name == arg:
                    if j + 1 < len(self.argument):
                        self.flag_value[i] = self.argument[j + 1]
                    else:
                        self.flag_value[i] = "EMPTY"

    def print_flag_values(self) -> None:
        for i, name in enumerate(self.flag_name):
            print(f"{i}   {name:>12}\t{self.flag_value[i]}")

    def return_flag_value(self, flag_index: int) -> str:
        if flag_index < len(self.flag_value):
            return self.flag_value[flag_index]
        return "ERROR"

    def return_float_flag_value(self, flag_index: int) -> float:
        if self.flag_value_is_number(flag_index):
            # C++ uses stof (float32) before assigning to double fields.
            value = float(self.flag_value[flag_index])
            return struct.unpack("f", struct.pack("f", value))[0]
        print("ERROR: flag value is not a number.")
        return 0.0

    def flag_value_is_empty(self, flag_index: int) -> bool:
        return flag_index < len(self.flag_value) and self.flag_value[flag_index] == "EMPTY"

    def flag_value_is_number(self, flag_index: int) -> bool:
        if flag_index < len(self.flag_value):
            try:
                float(self.flag_value[flag_index])
                return True
            except ValueError:
                return False
        print("ERROR: flag index exceeds limits.")
        return False
