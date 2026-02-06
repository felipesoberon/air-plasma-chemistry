"""Species model for the airGM Python port.

This mirrors src/species.h and src/species.cpp behavior.
"""

from __future__ import annotations

from .constants import NO_REACTIONS, SPECIES_FORMULAS


class Species:
    """Holds species state and source/loss reaction membership."""

    MINIMUM_DENSITY = 1.0e-3
    _MAX_NO_REACTIONS = NO_REACTIONS

    def __init__(self, index: int = 0) -> None:
        self.formula = "X"
        self.density = 0.0
        self.loss = 0.0
        self.source = 0.0

        # Index 0 stores the current count, matching the C++ layout.
        size = self._MAX_NO_REACTIONS + 1
        self.list_of_reaction_no_sources = [0] * size
        self.list_of_reaction_no_losses = [0] * size
        self.list_of_reaction_no_sources_multiplier = [0] * size
        self.list_of_reaction_no_losses_multiplier = [0] * size

        self.set_formula(index)

    def set_formula(self, species_index: int) -> None:
        self.formula = SPECIES_FORMULAS.get(species_index, "X")

    def set_density(self, density_value: float) -> None:
        if density_value < self.MINIMUM_DENSITY:
            self.density = 0.0
        else:
            self.density = float(density_value)

    def set_loss(self, loss_value: float) -> None:
        self.loss = float(loss_value)

    def set_source(self, source_value: float) -> None:
        self.source = float(source_value)

    def increment_number_of_reaction_no_sources(self) -> None:
        self.list_of_reaction_no_sources[0] += 1

    def increment_number_of_reaction_no_losses(self) -> None:
        self.list_of_reaction_no_losses[0] += 1

    def set_reaction_no_sources_item(self, index: int, value: int) -> None:
        self.list_of_reaction_no_sources[index] = value

    def set_reaction_no_losses_item(self, index: int, value: int) -> None:
        self.list_of_reaction_no_losses[index] = value

    def set_reaction_no_sources_multiplier_item(self, index: int, value: int) -> None:
        self.list_of_reaction_no_sources_multiplier[index] = value

    def set_reaction_no_losses_multiplier_item(self, index: int, value: int) -> None:
        self.list_of_reaction_no_losses_multiplier[index] = value

    def return_number_of_reaction_no_sources(self) -> int:
        return self.list_of_reaction_no_sources[0]

    def return_number_of_reaction_no_losses(self) -> int:
        return self.list_of_reaction_no_losses[0]

    def return_reaction_no_sources(self, index: int) -> int:
        return self.list_of_reaction_no_sources[index]

    def return_reaction_no_losses(self, index: int) -> int:
        return self.list_of_reaction_no_losses[index]

    def return_reaction_no_sources_multiplier(self, index: int) -> int:
        return self.list_of_reaction_no_sources_multiplier[index]

    def return_reaction_no_losses_multiplier(self, index: int) -> int:
        return self.list_of_reaction_no_losses_multiplier[index]

    def return_formula(self) -> str:
        return self.formula

    def return_density(self) -> float:
        return self.density

    def return_loss(self) -> float:
        return self.loss

    def return_source(self) -> float:
        return self.source

    def process_reduce_lists_of_reaction_no(self) -> None:
        no_reactions_sources = self.list_of_reaction_no_sources[0]
        no_reactions_losses = self.list_of_reaction_no_losses[0]

        for is_idx in range(1, no_reactions_sources + 1):
            source_reaction_no = self.list_of_reaction_no_sources[is_idx]

            for il_idx in range(1, no_reactions_losses + 1):
                loss_reaction_no = self.list_of_reaction_no_losses[il_idx]

                if source_reaction_no == loss_reaction_no:
                    source_multiplier = self.list_of_reaction_no_sources_multiplier[is_idx]
                    loss_multiplier = self.list_of_reaction_no_losses_multiplier[il_idx]

                    if source_multiplier == loss_multiplier:
                        self.list_of_reaction_no_sources_multiplier[is_idx] = 0
                        self.list_of_reaction_no_losses_multiplier[il_idx] = 0

                    if source_multiplier > loss_multiplier:
                        self.list_of_reaction_no_sources_multiplier[is_idx] = (
                            source_multiplier - loss_multiplier
                        )
                        self.list_of_reaction_no_losses_multiplier[il_idx] = 0

                    if source_multiplier < loss_multiplier:
                        self.list_of_reaction_no_sources_multiplier[is_idx] = 0
                        self.list_of_reaction_no_losses_multiplier[il_idx] = (
                            loss_multiplier - source_multiplier
                        )

        size = self._MAX_NO_REACTIONS + 1
        list_reaction_no = [0] * size
        list_multiplier = [0] * size

        # Tidy up sources.
        index = 0
        for is_idx in range(1, no_reactions_sources + 1):
            source_multiplier = self.list_of_reaction_no_sources_multiplier[is_idx]
            if source_multiplier > 0:
                index += 1
                list_reaction_no[index] = self.list_of_reaction_no_sources[is_idx]
                list_multiplier[index] = source_multiplier
        list_reaction_no[0] = index

        for is_idx in range(0, size):
            self.list_of_reaction_no_sources[is_idx] = list_reaction_no[is_idx]
            self.list_of_reaction_no_sources_multiplier[is_idx] = list_multiplier[is_idx]
            list_reaction_no[is_idx] = 0
            list_multiplier[is_idx] = 0

        # Tidy up losses.
        index = 0
        for il_idx in range(1, no_reactions_losses + 1):
            loss_multiplier = self.list_of_reaction_no_losses_multiplier[il_idx]
            if loss_multiplier > 0:
                index += 1
                list_reaction_no[index] = self.list_of_reaction_no_losses[il_idx]
                list_multiplier[index] = loss_multiplier
        list_reaction_no[0] = index

        for il_idx in range(0, size):
            self.list_of_reaction_no_losses[il_idx] = list_reaction_no[il_idx]
            self.list_of_reaction_no_losses_multiplier[il_idx] = list_multiplier[il_idx]
