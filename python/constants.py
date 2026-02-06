"""Core constants for the airGM Python port.

Values mirror the C++ implementation in src/globalmodel.h and related files.
"""

EV_TO_KELVIN = 11605.0

NO_SPECIES = 53
NO_REACTIONS = 624

# Species indices that contain hydrogen, used for [H2O] == 0 shortcut logic.
SPECIES_CONTAINING_H = (
    11, 12, 13, 14, 15, 16, 26, 27, 32, 44, 45, 46, 47, 48, 49, 50,
)

SPECIES_INDEX = {
    "M": 0,
    "N+": 1,
    "N2+": 2,
    "N3+": 3,
    "N4+": 4,
    "O+": 5,
    "O2+": 6,
    "O4+": 7,
    "NO+": 8,
    "N2O+": 9,
    "NO2+": 10,
    "H+": 11,
    "H2+": 12,
    "H3+": 13,
    "OH+": 14,
    "H2O+": 15,
    "H3O+": 16,
    "e": 17,
    "O-": 18,
    "O2-": 19,
    "O3-": 20,
    "O4-": 21,
    "NO-": 22,
    "N2O-": 23,
    "NO2-": 24,
    "NO3-": 25,
    "H-": 26,
    "OH-": 27,
    "N(2_D)": 28,
    "N2(A_3_Sigma)": 29,
    "N2(B_3_Pi)": 30,
    "O(1_D)": 31,
    "H": 32,
    "N": 33,
    "O": 34,
    "O2(a_1_Delta)": 35,
    "O3": 36,
    "NO": 37,
    "N2O": 38,
    "NO2": 39,
    "NO3": 40,
    "N2O3": 41,
    "N2O4": 42,
    "N2O5": 43,
    "H2": 44,
    "OH": 45,
    "HO2": 46,
    "H2O2": 47,
    "HNO": 48,
    "HNO2": 49,
    "HNO3": 50,
    "N2": 51,
    "O2": 52,
    "H2O": 53,
}

SPECIES_FORMULAS = {index: formula for formula, index in SPECIES_INDEX.items()}
