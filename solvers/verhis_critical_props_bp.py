# D:\ChemE_Calc_New\solvers\critical_props_bp.py
import os
import math
import datetime
from flask import Blueprint, render_template, request, jsonify, url_for

# --- RDKit Import (Keep it here as it's specific to Joback in this BP) ---
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    print("WARNING: RDKit not found (critical_props_bp). Joback SMILES processing will not work.")
    RDKIT_AVAILABLE = False
    Chem = None

# --- Blueprint Definition ---
# The first argument 'critical_props' is the Blueprint's name,
# used in url_for() calls (e.g., url_for('critical_props.calculator_joback'))
critical_props_bp = Blueprint(
    'critical_props', 
    __name__,
    template_folder='../templates', # Point to the main templates folder
    static_folder='../static'       # Point to the main static folder (if needed by BP-specific templates)
)

# --- JOBACK METHOD DATA (Moved from app.py) ---
JOBACK_GROUPS_DATA = {
    # Group Name: {"tfpk":..., "tbk":..., "tck": ..., "pck": ..., "vck": ..., "hfk":..., "gfk":..., ... "smarts": "SMARTS_PATTERN"}
    # Values from Table C-1. SMARTS are best-effort and will likely need refinement and testing.
    # Note on units: hfk (kJ/mol), gfk (cal/mol), hvk (cal/mol), hmk (cal/mol). Cp coefficients for J/mol-K.

    # --- Aliphatic Hydrocarbons ---
    "-CH3": { # CH3 (1)
        "tfpk": -5.10, "tbk": 23.58, "tck": 0.0141, "pck": -0.0012, "vck": 65.0,
        "hfk": -76.45, "gfk": -43.96, "hvk": 567.0, "hmk": 217.0,
        "CpAk": 19.500, "CpBk": -8.08E-03, "CpCk": 1.53E-04, "CpDk": -9.67E-08,
        "smarts": "[CH3X4;!$(*~[#7,#8,#9,#16,#17,#35,#53]);!R]" # Aliphatic CH3
    },
    "-CH2-": { # CH2 (2)
        "tfpk": 11.27, "tbk": 22.88, "tck": 0.0189, "pck": 0.0000, "vck": 56.0,
        "hfk": -20.64, "gfk": 8.42, "hvk": 532.0, "hmk": 619.0,
        "CpAk": -0.909, "CpBk": 9.50E-02, "CpCk": -5.44E-05, "CpDk": 1.19E-08,
        "smarts": "[CH2X4;!R;!$(C=[O,S,N])]" # Aliphatic CH2, not in ring, not adjacent to C=X
    },
    ">CH-": { # CH (3)
        "tfpk": 12.64, "tbk": 21.74, "tck": 0.0164, "pck": 0.0020, "vck": 41.0,
        "hfk": 29.89, "gfk": 58.36, "hvk": 404.0, "hmk": 179.0,
        "CpAk": -23.000, "CpBk": 2.04E-01, "CpCk": -2.65E-04, "CpDk": 1.20E-07,
        "smarts": "[CHX4;!R](~[!#1])(~[!#1])~[!#1]" # Aliphatic CH, not in ring
    },
    ">C<": { # C (4)
        "tfpk": 46.43, "tbk": 18.25, "tck": 0.0067, "pck": 0.0043, "vck": 27.0,
        "hfk": 82.23, "gfk": 116.02, "hvk": 152.0, "hmk": -349.0,
        "CpAk": -66.200, "CpBk": 4.27E-01, "CpCk": -6.41E-04, "CpDk": 3.01E-07,
        "smarts": "[CX4;!R](~[!#1])(~[!#1])(~[!#1])~[!#1]" # Aliphatic C, not in ring
    },

    # --- Olefinic Hydrocarbons ---
    "=CH2": { # =CH2 (1)
        "tfpk": -4.32, "tbk": 18.18, "tck": 0.0113, "pck": -0.0028, "vck": 56.0,
        "hfk": -9.63, "gfk": 3.77, "hvk": 412.0, "hmk": -113.0,
        "CpAk": -23.600, "CpBk": -3.81E-02, "CpCk": 1.72E-04, "CpDk": -1.03E-07,
        "smarts": "[CH2X3R0]=[#6]" # Terminal vinyl CH2=C
    },
    "=CH- (olefinic)": { # =CH (2)
        "tfpk": 8.73, "tbk": 24.96, "tck": 0.0129, "pck": -0.0006, "vck": 46.0,
        "hfk": 37.97, "gfk": 48.53, "hvk": 527.0, "hmk": 643.0,
        "CpAk": -8.000, "CpBk": 1.05E-01, "CpCk": -9.63E-05, "CpDk": 3.56E-08,
        "smarts": "[CHX3R0](=[#6])-[#6;!$(C=O)]" # Internal olefinic =CH-
    },
    "=C< (olefinic)": { # =C (3)
        "tfpk": 11.14, "tbk": 24.14, "tck": 0.0117, "pck": 0.0011, "vck": 38.0,
        "hfk": 83.99, "gfk": 92.36, "hvk": 511.0, "hmk": 732.0,
        "CpAk": -28.100, "CpBk": 2.08E-01, "CpCk": -3.06E-04, "CpDk": 1.46E-07,
        "smarts": "[CR0X3](=[#6])([#6;!$(C=O)])-[#6;!$(C=O)]" # Trisubstituted olefinic =C<
    },
    "=C= (allenic)": { # =C(d) (2)
        "tfpk": 17.78, "tbk": 26.15, "tck": 0.0026, "pck": 0.0028, "vck": 36.0,
        "hfk": 142.14, "gfk": 136.70, "hvk": 636.0, "hmk": 1128.0,
        "CpAk": 27.400, "CpBk": -5.57E-02, "CpCk": 1.01E-04, "CpDk": -5.02E-08,
        "smarts": "[CX2R0](=[#6])=[#6]" # Central carbon of an allene >C=C=C<
    },

    # --- Acetylenic Hydrocarbons ---
    "#CH": { # =-CH (1) (Triple bond CH)
        "tfpk": -11.18, "tbk": 9.20, "tck": 0.0027, "pck": -0.0008, "vck": 46.0,
        "hfk": 79.30, "gfk": 77.71, "hvk": 276.0, "hmk": 555.0,
        "CpAk": 24.500, "CpBk": -2.71E-02, "CpCk": 1.11E-04, "CpDk": -6.78E-08,
        "smarts": "[CHX2R0]#[#6]" # Terminal acetylenic CH
    },
    "#C-": { # =-C- (2) (Triple bond C)
        "tfpk": 64.32, "tbk": 27.38, "tck": 0.0020, "pck": 0.0016, "vck": 37.0,
        "hfk": 115.51, "gfk": 109.82, "hvk": 789.0, "hmk": 992.0,
        "CpAk": 7.870, "CpBk": 2.01E-02, "CpCk": -8.33E-06, "CpDk": 1.39E-09,
        "smarts": "[CX2R0](#[#6])-[#6]" # Internal acetylenic C
    },

    # --- Non-Aromatic Ring Carbons ---
    "-CH2- (ring)": { # -CH2(ss) (2)
        "tfpk": 7.75, "tbk": 27.15, "tck": 0.0100, "pck": 0.0025, "vck": 48.0,
        "hfk": -26.80, "gfk": -3.68, "hvk": 573.0, "hmk": 117.0,
        "CpAk": -6.030, "CpBk": 8.54E-02, "CpCk": -8.00E-06, "CpDk": 1.80E-08, # Corrected CpCk sign
        "smarts": "[CH2X4R1]" # CH2 in a non-aromatic ring
    },
    ">CH- (ring)": { # -CH(ss) (3)
        "tfpk": 19.88, "tbk": 21.78, "tck": 0.0122, "pck": 0.0004, "vck": 38.0,
        "hfk": 8.67, "gfk": 40.99, "hvk": 464.0, "hmk": 775.0,
        "CpAk": 8.670, "CpBk": 1.62E-01, "CpCk": -1.60E-04, "CpDk": 6.24E-08,
        "smarts": "[CHX4R1](~[#6R1])~[#6R1]" # CH in a non-aromatic ring, bonded to two other ring atoms
    },
    ">C< (ring)": { # -C(ss) (4)
        "tfpk": 60.15, "tbk": 21.32, "tck": 0.0042, "pck": 0.0061, "vck": 27.0,
        "hfk": 79.72, "gfk": 87.88, "hvk": 154.0, "hmk": -328.0,
        "CpAk": -90.900, "CpBk": 5.57E-01, "CpCk": -9.00E-04, "CpDk": 4.69E-07,
        "smarts": "[CX4R1](~[#6R1])(~[#6R1])~[#6R1]" # Quaternary C in a non-aromatic ring
    },
    "=CH- (ring)": { # =CH(ds) (2) ; 'ds' in Joback implies aromatic/conjugated-like
        "tfpk": 8.13, "tbk": 26.73, "tck": 0.0082, "pck": 0.0011, "vck": 41.0,
        "hfk": 2.09, "gfk": 11.30, "hvk": 608.0, "hmk": 263.0,
        "CpAk": -2.140, "CpBk": 5.74E-02, "CpCk": -1.64E-06, "CpDk": 1.59E-08, # Corrected CpCk sign
        "smarts": "[cH]" # Aromatic CH or CH in conjugated ring system
    },
    "=C< (ring)": { # =C(ds) (3)
        "tfpk": 37.02, "tbk": 31.01, "tck": 0.0143, "pck": 0.0008, "vck": 32.0,
        "hfk": 46.43, "gfk": 54.05, "hvk": 731.0, "hmk": 572.0,
        "CpAk": -8.250, "CpBk": 1.01E-01, "CpCk": -1.42E-04, "CpDk": 6.78E-08,
        "smarts": "[c](:[#6]):[#6]" # Aromatic C bonded to two other aromatic atoms
    },

    # --- Halogens ---
    "-F": { # F (1)
        "tfpk": -15.78, "tbk": -0.03, "tck": 0.0111, "pck": -0.0057, "vck": 27.0,
        "hfk": -251.92, "gfk": -247.19, "hvk": -160.0, "hmk": 334.0,
        "CpAk": 26.500, "CpBk": -9.13E-02, "CpCk": 1.91E-04, "CpDk": -1.03E-07,
        "smarts": "[F!H0X1]"
    },
    "-Cl": { # Cl (1)
        "tfpk": 13.55, "tbk": 38.13, "tck": 0.0105, "pck": -0.0049, "vck": 58.0,
        "hfk": -71.55, "gfk": -64.31, "hvk": 1083.0, "hmk": 601.0,
        "CpAk": 33.300, "CpBk": -9.63E-02, "CpCk": 1.87E-04, "CpDk": -9.96E-08,
        "smarts": "[Cl!H0X1]"
    },
    "-Br": { # Br (1)
        "tfpk": 43.43, "tbk": 66.86, "tck": 0.0133, "pck": 0.0057, "vck": 71.0, # pck is 0.0057 not -0.0057 from image
        "hfk": -29.48, "gfk": -38.06, "hvk": 1573.0, "hmk": 861.0,
        "CpAk": 28.600, "CpBk": -6.49E-02, "CpCk": 1.36E-04, "CpDk": -7.45E-08,
        "smarts": "[Br!H0X1]"
    },
    "-I": { # I (1)*
        "tfpk": 41.69, "tbk": 93.84, "tck": 0.0068, "pck": -0.0034, "vck": 97.0,
        "hfk": 21.06, "gfk": 5.74, "hvk": 2275.0, "hmk": 651.0,
        "CpAk": 32.100, "CpBk": -6.41E-02, "CpCk": 1.26E-04, "CpDk": -6.87E-08,
        "smarts": "[I!H0X1]"
    },

    # --- Oxygen Compounds ---
    "-OH (alcohol)": { # -OH (1)
        "tfpk": 44.45, "tbk": 92.88, "tck": 0.0741, "pck": 0.0112, "vck": 28.0,
        "hfk": -208.04, "gfk": -189.20, "hvk": 4021.0, "hmk": 575.0,
        "CpAk": 25.700, "CpBk": -6.91E-02, "CpCk": 1.77E-04, "CpDk": -9.88E-08,
        "smarts": "[OX2H1][#6X4;!$(C=O);!c]" # Alcohol OH on sp3 C, not carbonyl, not aromatic
    },
    "-OH (phenol)": { # ACOH (1) - Note: Table ACOH is likely Phenol for Joback context
        "tfpk": 82.83, "tbk": 76.34, "tck": 0.0240, "pck": 0.0184, "vck": -25.0, # vck is -25
        "hfk": -221.65, "gfk": -197.37, "hvk": 2987.0, "hmk": 1073.0,
        "CpAk": -2.810, "CpBk": 1.11E-01, "CpCk": -1.16E-04, "CpDk": 4.94E-08,
        "smarts": "[OX2H1]c" # OH directly on an aromatic ring
    },
    "-O- (non-ring ether)": { # -O- (2)
        "tfpk": 22.23, "tbk": 22.42, "tck": 0.0168, "pck": 0.0015, "vck": 18.0,
        "hfk": -132.22, "gfk": -105.00, "hvk": 576.0, "hmk": 284.0,
        "CpAk": 25.500, "CpBk": -6.32E-02, "CpCk": 1.11E-04, "CpDk": -5.48E-08,
        "smarts": "[OX2R0]([#6;!$(C=O);!c])([#6;!$(C=O);!c])" # Non-ring ether, not next to C=O or aromatic
    },
    "-O- (ring ether)": { # -O(ss) (2)
        "tfpk": 23.05, "tbk": 31.22, "tck": 0.0098, "pck": 0.0048, "vck": 13.0,
        "hfk": -138.16, "gfk": -98.22, "hvk": 1119.0, "hmk": 1405.0,
        "CpAk": 12.200, "CpBk": -1.26E-02, "CpCk": 6.03E-05, "CpDk": -3.86E-08,
        "smarts": "[OX2R1]" # Oxygen in a non-aromatic ring
    },
    ">C=O (ketone)": { # >C=O (2)
        "tfpk": 61.20, "tbk": 76.75, "tck": 0.0380, "pck": 0.0031, "vck": 62.0,
        "hfk": -133.22, "gfk": -120.50, "hvk": 2144.0, "hmk": 1001.0,
        "CpAk": 6.450, "CpBk": 6.70E-02, "CpCk": -3.57E-05, "CpDk": 2.86E-09,
        "smarts": "[CX3R0](=O)([#6;!c])([#6;!c])" # Ketone C=O, non-ring
    },
    ">C=O (ring ketone)": { # >C=O(ss) (2)
        "tfpk": 75.97, "tbk": 94.97, "tck": 0.0284, "pck": 0.0028, "vck": 55.0,
        "hfk": -164.50, "gfk": -126.27, "hvk": 1588.0, "hmk": None, # 'X'
        "CpAk": 30.400, "CpBk": -8.29E-02, "CpCk": 2.36E-04, "CpDk": -1.31E-07,
        "smarts": "[CX3R1](=O)" # Ketone C=O in a non-aromatic ring
    },
    "-CHO (aldehyde)": { # -CH=O (1)*
        "tfpk": 36.90, "tbk": 72.20, "tck": 0.0379, "pck": 0.0030, "vck": 82.0,
        "hfk": -162.03, "gfk": -143.48, "hvk": 2173.0, "hmk": 764.0,
        "CpAk": 30.900, "CpBk": -3.36E-02, "CpCk": 1.60E-04, "CpDk": -9.88E-08,
        "smarts": "[CX3H1R0](=O)[#6;!c]" # Aldehyde CHO
    },
    "-COOH (acid)": { # -COOH (1)
        "tfpk": 155.50, "tbk": 169.09, "tck": 0.0791, "pck": 0.0077, "vck": 89.0,
        "hfk": -426.72, "gfk": -387.87, "hvk": 4669.0, "hmk": 2641.0,
        "CpAk": 24.100, "CpBk": 4.27E-02, "CpCk": 8.04E-05, "CpDk": -6.87E-08,
        "smarts": "[CX3](=[OX1])[OX2H1]" # Carboxylic acid group
    },
    "-COO- (ester)": { # -COO- (2)
        "tfpk": 53.60, "tbk": 81.10, "tck": 0.0481, "pck": 0.0005, "vck": 82.0, # vck is 82
        "hfk": -337.92, "gfk": -301.95, "hvk": 2302.0, "hmk": 1663.0,
        "CpAk": 24.500, "CpBk": 4.02E-02, "CpCk": -4.02E-05, "CpDk": 4.52E-08, # Corrected CpCk sign
        "smarts": "[CX3R0](=[OX1])[OX2R0][#6]" # Ester group -C(=O)O-C
    },
    "=O (other than above)": { # =O (1)* (e.g. from SO2, NO2 but Joback usually has specific groups for those)
        "tfpk": 2.08, "tbk": -10.50, "tck": 0.0143, "pck": 0.0101, "vck": 36.0,
        "hfk": -247.61, "gfk": -250.83, "hvk": 1412.0, "hmk": 866.0,
        "CpAk": 6.820, "CpBk": 1.96E-02, "CpCk": -1.27E-05, "CpDk": 1.78E-08, # Corrected CpCk sign
        "smarts": "[OX1]=[#6;!$(C=O)]" # Double bonded O, not part of standard carbonyl or carboxyl. This is very general.
    },

    # --- Nitrogen Compounds ---
    "-NH2": { # -NH2 (1)
        "tfpk": 66.89, "tbk": 73.23, "tck": 0.0243, "pck": 0.0109, "vck": 38.0,
        "hfk": -22.02, "gfk": 14.07, "hvk": 2578.0, "hmk": 840.0,
        "CpAk": 26.900, "CpBk": 4.12E-02, "CpCk": 1.64E-04, "CpDk": -9.76E-08,
        "smarts": "[NH2X2][#6;!$(C=O);!c]" # Primary amine -NH2
    },
    "-NH- (non-ring)": { # -NH- (2)
        "tfpk": 52.66, "tbk": 50.17, "tck": 0.0295, "pck": 0.0077, "vck": 35.0,
        "hfk": 53.47, "gfk": 89.39, "hvk": 1538.0, "hmk": 1197.0,
        "CpAk": -1.210, "CpBk": 7.62E-02, "CpCk": -4.86E-05, "CpDk": 1.05E-08,
        "smarts": "[NHX3R0]([#6;!$(C=O);!c])([#6;!$(C=O);!c])" # Secondary non-ring amine
    },
    "-NH- (ring)": { # -NH(ss)- (2)
        "tfpk": 101.51, "tbk": 52.82, "tck": 0.0130, "pck": 0.0114, "vck": 29.0,
        "hfk": 31.65, "gfk": 75.61, "hvk": 1656.0, "hmk": 1790.0,
        "CpAk": 11.800, "CpBk": -2.30E-02, "CpCk": 1.07E-04, "CpDk": -6.28E-08,
        "smarts": "[NHX3R1]" # Secondary ring amine
    },
    ">N- (non-ring)": { # >N- (3)
        "tfpk": 48.84, "tbk": 11.74, "tck": 0.0169, "pck": 0.0074, "vck": 9.0, # vck is 9
        "hfk": 123.34, "gfk": 163.16, "hvk": 453.0, "hmk": 1124.0,
        "CpAk": -31.100, "CpBk": 2.27E-01, "CpCk": -3.20E-04, "CpDk": 1.46E-07,
        "smarts": "[NX4R0]([#6;!$(C=O);!c])([#6;!$(C=O);!c])([#6;!$(C=O);!c])" # Tertiary non-ring amine
    },
    "=N- (imine/schiff base)": { # =N- (2)
        "tfpk": None, "tbk": 74.60, "tck": 0.0255, "pck": -0.0099, "vck": None, # tfpk, vck are X
        "hfk": 23.61, "gfk": None, "hvk": 797.0, "hmk": None, # gfk, hmk are X
        "CpAk": None, "CpBk": None, "CpCk": None, "CpDk": None, # All Cp are X
        "smarts": "[NX2R0]=[#6]" # Imine =N- (bonded to carbon)
    },
    "=N- (aromatic N)": { # =N(ds)= (2) ; ds here implies part of aromatic system like pyridine
        "tfpk": 68.40, "tbk": 57.55, "tck": 0.0085, "pck": 0.0076, "vck": 34.0,
        "hfk": 55.52, "gfk": 79.93, "hvk": 1560.0, "hmk": 872.0,
        "CpAk": 8.830, "CpBk": -3.84E-03, "CpCk": 4.35E-05, "CpDk": -2.60E-08,
        "smarts": "n" # Aromatic nitrogen (e.g., in pyridine)
    },
    "-CN (nitrile)": { # -CN (1)*
        "tfpk": 59.89, "tbk": 125.66, "tck": 0.0496, "pck": -0.0101, "vck": 91.0,
        "hfk": 88.43, "gfk": 89.22, "hvk": 3071.0, "hmk": 577.0,
        "CpAk": 36.500, "CpBk": -7.33E-02, "CpCk": 1.84E-04, "CpDk": -1.03E-07,
        "smarts": "[CX1]#[NX1]" # Cyano group -C#N
    },
    "-NO2 (nitro)": { # -NO2 (1)
        "tfpk": 127.24, "tbk": 152.54, "tck": 0.0437, "pck": 0.0064, "vck": 91.0,
        "hfk": -66.57, "gfk": -16.83, "hvk": 4000.0, "hmk": 2313.0,
        "CpAk": 25.900, "CpBk": -3.74E-03, "CpCk": 1.29E-04, "CpDk": -8.88E-08,
        "smarts": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]" # Nitro group -NO2
    },
    ">NH (special single-bonded N)": { # -NH (1)* in table
        "tfpk": None, "tbk": None, "tck": None, "pck": None, "vck": None,
        "hfk": 93.70, "gfk": 119.66, "hvk": 2908.0, "hmk": None,
        "CpAk": 5.690, "CpBk": -4.12E-03, "CpCk": 1.28E-04, "CpDk": -8.88E-08,
        "smarts": "[NH1X2]([#6;!$(C=O)])" # NH attached to one non-carbonyl carbon.
    },

    # --- Sulfur Compounds ---
    "-SH (thiol)": { # -SH (1)
        "tfpk": 20.09, "tbk": 63.56, "tck": 0.0031, "pck": 0.0084, "vck": 63.0,
        "hfk": -17.33, "gfk": -22.99, "hvk": 1645.0, "hmk": 564.0,
        "CpAk": 35.300, "CpBk": -7.58E-02, "CpCk": 1.85E-04, "CpDk": -1.03E-07,
        "smarts": "[SH1X2][#6;!$(C=O);!c]" # Thiol -SH
    },
    "-S- (non-ring thioether)": { # -S- (2)
        "tfpk": 34.40, "tbk": 68.78, "tck": 0.0119, "pck": 0.0049, "vck": 54.0,
        "hfk": 41.87, "gfk": 33.12, "hvk": 1629.0, "hmk": 987.0,
        "CpAk": 19.600, "CpBk": -5.61E-03, "CpCk": 4.02E-05, "CpDk": -2.76E-08,
        "smarts": "[SX2R0]([#6;!$(C=O);!c])([#6;!$(C=O);!c])" # Non-ring thioether -S-
    },
    "-S- (ring thioether)": { # -S(ss)- (2)
        "tfpk": 79.93, "tbk": 52.10, "tck": 0.0019, "pck": 0.0051, "vck": 38.0,
        "hfk": 39.10, "gfk": 27.76, "hvk": 1430.0, "hmk": 372.0,
        "CpAk": 16.700, "CpBk": -4.81E-03, "CpCk": 2.77E-05, "CpDk": -2.11E-08,
        "smarts": "[SX2R1]" # Ring thioether -S-
    }
}
JOBACK_GROUP_NAMES_LIST = sorted(list(JOBACK_GROUPS_DATA.keys())) # Keep this updated

# --- JOBACK METHOD ROUTES (within critical_props_bp) ---
@critical_props_bp.route('/critical-properties/joback')
def calculator_joback_page(): # Renamed to avoid conflict with any app-level 'calculator_joback'
    return render_template('calculator_joback.html', JOBACK_GROUPS_LIST=JOBACK_GROUP_NAMES_LIST)

@critical_props_bp.route('/calculate_joback_critical', methods=['GET', 'POST']) # Allow GET for initial data
def calculate_joback_critical_api():
    # --- CITATION DATA (defined within the route to use url_for correctly after BP init) ---
    paper_citation_joback_data = {
        'authors': 'Joback, K. G., & Reid, R. C.', 'year': 1987,
        'title': 'Estimation of pure-component properties from group-contributions',
        'journal': 'Chemical Engineering Communications', 'volume': '57', 
        'issue': '1-6', 'pages': '233-243', 'doi': '10.1080/00986448708960487'
    }
    today = datetime.date.today()
    site_url = request.host_url.rstrip('/')
    # Use blueprint name in url_for: 'critical_props.calculator_joback_page'
    joback_page_url = site_url + url_for('critical_props.calculator_joback_page') 
    
    tool_citation_joback_data = {
        'author': 'ChemE Calc', 'year': today.year,
        'title': 'Joback Method Critical Properties Calculator [Web Application]',
        'retrieved_date': today.strftime('%Y-%m-%d'),
        'url': joback_page_url
    }
    applicability_joback = "Valid for organic compounds. Accuracy varies; generally better for non-polar molecules. Not recommended for very small molecules or complex structures without careful group definition. See Joback (1987)."

    if request.method == 'GET':
        # This part is for the JavaScript to fetch initial citation/applicability data
        return jsonify({
            'paper_citation_data': paper_citation_joback_data,
            'tool_citation_data': tool_citation_joback_data,
            'applicability_notes': applicability_joback
        })

    # --- POST request handling (actual calculation) ---
    print("--- JOBACK API (critical_props_bp) POST START ---")
    if not RDKIT_AVAILABLE:
        return jsonify({
            "error": "RDKit library not found on server. SMILES processing disabled.",
            'paper_citation_data': paper_citation_joback_data,
            'tool_citation_data': tool_citation_joback_data,
            'applicability_notes': applicability_joback
        }), 500

    data = request.get_json() or {}
    tb_str = data.get('Tb')
    smiles_str = data.get('smiles')

    error_payload = { # For returning consistent error structure
        'paper_citation_data': paper_citation_joback_data,
        'tool_citation_data': tool_citation_joback_data,
        'applicability_notes': applicability_joback
    }

    if not tb_str or not smiles_str:
        error_payload["error"] = "Normal Boiling Point (Tb) and SMILES string are required."
        return jsonify(error_payload), 400
    try:
        tb = float(tb_str)
        if tb <= 0: raise ValueError("Tb must be positive.")
    except ValueError:
        error_payload["error"] = "Normal Boiling Point (Tb) must be a valid positive number."
        return jsonify(error_payload), 400

    mol = Chem.MolFromSmiles(smiles_str.strip())
    if mol is None:
        error_payload["error"] = f"Invalid SMILES string: '{smiles_str}'. Please check the format."
        return jsonify(error_payload), 400
    
    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()

    group_counts = {}
    identified_groups_latex_parts = []

    for group_name, group_data in JOBACK_GROUPS_DATA.items():
        if "smarts" in group_data and group_data["smarts"]:
            try:
                pattern = Chem.MolFromSmarts(group_data["smarts"])
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    count = len(matches)
                    if count > 0:
                        group_counts[group_name] = count
                        latex_group_name = group_name.replace('-', '\\text{-}').replace('>', '\\text{>}')
                        identified_groups_latex_parts.append(f"\\text{{{latex_group_name}}}: {count}")
            except Exception as e:
                print(f"  Error processing SMARTS for group '{group_name}' ({group_data['smarts']}): {e}")
    
    sum_tck = sum(JOBACK_GROUPS_DATA[k]["tck"] * v for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]["tck"] is not None)
    sum_pck = sum(JOBACK_GROUPS_DATA[k]["pck"] * v for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]["pck"] is not None)
    sum_vck = sum(JOBACK_GROUPS_DATA[k]["vck"] * v for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]["vck"] is not None)

    tc, pc, vc = float('nan'), float('nan'), float('nan')
    try:
        tc_denominator = 0.584 + 0.965 * sum_tck - (sum_tck**2)
        tc = tb / tc_denominator if tc_denominator != 0 else float('inf')

        pc_base_val = (0.113 + 0.0032 * n_atoms - sum_pck)
        if pc_base_val == 0: pc = float('inf')
        elif isinstance(pc_base_val, complex) or (pc_base_val < 0 and (-2 % 1 != 0)): pc = float('nan')
        else: pc = pc_base_val**-2
        
        vc = 17.5 + sum_vck
    except Exception as e:
        error_payload["error"] = f"Calculation error: {str(e)}"
        return jsonify(error_payload), 400 # Or 500 for server-side calc error

    # --- Build LaTeX String (Same as before) ---
    latex_steps = rf"$$\text{{Input SMILES: \texttt{{{smiles_str.replace(' ', '').replace('\\', '\\\\')}}}}}$$"
    latex_steps += rf"$$\text{{Total Atoms (N}}_{{\text{{atoms}}}}\text{{) (after adding Hydrogens): {n_atoms}}}$$"

    if identified_groups_latex_parts:
        latex_steps += r"$$\text{\textbf{Identified Joback Groups:}}$$"
        group_list_str = r" \\ ".join(identified_groups_latex_parts) 
        latex_steps += f"$${group_list_str}$$"
    else:
        latex_steps += r"$$\text{(No specific Joback groups identified.)}$$"

    sum_tck_terms_latex = ' + '.join([f"({v} \\times {JOBACK_GROUPS_DATA[k]['tck']:.4f})" for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]['tck'] is not None]) if group_counts else "0"
    sum_pck_terms_latex = ' + '.join([f"({v} \\times {JOBACK_GROUPS_DATA[k]['pck']:.4f})" for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]['pck'] is not None]) if group_counts else "0"
    sum_vck_terms_latex = ' + '.join([f"({v} \\times {JOBACK_GROUPS_DATA[k]['vck']:.1f})" for k, v in group_counts.items() if k in JOBACK_GROUPS_DATA and JOBACK_GROUPS_DATA[k]['vck'] is not None]) if group_counts else "0"

    latex_steps += r"$$\text{\textbf{1. Joback Group Contributions Summation:}}$$"
    latex_steps += r"$$\begin{align*}"
    latex_steps += rf"\sum (N_k \cdot t_{{ck}}) &= {sum_tck_terms_latex} &&= {sum_tck:.4f} \\"
    latex_steps += rf"\sum (N_k \cdot p_{{ck}}) &= {sum_pck_terms_latex} &&= {sum_pck:.4f} \\"
    latex_steps += rf"\sum (N_k \cdot v_{{ck}}) &= {sum_vck_terms_latex} &&= {sum_vck:.1f}"
    latex_steps += r"\end{align*}$$"

    latex_steps += r"$$\text{\textbf{2. Critical Temperature (T}_c\text{):}}$$"
    latex_steps += r"$$T_c = T_b \left[ 0.584 + 0.965 \left( \sum N_k t_{ck} \right) - \left( \sum N_k t_{ck} \right)^2 \right]^{-1}$$"
    latex_steps += rf"$$T_c = {tb:.2f} \, \text{{K}} \left[ 0.584 + 0.965 \cdot ({sum_tck:.4f}) - ({sum_tck:.4f})^2 \right]^{{-1}} = {tc:.2f} \, \text{{K}}$$ "

    latex_steps += r"$$\text{\textbf{3. Critical Pressure (P}_c\text{):}}$$"
    latex_steps += r"$$P_c = \left[ 0.113 + 0.0032 N_{{\text{{atoms}}}} - \sum N_k p_{ck} \right]^{-2}$$"
    latex_steps += rf"$$P_c = \left[ 0.113 + 0.0032 \cdot {n_atoms} - ({sum_pck:.4f}) \right]^{{-2}} = {pc:.2f} \, \text{{bar}}$$ "

    latex_steps += r"$$\text{\textbf{4. Critical Volume (V}_c\text{):}}$$"
    latex_steps += r"$$V_c = 17.5 + \sum N_k v_{ck}$$"
    latex_steps += rf"$$V_c = 17.5 + ({sum_vck:.1f}) = {vc:.1f} \, \text{{cm}}^3/\text{{mol}}$$ "
    
    final_latex_string = latex_steps

    print("--- JOBACK API (critical_props_bp) END ---")
    return jsonify({
        'Tc': f"{tc:.2f}" if not (math.isnan(tc) or math.isinf(tc)) else "Error",
        'Pc': f"{pc:.2f}" if not (math.isnan(pc) or math.isinf(pc)) else "Error",
        'Vc': f"{vc:.1f}" if not (math.isnan(vc) or math.isinf(vc)) else "Error",
        'latex_steps': final_latex_string,
        'paper_citation_data': paper_citation_joback_data,
        'tool_citation_data': tool_citation_joback_data,
        'applicability_notes': applicability_joback
    })

    # Add this inside solvers/critical_props_bp.py
# You can place it after the calculate_joback_critical_api function

@critical_props_bp.route('/critical-properties/joback/estimate_tb', methods=['POST'])
def estimate_joback_tb_api():
    print("--- JOBACK TB ESTIMATION API (critical_props_bp) START ---")
    # ... (RDKit check, data retrieval, SMILES parsing - same as before) ...
    if not RDKIT_AVAILABLE: # Ensure this check is early
        return jsonify({"error": "RDKit library not found on server. SMILES processing disabled."}), 500

    data = request.get_json() or {}
    smiles_str = data.get('smiles')

    if not smiles_str:
        return jsonify({"error": "SMILES string is required for Tb estimation."}), 400

    mol = Chem.MolFromSmiles(smiles_str.strip())
    if mol is None:
        return jsonify({"error": f"Invalid SMILES string: '{smiles_str}'. Please check the format."}), 400
    
    mol = Chem.AddHs(mol)

    group_counts = {}
    identified_groups_for_tb_display = {} # For LaTeX and JSON response
    tbk_contributions_latex_parts = []

    print("Estimating Tb: Identifying Joback groups using SMARTS...")
    for group_name, group_data in JOBACK_GROUPS_DATA.items():
        if "smarts" in group_data and group_data["smarts"] and group_data.get("tbk") is not None:
            try:
                pattern = Chem.MolFromSmarts(group_data["smarts"])
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    count = len(matches)
                    if count > 0:
                        group_counts[group_name] = count
                        # For LaTeX and JSON response, store group name, count, and tbk value
                        identified_groups_for_tb_display[group_name] = {
                            "count": count,
                            "tbk_value": group_data["tbk"]
                        }
                        # For LaTeX sum terms
                        tbk_contributions_latex_parts.append(f"({count} \\times {group_data['tbk']:.2f})")
                        print(f"  Found group for Tb: {group_name}, Count: {count}, tbk: {group_data['tbk']}")
            except Exception as e:
                print(f"  Error processing SMARTS for Tb estimation, group '{group_name}': {e}")
    
    if not group_counts:
        return jsonify({"error": "No relevant Joback groups for Tb estimation found for the given SMILES."}), 400

    sum_tbk = sum(JOBACK_GROUPS_DATA[k]["tbk"] * v for k, v in group_counts.items())
    estimated_tb_k = 198 + sum_tbk
    
    # --- Build LaTeX String for Tb Estimation ---
    tb_latex_steps = rf"$$\text{{\textbf{{Joback T_b Estimation for SMILES: \texttt{{{smiles_str.replace(' ', '').replace('\\', '\\\\')}}}}}}}$$"
    
    if identified_groups_for_tb_display:
        tb_latex_steps += r"$$\text{\textbf{Identified Groups & T_{bk} Contributions:}}$$"
        group_detail_latex = []
        for name, info in identified_groups_for_tb_display.items():
            escaped_name = name.replace('-', '\\text{-}').replace('>', '\\text{>}')
            group_detail_latex.append(rf"\text{{{escaped_name}}}: {info['count']} \text{{ groups, T}}_{{bk}} = {info['tbk_value']:.2f}")
        tb_latex_steps += f"$${' \\\\ '.join(group_detail_latex)}$$" # Each group on a new line
    
    sum_tbk_terms_latex = ' + '.join(tbk_contributions_latex_parts) if tbk_contributions_latex_parts else "0"

    tb_latex_steps += r"$$\text{\textbf{Calculation:}}$$"
    tb_latex_steps += r"$$T_b = 198 + \sum (N_k \cdot T_{bk})$$"
    tb_latex_steps += rf"$$T_b = 198 + ({sum_tbk_terms_latex})$$"
    tb_latex_steps += rf"$$T_b = 198 + {sum_tbk:.2f} = {estimated_tb_k:.2f} \, \text{{K}}$$ "
    
    print(f"Estimated Tb (K) = {estimated_tb_k:.2f}")
    print("--- JOBACK TB ESTIMATION API (critical_props_bp) END ---")
    
    return jsonify({
        'estimated_tb_k': round(estimated_tb_k, 2),
        'groups_found_for_tb': identified_groups_for_tb_display, # Send detailed info
        'sum_tbk_contributions': round(sum_tbk, 4),
        'tb_estimation_latex_steps': tb_latex_steps # NEW: Send LaTeX for Tb steps
    })

    # In solvers/critical_props_bp.py

# ... (existing imports, Blueprint definition, RDKit, JOBACK_GROUPS_DATA, etc.) ...

# Helper patterns for CG_FIRST_ORDER_GROUPS_DATA in critical_props_bp.py
C_QUAT_ALL_ALKYL_PATTERN_REF = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
CONHCH3_PATTERN_REF = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
CONME_ET_N_ME_PART_PATTERN_REF = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
CONME2_PATTERN_REF = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
CH2F_PATTERN_REF = "[CH2X4H2](F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,Cl,Br,I])]"
FCH2O_PATTERN_REF = "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]"
AEO_FOR_CH3_EXCL_SMARTS_REF = "[OX2H0;!R;D2]-[#6;!R;!c;!$([CX3](=[OX1]))]"
AEO_FOR_CH2_CH_EXCL_SMARTS_REF = "[OX2H0;!R;D2;!$(*-[c]);!$(*-[CX3](=[OX1]))]"
CH2NH2_EXCL_SIMPLE_PATTERN_REF = "[CH2X4H2]([NH2X3])"
CHNH2_EXCL_SIMPLE_PATTERN_REF = "[CH1X4H1]([NH2X3])"
CH2NH_SEC_EXCL_SIMPLE_PATTERN_REF = "[CH2X4H2]([NH1X3H1])"
CH_NH_SEC_EXCL_SIMPLE_PATTERN_REF = "[CH1X4H1]([NH1X3H1])"
CH2N_TERT_EXCL_SIMPLE_PATTERN_REF = "[CH2X4H2]([NX3H0])"
CH_N_TERT_EXCL_SIMPLE_PATTERN_REF = "[CH1X4H1]([NX3H0])"
CH3_N_TERT_EXCL_SIMPLE_PATTERN_REF= "[CH3X4]([NX3H0])"
CH3_NH_SEC_EXCL_SIMPLE_PATTERN_REF= "[CH3X4]([NH1X3H1])"

# Final General Aliphatic SMARTS strings
FINAL_CH3_GEN_SMARTS = f"[CH3X4;!R;!$(*~[a]);!$(*~[#8X2H0;R]);!$(*~[#16]);!$(*~[F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_PATTERN_REF});!$({CONME_ET_N_ME_PART_PATTERN_REF});!$({CONME2_PATTERN_REF});!$(*-{AEO_FOR_CH3_EXCL_SMARTS_REF});!$({CH3_NH_SEC_EXCL_SIMPLE_PATTERN_REF});!$({CH3_N_TERT_EXCL_SIMPLE_PATTERN_REF});!$(*~[#7])]"
FINAL_CH2_GEN_SMARTS = f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*-[c]);!$(*-[CX3;!H1]=[OX1]);!$(*-{AEO_FOR_CH2_CH_EXCL_SMARTS_REF});!$(*~[F,Cl,Br,I,Si,P]);!$({CH2F_PATTERN_REF});!$({FCH2O_PATTERN_REF});!$({CH2NH2_EXCL_SIMPLE_PATTERN_REF});!$({CH2NH_SEC_EXCL_SIMPLE_PATTERN_REF});!$({CH2N_TERT_EXCL_SIMPLE_PATTERN_REF});!$(*-[C]#[N])]"
FINAL_CH_GEN_SMARTS  = f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*-[c]);!$(*-[CX3;!H1]=[OX1]);!$(*-{AEO_FOR_CH2_CH_EXCL_SMARTS_REF});!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$({CHNH2_EXCL_SIMPLE_PATTERN_REF});!$({CH_NH_SEC_EXCL_SIMPLE_PATTERN_REF});!$({CH_N_TERT_EXCL_SIMPLE_PATTERN_REF})]"
FINAL_C_GEN_SMARTS   = f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*-[c]);!$(*-[CX3;!H1]=[OX1]);!$(*-{AEO_FOR_CH2_CH_EXCL_SMARTS_REF});!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$({C_QUAT_ALL_ALKYL_PATTERN_REF})]"

# --- CONSTANTINOU-GANI (CG) METHOD DATA ---

CG_FIRST_ORDER_GROUPS_DATA = {
    # From Page C.6
    "CH3": { # CH3 (1) from Table C-2
        "tfp1k": 0.4640, "tb1k": 0.8894, "tc1k": 1.6781, "pc1k": 0.0199, "vc1k": 0.0750, "w1k": 0.296,
        "hf1k": -45.947, "gf1k": -8.030, "hv1k": 4.116, "vliq1k": 0.0261,
        "CpA1k": 35.1152, "CpB1k": 39.5923, "CpC1k": -9.9232,
        "smarts": FINAL_CH3_GEN_SMARTS
    },
    "CH2": { # CH2 (2) from Table C-2
        "tfp1k": 0.9246, "tb1k": 0.9225, "tc1k": 3.4920, "pc1k": 0.0106, "vc1k": 0.0558, "w1k": 0.147,
        "hf1k": -20.763, "gf1k": 8.231, "hv1k": 4.650, "vliq1k": 0.0164,
        "CpA1k": 22.6346, "CpB1k": 45.0933, "CpC1k": -15.7033,
        "smarts": FINAL_CH2_GEN_SMARTS
    },
    "OH": { # OH (1) (Alcohol) from Table C-2
        "tfp1k": 3.5979, "tb1k": 3.2152, "tc1k": 9.7292, "pc1k": 0.0051, "vc1k": 0.0390, "w1k": 0.737,
        "hf1k": -181.422, "gf1k": -158.589, "hv1k": 24.529, "vliq1k": 0.0055,
        "CpA1k": 27.2107, "CpB1k": 2.7609, "CpC1k": 1.3060,
        "smarts": "[OH1X2][CX4;!c;!$(C=O)]" # Mature
    },
    "C": { # C (4) - Quaternary aliphatic carbon
        "tfp1k": 1.6479, "tb1k": 0.2878, "tc1k": 4.8823, "pc1k": -0.0104, "vc1k": -0.0003, "w1k": -0.351,
        "hf1k": 17.119, "gf1k": 37.977, "hv1k": 1.284, "vliq1k": -0.0038,
        "CpA1k": 0.3456, "CpB1k": 74.0368, "CpC1k": -45.7878,
        "smarts": FINAL_C_GEN_SMARTS # This will be the fallback if C_quat_all_alkyl is not matched
    },
    "CH2=CH": {
        "tfp1k": 1.6472, "tb1k": 1.7827, "tc1k": 5.0146, "pc1k": 0.0250, "vc1k": 0.1165, "w1k": 0.408,
        "hf1k": 53.712, "gf1k": 84.926, "hv1k": 6.714, "vliq1k": 0.0373,
        "CpA1k": 49.2506, "CpB1k": 59.3840, "CpC1k": -21.7908,
        "smarts": "[CH2X3H2R0]=[CH1X3H1R0]"
    },
    "CH=CH": {
        "tfp1k": 1.6322, "tb1k": 1.8433, "tc1k": 7.3691, "pc1k": 0.0179, "vc1k": 0.0954, "w1k": 0.252,
        "hf1k": 69.939, "gf1k": 92.900, "hv1k": 7.370, "vliq1k": 0.0269,
        "CpA1k": 35.2248, "CpB1k": 62.1924, "CpC1k": -24.8156,
        "smarts": "[CH1X3H1R0;$(*-[#6;!a;!$(C=O)])]=[CH1X3H1R0;$(*-[#6;!a;!$(C=O)])]"
    },
    "CH2=C": {
        "tfp1k": 1.7899, "tb1k": 1.7117, "tc1k": 6.5081, "pc1k": 0.0223, "vc1k": 0.0918, "w1k": 0.223,
        "hf1k": 64.145, "gf1k": 88.402, "hv1k": 6.797, "vliq1k": 0.0270,
        "CpA1k": 37.6299, "CpB1k": 62.1285, "CpC1k": -26.0637,
        "smarts": "[CH2X3H2R0]=[#6X3H0R0]([#6;!a;!$(C=O)])([#6;!a;!$(C=O)])"
    },
    "CH=C": { # CH=C (3) from table, implying C is trisubstituted, CH is disubstituted (part of C=C and one other bond)
        "tfp1k": 2.0018, "tb1k": 1.7957, "tc1k": 8.9582, "pc1k": 0.0126, "vc1k": 0.0733, "w1k": 0.235,
        "hf1k": 82.528, "gf1k": 93.745, "hv1k": 8.178, "vliq1k": 0.0161,
        "CpA1k": 21.3528, "CpB1k": 66.3947, "CpC1k": -29.3703,
        "smarts": "[#6;!a;!$(C=O)][CH1X3H1R0]=[#6X3H0R0]([#6;!a;!$(C=O)])([#6;!a;!$(C=O)])"
    },
    "C=C": { # C=C (4) from table, tetrasubstituted alkene
        "tfp1k": 5.1175, "tb1k": 1.8881, "tc1k": 11.3764, "pc1k": 0.0020, "vc1k": 0.0762, "w1k": -0.210,
        "hf1k": 104.293, "gf1k": 116.613, "hv1k": 9.342, "vliq1k": 0.0030,
        "CpA1k": 10.2797, "CpB1k": 65.5372, "CpC1k": -30.6057,
        "smarts": "[#6X3D3;!R]([#6;!a;!$(C=[O,S])])([#6;!a;!$(C=[O,S])])=[#6X3D3;!R]([#6;!a;!$(C=[O,S])])([#6;!a;!$(C=[O,S])])"
    },
    "CH": { # CH (3)
        "tfp1k": 0.3557, "tb1k": 0.6033, "tc1k": 4.0330, "pc1k": 0.0013, "vc1k": 0.0315, "w1k": -0.071,
        "hf1k": -3.766, "gf1k": 19.848, "hv1k": 2.771, "vliq1k": 0.0071,
        "CpA1k": 8.9272, "CpB1k": 59.9786, "CpC1k": -29.5143,
        "smarts": FINAL_CH_GEN_SMARTS
    },
    "CH2=C=CH": {
        "tfp1k": 3.3439, "tb1k": 3.1243, "tc1k": 9.9318, "pc1k": 0.0313, "vc1k": 0.1483, "w1k": 0.152,
        "hf1k": 197.322, "gf1k": 221.308, "hv1k": 12.318, "vliq1k": 0.0434,
        "CpA1k": 66.0574, "CpB1k": 69.3936, "CpC1k": -25.1081,
        "smarts": "[CH2X3H2R0]=[#6X2H0R0]=[CH1X3H1R0][#6;!a;!$(C=O)]"
    },
    "ACH": {
        "tfp1k": 1.4669, "tb1k": 0.9297, "tc1k": 3.7337, "pc1k": 0.0075, "vc1k": 0.0422, "w1k": 0.027,
        "hf1k": 11.189, "gf1k": 22.533, "hv1k": 4.098, "vliq1k": 0.0132,
        "CpA1k": 16.3794, "CpB1k": 32.7433, "CpC1k": -13.1692,
        "smarts": "[cH1]"
    },
    "AC": {
        "tfp1k": 0.2098, "tb1k": 1.6254, "tc1k": 14.6409, "pc1k": 0.0021, "vc1k": 0.0398, "w1k": 0.334,
        "hf1k": 27.016, "gf1k": 30.485, "hv1k": 12.552, "vliq1k": 0.0044,
        "CpA1k": 10.4283, "CpB1k": 25.3634, "CpC1k": -12.7283,
        "smarts": "[c;H0;D3]"
    },
    "ACCH3": {
        "tfp1k": 1.8635, "tb1k": 1.9669, "tc1k": 8.2130, "pc1k": 0.0194, "vc1k": 0.1036, "w1k": 0.146,
        "hf1k": -19.243, "gf1k": 22.505, "hv1k": 9.776, "vliq1k": 0.0289,
        "CpA1k": 42.8569, "CpB1k": 65.6464, "CpC1k": -21.0670,
        "smarts": "[CH3X4;$(*-[c])]"
    },
    "ACCH2": {
        "tfp1k": 0.4177, "tb1k": 1.9478, "tc1k": 10.3239, "pc1k": 0.0122, "vc1k": 0.1010, "w1k": -0.088,
        "hf1k": 9.404, "gf1k": 41.228, "hv1k": 10.185, "vliq1k": 0.0192,
        "CpA1k": 32.8206, "CpB1k": 70.4153, "CpC1k": -28.9361,
        "smarts": "[CH2X4H2;!R;$(*-[c]);$(*-[#6;!a;!$(C=O)])]"
    },
    "ACCH": {
        "tfp1k": -1.7567, "tb1k": 1.7444, "tc1k": 10.4664, "pc1k": 0.0028, "vc1k": 0.0712, "w1k": 1.524,
        "hf1k": 27.671, "gf1k": 52.948, "hv1k": 8.834, "vliq1k": 0.0099,
        "CpA1k": 19.9504, "CpB1k": 81.8764, "CpC1k": -40.2864,
        "smarts": "[CH1X4H1;!R;$(*-[c]);$(*-[#6;!a;!$(C=O)]);$(*-[#6;!a;!$(C=O)])]"
    },
    "ACOH (2)": {
        "tfp1k": 13.7349, "tb1k": 4.4014, "tc1k": 25.9145, "pc1k": -0.0074, "vc1k": 0.0316, "w1k": 1.015,
        "hf1k": -164.609, "gf1k": -132.097, "hv1k": 40.246, "vliq1k": 0.0113,
        "CpA1k": 39.7712, "CpB1k": 35.5676, "CpC1k": -15.5875,
        "smarts": "[OH1X2]-[c]"
    },
    "CH3CO (1)": {
        "tfp1k": 4.8776, "tb1k": 3.5668, "tc1k": 13.2896, "pc1k": 0.0251, "vc1k": 0.1340, "w1k": 0.633,
        "hf1k": -182.329, "gf1k": -131.366, "hv1k": 18.999, "vliq1k": 0.0365,
        "CpA1k": 59.3032, "CpB1k": 67.8149, "CpC1k": -20.9948,
        "smarts": "[CH3X4;!R;$(*-[CX3R0;!c;!H1;!$(*-[O;X2H0;!R]-[#6])]=[OX1D1])]"
    },
    "CH2CO (2)": {
        "tfp1k": 5.6622, "tb1k": 3.8967, "tc1k": 14.6273, "pc1k": 0.0178, "vc1k": 0.1119, "w1k": 0.963,
        "hf1k": -164.410, "gf1k": -132.386, "hv1k": 20.041, "vliq1k": 0.0282,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH2X4H2;!R;$(*-[CX3R0;!c;!H1;!$(*-[O;X2H0;!R]-[#6])]=[OX1D1]);$(*-[#6;!a;!$(C=O)])]"
    },
    "CHO (1)*": {
        "tfp1k": 4.2927, "tb1k": 2.8526, "tc1k": 10.1986, "pc1k": 0.0141, "vc1k": 0.0863, "w1k": 1.133,
        "hf1k": -129.200, "gf1k": -107.858, "hv1k": 12.909, "vliq1k": 0.0200,
        "CpA1k": 40.7501, "CpB1k": 19.6990, "CpC1k": -5.4360,
        "smarts": "[CH1X3H1R0;!$(*~[OX2H0])](=[OX1D1])"
    },
    "CH3COO (1)": {
        "tfp1k": 4.0823, "tb1k": 3.6360, "tc1k": 12.5965, "pc1k": 0.0290, "vc1k": 0.1589, "w1k": 0.756,
        "hf1k": -389.737, "gf1k": -318.616, "hv1k": 22.709, "vliq1k": 0.0450,
        "CpA1k": 66.8423, "CpB1k": 102.4553, "CpC1k": -43.3306,
        "smarts": "[CH3X4][CX3;!c](=[O])[O;X2H0;!R]"
    },
    "CH2COO (2)": {
        "tfp1k": 3.5572, "tb1k": 3.3953, "tc1k": 13.8116, "pc1k": 0.0218, "vc1k": 0.1365, "w1k": 0.765,
        "hf1k": -359.258, "gf1k": -291.188, "hv1k": 17.759, "vliq1k": 0.0357,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH2X4;!$(C=O);!R;$(*-[#6;!a])][CX3;!c](=[O])[O;X2H0;!R]"
    },
    "HCOO (1)": {
        "tfp1k": 4.2250, "tb1k": 3.1459, "tc1k": 11.6057, "pc1k": 0.0138, "vc1k": 0.1056, "w1k": 0.526,
        "hf1k": -332.822, "gf1k": -288.902, "hv1k": None, "vliq1k": 0.0267,
        "CpA1k": 51.5048, "CpB1k": 44.4133, "CpC1k": -19.6155,
        "smarts": "[CH1X3H1;!R](=[O])[O;X2H0;!R]"

    },
    "CH3O (1)": {
        "tfp1k": 2.9248, "tb1k": 2.2536, "tc1k": 6.4737, "pc1k": 0.0204, "vc1k": 0.0875, "w1k": 0.442,
        "hf1k": -163.569, "gf1k": -105.767, "hv1k": 10.919, "vliq1k": 0.0327,
        "CpA1k": 50.5604, "CpB1k": 38.9681, "CpC1k": -4.7799,
        "smarts": "[CH3X4;!R;$(*-[OX2H0;!R;D2]-[#6;!R;!c;!$([CX3](=[OX1]))])])]"
    },
    "CH2O (2)": {
        "tfp1k": 2.0695, "tb1k": 1.6249, "tc1k": 6.0723, "pc1k": 0.0151, "vc1k": 0.0729, "w1k": 0.218,
        "hf1k": -151.143, "gf1k": -101.563, "hv1k": 7.478, "vliq1k": 0.0231,
        "CpA1k": 39.5784, "CpB1k": 41.8177, "CpC1k": -11.0837,
        "smarts": "[CH2X4;!R;$(*-[OX2H0;!R;D2;!$(*-[c]);!$(*-[CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))])])]"
    },
    "CH=O (3)": {
        "tfp1k": 4.0352, "tb1k": 1.1557, "tc1k": 5.0663, "pc1k": 0.0099, "vc1k": 0.0587, "w1k": 0.509,
        "hf1k": -129.488, "gf1k": -92.099, "hv1k": 5.708, "vliq1k": 0.0180,
        "CpA1k": 25.6750, "CpB1k": 24.7281, "CpC1k": 4.2419,
        "smarts": "[CH1X4;!R;$(*-[OX2H0;!R;D2;!$(*-[c]);!$(*-[CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))])])]"
    },
    "FCH2O (1)*": {
        "tfp1k": 4.5047, "tb1k": 2.5892, "tc1k": 9.5059, "pc1k": 0.0090, "vc1k": 0.0686, "w1k": 0.800,
        "hf1k": -140.313, "gf1k": -90.883, "hv1k": 11.227, "vliq1k": 0.0206,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]"
    },
    "CH2NH2 (1)": {
        "tfp1k": 6.7684, "tb1k": 3.1656, "tc1k": 12.1726, "pc1k": 0.0126, "vc1k": 0.1313, "w1k": None,
        "hf1k": -15.505, "gf1k": 58.085, "hv1k": 14.599, "vliq1k": 0.0265,
        "CpA1k": 57.6861, "CpB1k": 64.0768, "CpC1k": -21.0480,
        "smarts": "[CH2X4H2]([NH2X3])-[#6;!R;!c]"
    },
    "CHNH2 (2)": {
        "tfp1k": 4.1187, "tb1k": 2.5983, "tc1k": 10.2075, "pc1k": 0.0107, "vc1k": 0.0753, "w1k": 0.953,
        "hf1k": 3.320, "gf1k": 63.051, "hv1k": 11.876, "vliq1k": 0.0195,
        "CpA1k": 44.1122, "CpB1k": 77.2155, "CpC1k": -33.5086,
        "smarts": "[CH1X4H1]([NH2X3])(-[#6;!R;!c])-[#6;!R;!c]"
    },
    "CH3NH (2)": {
        "tfp1k": 4.5341, "tb1k": 3.1376, "tc1k": 9.8544, "pc1k": 0.0126, "vc1k": 0.1215, "w1k": 0.550,
        "hf1k": 5.432, "gf1k": 82.471, "hv1k": 14.452, "vliq1k": 0.0267,
        "CpA1k": 53.7012, "CpB1k": 71.7948, "CpC1k": -22.9685,
        "smarts": "[CH3X4;!R]([NH1X3H1]-[#6;!R;!c])"
    },
    "CH2NH (3)": {
        "tfp1k": 6.0609, "tb1k": 2.6127, "tc1k": 10.4677, "pc1k": 0.0104, "vc1k": 0.0996, "w1k": 0.386,
        "hf1k": 23.101, "gf1k": 95.888, "hv1k": 14.481, "vliq1k": 0.0232,
        "CpA1k": 44.6388, "CpB1k": 68.5041, "CpC1k": -26.7106,
        "smarts": "[CH2X4H2;!R]([NH1X3H1]-[#6;!R;!c])-[#6;!R;!c]"
    },
    "CHNH (4)*": {
        "tfp1k": 3.4100, "tb1k": 1.5780, "tc1k": 7.2121, "pc1k": -0.0005, "vc1k": 0.0916, "w1k": 0.384,
        "hf1k": 26.718, "gf1k": 85.001, "hv1k": None, "vliq1k": 0.0181,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH1X4H1;!R]([NH1X3H1]-[#6;!R;!c])(-[#6;!R;!c])-[#6;!R;!c]"
    },
    "CH3N (2)": {
        "tfp1k": 4.0580, "tb1k": 2.1647, "tc1k": 7.6924, "pc1k": 0.0159, "vc1k": 0.1260, "w1k": 0.075,
        "hf1k": 54.929, "gf1k": 128.602, "hv1k": 6.947, "vliq1k": 0.0191,
        "CpA1k": 41.4064, "CpB1k": 85.0996, "CpC1k": -35.6318,
        "smarts": "[CH3X4;!R;$(*-N(-[#6;!R;!c])-[#6;!R;!c])]"
    },
    "CH2N (3)": {
        "tfp1k": 0.9544, "tb1k": 1.2171, "tc1k": 5.5172, "pc1k": 0.0049, "vc1k": 0.0670, "w1k": 0.793,
        "hf1k": 69.885, "gf1k": 132.756, "hv1k": 6.918, "vliq1k": 0.0168,
        "CpA1k": 30.1561, "CpB1k": 81.6814, "CpC1k": -36.1441,
        "smarts": "[CH2X4H2;!R](N(-[#6;!R;!c])-[#6;!R;!c])-[#6;!R;!c]"
    },
    "ACNH2 (2)": {
        "tfp1k": 10.1031, "tb1k": 5.4736, "tc1k": 28.7570, "pc1k": 0.0011, "vc1k": 0.0636, "w1k": None,
        "hf1k": 20.079, "gf1k": 68.861, "hv1k": 28.453, "vliq1k": 0.0137,
        "CpA1k": 47.1311, "CpB1k": 51.3326, "CpC1k": -25.0276,
        "smarts": "[NH2X3]-[c]"
    },
    "C5H4N (1)": {
        "tfp1k": None, "tb1k": 6.2800, "tc1k": 29.1528, "pc1k": 0.0296, "vc1k": 0.2483, "w1k": None,
        "hf1k": 134.062, "gf1k": 199.958, "hv1k": 31.523, "vliq1k": 0.0608,
        "CpA1k": 84.7602, "CpB1k": 177.2513, "CpC1k": -72.3213,
        "smarts": "[n;r6;H0X2](:c):c"
    },
    "C5H3N (2)": { # Same as C5H4N (1) if it refers to the N atom itself in a substituted pyridine or fused system.
        "tfp1k": 12.6275, "tb1k": 5.9234, "tc1k": 27.9464, "pc1k": 0.0257, "vc1k": 0.1703, "w1k": None,
        "hf1k": 139.758, "gf1k": 199.288, "hv1k": 31.005, "vliq1k": 0.0524,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[n;r6;H0X2](:c):c"
    },
    "CH2CN (1)*": {
        "tfp1k": 4.1859, "tb1k": 5.0525, "tc1k": 20.3781, "pc1k": 0.0361, "vc1k": 0.1583, "w1k": 1.670,
        "hf1k": 88.298, "gf1k": 121.544, "hv1k": 23.340, "vliq1k": 0.0331,
        "CpA1k": 58.2837, "CpB1k": 49.6388, "CpC1k": -15.6291,
        "smarts": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "COOH (1)": {
        "tfp1k": 11.5630, "tb1k": 5.8337, "tc1k": 23.7593, "pc1k": 0.0115, "vc1k": 0.1019, "w1k": 0.570,
        "hf1k": -396.242, "gf1k": -349.439, "hv1k": 43.046, "vliq1k": 0.0223,
        "CpA1k": 46.5577, "CpB1k": 48.2322, "CpC1k": -20.4868,
        "smarts": "[#6X3R0](=[OX1R0])[OH1X2R0]"
    },
    "CH2Cl (1)": {
        "tfp1k": 3.3376, "tb1k": 2.9637, "tc1k": 11.0752, "pc1k": 0.0198, "vc1k": 0.1156, "w1k": None,
        "hf1k": -73.568, "gf1k": -33.373, "hv1k": 13.780, "vliq1k": 0.0337,
        "CpA1k": 48.4648, "CpB1k": 37.2370, "CpC1k": -13.0635,
        "smarts": "[CH2X4H2;!R]([ClX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CHCl (2)": {
        "tfp1k": 2.9933, "tb1k": 2.6948, "tc1k": 10.8632, "pc1k": 0.0114, "vc1k": 0.1035, "w1k": None,
        "hf1k": -63.795, "gf1k": -31.502, "hv1k": 11.985, "vliq1k": 0.0266,
        "CpA1k": 36.5885, "CpB1k": 47.6004, "CpC1k": -22.8148,
        "smarts": "[CH1X4H1;!R]([ClX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CCl (3)": {
        "tfp1k": 9.8409, "tb1k": 2.2073, "tc1k": 11.3959, "pc1k": 0.0031, "vc1k": 0.0792, "w1k": 0.716,
        "hf1k": -57.795, "gf1k": -25.261, "hv1k": 9.818, "vliq1k": 0.0202,
        "CpA1k": 29.1848, "CpB1k": 52.3817, "CpC1k": -30.8526,
        "smarts": "[CX4H0;!R]([ClX1])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CHCl2 (1)*": {
        "tfp1k": 5.1638, "tb1k": 3.9300, "tc1k": 16.3945, "pc1k": 0.0268, "vc1k": 0.1695, "w1k": None,
        "hf1k": -82.921, "gf1k": -35.814, "hv1k": 19.208, "vliq1k": 0.0468,
        "CpA1k": 60.8262, "CpB1k": 41.9908, "CpC1k": -20.4091,
        "smarts": "[CH1X4H1](Cl)(Cl)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CCl3 (1)": {
        "tfp1k": None, "tb1k": 3.5600, "tc1k": None, "pc1k": None, "vc1k": 0.0617, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": 17.574, "vliq1k": 0.0620,
        "CpA1k": 56.1685, "CpB1k": 46.9337, "CpC1k": -31.3325,
        "smarts": "[CX4H0](Cl)(Cl)(Cl)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CCl2 (2)": {
        "tfp1k": 10.2337, "tb1k": 4.5797, "tc1k": 18.5875, "pc1k": 0.0349, "vc1k": 0.2103, "w1k": None,
        "hf1k": -107.188, "gf1k": -53.332, "hv1k": None, "vliq1k": None,
        "CpA1k": 78.6054, "CpB1k": 32.1318, "CpC1k": -19.4033,
        "smarts": "[CX4H0](Cl)(Cl)(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "ACCl (2)": {
        "tfp1k": 2.7336, "tb1k": 2.6293, "tc1k": 14.1565, "pc1k": 0.0131, "vc1k": 0.1016, "w1k": 0.296,
        "hf1k": -16.752, "gf1k": -0.596, "hv1k": 11.883, "vliq1k": 0.0241,
        "CpA1k": 33.6450, "CpB1k": 23.2759, "CpC1k": -12.2406,
        "smarts": "[Cl]-[c]"
    },
    "CH2NO2 (1)*": {
        "tfp1k": 5.5424, "tb1k": 5.7619, "tc1k": 24.7369, "pc1k": 0.0210, "vc1k": 0.1653, "w1k": None,
        "hf1k": -66.138, "gf1k": 17.963, "hv1k": 30.644, "vliq1k": 0.0338,
        "CpA1k": 63.7851, "CpB1k": 83.4744, "CpC1k": -35.1171,
        "smarts": "[CH2X4H2]([$(N(=O)=O),$(N(=[O])-[O-])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CHNO2 (2)*": {
        "tfp1k": 4.9738, "tb1k": 5.0767, "tc1k": 23.2050, "pc1k": 0.0122, "vc1k": 0.1423, "w1k": None,
        "hf1k": -59.142, "gf1k": 18.088, "hv1k": 26.277, "vliq1k": 0.0262,
        "CpA1k": 51.1442, "CpB1k": 94.2934, "CpC1k": -45.2029,
        "smarts": "[CH1X4H1]([$(N(=O)=O),$(N(=[O])-[O-])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "ACNO2 (2)*": {
        "tfp1k": 8.4724, "tb1k": 6.0837, "tc1k": 34.5870, "pc1k": 0.0150, "vc1k": 0.1426, "w1k": None,
        "hf1k": -7.365, "gf1k": 60.161, "hv1k": None, "vliq1k": 0.0250,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[c]-[$(N(=O)=O),$(N(=[O])-[O-])]"
    },
    "CH2SH (1)": {
        "tfp1k": 3.0044, "tb1k": 3.2914, "tc1k": 13.8058, "pc1k": 0.0136, "vc1k": 0.1025, "w1k": None,
        "hf1k": -8.253, "gf1k": 16.731, "hv1k": 14.931, "vliq1k": 0.0345,
        "CpA1k": 58.2445, "CpB1k": 46.9958, "CpC1k": -10.5106,
        "smarts": "[CH2X4H2]([SX2H1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "I (1)*": {
        "tfp1k": 4.6089, "tb1k": 3.6650, "tc1k": 17.3947, "pc1k": 0.0028, "vc1k": 0.1081, "w1k": 0.233,
        "hf1k": 57.546, "gf1k": 46.945, "hv1k": 14.364, "vliq1k": 0.0279,
        "CpA1k": 29.1815, "CpB1k": -9.7846, "CpC1k": 3.4554,
        "smarts": "[#53X1]-[#6]"
    },
    "Br (1)": {
        "tfp1k": 3.7442, "tb1k": 2.6495, "tc1k": 10.5371, "pc1k": -0.0018, "vc1k": 0.0828, "w1k": 0.278,
        "hf1k": 1.834, "gf1k": -1.721, "hv1k": 11.423, "vliq1k": 0.0214,
        "CpA1k": 28.0260, "CpB1k": -7.1651, "CpC1k": 2.4332,
        "smarts": "[#35X1]-[#6]"
    },
    "CH#C (1)": {
        "tfp1k": 3.9106, "tb1k": 2.3678, "tc1k": 7.5433, "pc1k": 0.0148, "vc1k": 0.0933, "w1k": 0.618,
        "hf1k": 220.803, "gf1k": 217.003, "hv1k": 7.751, "vliq1k": None,
        "CpA1k": 45.9768, "CpB1k": 20.6417, "CpC1k": -8.3297,
        "smarts": "[CH1X2H1]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
    },
    "C#C (2)*": {
        "tfp1k": 9.5793, "tb1k": 2.5645, "tc1k": 11.4501, "pc1k": 0.0041, "vc1k": 0.0763, "w1k": None,
        "hf1k": 227.368, "gf1k": 216.328, "hv1k": 11.549, "vliq1k": 0.0145,
        "CpA1k": 26.7371, "CpB1k": 21.7676, "CpC1k": -6.4481,
        "smarts": "[#6;!a;!c;!R;!$([CX3](=[OX1]))]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
    },
    "Cl-C#C (3)*": {
        "tfp1k": 1.5598, "tb1k": 1.7824, "tc1k": 5.4334, "pc1k": 0.0160, "vc1k": 0.0569, "w1k": None,
        "hf1k": -36.097, "gf1k": -28.148, "hv1k": None, "vliq1k": 0.0153,
        "CpA1k": 25.8094, "CpB1k": -5.2241, "CpC1k": 1.4542,
        "smarts": "[ClX1]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
    },
    "ACF (2)": {
        "tfp1k": 2.5015, "tb1k": 0.9442, "tc1k": 2.8977, "pc1k": 0.0130, "vc1k": 0.0567, "w1k": 0.263,
        "hf1k": -161.740, "gf1k": -144.549, "hv1k": 4.877, "vliq1k": 0.0173,
        "CpA1k": 30.1696, "CpB1k": 26.9738, "CpC1k": -13.3722,
        "smarts": "[F]-[c]"
    },
    "HCON(CH3)2 (2)*": {
        "tfp1k": None, "tb1k": 7.2644, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": 0.500,
        "hf1k": None, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
    },
    "CF3 (1)": {
        "tfp1k": 3.2411, "tb1k": 1.2880, "tc1k": 2.4778, "pc1k": 0.0442, "vc1k": 0.1148, "w1k": None,
        "hf1k": -679.195, "gf1k": -626.580, "hv1k": 8.901, "vliq1k": None,
        "CpA1k": 63.2024, "CpB1k": 51.9366, "CpC1k": -28.6308,
        "smarts": "[CX4H0](F)(F)(F)"
    },
    "CF2 (2)": {
        "tfp1k": None, "tb1k": 0.6115, "tc1k": 1.7399, "pc1k": 0.0129, "vc1k": 0.0952, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": 1.860, "vliq1k": None,
        "CpA1k": 44.3567, "CpB1k": 44.5875, "CpC1k": -23.2820,
        "smarts": "[CX4;!$(C(F)(F)F)](F)(F)"
    },
    "CF (3)": {
        "tfp1k": None, "tb1k": 1.1739, "tc1k": 3.5192, "pc1k": 0.0047, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": 8.901, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CX4;!$(C(F)(F)F);!$(C(F)F)](F)"
    },
    "COO (2) Ester Carbonyl?": {
        "tfp1k": 3.4448, "tb1k": 2.6446, "tc1k": 12.1084, "pc1k": 0.0113, "vc1k": 0.0859, "w1k": None,
        "hf1k": -313.545, "gf1k": -281.495, "hv1k": None, "vliq1k": 0.0192,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "SMARTS_PATTERN_TO_BE_DEFINED_LATER"
    },
    "CCl2F (1)": {
        "tfp1k": 7.4756, "tb1k": 2.8881, "tc1k": 9.8408, "pc1k": 0.0354, "vc1k": 0.1821, "w1k": 0.503,
        "hf1k": -258.960, "gf1k": -209.337, "hv1k": 13.322, "vliq1k": 0.0538,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CX4H0](Cl)(Cl)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "HCClF (1)": {
        "tfp1k": None, "tb1k": 2.3086, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CH1X4H1](Cl)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "CClF2 (1)": {
        "tfp1k": 2.7523, "tb1k": 1.9163, "tc1k": 4.8923, "pc1k": 0.0390, "vc1k": 0.1475, "w1k": 0.547,
        "hf1k": -446.835, "gf1k": -392.975, "hv1k": 8.301, "vliq1k": 0.0538,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[CX4H0](Cl)(F)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
    },
    "FSpecial (1)": {
        "tfp1k": 1.9623, "tb1k": 1.0081, "tc1k": 1.5974, "pc1k": 0.0144, "vc1k": 0.0378, "w1k": None,
        "hf1k": -223.398, "gf1k": -212.718, "hv1k": None, "vliq1k": None,
        "CpA1k": 22.2082, "CpB1k": -2.8385, "CpC1k": 1.2679,
        "smarts": f"[F]-[#6;!c;!$({CH2F_PATTERN_REF});!$({FCH2O_PATTERN_REF});!$(C(F)(F));!$(C(F)(F)F)]"
    },
    "CONH2 (1)*": {
        "tfp1k": 31.2786, "tb1k": 10.3428, "tc1k": 65.1053, "pc1k": 0.0043, "vc1k": 0.1443, "w1k": None,
        "hf1k": -203.188, "gf1k": -136.742, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]"
    },
    "CONHCH3 (1)*": {
        "tfp1k": None, "tb1k": None, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": -67.778, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
    },
    "CONHCH2 (1)*": {
        "tfp1k": None, "tb1k": None, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": -182.096, "gf1k": None, "hv1k": 51.787, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]"
    },
    "CON(CH3)2 (1)*": {
        "tfp1k": 11.3770, "tb1k": 7.6904, "tc1k": 36.1403, "pc1k": 0.0401, "vc1k": 0.2503, "w1k": None,
        "hf1k": -189.888, "gf1k": -65.642, "hv1k": None, "vliq1k": 0.0548,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
    },
    "CONCH2CH3 (3)*": {
        "tfp1k": None, "tb1k": None, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": -46.562, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]"
    },
    "CON(CH2)2 (3)*": {
        "tfp1k": None, "tb1k": 6.7822, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "SMARTS_PATTERN_TO_BE_DEFINED_LATER" # Lactam N - still elusive for a general fragment
    },
    "C2H5O2 (1)*": {
        "tfp1k": None, "tb1k": 5.5566, "tc1k": 17.9668, "pc1k": 0.0254, "vc1k": 0.1675, "w1k": 0.428,
        "hf1k": -344.125, "gf1k": -241.373, "hv1k": None, "vliq1k": 0.0410,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "SMARTS_PATTERN_TO_BE_DEFINED_LATER"
    },
    "C2H4O2 (2)": {
        "tfp1k": None, "tb1k": 5.4248, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "SMARTS_PATTERN_TO_BE_DEFINED_LATER"
    },
    "CH3S (1)": {
        "tfp1k": 5.0506, "tb1k": 3.6796, "tc1k": 14.3969, "pc1k": 0.0160, "vc1k": 0.1302, "w1k": None,
        "hf1k": -2.084, "gf1k": 30.222, "hv1k": 16.921, "vliq1k": 0.0348,
        "CpA1k": 57.7670, "CpB1k": 44.1238, "CpC1k": -9.5565,
        "smarts": "[CH3X4;!R;$(*-S-[#6;!R;!c;!$(C=O)])]"
    },
    "CH2S (2)": {
        "tfp1k": 3.1468, "tb1k": 3.6763, "tc1k": 17.7916, "pc1k": 0.0111, "vc1k": 0.1165, "w1k": 0.438,
        "hf1k": 18.022, "gf1k": 38.346, "hv1k": 17.117, "vliq1k": 0.0273,
        "CpA1k": 45.0314, "CpB1k": 55.1432, "CpC1k": -18.7776,
        "smarts": "[CH2X4H2;!R]([SX2H0]-[#6;!R;!c;!$(C=O)])-[#6;!R;!c;!$(C=O)]"
    },
    "CHS (3)*": {
        "tfp1k": None, "tb1k": 2.6812, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": 0.739,
        "hf1k": None, "gf1k": None, "hv1k": 13.265, "vliq1k": None,
        "CpA1k": 40.5275, "CpB1k": 55.0141, "CpC1k": -31.7190,
        "smarts": "[CH1X4H1;!R]([SX2H0]-[#6;!R;!c;!$(C=O)])(-[#6;!R;!c;!$(C=O)])-[#6;!R;!c;!$(C=O)]"
    },
    "C4H3S (1)": {
        "tfp1k": None, "tb1k": 5.7093, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": 27.966, "vliq1k": None,
        "CpA1k": 80.3010, "CpB1k": 132.7786, "CpC1k": -58.3241,
        "smarts": "[s;r5;H0X2](:c):c"
    },
    "C4H2S (2)*": {
        "tfp1k": None, "tb1k": 5.8260, "tc1k": None, "pc1k": None, "vc1k": None, "w1k": None,
        "hf1k": None, "gf1k": None, "hv1k": None, "vliq1k": None,
        "CpA1k": None, "CpB1k": None, "CpC1k": None,
        "smarts": "[s;r5;H0X2](:c):c"
    }
}

CG_SECOND_ORDER_GROUPS_DATA = {
    # From Page C.9
    "(CH3)2CH": {
        "tfp2j": 0.0381, "tb2j": -0.1157, "tc2j": -0.5334, "pc2j": 0.000488, "vc2j": 0.00400, "w2j": 0.01740,
        "hf2j": -0.860, "gf2j": 0.297, "hv2j": 0.292, "Vliq2j": 0.00133,
        "CpA2j": 0.5830, "CpB2j": -1.2002, "CpC2j": -0.0584,
        "smarts": "[CH1X4D3H1]([CH3X4D1H3])([CH3X4D1H3])[#6]" # Isopropyl attached to another carbon (not H)
    },
    "(CH3)3C": {
        "tfp2j": -0.2355, "tb2j": -0.0489, "tc2j": -0.5143, "pc2j": 0.001410, "vc2j": 0.00572, "w2j": 0.01922,
        "hf2j": -1.338, "gf2j": -0.399, "hv2j": -0.720, "Vliq2j": 0.00179,
        "CpA2j": 0.3226, "CpB2j": 2.1309, "CpC2j": -1.5728,
        "smarts": "[#6X4D4H0]([CH3X4D1H3])([CH3X4D1H3])([CH3X4D1H3])[#6]" # Tert-butyl attached to another carbon
    },
    "CH(CH3)CH(CH3)": { # Two adjacent sec-methyl groups, e.g., internal part of 2,3-dimethylbutane.
        "tfp2j": 0.4401, "tb2j": 0.1798, "tc2j": 1.0699, "pc2j": -0.001850, "vc2j": -0.00398, "w2j": -0.00475,
        "hf2j": 6.771, "gf2j": 6.342, "hv2j": 0.868, "Vliq2j": -0.00203,
        "CpA2j": 0.9668, "CpB2j": -2.0762, "CpC2j": 0.3148,
        "smarts": "[CH1X4D3H1]([CH3X4D1H3])[CH1X4D3H1]([CH3X4D1H3])" # -(CH(CH3))-(CH(CH3))-
    },
    "CH(CH3)C(CH3)2": { # sec-methyl next to a gem-dimethyl. e.g. in 2,2,3-trimethylpentane (C(CH3)2 from C2, CH(CH3) from C3)
        "tfp2j": -0.4923, "tb2j": 0.3189, "tc2j": 1.9086, "pc2j": -0.005200, "vc2j": -0.01081, "w2j": -0.01828,
        "hf2j": 7.205, "gf2j": 7.466, "hv2j": 1.027, "Vliq2j": -0.00243,
        "CpA2j": -0.3082, "CpB2j": 1.8909, "CpC2j": -1.6454,
        "smarts": "[CH1X4D3H1]([CH3X4D1H3])[#6X4D4H0]([CH3X4D1H3])([CH3X4D1H3])" # -(CH(CH3))-(C(CH3)2)-
    },
    "C(CH3)2C(CH3)2": { # Two adjacent gem-dimethyl groups, e.g. in 2,2,3,3-tetramethylbutane
        "tfp2j": 6.0650, "tb2j": 0.7273, "tc2j": 5.8254, "pc2j": -0.013230, "vc2j": -0.02300, "w2j": -0.08632,
        "hf2j": 14.271, "gf2j": 16.224, "hv2j": 2.426, "Vliq2j": -0.00744,
        "CpA2j": -0.1201, "CpB2j": 4.2846, "CpC2j": -2.0262,
        "smarts": "[#6X4D4H0]([CH3X4D1H3])([CH3X4D1H3])[#6X4D4H0]([CH3X4D1H3])([CH3X4D1H3])" # -(C(CH3)2)-(C(CH3)2)-
    },
    "3 membered ring": {
        "tfp2j": 1.3772, "tb2j": 0.4745, "tc2j": -2.3305, "pc2j": 0.003714, "vc2j": -0.00014, "w2j": 0.17563,
        "hf2j": 104.800, "gf2j": 94.564, "hv2j": None, "Vliq2j": None,
        "CpA2j": 8.5546, "CpB2j": -22.9771, "CpC2j": 10.7278,
        "smarts": "[*;R1;r3]" # Any atom in a 3-membered ring. Counted per atom in the ring. Check C-4 for counting convention.
                                 # Table C-4: cyclopropane (1). This means the group is counted ONCE for the entire ring structure.
                                 # SMARTS should ideally match the ring structure itself, not per atom.
                                 # If count is 1 for cyclopropane: [*]1~[*]~[*]~1
    },
    "4 membered ring": {
        "tfp2j": None, "tb2j": 0.3563, "tc2j": -1.2978, "pc2j": 0.001171, "vc2j": -0.00851, "w2j": 0.22216,
        "hf2j": 99.455, "gf2j": 92.573, "hv2j": None, "Vliq2j": None,
        "CpA2j": 3.1721, "CpB2j": -10.0834, "CpC2j": 4.9674,
        "smarts": "[*]1~[*]~[*]~[*]~1" # One count per 4-membered ring structure
    },
    "5 membered ring": {
        "tfp2j": 0.6824, "tb2j": 0.1919, "tc2j": -0.6785, "pc2j": 0.000424, "vc2j": -0.00866, "w2j": 0.16284,
        "hf2j": 13.782, "gf2j": 5.733, "hv2j": -0.568, "Vliq2j": 0.00213,
        "CpA2j": -5.9060, "CpB2j": -1.8710, "CpC2j": 4.2945,
        "smarts": "[*]1~[*]~[*]~[*]~[*]~1" # One count per 5-membered ring structure
    },
    "6 membered ring": {
        "tfp2j": 1.5656, "tb2j": 0.1957, "tc2j": 0.8479, "pc2j": 0.002257, "vc2j": 0.01636, "w2j": -0.03065,
        "hf2j": -9.660, "gf2j": -8.180, "hv2j": -0.905, "Vliq2j": 0.00063,
        "CpA2j": -3.9682, "CpB2j": 17.7889, "CpC2j": -3.3639,
        "smarts": "[*]1~[*]~[*]~[*]~[*]~[*]~1" # One count per 6-membered ring structure
    },
    "7 membered ring": {
        "tfp2j": 6.9709, "tb2j": 0.3489, "tc2j": 3.6714, "pc2j": -0.004800, "vc2j": -0.02700, "w2j": -0.02094,
        "hf2j": 15.465, "gf2j": 20.597, "hv2j": -0.847, "Vliq2j": -0.00519,
        "CpA2j": -3.2746, "CpB2j": 32.1670, "CpC2j": -17.8246,
        "smarts": "[*]1~[*]~[*]~[*]~[*]~[*]~[*]~1" # One count per 7-membered ring structure
    },
    "CHn=CHm-CHp=CHk m,p E (0,1), k,n E (0,2)": { # Conjugated diene. Ex: 1,3-butadiene (1 group).
        "tfp2j": 1.9913, "tb2j": 0.1589, "tc2j": 0.4402, "pc2j": 0.004186, "vc2j": -0.00781, "w2j": 0.01648,
        "hf2j": -8.392, "gf2j": -5.505, "hv2j": 2.057, "Vliq2j": -0.00188,
        "CpA2j": 2.6142, "CpB2j": 4.4511, "CpC2j": -5.9808,
        "smarts": "[#6]=[#6]-[#6]=[#6]" # General conjugated diene structure
    },
    "CH3-CHm=CHn m E (0,1), n E (0,2)": { # Methyl group adjacent to a double bond (allylic methyl). Ex: 2-Methyl-2-Butene (3 groups). C(CH3)=C(CH3)CH3 has 3 such.
        "tfp2j": 0.2476, "tb2j": 0.0668, "tc2j": 0.0167, "pc2j": -0.000180, "vc2j": -0.00098, "w2j": 0.00619,
        "hf2j": 0.474, "gf2j": 0.950, "hv2j": -0.073, "Vliq2j": 0.00009,
        "CpA2j": -1.3913, "CpB2j": -1.5496, "CpC2j": 2.5899,
        "smarts": "[CH3X4D1][#6X3D2]=[#6X3]" # CH3-C=C
    },
    "CH2-CHm=CHn m E (0,1), n E (0,2)": { # Allylic methylene. Ex: 1,4-Pentadiene (2 groups). CH2=CH-CH2-CH=CH2 has two.
        "tfp2j": -0.5870, "tb2j": -0.1406, "tc2j": -0.5231, "pc2j": 0.003538, "vc2j": 0.00281, "w2j": -0.01150,
        "hf2j": 1.472, "gf2j": 0.699, "hv2j": -0.369, "Vliq2j": 0.00012,
        "CpA2j": 0.2630, "CpB2j": -2.3428, "CpC2j": 0.8975,
        "smarts": "[CH2X4D2]([#6X3D2]=[#6X3])[#6;!$(C=C)]" # -CH2-C=C (CH2 is not part of another double bond)
    },
    "CH-CHm=CHn or C-CHm=CHn* m E (0,1), n E (0,2)": { # Allylic methine or quaternary carbon. Ex: 4-Methyl-2-Pentene (2 groups). CH3-CH=CH-CH(CH3)2. One is CH, one is C (of CH(CH3)2).
        "tfp2j": -0.2361, "tb2j": -0.0900, "tc2j": -0.3850, "pc2j": 0.005675, "vc2j": 0.00826, "w2j": 0.02778,
        "hf2j": 4.504, "gf2j": 1.013, "hv2j": 0.345, "Vliq2j": 0.00142,
        "CpA2j": 6.5145, "CpB2j": -17.5541, "CpC2j": 10.6977,
        "smarts": "[CH1X4D3,CX4D4]([#6X3]=[#6X3])([#6;!$(C=C)])([#6;!$(C=C)])" # (CH or C)-C=C
    },
    "Alicyclic side-chain CcyclicCm m>1": { # C_cyclic attached to C_aliphatic_chain of length > 1.
        "tfp2j": -2.8298, "tb2j": 0.0511, "tc2j": 2.1160, "pc2j": -0.002550, "vc2j": -0.01755, "w2j": -0.11024,
        "hf2j": 1.252, "gf2j": 1.041, "hv2j": -0.114, "Vliq2j": -0.00107,
        "CpA2j": 4.1707, "CpB2j": -3.1964, "CpC2j": -1.1997,
        "smarts": "[#6R1][CH2X4,CH1X4,CX4][CH2X4,CH1X4,CX4;!R]" # Ring carbon attached to at least two non-ring carbons in a chain.
    },
    "CH3CH3": { # Ethane only
        "tfp2j": 1.4880, "tb2j": 0.6884, "tc2j": 2.0427, "pc2j": 0.002175, "vc2j": 0.00227, "w2j": -0.11240,
        "hf2j": -2.792, "gf2j": -1.062, "hv2j": None, "Vliq2j": None,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[CH3X4D1][CH3X4D1]" # Ethane
    },
    "CH3CHO or CCHO*": { # Acetaldehyde fragment or substituted aldehyde. Ex: 2-Methylbutylaldehyde (1).
        "tfp2j": 2.0547, "tb2j": -0.1074, "tc2j": -1.5826, "pc2j": 0.003659, "vc2j": -0.00664, "w2j": None,
        "hf2j": -2.092, "gf2j": -1.359, "hv2j": 0.207, "Vliq2j": -0.00009,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[CH3X4D1,CH2X4D2,CH1X4D3,CX4D4][CH1X3D1]=[OX1D1]" # R-CHO where R is CH3, CH2, CH, or C.
    },
    "CH3COCH2": { # Methyl ketone with an alpha-CH2. CH3-CO-CH2-R. Ex: 2-Pentanone (1).
        "tfp2j": -0.2951, "tb2j": 0.0224, "tc2j": 0.2996, "pc2j": 0.000474, "vc2j": -0.00051, "w2j": -0.00789,
        "hf2j": 0.975, "gf2j": 0.075, "hv2j": -0.668, "Vliq2j": -0.00030,
        "CpA2j": 3.7978, "CpB2j": -7.3251, "CpC2j": 2.5312,
        "smarts": "[CH3X4D1][#6X3D2](=[OX1D1])[CH2X4D2][#6;!c;!$(C=[O,S])]" # CH3-CO-CH2-R
    },
    "CH3COCH or CH3COC": { # Methyl ketone with alpha-CH or alpha-C. Ex: 3-Methyl-2-Pentanone (1).
        "tfp2j": -0.2986, "tb2j": 0.0920, "tc2j": 0.5018, "pc2j": -0.002300, "vc2j": -0.00122, "w2j": -0.16571,
        "hf2j": 4.753, "gf2j": None, "hv2j": 0.071, "Vliq2j": -0.00108,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[CH3X4D1][#6X3D2](=[OX1D1])[CH1X4D3,CX4D4]" # CH3-CO-(CH or C)
    },
    "Ccyclic=O": { # Carbonyl in a non-aromatic ring. Ex: Cyclohexanone (1).
        "tfp2j": 0.7143, "tb2j": 0.5580, "tc2j": 2.9571, "pc2j": 0.003818, "vc2j": -0.01066, "w2j": None,
        "hf2j": 14.145, "gf2j": 23.539, "hv2j": 0.744, "Vliq2j": -0.00111,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6R1X3D2]=[OX1D1]" # C=O where C is in a ring
    },
    "ACCHO*": { # Aromatic aldehyde. Ex: Benzaldehyde (1).
        "tfp2j": -0.6697, "tb2j": 0.0735, "tc2j": 1.1696, "pc2j": -0.002480, "vc2j": 0.00664, "w2j": None,
        "hf2j": -3.173, "gf2j": -2.602, "hv2j": -3.410, "Vliq2j": -0.00036,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6a][CH1X3D1]=[OX1D1]" # Ar-CHO
    },
    "CH3COOH or CCOOH*": { # Aliphatic carboxylic acid (not formic). Ex: 2-MethylButanoic acid (1).
        "tfp2j": -3.1034, "tb2j": -0.1552, "tc2j": -1.7493, "pc2j": 0.004920, "vc2j": 0.00559, "w2j": 0.08774,
        "hf2j": 1.279, "gf2j": 2.149, "hv2j": None, "Vliq2j": -0.00059,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[CH3X4,CH2X4,CH1X4,CX4;!c][#6X3D1](=[OX1D1])[OH1X2D1]" # R-COOH, R is aliphatic CH3, CH2, CH, C
    },
    "ACCOOH*": { # Aromatic carboxylic acid. Ex: Benzoic acid (1).
        "tfp2j": 28.4324, "tb2j": 0.7801, "tc2j": 6.1279, "pc2j": 0.000344, "vc2j": -0.00415, "w2j": None,
        "hf2j": 12.245, "gf2j": 10.715, "hv2j": 8.502, "Vliq2j": 0.00777,
        "CpA2j": -15.7667, "CpB2j": -0.1174, "CpC2j": 6.1191,
        "smarts": "[#6a][#6X3D1](=[OX1D1])[OH1X2D1]" # Ar-COOH
    },
    "CH3COOCH or CH3COOC": { # Acetate ester CH3-CO-O-R. Ex: 2-Methylethyl Ethanoate (1) -> Isopropyl acetate.
        "tfp2j": 0.4838, "tb2j": -0.2383, "tc2j": -1.3406, "pc2j": 0.000659, "vc2j": -0.00293, "w2j": -0.26623,
        "hf2j": -7.807, "gf2j": -6.208, "hv2j": -3.345, "Vliq2j": 0.00083,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[CH3X4D1][#6X3D2](=[OX1D1])[OX2D2H0][#6;!c;!$(C=[O,S])]" # CH3-CO-O-R (R is CH or C)
    },
    "COCH2COO or COCHCOO or COCCOO*": { # Beta-keto ester or malonate type. Ex: Ethylacetoethanoate (1).
        "tfp2j": 0.0127, "tb2j": 0.4456, "tc2j": 2.5413, "pc2j": 0.001067, "vc2j": -0.00591, "w2j": None,
        "hf2j": 37.462, "gf2j": 29.181, "hv2j": None, "Vliq2j": -0.00036,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X3D2](=[OX1D1])[CH2X4D2,CH1X4D3,CX4D4][#6X3D2](=[OX1D1])[OX2H0]" # R-CO-(C)-CO-OR'
    },
    "CO-O-CO": { # Acid anhydride. Ex: Acetic anhydride (1).
        "tfp2j": -2.3598, "tb2j": -0.1977, "tc2j": -2.7617, "pc2j": -0.004880, "vc2j": -0.00144, "w2j": 0.91939,
        "hf2j": -16.097, "gf2j": -11.809, "hv2j": 1.517, "Vliq2j": 0.00198,
        "CpA2j": -6.4072, "CpB2j": 15.2583, "CpC2j": -8.3149,
        "smarts": "[#6X3D1](=[OX1D1])[OX2D2H0][#6X3D1](=[OX1D1])" # R-CO-O-CO-R'
    },
    # From Page C.10
    "ACCOO*": { # Aromatic ester Ar-CO-O-R. Ex: Ethyl benzoate (1).
        "tfp2j": -2.0198, "tb2j": 0.0835, "tc2j": -3.4235, "pc2j": -0.000540, "vc2j": 0.02605, "w2j": None,
        "hf2j": -9.874, "gf2j": -7.415, "hv2j": None, "Vliq2j": 0.00001,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6a][#6X3D2](=[OX1D1])[OX2D2H0][#6;!c;!$(C=[O,S])]" # Ar-CO-O-R
    },
    "CHOH": { # Secondary alcohol. Ex: 2-Butanol (1).
        "tfp2j": -0.5480, "tb2j": -0.5385, "tc2j": -2.8035, "pc2j": -0.004390, "vc2j": -0.00777, "w2j": 0.03654,
        "hf2j": -3.887, "gf2j": -6.770, "hv2j": -1.398, "Vliq2j": -0.00092,
        "CpA2j": 2.4484, "CpB2j": -0.0765, "CpC2j": 0.1460,
        "smarts": "[CH1X4D2H1]([OH1X2D1])([#6X4;!c;!$(C=[O,N,S])])([#6X4;!c;!$(C=[O,N,S])])" # R1-CH(OH)-R2
    },
    "COH": { # Tertiary alcohol. Ex: 2-Methyl-2-Butanol (1).
        "tfp2j": 0.3189, "tb2j": -0.6331, "tc2j": -3.5442, "pc2j": 0.000178, "vc2j": 0.01511, "w2j": 0.21106,
        "hf2j": -24.125, "gf2j": -20.770, "hv2j": 0.320, "Vliq2j": 0.00175,
        "CpA2j": -1.5252, "CpB2j": -7.6380, "CpC2j": 8.1795,
        "smarts": "[#6X4D3H0]([OH1X2D1])([#6X4;!c;!$(C=[O,N,S])])([#6X4;!c;!$(C=[O,N,S])])([#6X4;!c;!$(C=[O,N,S])])" # R1R2R3C-OH
    },
    "CHm(OH)CHn(OH)*": { # Diol on adjacent carbons (m,n = 0,1,2 for C, CH, CH2). Ex: 1,2,3-Propantriol (1 count for CH(OH)CH(OH)).
        "tfp2j": 0.9124, "tb2j": 1.4108, "tc2j": 5.4941, "pc2j": 0.005052, "vc2j": 0.00397, "w2j": None,
        "hf2j": 0.366, "gf2j": 3.805, "hv2j": -3.661, "Vliq2j": 0.00235,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X4]([OH1X2D1])[#6X4]([OH1X2D1])" # Vicinal diol -C(OH)-C(OH)-
    },
    "CHm cyclic-OH m E (0,1)": { # OH on a cyclic carbon (CH or C). Ex: Cyclopentanol (1).
        "tfp2j": 9.5209, "tb2j": -0.0690, "tc2j": 0.3233, "pc2j": 0.006917, "vc2j": -0.02297, "w2j": None,
        "hf2j": -16.333, "gf2j": -5.487, "hv2j": 4.626, "Vliq2j": -0.00250,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6R1X4;CH1,CH0]([OH1X2D1])" # OH on cyclic CH or C
    },
    "CHn(OH)CHm(NHp)*": { # Amino alcohol on adjacent carbons. Ex: 1-Amino-2-Butanol (1).
        "tfp2j": 2.7826, "tb2j": 1.0682, "tc2j": 5.4864, "pc2j": 0.001408, "vc2j": 0.00433, "w2j": None,
        "hf2j": -2.992, "gf2j": -1.600, "hv2j": None, "Vliq2j": 0.00046,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X4]([OH1X2D1])[#6X4]([NX3H2,NH1,N])" # -C(OH)-C(N)-
    },
    "CHm(NH2)CHn(NH2)*": { # Diamine on adjacent carbons. Ex: 1,2-Diaminopropane (1).
        "tfp2j": 2.5114, "tb2j": 0.4247, "tc2j": 2.0699, "pc2j": 0.002148, "vc2j": 0.00580, "w2j": None,
        "hf2j": 2.855, "gf2j": 1.858, "hv2j": None, "Vliq2j": None,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X4]([NH2X3D1])[#6X4]([NH2X3D1])" # -C(NH2)-C(NH2)-
    },
    "CHm cyclic-NHp-CHn cyclic*": { # Cyclic secondary amine. Ex: Pyrrolidine (1).
        "tfp2j": 1.0729, "tb2j": 0.2499, "tc2j": 2.1345, "pc2j": -0.005950, "vc2j": -0.01380, "w2j": -0.13106,
        "hf2j": 0.351, "gf2j": 8.846, "hv2j": 2.311, "Vliq2j": -0.00179,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[NH1X3D2R1]([#6R1])([#6R1])" # N(H) in a ring, bonded to two ring carbons
    },
    "CHn-O-CHm-CHp*": { # Vinyl ether type structure R-O-CH=CH-R'. Ex: Ethylvinylether (1). CH2=CH-O-CH2CH3.
        "tfp2j": 0.2476, "tb2j": 0.1134, "tc2j": 1.0159, "pc2j": -0.000880, "vc2j": 0.00297, "w2j": None,
        "hf2j": -8.644, "gf2j": -13.167, "hv2j": None, "Vliq2j": -0.00206,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X4,CH3,CH2,CH,C][OX2H0][#6X3]=[#6X3]" # R-O-C=C
    },
    "AC-O-CHm* m E (0,3)": { # Aromatic ether Ar-O-R. Ex: Ethylphenylether (1).
        "tfp2j": 0.1175, "tb2j": -0.2596, "tc2j": -5.3307, "pc2j": -0.002250, "vc2j": -0.00045, "w2j": None,
        "hf2j": 1.532, "gf2j": -0.654, "hv2j": None, "Vliq2j": 0.01203,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6a][OX2H0][#6X4;!c;!$(C=O)]" # Ar-O-Alkyl
    },
    "CHm cyclic-S-CHn cyclic": { # Cyclic thioether. Ex: Tetrahydrothiophene (1).
        "tfp2j": -0.2914, "tb2j": 0.4408, "tc2j": 4.4847, "pc2j": None, "vc2j": None, "w2j": -0.01509,
        "hf2j": -0.329, "gf2j": -2.091, "hv2j": 0.972, "Vliq2j": -0.00023,
        "CpA2j": -2.7407, "CpB2j": 11.1033, "CpC2j": -11.0878,
        "smarts": "[SX2D2R1]([#6R1])([#6R1])" # S in a ring, bonded to two ring carbons
    },
    "CHn=CHm-F m E (0,1), n E (0,2)": { # Fluoroalkene. Ex: 1-Fluoro-1-propene (1). CH3-CH=CHF.
        "tfp2j": -0.0514, "tb2j": -0.1168, "tc2j": -0.4996, "pc2j": 0.000319, "vc2j": -0.00596, "w2j": None,
        "hf2j": None, "gf2j": None, "hv2j": None, "Vliq2j": None,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X3]=[#6X3][FX1D1]" # C=C-F
    },
    "CHn=CHm-Br* m E (0,1), n E (0,2)": { # Bromoalkene. Ex: 1-Bromo-1-propene (1). CH3-CH=CHBr.
        "tfp2j": -1.6425, "tb2j": -0.3201, "tc2j": -1.9334, "pc2j": 0.0, "vc2j": 0.00510, "w2j": None,
        "hf2j": 11.989, "gf2j": 12.373, "hv2j": None, "Vliq2j": 0.0,
        "CpA2j": -1.6978, "CpB2j": 1.0477, "CpC2j": 0.2002,
        "smarts": "[#6X3]=[#6X3][BrX1D1]" # C=C-Br
    },
    "CHn=CHm-I* m E (0,1), n E (0,2)": { # Iodoalkene. Ex: 1-Iodo-1-propene (1). CH3-CH=CHI.
        "tfp2j": None, "tb2j": -0.4453, "tc2j": None, "pc2j": None, "vc2j": None, "w2j": None,
        "hf2j": None, "gf2j": None, "hv2j": None, "Vliq2j": None,
        "CpA2j": None, "CpB2j": None, "CpC2j": None,
        "smarts": "[#6X3]=[#6X3][IX1D1]" # C=C-I
    },
    "ACBr*": { # Bromobenzene type. Ex: Bromobenzene (1).
        "tfp2j": 2.5832, "tb2j": -0.6776, "tc2j": -2.2974, "pc2j": 0.009027, "vc2j": -0.00832, "w2j": -0.03078,
        "hf2j": 12.285, "gf2j": 14.161, "hv2j": -7.488, "Vliq2j": 0.00178,
        "CpA2j": -2.2923, "CpB2j": 3.1142, "CpC2j": -1.4995,
        "smarts": "[#6a][BrX1D1]" # Ar-Br
    },
    "ACI*": { # Iodobenzene type. Ex: Iodobenzene (1).
        "tfp2j": -1.5511, "tb2j": -0.3678, "tc2j": 2.8907, "pc2j": 0.008247, "vc2j": -0.00341, "w2j": 0.00001,
        "hf2j": 11.207, "gf2j": 12.530, "hv2j": -4.864, "Vliq2j": 0.00171,
        "CpA2j": -0.3162, "CpB2j": 2.3711, "CpC2j": -1.4825,
        "smarts": "[#6a][IX1D1]" # Ar-I
    },
    "CHm(NH2)-COOH* m E (0,2)": { # Alpha-amino acid. Ex: 2-Aminohexanoic acid (1).
        "tfp2j": None, "tb2j": None, "tc2j": None, "pc2j": None, "vc2j": None, "w2j": None,
        "hf2j": 11.740, "gf2j": None, "hv2j": None, "Vliq2j": None,
        "CpA2j": None, "CpB2j": None, "CpC2j": -0.0584,
        "smarts": "[#6X4]([NH2X3D1])[#6X3D1](=[OX1D1])[OH1X2D1]" # C(NH2)-COOH
    }
}
# Merge the AI transcribed additions into the main dictionaries
#CG_FIRST_ORDER_GROUPS_DATA.update(CG_FIRST_ORDER_GROUPS_DATA_ADDITIONS)
#CG_SECOND_ORDER_GROUPS_DATA.update(CG_SECOND_ORDER_GROUPS_DATA_ADDITIONS)

# Add this new route function to solvers/critical_props_bp.py
# You can place it after the Joback API functions.

@critical_props_bp.route('/critical-properties/constantinou-gani')
def calculator_cg_page_render():
    return render_template('calculator_constantinou_gani.html')

@critical_props_bp.route('/api/calculate_cg_critical', methods=['GET', 'POST'])
def calculate_cg_critical_api():
    # --- DEFINE CG CITATION AND APPLICABILITY DATA ---
    cg_paper_citation_data = {
        "authors": "Constantinou, L., & Gani, R.", "year": 1994, 
        "title": "New group contribution method for estimating properties of pure compounds", 
        "journal": "AIChE Journal", "volume": "40", "issue": "10",
        "pages": "1697-1710", "doi": "10.1002/aic.690401011"
    }
    today = datetime.date.today()
    site_url = request.host_url.rstrip('/')
    cg_page_url_for_citation = site_url + url_for('critical_props.calculator_cg_page_render', _external=True)

    cg_tool_citation_data = {
        'author': 'ChemE Calc', 'year': today.year,
        'title': 'Constantinou-Gani Critical Properties Calculator [Web Application]',
        'retrieved_date': today.strftime('%Y-%m-%d'),
        'url': cg_page_url_for_citation
    }
    cg_applicability_notes = (
        "Constantinou-Gani method. Generally reliable for critical properties. "
        "May have errors for very small (<3 carbons) or very large molecules, "
        "especially highly fluorinated compounds and larger ring structures. "
        "Second-order contributions generally offer marginal improvements but can be important for "
        "specific structures. Consult original papers (Constantinou & Gani, 1994; 1995)."
    )

    if request.method == 'GET':
        print("--- CG CRITICAL PROPS API GET (Initial Data) ---")
        return jsonify({
            'paper_citation_data': cg_paper_citation_data,
            'tool_citation_data': cg_tool_citation_data, # Ensure this is being sent
            'applicability_notes': cg_applicability_notes
        })

    # --- POST request handling for CG ---
    print("--- CG CRITICAL PROPS API POST START ---")
    if not RDKIT_AVAILABLE:
        return jsonify({"error": "RDKit library not found...", 'paper_citation_data': cg_paper_citation_data, 'tool_citation_data': cg_tool_citation_data, 'applicability_notes': cg_applicability_notes}), 500

    data = request.get_json()
    if not data:
        return jsonify({"error": "No input data provided.", 'paper_citation_data': cg_paper_citation_data, 'tool_citation_data': cg_tool_citation_data, 'applicability_notes': cg_applicability_notes}), 400

    smiles_str = data.get('smiles')
    calculation_order_str = data.get('calculation_order', '0')

    error_payload_cg = { # For returning consistent error structure with CG data
        'paper_citation_data': cg_paper_citation_data,
        'tool_citation_data': cg_tool_citation_data,
        'applicability_notes': cg_applicability_notes
    }

    if not smiles_str:
        error_payload_cg["error"] = "SMILES string is required."
        return jsonify(error_payload_cg), 400
    try:
        W = int(calculation_order_str)
        if W not in [0, 1]: raise ValueError("Calculation order must be 0 or 1.")
    except ValueError as e:
        error_payload_cg["error"] = f"Invalid calculation order: {e}"
        return jsonify(error_payload_cg), 400

    mol = Chem.MolFromSmiles(smiles_str.strip())
    if mol is None:
        error_payload_cg["error"] = f"Invalid SMILES string: '{smiles_str}'."
        return jsonify(error_payload_cg), 400
    
    mol = Chem.AddHs(mol)
    
    sum_tc1k, sum_pc1k, sum_vc1k = 0.0, 0.0, 0.0
    sum_tc2j, sum_pc2j, sum_vc2j = 0.0, 0.0, 0.0
    identified_first_order_groups = {}
    identified_second_order_groups = {}
    
    # --- 1. First-Order Group Contributions ---
    for group_name, group_data in CG_FIRST_ORDER_GROUPS_DATA.items():
        smarts = group_data.get("smarts", "")
        if smarts == "SMARTS_PATTERN_TO_BE_DEFINED" or not smarts: continue 
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                Nk = len(matches)
                if Nk > 0:
                    identified_first_order_groups[group_name] = Nk
                    if group_data.get("tc1k") is not None: sum_tc1k += Nk * group_data["tc1k"]
                    if group_data.get("pc1k") is not None: sum_pc1k += Nk * group_data["pc1k"]
                    if group_data.get("vc1k") is not None: sum_vc1k += Nk * group_data["vc1k"]
        except Exception as e: print(f"Error (1st order CG SMARTS) '{group_name}': {e}")

    if W == 1:
        for group_name, group_data in CG_SECOND_ORDER_GROUPS_DATA.items():
            smarts = group_data.get("smarts", "")
            if smarts == "SMARTS_PATTERN_TO_BE_DEFINED" or not smarts: continue
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    Mj = len(matches)
                    if Mj > 0:
                        identified_second_order_groups[group_name] = Mj
                        if group_data.get("tc2j") is not None: sum_tc2j += Mj * group_data["tc2j"]
                        if group_data.get("pc2j") is not None: sum_pc2j += Mj * group_data["pc2j"]
                        if group_data.get("vc2j") is not None: sum_vc2j += Mj * group_data["vc2j"]
            except Exception as e: print(f"Error (2nd order CG SMARTS) '{group_name}': {e}")
    
    Tc, Pc, Vc = float('nan'), float('nan'), float('nan')
    latex_steps = rf"$$\text{{Input SMILES: \texttt{{{smiles_str.replace(' ', '').replace('\\', '\\\\')}}}}}$$"
    latex_steps += rf"$$\text{{Calculation Order (W): {W}}}$$"

    try:
        tc_sum_total = sum_tc1k + (W * sum_tc2j)
        Tc = 181.128 * math.log(tc_sum_total) if tc_sum_total > 0 else float('nan')
        pc_sum_total = sum_pc1k + (W * sum_pc2j)
        pc_base = pc_sum_total + 0.10022
        Pc = (pc_base ** -2) + 1.3705 if pc_base != 0 else float('inf')
        if isinstance(pc_base, complex) or (pc_base < 0 and (-2 % 1 != 0)): Pc = float('nan')
        vc_sum_total = sum_vc1k + (W * sum_vc2j)
        Vc = (-0.00435 + vc_sum_total) * 1000.0

        latex_steps += rf"$$\sum N_k(tc1k) = {sum_tc1k:.4f}; \quad W \cdot \sum M_j(tc2j) = {W * sum_tc2j:.4f}$$"
        latex_steps += rf"$$T_c = 181.128 \cdot \ln({tc_sum_total:.4f}) = {Tc:.2f} \, K$$"
        latex_steps += rf"$$\sum N_k(pc1k) = {sum_pc1k:.4f}; \quad W \cdot \sum M_j(pc2j) = {W * sum_pc2j:.4f}$$"
        latex_steps += rf"$$P_c = ({pc_sum_total:.4f} + 0.10022)^{{-2}} + 1.3705 = {Pc:.2f} \, bar$$"
        latex_steps += rf"$$\sum N_k(vc1k) = {sum_vc1k:.6f}; \quad W \cdot \sum M_j(vc2j) = {W * sum_vc2j:.6f}$$"
        latex_steps += rf"$$V_c = (-0.00435 + {vc_sum_total:.6f}) \cdot 1000 = {Vc:.1f} \, cm^3/mol$$"

    except Exception as e:
        print(f"Error in CG formula calculation: {e}")
        error_payload_cg["error"] = f"Calculation error: {e}"
        return jsonify(error_payload_cg), 500

    print(f"CG Results: Tc={Tc:.2f} K, Pc={Pc:.2f} bar, Vc={Vc:.1f} cm3/mol")
    print("--- CG CRITICAL PROPS API END ---")

    return jsonify({
        'Tc': f"{Tc:.2f}" if not (math.isnan(Tc) or math.isinf(Tc)) else "Error",
        'Pc': f"{Pc:.2f}" if not (math.isnan(Pc) or math.isinf(Pc)) else "Error",
        'Vc': f"{Vc:.1f}" if not (math.isnan(Vc) or math.isinf(Vc)) else "Error",
        'latex_steps': latex_steps,
        'identified_first_order_groups': identified_first_order_groups,
        'identified_second_order_groups': identified_second_order_groups,
        'paper_citation_data': cg_paper_citation_data,    
        'tool_citation_data': cg_tool_citation_data,      
        'applicability_notes': cg_applicability_notes     
    })