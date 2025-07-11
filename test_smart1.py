# test_smart1.py
from rdkit import Chem

def test_single_smarts(molecule_smiles, group_name, smarts_pattern, expected_mol_for_log=""):
    """Helper function to test a single SMARTS pattern against a molecule."""
    if not expected_mol_for_log:
        expected_mol_for_log = molecule_smiles
    mol = Chem.MolFromSmiles(molecule_smiles)
    if not mol:
        print(f"  ERROR: Could not parse SMILES '{molecule_smiles}' for group '{group_name}' test.")
        return 0, [] 
    mol = Chem.AddHs(mol)
    print(f"  Testing '{group_name}' on {expected_mol_for_log} ('{molecule_smiles}') with SMARTS: '{smarts_pattern}'")
    count = 0
    matches_found = []
    try:
        pattern = Chem.MolFromSmarts(smarts_pattern)
        if pattern:
            matches_found = mol.GetSubstructMatches(pattern)
            count = len(matches_found)
            print(f"    Matches: {matches_found} (Count: {count})")
        else:
            print(f"    INVALID SMARTS pattern for '{group_name}' (pattern is None)")
    except Exception as e:
        print(f"    ERROR processing SMARTS for '{group_name}': {e}")
    return count, matches_found

# ==============================================================================
# --- GOLDEN/MATURE SMARTS (Last Update: After "Halogen Cleanup & Nitrogen Intro (Round 1)" successful specific group tests) ---
# ==============================================================================
# Esters:
#   "CH3COO (1)": "[CH3X4][CX3;!c](=[O])[O;X2H0;!R]"
#   "CH2COO (2)": "[CH2X4;!$(C=O);!R;$(*-[#6;!a])][CX3;!c](=[O])[O;X2H0;!R]"
#   "HCOO (1)":   "[CH1X3H1;!R](=[O])[O;X2H0;!R]"
# Ethers:
#   "CH3O (1)":   "[CH3X4;!R;$(*-[OX2H0;!R;D2]-[#6;!R;!c;!$([CX3](=[OX1]))])]"
#   "CH2O (2)":   "[CH2X4;!R;$(*-[OX2H0;!R;D2;!$(*-[c]);!$(*-[CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))])]"
#   "CH=O (3)":   "[CH1X4;!R;$(*-[OX2H0;!R;D2;!$(*-[c]);!$(*-[CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))]);$(*-[#6;!R;!c;!$([CX3](=[OX1]))])]"
#   "FCH2O (1)*": "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]" 
# Ketones:
#   "CH3CO (1)":  "[CH3X4;!R;$(*-[CX3R0;!c;!H1;!$(*-[O;X2H0;!R]-[#6])]=[OX1D1])]"
#   "CH2CO (2)":  "[CH2X4H2;!R;$(*-[CX3R0;!c;!H1;!$(*-[O;X2H0;!R]-[#6])]=[OX1D1]);$(*-[#6;!a;!$(C=O)])]"
# Aldehyde:
#   "CHO (1)*":   "[CH1X3H1R0;!$(*~[OX2H0])](=[OX1D1])"
# Alkyl Chlorides (Simple & R-Attached):
#   "CH2Cl (1)":  "[CH2X4H2;!R]([ClX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHCl (2)":   "[CH1X4H1;!R]([ClX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CCl (3)":    "[CX4H0;!R]([ClX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Alkyl Bromides:
#   "CH2Br (1)": "[CH2X4H2;!R]([BrX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHBr (2)":  "[CH1X4H1;!R]([BrX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CBr (3)":   "[CX4H0;!R]([BrX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Alkyl Iodides:
#   "CH2I (1)":  "[CH2X4H2;!R]([IX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHI (2)":   "[CH1X4H1;!R]([IX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CI (3)":    "[CX4H0;!R]([IX1])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Specific Polyhalogenated (Terminal):
#   "CH2Cl2_explicit": "[CH2X4H2](Cl)(Cl)" 
#   "CHCl3_explicit":  "[CH1X4H1](Cl)(Cl)(Cl)"
#   "CCl4_explicit":   "[CX4H0](Cl)(Cl)(Cl)(Cl)"
# General CCl3 unit:
#   "CCl3_unit":       "[CX4H0](Cl)(Cl)(Cl)" 
# Aromatic Halides:
#   "ACF":   "[F]-[c]"
#   "ACCl":  "[Cl]-[c]" 
#   "ACBr":  "[Br]-[c]"
#   "ACI":   "[I]-[c]"
# Mixed Cl/F (Explicit Terminal):
#   "CHFCl2_explicit": "[CH1X4H1](F)(Cl)(Cl)" 
#   "CFCl3_explicit":  "[CX4H0](F)(Cl)(Cl)(Cl)" 
#   "CF2Cl2_explicit": "[CX4H0](F)(F)(Cl)(Cl)" 
# Mixed Cl/F (Alkyl-Attached):
#   "CCl2F-R":   "[CX4H0](Cl)(Cl)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHFCl-R":   "[CH1X4H1](Cl)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CClF2-R":   "[CX4H0](Cl)(F)(F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Polyfluorinated Aliphatics:
#   "CF3":   "[CX4H0](F)(F)(F)"
#   "CF2":   "[CX4;!$(C(F)(F)F)](F)(F)"
#   "CF":    "[CX4;!$(C(F)(F)F);!$(C(F)F)](F)"
# Primary Fluoride:
#   "CH2F": "[CH2X4H2](F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,Cl,Br,I])]"
# Primary Aliphatic Amines:
#   "CH2NH2 (1)": "[CH2X4H2]([NH2X3])-[#6;!R;!c]"
#   "CHNH2 (2)":  "[CH1X4H1]([NH2X3])(-[#6;!R;!c])-[#6;!R;!c]"
# Secondary Aliphatic Amines:
#   "CH3NH (2)": "[CH3X4;!R]([NH1X3H1]-[#6;!R;!c])"
#   "CH2NH (3)": "[CH2X4H2;!R]([NH1X3H1]-[#6;!R;!c])-[#6;!R;!c]" 
#   "CHNH (4)*": "[CH1X4H1;!R]([NH1X3H1]-[#6;!R;!c])(-[#6;!R;!c])-[#6;!R;!c]"
# Tertiary Aliphatic Amines:
#   "CH3N (2)": "[CH3X4;!R;$(*-N(-[#6;!R;!c])-[#6;!R;!c])]" 
#   "CH2N (3)": "[CH2X4H2;!R](N(-[#6;!R;!c])-[#6;!R;!c])-[#6;!R;!c]" 
#   "CHN_tert": "[CH1X4H1;!R](N(-[#6;!R;!c])-[#6;!R;!c])(-[#6;!R;!c])-[#6;!R;!c]"
# Aromatic Amine:
#   "ACNH2 (2)": "[NH2X3]-[c]"
# --- NEWLY MATURED (After Sulfur & Quaternary C Cleanup) ---
# Carboxylic Acid:
#   "COOH (1)": "[CX3;!c;!R](=[O])[OH1X2]"
# Quaternary Aliphatic Carbon (specific definition for all-alkyl non-functionalized):
#   "C_quat_all_alkyl":   "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]" (was C_quat_final)
# Thiols:
#   "CH2SH (1)": "[CH2X4H2]([SX2H1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Thioethers:
#   "CH3S (1)":  "[CH3X4;!R;$(*-S-[#6;!R;!c;!$(C=O)])]" 
#   "CH2S (2)":  "[CH2X4H2;!R]([SX2H0]-[#6;!R;!c;!$(C=O)])-[#6;!R;!c;!$(C=O)]" 
#   "CHS (3)*":  "[CH1X4H1;!R]([SX2H0]-[#6;!R;!c;!$(C=O)])(-[#6;!R;!c;!$(C=O)])-[#6;!R;!c;!$(C=O)]" 
# Thiophenes:
#   "C4H3S (1)" (Thiophene S atom): "[s;r5;H0X2](:c):c"
# --- NEWLY MATURED (After Phenol, Alkynes, Aromatic N-Rings Iteration) ---
# Phenol:
#   "ACOH (2)": "[OH1X2]-[c]"
# Alkynes:
#   "CH#C (1)": "[CH1X2H1]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
#   "C#C (2)*": "[#6;!a;!c;!R;!$([CX3](=[OX1]))]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
#   "Cl-C#C (3)*": "[ClX1]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
# Aromatic Nitrogen in Ring (Pyridine-type):
#   "C5H4N (1)" (and C5H3N (2)): "[n;r6;H0X2](:c):c"

# --- NEWLY MATURED (After Quaternary C, COOH, Nitriles, Nitro Iteration) ---
# Carboxylic Acid:
#   "COOH (1)": "[#6X3R0](=[OX1R0])[OH1X2R0]" (matches user's bp.py and works)
# Quaternary Aliphatic Carbon (specific all-alkyl):
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitriles (Alkyl):
#   "CH2CN (1)*": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitro Compounds:
#   "CH2NO2 (1)*" (Alkyl): "[CH2X4H2]([$(N(=O)=O),$(N(=[O])-[O-])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHNO2 (2)*" (Alkyl):  "[CH1X4H1]([$(N(=O)=O),$(N(=[O])-[O-])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "ACNO2 (2)*" (Aromatic): "[c]-[$(N(=O)=O),$(N(=[O])-[O-])]"

# --- NEWLY MATURED (After Quaternary C, COOH, Nitriles, Nitro Iteration) ---
# Carboxylic Acid:
#   "COOH (1)": "[#6X3R0](=[OX1R0])[OH1X2R0]"
# Quaternary Aliphatic Carbon (specific all-alkyl):
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitriles (Alkyl R-CH2-CN type):
#   "CH2CN (1)*": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitro Compounds:
#   "CH2NO2 (1)*" (Alkyl R-CH2-NO2): "[CH2X4H2]([$(N(=O)=O),$(N(=[O])-[O-])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "CHNO2 (2)*" (Alkyl R1R2CH-NO2):  "[CH1X4H1]([$(N(=O)=O),$(N(=[O])-[O-])])(-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   "ACNO2 (2)*" (Aromatic Ar-NO2): "[c]-[$(N(=O)=O),$(N(=[O])-[O-])]"

# --- NEWLY MATURED (After General C/CH2 Cleanup & Amide Intro - Round 1 - Specific Groups) ---
# Carboxylic Acid:
#   "COOH (1)": "[CX3;!c;!R](=[O])[OH1X2]" (was COOH_bp)
# Specific Quaternary Carbon:
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]" (was C_quat_final)
# Nitriles (Alkyl R-CH2-CN):
#   "CH2CN (1)*": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]" (was CH2CN_gold)
# Amides (Initial):
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]" (was CONH2_v1)
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]" (was CONHCH3_v1)

# --- NEWLY MATURED (After Amides Round 1 & General Aliphatic Simplification) ---
# Carboxylic Acid: (Confirmed from previous)
#   "COOH (1)": "[#6X3R0](=[OX1R0])[OH1X2R0]"
# Specific Quaternary Carbon: (Confirmed from previous)
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitriles (Alkyl R-CH2-CN type): (Confirmed from previous)
#   "CH2CN (1)*": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Amides (Initial positive matches):
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]"
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   "CONHCH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]" (was CONHCH2_v1)
#   "CON(CH3)2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]" (was CON_CH3_2_v1)
#   "CON(Et)(Me)-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]" (was CON_Et_Me_v1, for R-CO-N(Et)(Me))
#   "CON(CH2R)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[#6;!R;!c])-[CH2X4H2]-[#6;!R;!c]" (was CON_CH2R_2_v1 for R-CO-N(CH2R1)(CH2R2))

# --- NEWLY MATURED (After Amides Round 1 / General C & CH2 Cleanup attempt) ---
# Carboxylic Acid:
#   "COOH (1)": "[CX3;!c;!R](=[O])[OH1X2]"
# Specific Quaternary Carbon:
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Nitriles (Alkyl R-CH2-CN type):
#   "CH2CN (1)*": "[CH2X4H2]([CX2]#[NX1])-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Amides (Initial positive matches - these are the ones we want general aliphatics to yield to):
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]"
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   "CONHCH2-R": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]" (was CONHCH2_v1)
#   "CON(CH3)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]" (was CON_CH3_2_v1)
#   "CON(Et)(Me)-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]" (was CON_Et_Me_v1)
#   "CON(CH2R)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[#6;!R;!c])-[CH2X4H2]-[#6;!R;!c]" (was CON_CH2R_2_v1)

# --- NEWLY MATURED (After Amides Round 2) ---
# Amides (R-CO-NH2, R-CO-NHCH3, R-CO-NHCH2R', R-CO-N(CH3)2, R-CO-N(Et)(Me) ):
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]"
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   "CONHCH2-R": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]"
#   "CON(CH3)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CON(Et)(Me)-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]"
#   "CON(CH2R)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[#6;!R;!c])-[CH2X4H2]-[#6;!R;!c]" (This one for N,N-Diethylacetamide worked)

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish) ---
# General Aliphatics (pending full review against all functional groups, but very good progress):
#   "CH3_gen_final": "[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]))]"
#   "CH2_gen_final": "[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N])]"
#   "CH_gen_final":  "[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]))]"
#   "C_gen_final":   "[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$([CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])]" (The C_gen_final with its complex exclusion PARSED and WORKED!)
# Amides:
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]" (Re-confirmed)
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]" (Re-confirmed)
#   "CONHCH2-R": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]" (Re-confirmed)
#   "CON(CH3)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]" (Re-confirmed)
#   "CON(Et)(Me)-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]" (Re-confirmed)
#   "CON(CH2R)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[#6;!R;!c])-[CH2X4H2]-[#6;!R;!c]" (Re-confirmed)
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish) ---
# General Quaternary Carbon (All-Alkyl, specific):
#   "C_quat_all_alkyl": "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
# Dimethylformamide-type group:
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
# General Aliphatic Carbons (These are very robust now):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]))]" (Excludes N-CH2 from sec/tert amines/amides too)
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]" (Added general sec/tert amine CH exclusion)
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish) ---
# General Aliphatic Carbon C:
#   "C_gen": f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]" 
#   (where C_quat_all_alkyl_pattern_ref is the SMARTS for C_quat_all_alkyl)
# Dimethylformamide-type group:
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish) ---
# General Aliphatic Carbons (These are now considered highly robust):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]))]"
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"
# Dimethylformamide-type group:
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish) ---
# General Aliphatic Carbons (These are now considered highly robust after parsing success):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CH3_gen_final": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen_final": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]))]"
#   "CH_gen_final":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen_final":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"
# Amides (Specific examples that worked well):
#   "CONH2 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH2X3]"
#   "CONHCH3 (1)*": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   "CONHCH2-R": "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH2X4H2;!R]-[#6;!c;!R]"
#   "CON(CH3)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CON(Et)(Me)-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[CH3X4;!R])-[CH3X4;!R]"
#   "CON(CH2R)2-R": "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH2X4H2]-[#6;!R;!c])-[CH2X4H2]-[#6;!R;!c]"
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"

# --- NEWLY MATURED (After Phenol, Alkynes, Aromatic N-Rings, COOH Iteration) ---
# Phenol:
#   "ACOH (2)": "[OH1X2]-[c]"
# Carboxylic Acid:
#   "COOH (1)": "[#6X3R0](=[OX1R0])[OH1X2R0]"
# Alkynes:
#   "CH#C (1)": "[CH1X2H1]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
#   "C#C (2)*": "[#6;!a;!c;!R;!$([CX3](=[OX1]))]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
#   "Cl-C#C (3)*": "[ClX1]-[CX2H0]#[CX2H0]-[#6;!a;!c;!R;!$([CX3](=[OX1]))]"
# Aromatic Nitrogen in Ring (Pyridine-type):
#   "C5H4N (1)" (covers C5H3N (2) as well): "[n;r6;H0X2](:c):c"
# General Aliphatics (from last successful full parse):

# --- NEWLY MATURED (After Lactam Definition & Final General Aliphatic Review) ---
# General Aliphatic Carbons (Now considered FINAL and ROBUST):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]))]"
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"

# --- NEWLY MATURED (After Amide & General Aliphatic Final Polish - Iteration "Lactams & General CH3/C Polish" in script) ---
# General Aliphatic Carbons (Now considered FINAL and ROBUST after parsing success and exclusivity tests):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]" # N-Me part of R-CO-N(Me)(Et)
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]))]"
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"
# Dimethylformamide-type group:
#   "HCON(CH3)2 (2)*": "[CH1X3H1;!R](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"

# --- NEWLY MATURED (After Lactam Definition & Final General Aliphatic Review) ---
# General Aliphatic Carbons (FINALIZED - These are now considered very robust):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   CH2F_gold_pattern_ref = "[CH2X4H2](F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,Cl,Br,I])]"
#   FCH2O_gold_pattern_ref = "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]));!$({CH2F_gold_pattern_ref});!$({FCH2O_gold_pattern_ref})]"
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"
# Specific Fluorine group:
#   "FSpecial (1)": f"[F]-[#6;!c;!$({CH2F_gold_pattern_ref});!$({FCH2O_gold_pattern_ref});!$(C(F)(F));!$(C(F)(F)F)]"

# --- NEWLY MATURED (After Lactam Finalization & Pre-Mimi Review - Iteration: Gen Alph Final) ---
# General Aliphatic Carbons (FINALIZED - These are robust):
#   C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
#   CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
#   CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
#   CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
#   CH2F_gold_pattern_ref = "[CH2X4H2](F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,Cl,Br,I])]"
#   FCH2O_gold_pattern_ref = "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]"
#   "CH3_gen": f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
#   "CH2_gen": f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]));!$({CH2F_gold_pattern_ref});!$({FCH2O_gold_pattern_ref})]"
#   "CH_gen":  f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
#   "C_gen":   f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"
# Specific Fluorine group:
#   "FSpecial (1)": f"[F]-[#6;!c;!$({CH2F_gold_pattern_ref});!$({FCH2O_gold_pattern_ref});!$(C(F)(F));!$(C(F)(F)F)]"

## General Aliphatics (FINALIZED from previous run):
C_quat_all_alkyl_pattern_ref = "[CX4H0;!R](-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])(-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])])-[#6;!R;!c;!$([CX3]=[O]);!$(*~[#7,#8,#16,F,Cl,Br,I])]"
CONHCH3_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NH1X3H1]-[CH3X4;!R]"
CONMe_Et_pattern_N_Me_part_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH2X4H2;!R]"
CONMe2_pattern_ref = "[#6;!c;!R]-[CX3](=[O])-[NX3H0](-[CH3X4;!R])-[CH3X4;!R]"
CH2F_gold_pattern_ref = "[CH2X4H2](F)-[#6;!R;!c;!$([CX3](=[OX1]));!$(*~[#7,#8,#16,Cl,Br,I])]"
FCH2O_gold_pattern_ref = "[CH2X4;!R]([FX1])-[OX2H0;!R;!$(*~[CX3](=[OX1]))]"

# ==============================================================================

print("\n--- CG SMARTS DEVELOPMENT & INTEGRATION TEST (test_smart1.py) ---")
print("--- Iteration: Lactam Final Attempt v2 ---")

# Using the final general aliphatics from the previous successful run
CH3_gen_final = f"[CH3X4;!R;!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({CONHCH3_pattern_ref});!$({CONMe_Et_pattern_N_Me_part_ref});!$({CONMe2_pattern_ref})]"
CH2_gen_final = f"[CH2X4H2;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[C]#[N]);!$(*-[NH1X3H1]-[CH2X4H2;!R]);!$(*-[NX3H0](-[CH2X4H2;!R]));!$({CH2F_gold_pattern_ref});!$({FCH2O_gold_pattern_ref})]"
CH_gen_final  = f"[CH1X4H1;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$(*-[NH1X3H1]);!$(*-[NX3H0])]"
C_gen_final   = f"[CX4H0;!R;!$([CX3]=[O]);!$(*~[a]);!$(*~[#7,#8,#16,F,Cl,Br,I,Si,P]);!$(*-[CX3](=[O]));!$({C_quat_all_alkyl_pattern_ref})]"

smarts_under_test_lactam_final_v2 = {
    # --- General Aliphatics (for context) ---
    "CH3_gen_final": CH3_gen_final, # Using the string directly
    "CH2_gen_final": CH2_gen_final, # Using the string directly
    "CH_gen_final":  CH_gen_final,  # Using the string directly
    "C_gen_final":   C_gen_final,   # Using the string directly

    # --- LACTAM GROUP CON(CH2)2 (3)* ---
    # Nitrogen in a ring, part of an amide CONH, and connected to two CH2s also in the ring.
    "LactamN_v5": "[N;X3;H0;R;$(N([C;R](=[O;!R]))([CH2;R])[CH2;R])]", 
}

# --- Test Section 1: Lactam Group ---
print("\n\n--- SECTION 1: LACTAM GROUP ---")
print("\n--- Testing Pyrrolidinone (5-membered lactam): O=C1CCCN1 ---") 
test_single_smarts("O=C1CCCN1", "LactamN_v5", smarts_under_test_lactam_final_v2["LactamN_v5"]) # Exp 1 (the N atom)
test_single_smarts("O=C1CCCN1", "CH2_gen_final (ring CH2s)", smarts_under_test_lactam_final_v2["CH2_gen_final"]) 
# Expectation for CH2_gen_final in Pyrrolidinone (O=C-CH2-CH2-CH2-NH-):
# - CH2 alpha to C=O: Should be 0 (due to !$(*-[CX3](=[O])))
# - CH2 alpha to NH: Should be 0 (due to !$(*~[#7]))
# - CH2 beta to C=O and NH (the middle one): Should be 1
# So, total expected CH2_gen_final = 1

print("\n--- Testing Caprolactam (7-atom lactam ring, O=C1CCCCCN1) ---") 
test_single_smarts("O=C1CCCCCN1", "LactamN_v5", smarts_under_test_lactam_final_v2["LactamN_v5"]) # Exp 1 
test_single_smarts("O=C1CCCCCN1", "CH2_gen_final (ring CH2s)", smarts_under_test_lactam_final_v2["CH2_gen_final"])
# Expectation for CH2_gen_final in Caprolactam (O=C-CH2-CH2-CH2-CH2-CH2-NH-):
# - CH2 alpha to C=O: 0
# - CH2 alpha to NH: 0
# - The three CH2s in the middle: 3
# So, total expected CH2_gen_final = 3

print("\n--- Testing N-Methylpyrrolidinone: O=C1CCCN1C ---")
# LactamN_v5 expects H0 on N from ring bonds only. N-alkylation makes it not match.
test_single_smarts("O=C1CCCN1C", "LactamN_v5", smarts_under_test_lactam_final_v2["LactamN_v5"]) # Exp 0
test_single_smarts("O=C1CCCN1C", "CH3_gen_final (N-methyl on lactam)", smarts_under_test_lactam_final_v2["CH3_gen_final"]) # Exp 0 (due to !$(*~[#7]))

print("\n\n--- test_smart1.py finished ---")