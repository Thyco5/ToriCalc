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

print("--- Starting CG First-Order SMARTS Test Session (Esters Refined v5) ---")

smarts_db_esters_v5 = {
    # Ester Core Groups (Winners from last test)
    "CH3COO (1)": "[CH3X4][CX3;!c](=[O])[OX2H0;!R]", 
    "CH2COO (2)": "[CH2X4;!R;!$(C=O);$(*-[#6;!a])][CX3;!c](=[O])[OX2H0;!R]", 
    "HCOO (1)":   "[CH1X3H1;!R](=[O])[OX2H0;!R]",       
    
    # Revised General Aliphatics (from our "mostly YES!" state)
    "CH3_gen_final": "[CH3X4;!R;!$(*-[c]);!$(*-[CX3;!H1]=[OX1]);!$(*~[#7,#16,F,Cl,Br,I])]",
    "CH2_gen_final": "[CH2X4H2;!R;!$(C=O);!$(*-[c]);!$(*-[CX3;!H1]=[OX1])]"
}

# --- Test CH3COO (1) ---
print("\n--- Testing CH3COO (1) with Ethyl Acetate: CC(=O)OCC ---")
# Ethyl Acetate: CH3-C(=O)-O-CH2-CH3. 
# Expect: 1 "CH3COO (1)", 1 "CH2_gen_final", 1 "CH3_gen_final"
test_single_smarts("CC(=O)OCC", "CH3COO (Ethyl Acetate)", smarts_db_esters_v5["CH3COO (1)"])
print("  Context for Ethyl Acetate:")
test_single_smarts("CC(=O)OCC", "  CH2_gen_final part", smarts_db_esters_v5["CH2_gen_final"])
test_single_smarts("CC(=O)OCC", "  CH3_gen_final part", smarts_db_esters_v5["CH3_gen_final"])

print("\n--- Testing CH3COO (1) with Methyl Acetate: CC(=O)OC ---")
# Methyl acetate: CH3-C(=O)-O-CH3. 
# Expect: 1 "CH3COO (1)", 1 "CH3_gen_final" (for the O-CH3 part)
test_single_smarts("CC(=O)OC", "CH3COO (Methyl Acetate)", smarts_db_esters_v5["CH3COO (1)"])
print("  Context for Methyl Acetate:")
test_single_smarts("CC(=O)OC", "  CH3_gen_final (O-CH3 part)", smarts_db_esters_v5["CH3_gen_final"])

# --- Test CH2COO (2) ---
print("\n--- Testing CH2COO (2) with Ethyl Propanoate: CCC(=O)OCC ---")
# Ethyl propanoate: CH3-CH2-C(=O)-O-CH2-CH3
# Expect: 1 "CH2COO (2)", 2 "CH3_gen_final", 1 "CH2_gen_final" (the CH2 from O-CH2-CH3)
test_single_smarts("CCC(=O)OCC", "CH2COO (Ethyl Propanoate)", smarts_db_esters_v5["CH2COO (2)"])
print("  Context for Ethyl Propanoate:")
test_single_smarts("CCC(=O)OCC", "  CH3_gen_final parts", smarts_db_esters_v5["CH3_gen_final"])
test_single_smarts("CCC(=O)OCC", "  CH2_gen_final (O-CH2- part)", smarts_db_esters_v5["CH2_gen_final"])


# --- Test HCOO (1) ---
print("\n--- Testing HCOO (1) with Methyl Formate: COC=O ---")
# Methyl formate: H-C(=O)-O-CH3. 
# Expect: 1 "HCOO", 1 "CH3_gen_final"
test_single_smarts("COC=O", "HCOO (Methyl Formate)", smarts_db_esters_v5["HCOO (1)"])
print("  Context for Methyl Formate:")
test_single_smarts("COC=O", "  CH3_gen_final part", smarts_db_esters_v5["CH3_gen_final"])

# Add to test_smarts.py
print("\n--- Testing CH3CO on Acetone ---")
smarts_ch3co = "[CH3X4][CX3R0;!c;!H1](=[OX1D1])"
test_single_smarts("CC(=O)C", "CH3CO (Acetone)", smarts_ch3co) # Expect 2

print("\n--- Ester Test Session Finished ---")