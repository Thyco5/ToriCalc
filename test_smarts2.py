# test_smarts2.py - SECOND-ORDER DEBUGGING - CXnHmDk style
from rdkit import Chem
from rdkit.Chem import Draw
import os

# --- test_single_smarts function (keep as is) ---
def test_single_smarts(molecule_smiles, group_name_in_test, smarts_pattern, expected_match_count=None, expected_mol_for_log="", draw_matches=False):
    # ... (Full function from previous, with AddHs and atom info prints commented out for brevity unless needed) ...
    if not expected_mol_for_log:
        expected_mol_for_log = molecule_smiles
    mol = Chem.MolFromSmiles(molecule_smiles)
    if not mol:
        print(f"  ERROR: Could not parse SMILES '{molecule_smiles}' for group '{group_name_in_test}' test.")
        return 0, []
    mol_for_drawing = Chem.Mol(mol)
    mol = Chem.AddHs(mol) 

    print(f"\n  Testing '{group_name_in_test}' on {expected_mol_for_log} ('{molecule_smiles}') with SMARTS: '{smarts_pattern}'")
    
    count = 0
    matches_found_tuples = []
    try:
        pattern = Chem.MolFromSmarts(smarts_pattern)
        if pattern:
            matches_found_tuples = mol.GetSubstructMatches(pattern)
            count = len(matches_found_tuples)
            status = ""
            if expected_match_count is not None:
                if count == expected_match_count:
                    status = "OK"
                else:
                    status = f"FAIL (Expected {expected_match_count})"
            print(f"    Matches (atom index tuples): {matches_found_tuples} (Count: {count}) {status}")

            if draw_matches and matches_found_tuples: # Simplified drawing part
                output_dir = "test_outputs"
                if not os.path.exists(output_dir): os.makedirs(output_dir)
                filename_all = os.path.join(output_dir, f"{molecule_smiles.replace('/','_')}_{group_name_in_test}_all_matches.png")
                try:
                    all_match_atoms_flat = [atom_idx for tpl in matches_found_tuples for atom_idx in tpl]
                    img_all = Draw.MolToImage(mol, highlightAtoms=list(set(all_match_atoms_flat))) 
                    img_all.save(filename_all)
                    print(f"    Saved match image to {filename_all}")
                except Exception as e_draw: print(f"    Could not save image: {e_draw}")
        else:
            print(f"    INVALID SMARTS pattern for '{group_name_in_test}' (pattern is None): {smarts_pattern}")
    except Exception as e:
        print(f"    ERROR processing SMARTS for '{group_name_in_test}': {e}")
    return count, matches_found_tuples

print("\n--- CG SMARTS DEVELOPMENT & ISOLATED TEST (test_smarts2.py) ---")
print("--- Iteration: Second-Order SMARTS with no 'D' degree specifiers ---")

smarts_to_test_no_D = {
    "CH3CH3_noD": "[CH3X4]-[CH3X4]", # Known to work
    "(CH3)2CH_noD": "[CH1X4H1]([CH3X4])([CH3X4])[#6]", 
    "(CH3)3C_noD": "[CX4H0]([CH3X4])([CH3X4])([CH3X4])[#6]",
    "CH(CH3)CH(CH3)_noD": "[CH1X4H1]([CH3X4])-[CH1X4H1]([CH3X4])",
    "Allylic_CH3_noD": "[CH3X4]-[#6X3]=[#6X3]",
    "CHOH_noD": "[CH1X4H1]([OH1X2])([#6X4;!c;!$(C=[O,N,S])])([#6X4;!c;!$(C=[O,N,S])])",
}

print("\n\n--- Testing Second-Order SMARTS with no 'D' specifiers ---")
test_single_smarts("CC", "CH3CH3_noD", smarts_to_test_no_D["CH3CH3_noD"], expected_match_count=1, draw_matches=True)
test_single_smarts("CC(C)CCC", "(CH3)2CH_noD on 2-Methylpentane", smarts_to_test_no_D["(CH3)2CH_noD"], expected_match_count=1, draw_matches=True)
test_single_smarts("CC(C)(C)CC", "(CH3)3C_noD on Neohexane", smarts_to_test_no_D["(CH3)3C_noD"], expected_match_count=1, draw_matches=True)
test_single_smarts("CC(C)C(C)C", "CH(CH3)CH(CH3)_noD on 2,3-Dimethylbutane", smarts_to_test_no_D["CH(CH3)CH(CH3)_noD"], expected_match_count=1, draw_matches=True) # Hoping for 1 now
test_single_smarts("CC=C", "Allylic_CH3_noD on Propene", smarts_to_test_no_D["Allylic_CH3_noD"], expected_match_count=1, draw_matches=True)
test_single_smarts("CC=C(C)C", "Allylic_CH3_noD on 2-Methyl-2-butene", smarts_to_test_no_D["Allylic_CH3_noD"], expected_match_count=3, draw_matches=True)
test_single_smarts("CCC(O)C", "CHOH_noD on 2-Butanol", smarts_to_test_no_D["CHOH_noD"], expected_match_count=1, draw_matches=True)

print("\n\n--- test_smarts2.py (Second-Order no 'D') finished ---")