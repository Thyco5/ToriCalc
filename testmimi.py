# testmimi.py (Comprehensive for First-Order Groups)
import requests
import json 

def test_cg_api(smiles_string, chemical_name_for_log="", expected_groups=None):
    print(f"\n--- testmimi.py: Starting CG API call for: {chemical_name_for_log} ({smiles_string}) ---")
    if expected_groups:
        print(f"    Expecting: {expected_groups}")
    try:
        response = requests.post(
            "http://127.0.0.1:5000/api/calculate_cg_critical",
            json={
                "smiles": smiles_string,
                "calculation_order": "0"  # First-order only for these tests
            },
            timeout=30 # Increased timeout slightly for complex molecules if many SMARTS
        )
        print(f"Status Code for {chemical_name_for_log}: {response.status_code}")
        
        try:
            response_data = response.json()
            # Uncomment below to see full JSON for debugging individual cases
            # print(f"Parsed Response JSON for {chemical_name_for_log}:\n{json.dumps(response_data, indent=2)}") 
            
            identified_groups = response_data.get("identified_first_order_groups", {})
            # Sort the identified groups for consistent comparison if needed, though dict comparison works
            # identified_groups = dict(sorted(identified_groups.items()))
            # if expected_groups:
            #     expected_groups = dict(sorted(expected_groups.items()))
            
            print(f"==> IDENTIFIED FIRST ORDER GROUPS for {chemical_name_for_log}: {identified_groups}")
            
            if expected_groups:
                match = True
                # Check if all expected groups are present with the correct count
                for group, count in expected_groups.items():
                    if identified_groups.get(group) != count:
                        match = False
                        print(f"    MISMATCH: Expected {group}: {count}, Got: {identified_groups.get(group)}")
                        break
                # Check if there are any unexpected groups identified
                if match: # Only if previous check passed
                    for group, count in identified_groups.items():
                        if group not in expected_groups:
                            match = False
                            print(f"    UNEXPECTED GROUP: Identified {group}: {count}")
                            break
                # Final check on length if all individual items matched but lengths are different (e.g. an expected group was missing entirely)
                if match and len(identified_groups) != len(expected_groups):
                    match = False
                    print(f"    MISMATCH: Lengths differ. Expected {len(expected_groups)} groups, Got {len(identified_groups)} groups.")


                print(f"    MATCH EXPECTED: {'YES!' if match else 'NO! Review SMARTS.'}")

            if "Tc" in response_data and response_data["Tc"] != "Error" and response_data["Tc"] is not None:
                 print(f"    Calculated Tc: {response_data['Tc']} K, Pc: {response_data['Pc']} bar, Vc: {response_data['Vc']} cm3/mol")
            else:
                print(f"    Calculation resulted in Error or NaN for Tc.")

        except requests.exceptions.JSONDecodeError:
            print(f"Raw Response Text for {chemical_name_for_log} (Not valid JSON):", response.text)

    except requests.exceptions.ConnectionError as e:
        print(f"!!! CONNECTION ERROR for {chemical_name_for_log} - Is Flask app running? Error: {e}")
    except requests.exceptions.Timeout as e:
        print(f"!!! TIMEOUT ERROR for {chemical_name_for_log} - Error: {e}")
    except Exception as e:
        print(f"!!! UNEXPECTED SCRIPT ERROR for {chemical_name_for_log} - Error: {e}, {type(e)}")
    print(f"--- Test for {chemical_name_for_log} finished ---")

# --- Test Cases with REVISED EXPECTATIONS based on Table C-2 ---
print("======================================================================")
print("RUNNING REVISED CG FIRST-ORDER GROUP ID TESTS (After SMARTS Finalization)")
print("======================================================================")

# --- Cases that should NOT find groups from Table C-2 ---
print("\n--- Testing Small Molecules Not in Table C-2 ---")
test_cg_api("C", "Methane", {}) # Not in Table C-2
test_cg_api("C=O", "Formaldehyde", {}) # Not in Table C-2
test_cg_api("C=C=C", "Allene (Propadiene)", {}) # CH2=C=CH is for R-CH=C=CH2, not propadiene itself

# --- General Aliphatics & Simple Alcohol (Expectations based on Table C-2 & successful SMARTS) ---
print("\n--- Testing General Aliphatics & Simple Alcohol ---")
test_cg_api("CC", "Ethane", {"CH3": 2})
test_cg_api("CCC", "Propane", {"CH3": 2, "CH2": 1})
test_cg_api("CCCC", "Butane", {"CH3": 2, "CH2": 2})
test_cg_api("CC(C)C", "Isobutane", {"CH3": 3, "CH": 1})
test_cg_api("CC(C)(C)C", "Neopentane", {"CH3": 4, "C": 1}) # C_GEN_SMARTS should now find the C
test_cg_api("CCO", "Ethanol", {"CH3": 1, "CH2": 1, "OH": 1}) # OH recursive, CH2 general
test_cg_api("CC(C)O", "Isopropanol", {"CH3": 2, "CH": 1, "OH": 1}) # OH recursive, CH general
test_cg_api("CC(C)(C)O", "tert-Butanol", {"CH3": 3, "C": 1, "OH": 1}) # OH recursive, C general

# --- Olefins (Table C-2 groups) ---
print("\n--- Testing Olefins ---")
test_cg_api("CC=C", "Propene", {"CH3": 1, "CH2=CH": 1})      # CH2=CH (1)
test_cg_api("CC=CC", "2-Butene", {"CH3": 2, "CH=CH": 1})     # CH=CH (2)
test_cg_api("CC(=C)C", "Isobutene", {"CH3": 2, "CH2=C": 1})  # CH2=C (2)
test_cg_api("CC=C(C)C", "2-Methyl-2-butene", {"CH3": 3, "CH=C": 1}) # CH=C (3)
test_cg_api("CC=C=C", "1,2-Butadiene (Methylallene)", {"CH3": 1, "CH2=C=CH": 1}) # CH2=C=CH (1)

# --- Aromatics & Linkages (Table C-2 groups) ---
print("\n--- Testing Aromatics & Linkages ---")
test_cg_api("c1ccccc1", "Benzene", {"ACH": 6})               # ACH (2)
test_cg_api("Cc1ccccc1", "Toluene", {"ACCH3": 1, "ACH": 5, "AC": 1}) # ACCH3 (2), ACH (2), AC (3)
test_cg_api("CCc1ccccc1", "Ethylbenzene", {"CH3": 1, "ACCH2": 1, "ACH": 5, "AC": 1}) # ACCH2 (3)
test_cg_api("CC(C)c1ccccc1", "Cumene", {"CH3": 2, "ACCH": 1, "ACH": 5, "AC": 1})   # ACCH (4)
test_cg_api("Oc1ccccc1", "Phenol", {"ACOH (2)": 1, "ACH": 5, "AC": 1})             # ACOH (2)

# --- Aldehydes & Ketones (Table C-2 groups) ---
print("\n--- Testing Aldehydes & Ketones ---")
test_cg_api("CC=O", "Acetaldehyde", {"CH3": 1, "CHO (1)*": 1}) # CHO (1)*
test_cg_api("CCC=O", "Propanal", {"CH3": 1, "CH2": 1, "CHO (1)*": 1})
test_cg_api("O=Cc1ccccc1", "Benzaldehyde", {"ACH": 5, "AC": 1, "CHO (1)*": 1}) # CHO (1)*
test_cg_api("CC(=O)C", "Acetone", {"CH3CO (1)": 2})          # CH3CO (1)
test_cg_api("CCC(=O)C", "2-Butanone", {"CH3CO (1)": 1, "CH2CO (2)": 1, "CH3": 1}) # CH2CO (2)
test_cg_api("CCC(=O)CC", "3-Pentanone", {"CH2CO (2)": 2, "CH3": 2})

# --- Esters (Table C-2 groups) ---
print("\n--- Testing Esters ---")
test_cg_api("CC(=O)OCC", "Ethyl Acetate", {"CH3COO (1)": 1, "CH2": 1, "CH3": 1}) # CH3COO (1)
test_cg_api("CC(=O)OC", "Methyl Acetate", {"CH3COO (1)": 1, "CH3": 1})
test_cg_api("CCC(=O)OCC", "Ethyl Propanoate", {"CH2COO (2)": 1, "CH3": 2, "CH2": 1}) # CH2COO (2)
test_cg_api("COC=O", "Methyl Formate", {"HCOO (1)": 1, "CH3": 1})                 # HCOO (1)
test_cg_api("O=COCC", "Ethyl Formate", {"HCOO (1)": 1, "CH2": 1, "CH3": 1})

# --- Ethers (Table C-2 groups) ---
print("\n--- Testing Ethers ---")
test_cg_api("COC", "Dimethyl Ether", {"CH3O (1)": 2})       # CH3O (1)
test_cg_api("COCC", "Methyl Ethyl Ether", {"CH3O (1)": 1, "CH2O (2)": 1, "CH3": 1}) # CH2O (2)
test_cg_api("CCOCC", "Diethyl Ether", {"CH2O (2)": 2, "CH3": 2})
test_cg_api("CC(C)OC(C)C", "Diisopropyl Ether", {"CH=O (3)": 2, "CH3": 4}) # CH=O (3) - ether methine
test_cg_api("FCOC", "Fluoromethyl methyl ether", {"FCH2O (1)*": 1, "CH3O (1)": 1}) # FCH2O (1)*
test_cg_api("COc1ccccc1", "Anisole", {"CH3": 1, "AC": 1, "ACH": 5}) # No Ar-O-CH3 in Table C-2. Decomposes to CH3 + AC(for O-C<) + 5 ACH. Oxygen is implicit.

# --- Amines (Table C-2 groups) ---
print("\n--- Testing Amines ---")
test_cg_api("CCN", "Ethylamine", {"CH2NH2 (1)": 1, "CH3": 1})    # CH2NH2 (1)
test_cg_api("CC(N)C", "Isopropylamine", {"CHNH2 (2)": 1, "CH3": 2}) # CHNH2 (2)
test_cg_api("CCNCC", "Diethylamine", {"CH2NH (3)": 2, "CH3": 2})   # CH2NH (3)
test_cg_api("CNCCC", "N-Methylpropylamine", {"CH3NH (2)":1, "CH2NH (3)":1, "CH2":1, "CH3":1}) # CH3NH (2)
test_cg_api("CC(C)NC(C)C", "Diisopropylamine", {"CHNH (4)*": 2, "CH3": 4}) # CHNH (4)*
test_cg_api("CN(C)C", "Trimethylamine", {"CH3N (2)": 3})            # CH3N (2)
test_cg_api("CCN(C)C", "N,N-Dimethylethylamine", {"CH2N (3)":1, "CH3N (2)":2, "CH3":1}) # CH2N (3)
test_cg_api("Nc1ccccc1", "Aniline", {"ACNH2 (2)": 1, "ACH": 5, "AC": 1}) # ACNH2 (2)

# --- COOH, Nitriles, Nitro (Table C-2 groups) ---
print("\n--- Testing COOH, Nitriles, Nitro ---")
test_cg_api("CC(=O)O", "Acetic Acid", {"COOH (1)": 1, "CH3": 1}) # COOH (1)
test_cg_api("CCC#N", "Propionitrile", {"CH2CN (1)*": 1, "CH3": 1}) # CH2CN (1)* (Table C-2 uses CH2CN not R-CH2CN)
test_cg_api("CC[N+](=O)[O-]", "Nitroethane", {"CH2NO2 (1)*":1, "CH3":1}) # CH2NO2 (1)*
test_cg_api("CC([N+](=O)[O-])C", "2-Nitropropane", {"CHNO2 (2)*":1, "CH3":2}) # CHNO2 (2)*
test_cg_api("c1ccc(cc1)[N+](=O)[O-]", "Nitrobenzene", {"ACNO2 (2)*":1, "ACH":5, "AC":1}) # ACNO2 (2)*

# --- Halogenated Compounds (Table C-2 groups) ---
print("\n--- Testing Halogenated Compounds ---")
test_cg_api("CCCCl", "1-Chloropropane", {"CH2Cl (1)":1, "CH2":1, "CH3":1}) # CH2Cl (1)
test_cg_api("CC(Cl)C", "2-Chloropropane", {"CHCl (2)":1, "CH3":2})      # CHCl (2)
test_cg_api("CC(C)(Cl)C", "tert-Butyl chloride", {"CCl (3)":1, "CH3":3})# CCl (3)
test_cg_api("ClCCCl", "1,2-Dichloroethane", {"CH2Cl (1)":2}) # Should now work with refined CH2Cl (1) SMARTS
# Specific polyhalogenated not in Table C-2 as primary first order groups (CH2Cl2, CHCl3, CCl4)
# They are either second order, or their properties are built from C, H, Cl atoms if no other group fits.
# For now, expect them to decompose to simpler existing groups if possible, or not match.
# Dichloromethane, Chloroform, CCl4 are *not* in Table C-2 as simple first-order R-X groups.
# They might be covered by second-order effects or combinations of simpler groups if any apply.
# Let's assume they are not directly found as unique 1st order if not explicitly R-CXnYm in Table C-2.
test_cg_api("C(Cl)Cl", "Dichloromethane (CH2Cl2)", {}) # Expect CH2 + 2x Cl(1) IF Cl(1) existed as generic.
                                                      # Table C-2 has CH2Cl(1), CHCl(2), CCl(3), CHCl2(1)*, CCl2(2), CCl3(1).
                                                      # For CH2Cl2: it's two CH2Cl groups if ClCH2-Cl? No. It's one C.
                                                      # Best fit might be CH2 + specific halide interaction if we had Cl(1).
                                                      # Since we don't, and no CH2Cl2 group, expect {}.
test_cg_api("C(Cl)(Cl)Cl", "Chloroform (CHCl3)", {}) # Similarly, expect {}
test_cg_api("CC(Cl)(Cl)Cl", "1,1,1-Trichloroethane", {"CCl3 (1)":1, "CH3":1}) # CCl3 (1)
test_cg_api("Clc1ccccc1", "Chlorobenzene", {"ACCl (2)":1, "ACH":5, "AC":1}) # ACCl (2)

# Fluorinated (Table C-2 has CF3, CF2, CF, FCH2O, ACF, FSpecial)
test_cg_api("Fc1ccccc1", "Fluorobenzene", {"ACF (2)":1, "ACH":5, "AC":1}) # ACF (2)
test_cg_api("CCCF", "1-Fluoropropane", {"FSpecial (1)":1, "CH2":2, "CH3":1}) # Primary F on CH2 not FCH2O
test_cg_api("CC(F)C", "2-Fluoropropane", {"CF (3)":1, "CH3":2}) # CH-F should be CF(3)
test_cg_api("CC(F)(F)F", "1,1,1-Trifluoroethane", {"CF3 (1)":1, "CH3":1}) # CF3 (1)
test_cg_api("FC(F)(F)C(F)(F)C(F)(F)F", "Perfluoropropane (n-C3F8)", {"CF3 (1)":2, "CF2 (2)":1}) # CF2 (2)

# --- Amides (Table C-2 groups) ---
print("\n--- Testing Amides ---")
test_cg_api("CC(=O)N", "Acetamide", {"CONH2 (1)*":1, "CH3":1}) # CONH2 (1)*
test_cg_api("CC(=O)NC", "N-Methylacetamide", {"CONHCH3 (1)*":1, "CH3":1}) # CONHCH3 (1)*
test_cg_api("CC(=O)N(C)C", "N,N-Dimethylacetamide", {"CON(CH3)2 (1)*":1, "CH3":1}) # CON(CH3)2 (1)*
test_cg_api("O=CN(C)C", "Dimethylformamide", {"HCON(CH3)2 (2)*":1}) # HCON(CH3)2 (2)*
test_cg_api("O=C1CCCN1", "Pyrrolidinone (Lactam)", {"CON(CH2)2 (3)*": 1, "CH2": 3}) # CON(CH2)2 (3)*

print("\n--- testmimi.py: All REVISED first-order core tests finished ---")