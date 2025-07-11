# testmimi_2nd_order.py (or replace your existing testmimi.py)
import requests
import json 

# MODIFIED test_cg_api function to accept and use calculation_order
def test_cg_api(smiles_string, chemical_name_for_log="", expected_groups=None, calculation_order="0"): # Added calculation_order
    print(f"\n--- testmimi: Starting CG API call for: {chemical_name_for_log} ({smiles_string}) ---")
    print(f"    Calculation Order: {calculation_order}")
    if expected_groups:
        print(f"    Expecting: {expected_groups}")
    
    payload = {
        "smiles": smiles_string,
        "calculation_order": calculation_order # Use the passed argument
    }

    try:
        response = requests.post(
            "http://127.0.0.1:5000/api/calculate_cg_critical",
            json=payload, # Use the payload dictionary
            timeout=30 
        )
        print(f"Status Code for {chemical_name_for_log}: {response.status_code}")
        
        try:
            response_data = response.json()
            # Uncomment to see full JSON for debugging
            # print(f"Parsed Response JSON for {chemical_name_for_log}:\n{json.dumps(response_data, indent=2)}") 
            
            identified_first_order = response_data.get("identified_first_order_groups", {})
            identified_second_order = response_data.get("identified_second_order_groups", {})
            
            # Combine for easier comparison if needed, or compare separately
            all_identified_groups = {**identified_first_order, **identified_second_order}
            
            print(f"==> IDENTIFIED FIRST ORDER GROUPS for {chemical_name_for_log}: {identified_first_order}")
            if identified_second_order: # Only print if there are second order groups
                 print(f"==> IDENTIFIED SECOND ORDER GROUPS for {chemical_name_for_log}: {identified_second_order}")
            
            if expected_groups:
                match = True
                # Check if all expected groups are present with the correct count
                for group, count in expected_groups.items():
                    if all_identified_groups.get(group) != count:
                        match = False
                        print(f"    MISMATCH: Expected {group}: {count}, Got: {all_identified_groups.get(group)}")
                        break
                # Check if there are any unexpected groups identified
                if match: 
                    for group, count in all_identified_groups.items():
                        if group not in expected_groups:
                            match = False
                            print(f"    UNEXPECTED GROUP: Identified {group}: {count}")
                            break
                # Final check on length
                if match and len(all_identified_groups) != len(expected_groups):
                    match = False
                    print(f"    MISMATCH: Lengths differ. Expected {len(expected_groups)} groups, Got {len(all_identified_groups)} groups.")

                print(f"    MATCH EXPECTED: {'YES!' if match else 'NO! Review SMARTS.'}")

            if "Tc" in response_data and response_data["Tc"] != "Error" and response_data["Tc"] is not None:
                 print(f"    Calculated Tc: {response_data['Tc']} K, Pc: {response_data['Pc']} bar, Vc: {response_data['Vc']} cm3/mol, Omega: {response_data.get('omega', 'N/A')}")
            else:
                print(f"    Calculation resulted in Error or NaN for Tc.")

        except requests.exceptions.JSONDecodeError:
            print(f"Raw Response Text for {chemical_name_for_log} (Not valid JSON):", response.text)
        except Exception as e_json: # Catch other errors during JSON processing
            print(f"!!! ERROR processing response JSON for {chemical_name_for_log} - Error: {e_json}, {type(e_json)}")
            print(f"Raw Response Text: {response.text[:500]}...") # Print first 500 chars of raw response


    except requests.exceptions.ConnectionError as e:
        print(f"!!! CONNECTION ERROR for {chemical_name_for_log} - Is Flask app running? Error: {e}")
    except requests.exceptions.Timeout as e:
        print(f"!!! TIMEOUT ERROR for {chemical_name_for_log} - Error: {e}")
    except Exception as e:
        print(f"!!! UNEXPECTED SCRIPT ERROR for {chemical_name_for_log} - Error: {e}, {type(e)}")
    print(f"--- Test for {chemical_name_for_log} finished ---")

# --- Test Cases for First-Order Groups (from your previous successful testmimi.py) ---
# (You can copy-paste your full suite of first-order tests here if you want one consolidated script)
# For brevity here, I'll just put a few representative ones.
print("======================================================================")
print("RUNNING CG FIRST-ORDER GROUP IDENTIFICATION TESTS (CONFIRMATION)")
print("======================================================================")
test_cg_api("CC", "Ethane", {"CH3": 2}, calculation_order="0")
test_cg_api("CCO", "Ethanol", {"CH3": 1, "CH2": 1, "OH": 1}, calculation_order="0")
test_cg_api("O=Cc1ccccc1", "Benzaldehyde", {"ACH": 5, "AC": 1, "CHO (1)*": 1}, calculation_order="0")
test_cg_api("FC(F)(F)C(F)(F)C(F)(F)F", "Perfluoropropane (n-C3F8)", {"CF3 (1)":2, "CF2 (2)":1}, calculation_order="0")
test_cg_api("O=C1CCCN1", "Pyrrolidinone (Lactam)", {"CON(CH2)2 (3)*": 1, "CH2": 3}, calculation_order="0")


# --- Test Cases for Second-Order Groups ---
print("\n\n======================================================================")
print("RUNNING CONSTANTINOU-GANI SECOND-ORDER GROUP IDENTIFICATION TESTS")
print("======================================================================")

# Example from Table C-4: (CH3)2CH
test_cg_api(
    "CC(C)CCC", 
    "2-Methylpentane", 
    expected_groups={
        "CH3": 3, "CH2": 2, "CH": 1, 
        "(CH3)2CH": 1
    },
    calculation_order="1" # IMPORTANT
)

# Example from Table C-4: (CH3)3C
test_cg_api(
    "CC(C)(C)CC", 
    "Neohexane (2,2-Dimethylbutane)", 
    expected_groups={
        "CH3": 4, "CH2": 1, "C": 1,
        "(CH3)3C": 1 
    },
    calculation_order="1"
)

# Example from Table C-4: CH(CH3)CH(CH3)
test_cg_api(
    "CC(C)C(C)C", 
    "2,3-Dimethylbutane", 
    expected_groups={
        "CH3": 4, "CH": 2, 
        "(CH3)2CH": 2,  # ADDED THIS
        "CH(CH3)CH(CH3)": 1 # Keep this expectation at 1
    },
    calculation_order="1"
)

# Example for Rings
test_cg_api(
    "C1CC1", 
    "Cyclopropane", 
    expected_groups={
        "CH2": 3, 
        "3 membered ring": 1
    },
    calculation_order="1"
)

test_cg_api(
    "C1CCCCC1", 
    "Cyclohexane", 
    expected_groups={
        "CH2": 6, 
        "6 membered ring": 1
    },
    calculation_order="1"
)

# Example: CHn=CHm-CHp=CHk (Conjugated diene)
test_cg_api(
    "C=CC=C", 
    "1,3-Butadiene", 
    expected_groups={
        "CH2=CH": 2, 
        "CHn=CHm-CHp=CHk m,p E (0,1), k,n E (0,2)": 1
    },
    calculation_order="1"
)

# Example: CH3-CHm=CHn (Allylic methyl)
test_cg_api(
    "CC=C", 
    "Propene (for allylic CH3)", 
    expected_groups={
        "CH3": 1, "CH2=CH": 1,
        "CH3-CHm=CHn m E (0,1), n E (0,2)": 1
    },
    calculation_order="1"
)
test_cg_api(
    "CC=C(C)C", 
    "2-Methyl-2-butene (for allylic CH3)", 
    expected_groups={'CH=C': 1, 
        'CH3-CHm=CHn m E (0,1), n E (0,2)': 3
    },
    calculation_order="1"
)

# Example: CH3CH3 (Ethane only)
test_cg_api(
    "CC", 
    "Ethane (for CH3CH3 2nd order)", 
    expected_groups={
        "CH3": 2, 
        "CH3CH3": 1
    },
    calculation_order="1"
)

# Example: Ccyclic=O
# For Cyclohexanone O=C1CCCCC1
# 1st Order: CH2 (5) (C of C=O is not identified by general C due to !$(C=O) exclusion)
# 2nd Order: Ccyclic=O (1), 6 membered ring (1)
test_cg_api(
    "O=C1CCCCC1", 
    "Cyclohexanone (with ring)", 
    expected_groups={
        "CH2": 5, 
        "Ccyclic=O": 1,
        "6 membered ring": 1
    },
    calculation_order="1"
)

# Example: CHOH (Secondary alcohol)
test_cg_api(
    "CCC(O)C", 
    "2-Butanol", 
    expected_groups={
        "CH3": 2, "CH2": 1, "CH": 1, "OH": 1,
        "CHOH": 1
    },
    calculation_order="1"
)

print("\n--- testmimi_2nd_order.py: All tests finished ---")