import traceback
from flask import Blueprint, render_template, request, jsonify, url_for

# Import the necessary functions and dataframes for the older library version
from chemicals.identifiers import search_chemical, int_to_CAS
from chemicals.heat_capacity import Cp_data_Poling

ideal_gas_bp = Blueprint('ideal_gas', __name__, template_folder='../templates')

@ideal_gas_bp.route('/calculators/ideal-gas-properties')
def ideal_gas_page():
    return render_template('ideal_gas_calculator.html')

@ideal_gas_bp.route('/api/ideal_gas/get_cp_coeffs', methods=['POST'])
def get_cp_coeffs_api():
    data = request.get_json() or {}
    chemical_name = data.get('chemical_name')
    if not chemical_name: return jsonify({"error": "Chemical name is required."}), 400

    try:
        metadata = search_chemical(chemical_name)
        if not metadata: raise ValueError(f"Could not identify '{chemical_name}'.")

        cas_rn_str = int_to_CAS(metadata.CAS)

        # In v1.3.3, we access the coefficients from the loaded dataframe
        if cas_rn_str not in Cp_data_Poling.index:
            raise ValueError(f"No Poling coefficients found for '{metadata.name}'.")

        # Get the row of coefficients for the chemical
        coeffs_row = Cp_data_Poling.loc[cas_rn_str]
        
        # Extract the 5 Poling coefficients (a, b, c, d, e)
        coeffs = coeffs_row[['a', 'b', 'c', 'd', 'e']].tolist()

        return jsonify({
            "chemical_name": metadata.name,
            "T_low": coeffs_row['Tmin'],
            "T_high": coeffs_row['Tmax'],
            "coeffs": coeffs
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
    
@ideal_gas_bp.route('/api/ideal_gas/joback_cp', methods=['POST'])
def joback_cp_api():
    """
    Calculates Ideal Gas Heat Capacity coefficients (A, B, C, D) using the
    Joback group contribution method and provides a "Glass Box" explanation.
    """
    if not RDKIT_AVAILABLE:
        return jsonify({"error": "RDKit library is required for this calculation."}), 500

    data = request.get_json() or {}
    chemical_name = data.get('chemical_name')
    if not chemical_name:
        return jsonify({"error": "Chemical name is required."}), 400

    try:
        # --- Step 1: Identify Joback Groups ---
        metadata = search_chemical(chemical_name)
        if not metadata or not metadata.smiles:
            raise ValueError(f"Could not find a SMILES string for '{chemical_name}'.")

        mol = Chem.MolFromSmiles(metadata.smiles)
        mol_with_hs = Chem.AddHs(mol)
        
        group_counts = {}
        for name, group_data in JOBACK_GROUPS_DATA.items():
            if "smarts" in group_data and group_data["smarts"]:
                pattern = Chem.MolFromSmarts(group_data["smarts"])
                if pattern and mol_with_hs.HasSubstructMatch(pattern):
                    matches = mol_with_hs.GetSubstructMatches(pattern)
                    if matches:
                        group_counts[name] = len(matches)

        if not group_counts:
            raise ValueError("No Joback groups could be identified for this molecule.")

        # --- Step 2: Sum Group Contributions ---
        sum_a, sum_b, sum_c, sum_d = 0.0, 0.0, 0.0, 0.0
        contributions_table_rows = ""
        for name, count in sorted(group_counts.items()):
            g = JOBACK_GROUPS_DATA[name]
            cpa, cpb, cpc, cpd = g.get('CpAk', 0), g.get('CpBk', 0), g.get('CpCk', 0), g.get('CpDk', 0)
            sum_a += count * cpa
            sum_b += count * cpb
            sum_c += count * cpc
            sum_d += count * cpd
            contributions_table_rows += f"<tr><td>{name}</td><td>{count}</td><td>{cpa:.4f}</td><td>{cpb:.4f}</td><td>{cpc:g}</td><td>{cpd:g}</td></tr>"

        # --- Step 3: Apply Adjustment Constants to get Final A, B, C, D ---
        coeff_A = sum_a - 37.93
        coeff_B = sum_b + 0.210
        coeff_C = sum_c - 3.91e-4
        coeff_D = sum_d + 2.06e-7

        # --- Step 4: Generate the "Glass Box" HTML Explanation ---
        plain_steps = "<h6>Step 1: Group Contribution Summation (Σ)</h6>"
        plain_steps += f"""
            <table class="table table-sm table-bordered">
                <thead class="table-light"><tr><th>Group</th><th>Count (N)</th><th>Δ<sub>CpA</sub></th><th>Δ<sub>CpB</sub></th><th>Δ<sub>CpC</sub></th><th>Δ<sub>CpD</sub></th></tr></thead>
                <tbody>{contributions_table_rows}</tbody>
                <tfoot><tr><td><b>Σ =</b></td><td></td><td><b>{sum_a:.4f}</b></td><td><b>{sum_b:.4f}</b></td><td><b>{sum_c:g}</b></td><td><b>{sum_d:g}</b></td></tr></tfoot>
            </table>
        """
        plain_steps += "<h6>Step 2: Coefficient Adjustment</h6>"
        plain_steps += f"""
            <ul class="list-unstyled">
                <li><b>A</b> = (ΣΔ<sub>CpA</sub>) - 37.93 = {sum_a:.4f} - 37.93 = <strong>{coeff_A:.4f}</strong></li>
                <li><b>B</b> = (ΣΔ<sub>CpB</sub>) + 0.210 = {sum_b:.4f} + 0.210 = <strong>{coeff_B:.4f}</strong></li>
                <li><b>C</b> = (ΣΔ<sub>CpC</sub>) - 3.91e-4 = {sum_c:g} - 3.91e-4 = <strong>{coeff_C:g}</strong></li>
                <li><b>D</b> = (ΣΔ<sub>CpD</sub>) + 2.06e-7 = {sum_d:g} + 2.06e-7 = <strong>{coeff_D:g}</strong></li>
            </ul>
        """
        plain_steps += "<h6>Step 3: Final Polynomial</h6>"
        plain_steps += f"<p class='fs-6'><b>C<sub>p</sub> (J/mol·K) = {coeff_A:.2f} + {coeff_B:.3f}T + {coeff_C:.2e}T<sup>2</sup> + {coeff_D:.2e}T<sup>3</sup></b></p>"

        # --- Step 5: Return the Final Data ---
        return jsonify({
            'A': coeff_A,
            'B': coeff_B,
            'C': coeff_C,
            'D': coeff_D,
            'plain_steps': plain_steps
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500