import traceback
import datetime
from flask import Blueprint, render_template, request, jsonify

# Use imports that are compatible with chemicals v1.3.3
from chemicals.identifiers import search_chemical, int_to_CAS
from chemicals.heat_capacity import Cp_dict_PerryI
from solvers.critical_props_bp import JOBACK_GROUPS_DATA, RDKIT_AVAILABLE, Chem

# --- Blueprint Definition ---
ideal_gas_bp = Blueprint('ideal_gas', __name__, template_folder='../templates')

# --- Page Rendering Route ---
@ideal_gas_bp.route('/calculators/ideal-gas-properties')
def ideal_gas_page():
    return render_template('calculator_ideal_gas.html')

# --- The Single, Unified API Endpoint ---
@ideal_gas_bp.route('/api/ideal_gas/get_comparison_data', methods=['POST'])
def get_comparison_data_api():
    data = request.get_json() or {}
    chemical_name = data.get('chemical_name')
    if not chemical_name:
        return jsonify({"error": "Chemical identifier is required."}), 400

    # --- Helper function to format numbers into LaTeX scientific notation ---
    def _to_latex_sci(num, precision=3):
        if num == 0:
            return "0"
        # Use scientific notation for very small or large numbers
        if abs(num) < 0.001 or abs(num) > 10000:
            return f"{{:{precision}e}}".format(num).replace('e', ' \\times 10^{{') + '}}'
        else:
            return f"{{:.{precision}f}}".format(num)

    results = {}
    citations = []
    try:
        metadata = search_chemical(chemical_name)
        cas_rn_str = int_to_CAS(metadata.CAS)
        results['chemical_name'] = metadata.common_name or chemical_name
        
        # --- Part 1: Perry's Coefficients ---
        try:
            if cas_rn_str in Cp_dict_PerryI and 'g' in Cp_dict_PerryI[cas_rn_str]:
                coeffs_dict = Cp_dict_PerryI[cas_rn_str]['g']
                coeffs = { 'a': coeffs_dict.get('Const', 0), 'b': coeffs_dict.get('Lin', 0), 'c': coeffs_dict.get('Quadinv', 0), 'd': coeffs_dict.get('Quad', 0) }
                # Build the correct LaTeX string
                latex_steps = f"C_p(T) [cal/mol·K] = {_to_latex_sci(coeffs['a'])} + {_to_latex_sci(coeffs['b'])}T + \\frac{{{_to_latex_sci(coeffs['c'])}}}{{T^2}} + {_to_latex_sci(coeffs['d'])}T^2"
                results['perry_results'] = { "coeffs": coeffs, "T_low": coeffs_dict.get('Tmin'), "T_high": coeffs_dict.get('Tmax'), "latex_steps": latex_steps }
                citations.append({"id": "perry_1997", "type": "Book", "data": {"authors": ["Perry, R.H.", "Green, D.W."], "year": 1997, "title": "Perry's Chemical Engineers' Handbook, 7th Edition"}})
            else:
                results['perry_results'] = {"error": "No Perry's Handbook gas coefficients found."}
        except Exception as e:
            results['perry_results'] = {"error": str(e)}

        # --- Part 2: Joback Calculation ---
        try:
            if not RDKIT_AVAILABLE: raise ValueError("RDKit library not available.")
            if not metadata.smiles: raise ValueError("Could not find SMILES for Joback method.")
            
            mol = Chem.MolFromSmiles(metadata.smiles)
            mol_with_hs = Chem.AddHs(mol)
            
            group_counts = {name: len(mol_with_hs.GetSubstructMatches(Chem.MolFromSmarts(g["smarts"]))) for name, g in JOBACK_GROUPS_DATA.items() if "smarts" in g and g["smarts"] and mol_with_hs.HasSubstructMatch(Chem.MolFromSmarts(g["smarts"]))}
            if not group_counts: raise ValueError("No Joback groups identified (method is for organic compounds).")
            
            sums = {k: sum(count * (JOBACK_GROUPS_DATA.get(name, {}).get(f'Cp{k.upper()}k', 0) or 0) for name, count in group_counts.items()) for k in ['a', 'b', 'c', 'd']}
            coeffs = {'A': sums['a'] - 37.93, 'B': sums['b'] + 0.210, 'C': sums['c'] - 3.91e-4, 'D': sums['d'] + 2.06e-7}
            
            contributions_rows = "".join([f"<tr><td>{name}</td><td>{count}</td><td>{JOBACK_GROUPS_DATA.get(name, {}).get('CpAk', 0):.4f}</td><td>{JOBACK_GROUPS_DATA.get(name, {}).get('CpBk', 0):.4f}</td><td>{JOBACK_GROUPS_DATA.get(name, {}).get('CpCk', 0):.2e}</td><td>{JOBACK_GROUPS_DATA.get(name, {}).get('CpDk', 0):.2e}</td></tr>" for name, count in sorted(group_counts.items())])
            plain_steps = f"<h6>Step 1: Group Summation (Σ)</h6><div class='table-responsive'><table class='table table-sm table-bordered'><thead><tr><th>Group</th><th>N</th><th>$\\Delta_{{A}}$</th><th>$\\Delta_{{B}}$</th><th>$\\Delta_{{C}}$</th><th>$\\Delta_{{D}}$</th></tr></thead><tbody>{contributions_rows}</tbody><tfoot><tr><td colspan='2'><b>Σ</b></td><td><b>{sums['a']:.4f}</b></td><td><b>{sums['b']:.4f}</b></td><td><b>{sums['c']:.2e}</b></td><td><b>{sums['d']:.2e}</b></td></tr></tfoot></table></div>"
            plain_steps += f"<h6>Step 2: Coefficient Adjustment</h6><ul class='list-unstyled small'><li>$A = (\\Sigma\\Delta_{{A}}) - 37.93 = <strong>{coeffs['A']:.4f}</strong></li><li>$B = (\\Sigma\\Delta_{{B}}) + 0.210 = <strong>{coeffs['B']:.4f}</strong></li><li>$C = (\\Sigma\\Delta_{{C}}) - 3.91 \\times 10^{{-4}} = <strong>{coeffs['C']:.2e}</strong></li><li>$D = (\\Sigma\\Delta_{{D}}) + 2.06 \\times 10^{{-7}} = <strong>{coeffs['D']:.2e}</strong></li></ul>"
            plain_steps += f"<h6>Step 3: Final Polynomial</h6><p class='text-center'><b>$C_p(T) \\approx {coeffs['A']:.2f} + {coeffs['B']:.3f}T {coeffs['C']:+.2e}T^2 {coeffs['D']:+.2e}T^3$</b></p>"
            
            # Build the full, correct LaTeX string for the copy button
            latex_steps = (f"C_p(T) [J/mol·K] \\approx {_to_latex_sci(coeffs['A'])} "
                           f"{'+' if coeffs['B'] > 0 else ''}{_to_latex_sci(coeffs['B'])}T "
                           f"{'+' if coeffs['C'] > 0 else ''}{_to_latex_sci(coeffs['C'])}T^2 "
                           f"{'+' if coeffs['D'] > 0 else ''}{_to_latex_sci(coeffs['D'])}T^3")

            results['joback_results'] = {'coeffs': coeffs, 'plain_steps': plain_steps, 'latex_steps': latex_steps}
            citations.append({"id": "joback_1987", "type": "Journal Article", "data": {"authors": ["Joback, K.G.", "Reid, R.C."], "year": 1987, "title": "Estimation of Pure-Component Properties from Group-Contributions", "journal": "Chemical Engineering Communications"}})
        except Exception as e:
            results['joback_results'] = {"error": str(e)}

        results['citations'] = citations
        return jsonify(results)
    except Exception as e:
        return jsonify({"error": f"Could not identify '{chemical_name}'. Please try a different name."}), 400