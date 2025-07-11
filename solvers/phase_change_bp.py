import math
import traceback
import datetime
from flask import Blueprint, render_template, request, jsonify, url_for

# --- Project-Specific Imports ---
from .critical_props_bp import JOBACK_GROUPS_DATA, RDKIT_AVAILABLE, Chem

# --- `chemicals` Library Imports ---
from chemicals.identifiers import CAS_from_any, check_CAS, search_chemical
from chemicals.critical import Tc, Pc
from chemicals.acentric import omega
from chemicals.phase_change import (
    Tb, Tb_methods, 
    Tm, Tm_methods, 
    Hfus, Hfus_methods,
    Riedel, Chen, Vetere, Pitzer, Watson
)

# --- Blueprint Definition ---
phase_change_bp = Blueprint('phase_change', __name__)

# --- Page Rendering Route ---
@phase_change_bp.route('/calculators/phase-change-properties')
def phase_change_page():
    return render_template('calculator_phase_change_properties.html')

# --- API Endpoints ---

@phase_change_bp.route('/api/phase_change/find_methods', methods=['POST'])
def find_chemical_methods_api():
    # ... (This function is already correct from our last step, no changes needed) ...
    data = request.get_json() or {}
    chemical_name = data.get('name')
    if not chemical_name: return jsonify({"error": "Chemical name is required."}), 400
    try:
        cas_rn = CAS_from_any(chemical_name)
        if not cas_rn or not check_CAS(cas_rn):
             return jsonify({"error": f"Could not find a valid CAS number for '{chemical_name}'."}), 400
        
        metadata = search_chemical(chemical_name)
        canonical_name = metadata.common_name if metadata and metadata.common_name else chemical_name
        
        try: tc_val = Tc(cas_rn)
        except: tc_val = None
        try: pc_val = Pc(cas_rn)
        except: pc_val = None
        try: omega_val = omega(cas_rn)
        except: omega_val = None
        try: tb_val = Tb(cas_rn)
        except: tb_val = None
        
        hvap_correlations_data = [
            {"name": "Watson", "required_inputs": ["T", "Hvap_ref", "T_ref", "Tc"]},
            {"name": "Riedel", "required_inputs": ["Tb", "Tc", "Pc_pa"]},
            {"name": "Chen", "required_inputs": ["Tb", "Tc", "Pc_pa"]},
            {"name": "Pitzer", "required_inputs": ["T", "Tc", "omega"]},
        ]
        
        return jsonify({
            "cas": cas_rn,
            "chemical_name": canonical_name,
            "prerequisite_data": { "tc": tc_val, "pc_pa": pc_val, "omega": omega_val, "tb": tb_val },
            "available_methods": { 'tb': Tb_methods(cas_rn), 'tm': Tm_methods(cas_rn), 'hfus': Hfus_methods(cas_rn) },
            "hvap_correlations": hvap_correlations_data
        })
    except Exception as e:
        return jsonify({"error": f"Could not identify '{chemical_name}'. Please try a different name or synonym."}), 400

@phase_change_bp.route('/api/phase_change/get_value', methods=['POST'])
def get_property_value_api():
    # ... (This function is already correct, no changes needed) ...
    data = request.get_json() or {}
    cas_rn, prop_key, method = data.get('cas'), data.get('property'), data.get('method')
    if not all([cas_rn, prop_key, method]): return jsonify({"error": "CAS, property, and method are required."}), 400
    try:
        value = None
        if prop_key == 'tb': value = Tb(cas_rn, method=method)
        elif prop_key == 'tm': value = Tm(cas_rn, method=method)
        elif prop_key == 'hfus': value = Hfus(cas_rn, method=method)
        else: return jsonify({"error": f"Unknown property key: '{prop_key}'"}), 400
        if value is None: return jsonify({"error": f"No value found for '{prop_key}' with method '{method}'."}), 404
        return jsonify({ "value": value })
    except Exception as e: return jsonify({"error": f"An error occurred: {str(e)}"}), 500

@phase_change_bp.route('/api/phase_change/joback_calculate', methods=['POST'])
def joback_calculate_api():
    # ... (This function is already correct, but we'll add a citation to its response) ...
    data = request.get_json() or {}
    # ... (rest of the logic is the same)
    if not RDKIT_AVAILABLE: return jsonify({"error": "RDKit library is required."}), 500
    chemical_name, prop_key = data.get('chemical_name'), data.get('property')
    if not chemical_name or not prop_key: return jsonify({"error": "Chemical name and property key are required."}), 400
    try:
        mol = Chem.MolFromSmiles(search_chemical(chemical_name).smiles)
        mol_with_hs = Chem.AddHs(mol)
        group_counts = {name: len(mol_with_hs.GetSubstructMatches(Chem.MolFromSmarts(data["smarts"]))) for name, data in JOBACK_GROUPS_DATA.items() if "smarts" in data and data["smarts"] and mol_with_hs.HasSubstructMatch(Chem.MolFromSmarts(data["smarts"]))}
        if not group_counts: raise ValueError("No Joback groups identified.")
        prop_map = {'tb': 'tbk', 'tm': 'tfpk', 'hfus': 'hfk'}
        joback_prop_key = prop_map.get(prop_key)
        sum_contrib = sum(JOBACK_GROUPS_DATA[k][joback_prop_key] * v for k, v in group_counts.items() if JOBACK_GROUPS_DATA[k].get(joback_prop_key) is not None)
        formula_map = {'tb': (198.2, "T_b (K) = 198.2 + \\Sigma (N_k T_{bk})"), 'tm': (122.5, "T_m (K) = 122.5 + \\Sigma (N_k T_{fk})"), 'hfus': (-0.88, "H_{fus} (kJ/mol) = -0.88 + \\Sigma (N_k H_{fk})")}
        constant, formula = formula_map[prop_key]
        result_val = constant + sum_contrib
        html_steps = "<h6>1. Identified Joback Groups (N<sub>k</sub>)</h6><table class='table table-sm table-bordered'><thead><tr><th>Group</th><th>Count</th><th>Contribution</th><th>Total</th></tr></thead><tbody>"
        for name, count in sorted(group_counts.items()):
            contrib_val = JOBACK_GROUPS_DATA[name].get(joback_prop_key, 0)
            html_steps += f"<tr><td>{name}</td><td>{count}</td><td>{contrib_val:.2f}</td><td>{count * contrib_val:.2f}</td></tr>"
        html_steps += f"</tbody><tfoot><tr><td colspan='3' class='text-end'><b>Sum</b></td><td><b>{sum_contrib:.2f}</b></td></tr></tfoot></table><h6>2. Final Calculation</h6><p>{formula.replace('_','<sub>').replace('{','').replace('}','')} = {constant} + {sum_contrib:.2f} = <b>{result_val:.2f}</b></p>"
        if prop_key == 'hfus': result_val *= 1000
        citation = {"id": "joback_1987", "type": "Journal Article", "data": {"authors": ["Joback, K.G.", "Reid, R.C."], "year": 1987, "title": "Estimation of Pure-Component Properties from Group-Contributions", "journal": "Chemical Engineering Communications", "volume": "57", "issue": "1-6", "pages": "233-243"}}
        return jsonify({"value": result_val, "glass_box_html": html_steps, "citation": citation})
    except Exception as e: return jsonify({"error": str(e)}), 500



# --- UPGRADED: Hvap "Glass Box" Calculator ---
@phase_change_bp.route('/api/phase_change/calculate_hvap', methods=['POST'])
def calculate_hvap_api():
    data = request.get_json() or {}
    correlation, inputs = data.get('correlation'), data.get('inputs', {})
    try:
        float_inputs = {k: float(v) for k, v in inputs.items() if v is not None and v != ''}
        
        correlation_map = {
            'Riedel': (Riedel, "H_{vap} = \\frac{1.093RT_b(\\ln(P_c/101325) - 1.013)}{0.930 - T_{br}}"),
            'Chen': (Chen, "H_{vap} = \\frac{RT_b(3.978T_{br}-3.958+1.555\\ln P_c)}{1.07-T_{br}}"),
            'Pitzer': (Pitzer, "H_{vap} = RT_c(7.08(1-T_r)^{0.354} + 10.95\\omega(1-T_r)^{0.456})"),
            'Watson': (Watson, "H_{vap,T} = H_{vap,ref} \\left( \\frac{1 - T_{r,target}}{1 - T_{r,ref}} \\right)^{0.38}")
        }
        if correlation not in correlation_map: return jsonify({"error": "Invalid correlation."}), 400

        func, formula_str = correlation_map[correlation]
        arg_map = {'Pc_pa': 'Pc', 'Tb': 'Tb', 'Tc': 'Tc', 'T': 'T', 'omega': 'omega', 'Hvap_ref': 'Hvap_ref', 'T_ref': 'T_ref'}
        kwargs = {arg_map[k]: v for k, v in float_inputs.items() if k in arg_map}
        result = func(**kwargs)

        # --- Build Detailed, LaTeX-Consistent Explanation ---
        plain_steps = ""
        latex_steps = ""
        if correlation == 'Watson':
            t, hvap_ref, t_ref, tc = float_inputs['T'], float_inputs['Hvap_ref'], float_inputs['T_ref'], float_inputs['Tc']
            tr_ref, tr_target = t_ref / tc, t / tc
            plain_steps = f"<p><b>1. Assemble Inputs:</b></p><ul><li>$T_{{ref}}$: {t_ref:.2f} K</li><li>$H_{{vap,ref}}$: {hvap_ref:.2f} J/mol</li><li>Target $T$: {t:.2f} K</li><li>$T_c$: {tc:.2f} K</li></ul>"
            plain_steps += f"<p><b>2. Calculate Reduced Temperatures ($T_r$):</b></p><ul><li>$T_{{r,ref}} = T_{{ref}} / T_c = {t_ref:.2f} / {tc:.2f} = {tr_ref:.4f}$</li><li>$T_{{r,target}} = T / T_c = {t:.2f} / {tc:.2f} = {tr_target:.4f}$</li></ul>"
            plain_steps += f"<p><b>3. Apply Watson Correlation:</b></p><p>$${formula_str}$$</p>"
            plain_steps += f"<p>$$ = {hvap_ref:.2f} \\left( \\frac{{1 - {tr_target:.4f}}}{{1 - {tr_ref:.4f}}} \\right)^{{0.38}} \\approx \\mathbf{{{result:.2f}}}$$ J/mol</p>"
            latex_steps = rf"$$T_{{r,ref}} = {tr_ref:.4f} \quad T_{{r,target}} = {tr_target:.4f} \\ \Delta H_{{vap,T}} = {hvap_ref:.2f} \left( \frac{{1 - {tr_target:.4f}}}{{1 - {tr_ref:.4f}}} \right)^{{0.38}} \approx \mathbf{{{result:.2f}}}$$"

        elif correlation == 'Riedel':
            tb, tc, pc = float_inputs['Tb'], float_inputs['Tc'], float_inputs['Pc_pa']
            tbr = tb / tc
            plain_steps = f"<p><b>1. Assemble Inputs:</b></p><ul><li>$T_b$: {tb:.2f} K</li><li>$T_c$: {tc:.2f} K</li><li>$P_c$: {pc:.2f} Pa</li></ul>"
            plain_steps += f"<p><b>2. Calculate Reduced Boiling Temp ($T_{{br}}$):</b> $T_{{br}} = T_b / T_c = {tb:.2f} / {tc:.2f} = {tbr:.4f}$</p>"
            plain_steps += f"<p><b>3. Apply Riedel Correlation:</b></p><p>$${formula_str}$$</p><p>...evaluates to $\\approx \\mathbf{{{result:.2f}}}$ J/mol</p>"
            latex_steps = rf"$$T_{{br}} = {tbr:.4f} \\ \Delta H_{{vap}} \approx \mathbf{{{result:.2f}}}$$"

        # (Other correlations like Chen/Pitzer would follow the same pattern)
        else:
            plain_steps = f"Calculation for {correlation} completed. Result: <b>{result:.2f} J/mol</b>"
            latex_steps = f"Result \\approx {result:.2f}"

        citation = {"id": correlation.lower(), "type": "Correlation", "data": {"authors": [correlation], "year": "N/A", "title": f"{correlation} Correlation", "journal": "The Properties of Gases and Liquids", "volume": "N/A", "issue": "N/A", "pages": "N/A"}}
        return jsonify({"value": result, "plain_steps": plain_steps, "latex_steps": latex_steps, "citation": citation})
    except Exception as e:
        return jsonify({"error": f"Calculation failed: {str(e)}"}), 500