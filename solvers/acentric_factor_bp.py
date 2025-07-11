# D:\ChemE_Calc_New\solvers\acentric_factor_bp.py
import math
import traceback
import datetime
from flask import Blueprint, render_template, request, jsonify, url_for

# --- Corrected Imports ---
from chemicals.acentric import LK_omega
from chemicals.identifiers import CAS_from_any, check_CAS, search_chemical
from chemicals.critical import Tc_methods, Pc_methods, Tc, Pc
from chemicals.phase_change import Tb_methods, Tb

# --- Import our own work from the other blueprint! ---
from .critical_props_bp import identify_cg_groups, CG_FIRST_ORDER_GROUPS_DATA, CG_SECOND_ORDER_GROUPS_DATA, RDKIT_AVAILABLE, Chem

# --- Blueprint Definition ---
acentric_factor_bp = Blueprint(
    'acentric_factor',
    __name__,
    template_folder='../templates',
    static_folder='../static'
)

# --- API to find chemical properties and methods ---
@acentric_factor_bp.route('/api/find_chemical_methods', methods=['POST'])
def find_chemical_methods_api():
    data = request.get_json() or {}
    chemical_name = data.get('name')
    if not chemical_name:
        return jsonify({"error": "Chemical name is required."}), 400
    try:
        cas_rn = CAS_from_any(chemical_name)
        if not cas_rn or not check_CAS(cas_rn):
             return jsonify({"error": f"Could not find a valid CAS number for '{chemical_name}'."}), 404
        
        metadata = search_chemical(chemical_name)
        smiles = metadata.smiles if metadata else None

        available_methods = {
            'tc': Tc_methods(cas_rn),
            'pc': Pc_methods(cas_rn),
            'tb': Tb_methods(cas_rn)
        }
        return jsonify({
            "cas": cas_rn,
            "chemical_name": chemical_name,
            "smiles": smiles,
            "available_methods": available_methods
        })
    except Exception as e:
        return jsonify({"error": f"Could not identify '{chemical_name}'. Please try a different name or synonym."}), 404

# --- Helper for citations ---
def get_method_citation(method_key):
    citations = {
        'CRC_ORG': {"id": "crc_2014", "type": "Book", "source_type": "Data Source", "data": {"author": "Haynes, W.M., et al.", "year": 2014, "title": "CRC Handbook of Chemistry and Physics, 95E", "publisher": "CRC Press"}},
        'YAWS': {"id": "yaws_2014", "type": "Book", "source_type": "Data Source", "data": {"author": "Yaws, C. L.", "year": 2014, "title": "Thermophysical Properties of Chemicals and Hydrocarbons, Second Edition", "publisher": "Gulf Professional Publishing"}},
        'JOBACK': {"id": "joback_1987", "type": "Journal Article", "source_type": "Estimation Method", "data": {"authors": ["Joback, K.G.", "Reid, R.C."], "year": 1987, "title": "Estimation of Pure-Component Properties from Group-Contributions", "journal": "Chemical Engineering Communications", "volume": "57", "pages": "233-243"}},
        'IUPAC': {"id": "iupac_crit_review", "type": "Journal Series", "source_type": "Critically Evaluated Data", "data": {"author": "Ambrose, D., Tsonopoulos, C., et al.", "year": "1995-2015", "title": "Vapor-Liquid Critical Properties of Elements and Compounds (Series)", "journal": "Journal of Chemical & Engineering Data"}},
        'WILSON_JASPERSON': {"id": "wilson_1996", "type": "Conference Proceeding", "source_type": "Estimation Method", "data": {"authors": ["Wilson, G. M.", "Jasperson, L. V."], "year": 1996, "title": "Critical Constants Tc, Pc, Estimation Based on Zero, First and Second Order Methods", "booktitle": "Proceedings of the AIChE Spring Meeting"}}
    }
    default_citation = {"id": "unknown", "type": "Data", "source_type": "Unspecified Source", "data": {"title": f"Data retrieved using the '{method_key}' method from the 'chemicals' library."}}
    return citations.get(method_key, default_citation)

# --- API to get a specific property value ---
@acentric_factor_bp.route('/api/get_property_value', methods=['POST'])
def get_property_value_api():
    data = request.get_json() or {}
    cas_rn, prop_key, method = data.get('cas'), data.get('property'), data.get('method')
    if not all([cas_rn, prop_key, method]):
        return jsonify({"error": "CAS number, property, and method are required."}), 400
    try:
        value = None
        if prop_key == 'tc': value = Tc(cas_rn, method=method)
        elif prop_key == 'pc': value = Pc(cas_rn, method=method)
        elif prop_key == 'tb': value = Tb(cas_rn, method=method)
        else: return jsonify({"error": f"Unknown property key: '{prop_key}'"}), 400
        if value is None: return jsonify({"error": f"No value found for '{prop_key}' with method '{method}'."}), 404
        return jsonify({"property": prop_key, "value": value, "citation": get_method_citation(method)})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"An error occurred: {str(e)}"}), 500

# --- Page-Rendering Routes ---
@acentric_factor_bp.route('/calculators/acentric-factor')
def acentric_factor_page():
    return render_template('acentric_factor_solvers.html')

@acentric_factor_bp.route('/solvers/acentric-factor-cg')
def cg_acentric_factor_solver_page():
    return render_template('solvers/cg_acentric_factor_calculator.html')

# --- Lee-Kesler API Endpoint ---
@acentric_factor_bp.route('/api/calculate_lk_omega', methods=['POST'])
def calculate_lk_omega_api():
    data = request.get_json() or {}
    try:
        tb, tc, pc_bar = float(data.get('Tb')), float(data.get('Tc')), float(data.get('Pc_bar'))
    except (TypeError, ValueError):
        return jsonify({"error": "Valid numerical inputs for Tb, Tc, and Pc are required."}), 400
    try:
        omega, tbr = LK_omega(Tb=tb, Tc=tc, Pc=pc_bar*1e5), tb / tc
        tau = 1.0 - tbr
        f0 = -5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 - 1.06841*tau**5
        f1 = -5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 - 7.46628*tau**5
        plain_steps = f"<p><b>Inputs:</b><br>T<sub>b</sub> = {tb:.2f} K<br>T<sub>c</sub> = {tc:.2f} K<br>P<sub>c</sub> = {pc_bar:.2f} bar</p>"
        plain_steps += f"<p><b>1. Calc. Reduced Boiling Temp (T<sub>br</sub>):</b><br>T<sub>br</sub> = {tb:.2f} / {tc:.2f} = {tbr:.4f}</p>"
        plain_steps += f"<p><b>2. Calc. Lee-Kesler Functions:</b><br>f<sup>(0)</sup>(T<sub>br</sub>) &approx; {f0:.4f}<br>f<sup>(1)</sup>(T<sub>br</sub>) &approx; {f1:.4f}</p>"
        plain_steps += f"<p><b>3. Calc. Acentric Factor (ω):</b><br>ω &approx; [{f0:.4f} - ln({pc_bar:.2f}/1.01325)] / {f1:.4f} &approx; {omega:.4f}</p>"
        latex_steps = rf"$$T_{{br}} = \frac{{{tb:.2f}}}{{{tc:.2f}}} = {tbr:.4f}$$"
        latex_steps += rf"$$f^{{(0)}}(T_{{br}}) \approx {f0:.4f} \quad f^{{(1)}}(T_{{br}}) \approx {f1:.4f}$$"
        latex_steps += rf"$$\omega \approx \frac{{f^{{(0)}} - \ln(P_c/1.01325)}}{{f^{{(1)}}}} \approx \frac{{{f0:.4f} - \ln({pc_bar:.2f}/1.01325)}}{{{f1:.4f}}} \approx \mathbf{{{omega:.4f}}}$$"
        today, page_url = datetime.date.today(), request.url_root.rstrip('/') + url_for('acentric_factor.acentric_factor_page')
        citations = [
            {"id": "tori_calc", "type": "Software", "source_type": "Tool", "data": {"author": "ToriCalc", "year": today.year, "title": "Lee-Kesler Acentric Factor Calculator", "retrieved_date": today.strftime('%Y-%m-%d'), "url": page_url}},
            {"id": "lee_kesler_1975", "type": "Journal Article", "source_type": "Source Method", "data": {"authors": ["Lee, B. I.", "Kesler, M. G."], "year": 1975, "title": "A Generalized Thermodynamic Correlation Based on Three-Parameter Corresponding States", "journal": "AIChE Journal", "volume": "21", "issue": "3", "pages": "510-527"}},
            {"id": "chemicals_lib", "type": "Software", "source_type": "Calculation Engine", "data": {"author": "Caleb Bell, Yoel Rene Cortes-Pena, and Contributors", "year": "2016-2024", "title": "Chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)", "url": "https://github.com/CalebBell/chemicals"}}
        ]
        return jsonify({'omega': f"{omega:.4f}", 'latex_steps': latex_steps, 'plain_steps': plain_steps, 'citations': citations})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"An error occurred during calculation: {str(e)}"}), 500

# --- Constantinou-Gani API Endpoint ---
@acentric_factor_bp.route('/api/calculate_cg_omega', methods=['POST'])
def calculate_cg_omega_api():
    """
    Calculates acentric factor using CG method. NOW with highly detailed,
    table-based explanations for both HTML and LaTeX.
    """
    # --- Citation and Notes Data ---
    today = datetime.date.today()
    page_url = request.url_root.rstrip('/') + url_for('acentric_factor.cg_acentric_factor_solver_page')
    tool_citation = {"id": "toricalc_cg", "type": "Software", "source_type": "Tool", "data": { "author": "ToriCalc", "year": today.year, "title": "Constantinou-Gani Acentric Factor Calculator", "retrieved_date": today.strftime('%Y-%m-%d'), "url": page_url }}
    cg_paper_citation = {"id": "cg_1995", "type": "Journal Article", "source_type": "Source Method", "data": { "authors": ["Constantinou, L.", "Gani, R.", "O'Connell, J. P."], "year": 1995, "title": "Estimation of the acentric factor and the liquid molar volume at 298 K using a new group contribution method", "journal": "Fluid Phase Equilibria", "volume": "103", "issue": "1", "pages": "11-22"}}
    library_citation = {"id": "chemicals_lib", "type": "Software", "source_type": "Calculation Engine", "data": { "author": "Caleb Bell, Yoel Rene Cortes-Pena, and Contributors", "year": "2016-2024", "title": "Chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)", "url": "https://github.com/CalebBell/chemicals"}}
    applicability_notes = "Highly versatile group contribution method for diverse organic compounds. Accuracy is improved with second-order groups for complex structures. Always verify against experimental data where possible."

    # --- Request Handling & Calculation ---
    data = request.get_json() or {}
    chemical_name = data.get('chemical_name')
    calculation_order_str = data.get('calculation_order', '1')

    if not chemical_name: return jsonify({"error": "Chemical name is required."}), 400

    try:
        metadata = search_chemical(chemical_name)
        if not metadata or not metadata.smiles: raise ValueError(f"Could not identify '{chemical_name}'.")
        smiles_str = metadata.smiles
        mol = Chem.MolFromSmiles(smiles_str.strip())
        if mol is None: raise ValueError(f"Invalid SMILES: '{smiles_str}'.")

        first_order, second_order = identify_cg_groups(mol, CG_FIRST_ORDER_GROUPS_DATA, CG_SECOND_ORDER_GROUPS_DATA, calculation_order_str)
        if not first_order: raise ValueError("No first-order Constantinou-Gani groups were identified for this molecule.")

        sum_w1k = sum(CG_FIRST_ORDER_GROUPS_DATA[k]["w1k"] * v for k, v in first_order.items() if CG_FIRST_ORDER_GROUPS_DATA[k].get("w1k") is not None)
        sum_w2j = sum(CG_SECOND_ORDER_GROUPS_DATA[k]["w2j"] * v for k, v in second_order.items() if CG_SECOND_ORDER_GROUPS_DATA[k].get("w2j") is not None)

        omega_ln_arg = sum_w1k + sum_w2j + 1.1507
        if omega_ln_arg <= 0: raise ValueError("Logarithm argument for omega calculation is not positive.")
        omega = (1.0 / 0.5050) * (0.4085 * math.log(omega_ln_arg))

        # --- NEW: DETAILED "GLASS BOX" EXPLANATION GENERATION ---

        # --- Plain HTML Steps ---
        plain_steps = f"<h5>1. Identify First-Order Groups (N<sub>k</sub>)</h5>"
        plain_steps += "<p>Based on the structure, the following primary groups are identified:</p>"
        plain_steps += "<table class='table table-bordered table-sm'><thead><tr><th>Group Name</th><th>Count (N<sub>k</sub>)</th><th>ω<sub>1k</sub> Contribution</th><th>Total (N<sub>k</sub> * ω<sub>1k</sub>)</th></tr></thead><tbody>"
        for name, count in first_order.items():
            if count == 0: continue
            w1k_val = CG_FIRST_ORDER_GROUPS_DATA[name].get('w1k', 0)
            plain_steps += f"<tr><td>{name}</td><td>{count}</td><td>{w1k_val:.4f}</td><td>{count * w1k_val:.4f}</td></tr>"
        plain_steps += f"</tbody><tfoot><tr><td colspan='3' class='text-end'><b>Sum of First-Order Contributions (Σ N<sub>k</sub>ω<sub>1k</sub>)</b></td><td><b>{sum_w1k:.4f}</b></td></tr></tfoot></table>"

        if second_order and any(v > 0 for v in second_order.values()):
            plain_steps += f"<h5>2. Identify Second-Order Groups (M<sub>j</sub>)</h5>"
            plain_steps += "<table class='table table-bordered table-sm'><thead><tr><th>Group Name</th><th>Count (M<sub>j</sub>)</th><th>ω<sub>2j</sub> Contribution</th><th>Total (M<sub>j</sub> * ω<sub>2j</sub>)</th></tr></thead><tbody>"
            for name, count in second_order.items():
                if count == 0: continue
                w2j_val = CG_SECOND_ORDER_GROUPS_DATA[name].get('w2j', 0)
                plain_steps += f"<tr><td>{name}</td><td>{count}</td><td>{w2j_val:.5f}</td><td>{count * w2j_val:.5f}</td></tr>"
            plain_steps += f"</tbody><tfoot><tr><td colspan='3' class='text-end'><b>Sum of Second-Order Contributions (Σ M<sub>j</sub>ω<sub>2j</sub>)</b></td><td><b>{sum_w2j:.5f}</b></td></tr></tfoot></table>"

        plain_steps += f"<h5>3. Final Calculation</h5>"
        plain_steps += "<p>The acentric factor (ω) is calculated using the Constantinou & Gani formula:</p>"
        plain_steps += "<p>ω = (1/0.5050) * [0.4085 * ln(ΣN<sub>k</sub>ω<sub>1k</sub> + ΣM<sub>j</sub>ω<sub>2j</sub> + 1.1507)]</p>"
        plain_steps += f"<p>ω = 1.9802 * [0.4085 * ln({sum_w1k:.4f} + {sum_w2j:.5f} + 1.1507)]</p>"
        plain_steps += f"<p>ω = 1.9802 * [0.4085 * ln({omega_ln_arg:.4f})]</p>"
        plain_steps += f"<h6><b>ω ≈ {omega:.4f}</b></h6>"

        # --- LaTeX Steps ---
        latex_steps = r"$$\text{1. First-Order Group Contributions}$$"
        latex_steps += r"$$\begin{array}{|l|c|c|r|} \hline \textbf{Group} & \mathbf{N_k} & \mathbf{\omega_{1k}} & \mathbf{N_k \cdot \omega_{1k}} \\ \hline "
        for name, count in first_order.items():
            if count == 0: continue
            w1k_val = CG_FIRST_ORDER_GROUPS_DATA[name].get('w1k', 0)
            latex_steps += rf"\text{{{name.replace('-', ' - ')}}} & {count} & {w1k_val:.4f} & {count * w1k_val:.4f} \\\\ "
        latex_steps += rf"\hline \multicolumn{{3}}{{r|}}{{\textbf{{Sum}}}} & \mathbf{{{sum_w1k:.4f}}} \\\\ \hline \end{{array}}$$ "
        
        if second_order and any(v > 0 for v in second_order.values()):
            latex_steps += r"$$\text{2. Second-Order Group Contributions}$$"
            latex_steps += r"$$\begin{array}{|l|c|c|r|} \hline \textbf{Group} & \mathbf{M_j} & \mathbf{\omega_{2j}} & \mathbf{M_j \cdot \omega_{2j}} \\ \hline "
            for name, count in second_order.items():
                if count == 0: continue
                w2j_val = CG_SECOND_ORDER_GROUPS_DATA[name].get('w2j', 0)
                latex_steps += rf"\text{{{name.replace('-', ' - ')}}} & {count} & {w2j_val:.5f} & {count * w2j_val:.5f} \\\\ "
            latex_steps += rf"\hline \multicolumn{{3}}{{r|}}{{\textbf{{Sum}}}} & \mathbf{{{sum_w2j:.5f}}} \\\\ \hline \end{{array}}$$ "
        
        latex_steps += r"$$\text{3. Final Calculation}$$"
        latex_steps += r"$$\omega = \frac{1}{0.5050} \left[ 0.4085 \cdot \ln\left(\sum w_{1k} + \sum w_{2j} + 1.1507\right) \right]$$"
        latex_steps += rf"$$\omega \approx \frac{1}{0.5050} \left[ 0.4085 \cdot \ln\left({sum_w1k:.4f} + {sum_w2j:.5f} + 1.1507\right) \right] = \mathbf{{{omega:.4f}}}$$"
        
        return jsonify({
            'omega': f"{omega:.4f}",
            'latex_steps': latex_steps,
            'plain_steps': plain_steps,
            'found_smiles': smiles_str,
            'first_order_groups': first_order,
            'second_order_groups': second_order,
            'citations': [tool_citation, cg_paper_citation, library_citation],
            'applicability_notes': applicability_notes
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"An error occurred: {str(e)}"}), 500 