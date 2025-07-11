# D:\ChemE_Calc_New\solvers\diffusion_bp.py
from flask import Blueprint, render_template, request, jsonify, url_for # <<<< ADD url_for HERE
import math
import datetime

# --- Create the Blueprint ---
diffusion_bp = Blueprint('diffusion', __name__,
                         template_folder='../templates', # Relative to this blueprint file's location
                         static_folder='../static',      # If this blueprint had its own static files
                         url_prefix='/calculators/diffusion')

# --- Helper Calculation Function (Chapman-Enskog) ---
# (Copied from your original calculator.py for encapsulation within this blueprint)
def local_calculate_chapman_enskog_DAB(T, P_bar, M_A, M_B, sigma_A, sigma_B, epsilon_k_A, epsilon_k_B):
    """
    Calculate the binary gas diffusion coefficient (D_AB) at low pressures using the Chapman-Enskog theory.
    P_bar is pressure in bar.
    """
    if T <= 0: raise ValueError("Temperature must be positive.")
    if P_bar <= 0: raise ValueError("Pressure must be positive.")
    if sigma_A <= 0 or sigma_B <= 0: raise ValueError("LJ sigma values must be positive.")
    if epsilon_k_A <= 0 or epsilon_k_B <= 0: raise ValueError("LJ epsilon/k values must be positive.")
    if M_A <= 0 or M_B <= 0: raise ValueError("Molar masses must be positive.")

    P_atm = P_bar / 1.01325 # Convert P from bar to atm for the formula constant 0.00266

    # Handle potential division by zero for M_AB
    if M_A == 0 or M_B == 0 or (1/M_A + 1/M_B) == 0 :
        M_AB = float('inf')
    else:
        M_AB = 2 / (1/M_A + 1/M_B)

    sigma_AB = (sigma_A + sigma_B) / 2
    epsilon_AB_k = math.sqrt(epsilon_k_A * epsilon_k_B) if epsilon_k_A >=0 and epsilon_k_B >=0 else float('nan')
    
    if epsilon_AB_k == 0:
        T_star = float('inf')
    else:
        T_star = T / epsilon_AB_k
    
    A, B, C, D_const, E, F, G, H = 1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996, 1.76474, 3.89411
    try:
        term1 = (A / (T_star ** B)) if T_star > 0 else float('inf') # Avoid T_star=0 or negative if B is not integer
        term2 = (C / math.exp(D_const * T_star))
        term3 = (E / math.exp(F * T_star))
        term4 = (G / math.exp(H * T_star))
        Omega_D = term1 + term2 + term3 + term4
    except (OverflowError, ValueError): # Catch potential math errors with T_star
        Omega_D = float('inf')

    denominator = (P_atm * (M_AB ** 0.5) * (sigma_AB ** 2) * Omega_D)
    D_AB = (0.00266 * T ** 1.5) / denominator if denominator != 0 and not math.isinf(denominator) else float('inf')
    return D_AB

# --- CHAPMAN-ENSKOG ROUTES ---
@diffusion_bp.route('/chapman-enskog')
def chapman_enskog():
    return render_template('calculator_chapman_enskog.html')

@diffusion_bp.route('/calculate_chapman_enskog', methods=['GET','POST'])
def calculate_chapman_enskog_api():
    today = datetime.date.today()
    site_url = request.host_url.rstrip('/')
    # Note: url_for inside blueprint needs '.endpoint_name'
    ce_page_url = site_url + url_for('diffusion.chapman_enskog')


    paper_citation_data = {
        'authors': 'Chapman, S., & Cowling, T. G.', 'year': 1970,
        'title': 'The Mathematical Theory of Non-uniform Gases', 'edition': '3rd Edition',
        'publisher': 'Cambridge University Press', 'doi': ''
    }
    tool_citation_data = {
        'author': 'ChemE Calc', 'year': today.year,
        'title': 'Chapman–Enskog Diffusion Calculator', 'retrieved_date': today.strftime('%Y-%m-%d'),
        'url': ce_page_url
    }

    if request.method == 'GET':
        return jsonify({
            'paper_citation_data': paper_citation_data,
            'tool_citation_data': tool_citation_data,
            'applicability_notes': "Valid for dilute gas mixtures (~50–5000 K, <20 bar). For D_AB in cm²/s, use T in K, P in bar, M in g/mol, σ in Å for inputs to this calculator."
        })

    data = request.get_json() or {}
    try:
        # These keys 'T_ce', 'P_ce' etc. must match what JS sends
        T = float(data.get('T_ce', 0)) # CORRECTED
        P_bar = float(data.get('P_ce', 0)) # CORRECTED
        M_A = float(data.get('M_A_ce', 0)) # CORRECTED
        M_B = float(data.get('M_B_ce', 0)) # CORRECTED
        sigma_A = float(data.get('sigma_A_ce', 0)) # CORRECTED
        sigma_B = float(data.get('sigma_B_ce', 0)) # CORRECTED
        epsilon_k_A = float(data.get('epsilon_k_A_ce', 0)) # CORRECTED
        epsilon_k_B = float(data.get('epsilon_k_B_ce', 0)) # CORRECTED

        D_AB = local_calculate_chapman_enskog_DAB(T, P_bar, M_A, M_B, sigma_A, sigma_B, epsilon_k_A, epsilon_k_B)

        output_unit = data.get('output_unit', 'cm^2/s') # Assuming JS sends output_unit from output_unit_ce
        D_AB_out = D_AB
        unit_latex = 'cm^2/s'
        if output_unit == 'm^2/s':
            D_AB_out = D_AB * 1e-4
            unit_latex = 'm^2/s'
        # Add other output unit conversions if needed

        sigma_unit_latex = '\\text{\\AA}'
        P_atm_display = P_bar / 1.01325
        M_AB_display = 2 / (1/M_A + 1/M_B) if M_A > 0 and M_B > 0 and (1/M_A + 1/M_B) !=0 else float('inf')
        sigma_AB_display = (sigma_A + sigma_B) / 2
        epsilon_AB_k_display = math.sqrt(epsilon_k_A * epsilon_k_B) if epsilon_k_A >=0 and epsilon_k_B >=0 else float('nan')
        T_star_display = T / epsilon_AB_k_display if epsilon_AB_k_display != 0 else float('inf')
        A_c, B_c, C_c, D_c_const, E_c, F_c, G_c, H_c = 1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996, 1.76474, 3.89411
        Omega_D_display = (A_c / (T_star_display ** B_c)) + (C_c / math.exp(D_c_const * T_star_display)) + (E_c / math.exp(F_c * T_star_display)) + (G_c / math.exp(H_c * T_star_display))

        latex_steps = (
            r"$$\text{\textbf{1. Chapman-Enskog Formula:}}$$"
            rf"$$D_{{AB}} = \frac{{0.00266 \, T^{{1.5}}}}{{P_{{atm}} \, \sqrt{{M_{{AB}}}} \, \sigma_{{AB}}^2 \, \Omega_D}}$$"
            rf"$$\text{{[Units for formula constant: }}D_{{AB}} \text{{ in cm²s⁻¹, }}T\text{{ in K, }}P_{{atm}}\text{{ in atm, }}M\text{{ in g/mol, }}\sigma\text{{ in {sigma_unit_latex}]}}$$"
            r"$$\text{\textbf{2a. Combined Molar Mass }} (M_{{AB}}):$$"
            rf"$$M_{{AB}} = 2 \left( \frac{{1}}{{{M_A:.4f}}} + \frac{{1}}{{{M_B:.4f}}} \right)^{{-1}} = {M_AB_display:.4f} \, \text{{g/mol}}$$"
            r"$$\text{\textbf{2b. Combined Collision Diameter }} (\sigma_{{AB}}):$$"
            rf"$$\sigma_{{AB}} = \frac{{{sigma_A:.4f} + {sigma_B:.4f}}}{{2}} = {sigma_AB_display:.4f} \, {sigma_unit_latex}$$"
            r"$$\text{\textbf{2c. Combined Characteristic Energy }} ((\epsilon/\kappa)_{{AB}}):$$"
            rf"$$(\epsilon/\kappa)_{{AB}} = \sqrt{{{epsilon_k_A:.2f} \times {epsilon_k_B:.2f}}} = {epsilon_AB_k_display:.2f} \, \text{{K}}$$"
            r"$$\text{\textbf{2d. Reduced Temperature }} (T^*_{{AB}}):$$"
            rf"$$T^*_{{AB}} = \frac{{{T:.2f}}}{{{epsilon_AB_k_display:.2f}}} = {T_star_display:.4f}$$"
            r"$$\text{\textbf{2e. Collision Integral }} (\Omega_D) \text{{ [Using Neufeld et al. approximation]}}:$$"
            r"$$\Omega_D = \frac{{A}}{{(T^*)^B}} + \frac{{C}}{{\exp(D T^*)}} + \frac{{E}}{{\exp(F T^*)}} + \frac{{G}}{{\exp(H T^*)}}$$"
            rf"$$\Omega_D = {Omega_D_display:.5f}$$"
            r"$$\text{\textbf{3. Substituting Values (P converted to atm for calculation):}}$$"
            rf"$$D_{{AB}} = \frac{{0.00266 \times ({T:.2f})^{{1.5}}}}{{({P_atm_display:.4f}) \sqrt{{{M_AB_display:.4f}}} ({sigma_AB_display:.4f})^2 ({Omega_D_display:.5f})}}$$"
            r"$$\text{\textbf{4. Final Result:}}$$"
            rf"$$D_{{AB}} = {D_AB_out:.6g} \, \text{{{unit_latex.replace('^2','²')}}}$$"
        )
        applicability_notes = "Valid for dilute gas mixtures (~50–5000 K, <20 bar for this calculator's typical usage of the 0.00266 constant form). User inputs P in bar."

        return jsonify({
            'value_base': D_AB, 'unit_base': 'cm^2/s', # Send base value in cm^2/s
            'latex_steps': latex_steps,
            'paper_citation_data': paper_citation_data, 'tool_citation_data': tool_citation_data,
            'applicability_notes': applicability_notes
        })
    except Exception as e:
        return jsonify({'error': str(e), 'paper_citation_data': paper_citation_data, 'tool_citation_data': tool_citation_data}), 400


# --- FULLER-SCHETTLER-GIDDINGS ROUTES ---
@diffusion_bp.route('/fuller-schettler-giddings')
def calculator_fsg():
    return render_template('calculator_fsg.html')

@diffusion_bp.route('/calculate_fsg', methods=['POST'])
def calculate_fsg_api():
    # This is the exact FSG API logic you had in app.py
    data = request.get_json() or {}
    try:
        T = float(data.get('T', 0)) # Assuming T from JS is in K
        P_bar = float(data.get('P', 0)) # Assuming P from JS is in bar
        M_A = float(data.get('M_A', 0))
        M_B = float(data.get('M_B', 0))
        sum_v_A = float(data.get('sum_v_A', 0))
        sum_v_B = float(data.get('sum_v_B', 0))

        # Validation
        errors = []
        if T <= 0: errors.append("Temperature must be > 0 K.")
        if P_bar <= 0: errors.append("Pressure must be > 0 bar.")
        if M_A <= 0 or M_B <= 0: errors.append("Molar masses must be > 0.")
        if sum_v_A <= 0 or sum_v_B <= 0: errors.append("Sum of atomic diffusion volumes must be > 0.")
        if errors:
            return jsonify({'error': "; ".join(errors)}), 400
        
        M_AB = 2 / (1/M_A + 1/M_B) if M_A > 0 and M_B > 0 and (1/M_A + 1/M_B) !=0 else float('inf')
        
        # FSG Calculation (Eq. 11-4.4 - uses P in bar directly)
        term_v_sum_cubed = (sum_v_A**(1/3)) + (sum_v_B**(1/3))
        denominator = (P_bar * (M_AB ** 0.5) * (term_v_sum_cubed**2))
        D_AB_fsg = (0.00143 * T**1.75) / denominator if denominator != 0 else float('inf')

        output_unit = data.get('output_unit', 'cm^2/s') # Key from FSG JS
        D_AB_fsg_out = D_AB_fsg
        unit_latex_fsg = 'cm^2/s'
        if output_unit == 'm^2/s':
            D_AB_fsg_out = D_AB_fsg * 1e-4
            unit_latex_fsg = 'm^2/s'
        # Add other output unit conversions for FSG if needed

        latex_steps_fsg = (
            r"$$\text{\textbf{1. Fuller-Schettler-Giddings Formula:}}$$"
            rf"$$D_{{AB}} = \frac{{0.00143 \, T^{{1.75}}}}{{P \, \sqrt{{M_{{AB}}}} \, \left( (\Sigma v_A)^{{1/3}} + (\Sigma v_B)^{{1/3}} \right)^2}}$$"
            rf"$$\text{{[Units for formula constant: }}D_{{AB}} \text{{ in cm²s⁻¹, }}T\text{{ in K, }}P\text{{ in bar, }}M\text{{ in g/mol, }}\Sigma v\text{{ unitless]}}$$"
            r"$$\text{\textbf{2a. Combined Molar Mass }} (M_{{AB}}):$$"
            rf"$$M_{{AB}} = 2 \left( \frac{{1}}{{{M_A:.4f}}} + \frac{{1}}{{{M_B:.4f}}} \right)^{{-1}} = {M_AB:.4f} \, \text{{g/mol}}$$"
            r"$$\text{\textbf{2b. Sum Atomic Diffusion Volumes Term:}}$$"
            rf"$$(\Sigma v_A)^{{1/3}} + (\Sigma v_B)^{{1/3}} = ({sum_v_A:.2f})^{{1/3}} + ({sum_v_B:.2f})^{{1/3}} = {term_v_sum_cubed:.4f}$$"
            r"$$\text{\textbf{3. Substituting Values into FSG Formula:}}$$"
            rf"$$D_{{AB}} = \frac{{0.00143 \times ({T:.2f})^{{1.75}}}}{{({P_bar:.4f}) \sqrt{{{M_AB:.4f}}} ({term_v_sum_cubed:.4f})^2 }}$$"
            r"$$\text{\textbf{4. Final Result:}}$$"
            rf"$$D_{{AB}} = {D_AB_fsg_out:.6g} \, \text{{{unit_latex_fsg.replace('^2','²')}}}$$"
        )
        
        site_url = request.host_url.rstrip('/')
        fsg_page_url = site_url + url_for('diffusion.calculator_fsg')
        today = datetime.date.today()
        paper_citation_fsg = {
            'authors': 'Fuller, E.N., Schettler, P.D., & Giddings, J.C.', 'year': 1966,
            'title': 'A New Method for Prediction of Binary Gas-Phase Diffusion Coefficients',
            'journal': 'J. Phys. Chem.', 'volume': '70(11)', 'pages': '3676–3680', 'doi': '10.1021/j100883a041' # Example DOI
        }
        tool_citation_fsg = {
            'author': 'ChemE Calc', 'year': today.year,
            'title': 'Fuller-Schettler-Giddings Diffusion Calculator', 'retrieved_date': today.strftime('%Y-%m-%d'),
            'url': fsg_page_url
        }
        applicability_fsg = "Generally reliable for nonpolar/polar gas pairs at low pressure. Avg. error ~4%. Relies on atomic diffusion volumes (Σv)."

        return jsonify({
            'value_base': D_AB_fsg, 'unit_base': 'cm^2/s',
            'latex_steps': latex_steps_fsg,
            'paper_citation_data': paper_citation_fsg,
            'tool_citation_data': tool_citation_fsg,
            'applicability_notes': applicability_fsg
        })
    except Exception as e:
        # Fallback for errors
        return jsonify({'error': str(e)}), 400