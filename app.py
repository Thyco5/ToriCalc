# D:\ChemE_Calc_New\app.py
import os
import math
import datetime
from flask import Flask, render_template, request, jsonify # Added request and jsonify
import requests # NEW IMPORT for making HTTP requests to PubChem
from urllib.parse import quote as url_quote # NEW IMPORT for URL encoding the chemical name

# --- Import Blueprints ---
from solvers.diffusion_bp import diffusion_bp
from solvers.critical_props_bp import critical_props_bp
from solvers.acentric_factor_bp import acentric_factor_bp
from solvers.phase_change_bp import phase_change_bp
from solvers.dipole_moment_bp import dipole_moment_bp
from solvers.ideal_gas_bp import ideal_gas_bp

# --- Base Directory Setup ---
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# --- App Initialization ---
app = Flask(__name__,
            static_folder=os.path.join(BASE_DIR, 'static'),
            template_folder=os.path.join(BASE_DIR, 'templates'))

app.config['TEMPLATES_AUTO_RELOAD'] = True

# --- Context Processor ---
@app.context_processor
def inject_current_year():
    return {'current_year': datetime.date.today().year}

# --- Register Blueprints ---
app.register_blueprint(diffusion_bp)
app.register_blueprint(critical_props_bp)
app.register_blueprint(acentric_factor_bp)
app.register_blueprint(phase_change_bp)
app.register_blueprint(dipole_moment_bp)
app.register_blueprint(ideal_gas_bp)

# --- Main App Routes ---
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/cheme-tools')
def cheme_tools():
    return render_template('cheme_tools.html')

@app.route('/sponsors')
def sponsors():
    return render_template('sponsors.html')

@app.route('/diffusion-calculators')
def diffusion_landing():
    return render_template('diffusion_landing.html')

@app.route('/calculators')
def calculators_landing_page():
    return render_template('calculators_landing_page.html')

@app.route('/calculators/critical-properties')
def critical_properties_solvers_page():
    return render_template('critical_properties_solvers.html')

# --- NEW API Endpoint: Chemical Name to SMILES ---
@app.route('/api/name_to_smiles', methods=['POST'])
def name_to_smiles_api():
    data = request.get_json()
    if not data or 'chemical_name' not in data:
        return jsonify({'error': 'Chemical name not provided.'}), 400

    chemical_name = data['chemical_name'].strip()
    if not chemical_name:
        return jsonify({'error': 'Chemical name cannot be empty.'}), 400

    print(f"API: Received request to convert name to SMILES: '{chemical_name}'")

    # URL encode the chemical name to handle spaces and special characters
    encoded_name = url_quote(chemical_name)
    
    # PubChem PUG REST API URL
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/JSON"

    try:
        response = requests.get(pubchem_url, timeout=10) # Added timeout
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)
        
        pubchem_data = response.json()
        
        # Navigate the PubChem JSON structure
        if pubchem_data and "PropertyTable" in pubchem_data and "Properties" in pubchem_data["PropertyTable"]:
            properties = pubchem_data["PropertyTable"]["Properties"]
            if properties and len(properties) > 0 and "CanonicalSMILES" in properties[0]:
                smiles = properties[0]["CanonicalSMILES"]
                print(f"API: Found SMILES for '{chemical_name}': '{smiles}'")
                return jsonify({'smiles': smiles, 'name_searched': chemical_name})
            else:
                print(f"API: CanonicalSMILES not found in PubChem response for '{chemical_name}'. Response: {pubchem_data}")
                return jsonify({'error': f"Could not find SMILES for '{chemical_name}'. The substance might not be in PubChem or the name is ambiguous."}), 404
        else:
            print(f"API: Unexpected PubChem JSON structure for '{chemical_name}'. Response: {pubchem_data}")
            return jsonify({'error': 'Unexpected response structure from PubChem.'}), 500

    except requests.exceptions.HTTPError as http_err:
        if response.status_code == 404:
            print(f"API: PubChem returned 404 for '{chemical_name}': {http_err}")
            return jsonify({'error': f"'{chemical_name}' not found in PubChem."}), 404
        else:
            print(f"API: HTTP error occurred with PubChem for '{chemical_name}': {http_err}")
            return jsonify({'error': f"Error communicating with PubChem: {http_err}"}), 502 # Bad Gateway
    except requests.exceptions.Timeout:
        print(f"API: Timeout while connecting to PubChem for '{chemical_name}'")
        return jsonify({'error': 'Timeout connecting to PubChem name resolver service.'}), 504 # Gateway Timeout
    except requests.exceptions.RequestException as req_err:
        print(f"API: Request error with PubChem for '{chemical_name}': {req_err}")
        return jsonify({'error': f"Error during request to PubChem: {req_err}"}), 500
    except ValueError as json_err: # Includes JSONDecodeError
        print(f"API: Could not decode JSON response from PubChem for '{chemical_name}': {json_err}")
        return jsonify({'error': 'Invalid JSON response from PubChem name resolver service.'}), 500


# --- Run the App ---
if __name__ == '__main__':
    app.run(debug=True, use_reloader=True)