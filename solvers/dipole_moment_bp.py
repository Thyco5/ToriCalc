import traceback
import datetime
from flask import Blueprint, render_template, request, jsonify, url_for
from chemicals.identifiers import search_chemical, int_to_CAS
from chemicals.dipole import dipole_moment, dipole_moment_methods

dipole_moment_bp = Blueprint('dipole_moment', __name__, template_folder='../templates')

@dipole_moment_bp.route('/calculators/dipole-moment')
def dipole_moment_page():
    return render_template('calculator_dipole_moment.html')

@dipole_moment_bp.route('/api/dipole/get_all_data', methods=['POST'])
def get_dipole_data_api():
    data = request.get_json() or {}
    chemical_name = data.get('chemical_name')
    if not chemical_name: return jsonify({"error": "Chemical name required."}), 400

    try:
        metadata = search_chemical(chemical_name)
        if not metadata: raise ValueError(f"Could not identify '{chemical_name}'.")
        
        cas_rn_int = metadata.CAS
        cas_rn_str = int_to_CAS(cas_rn_int) # CRITICAL FIX: Convert int to dashed-string
        
        available_methods = dipole_moment_methods(cas_rn_str)
        if not available_methods: raise ValueError(f"No dipole moment data found for '{metadata.common_name}'.")

        results = []
        for method in available_methods:
            value = dipole_moment(cas_rn_str, method=method)
            if value is not None:
                results.append({"method": method, "value": value})

        if not results: raise ValueError(f"Found data sources for '{metadata.common_name}', but none returned a value.")
        
        # --- Build Full, Structured Citation List ---
        today = datetime.date.today()
        page_url = request.url_root.rstrip('/') + url_for('dipole_moment.dipole_moment_page')
        tool_citation = {"id": "toricalc_dipole", "type": "Software", "source_type": "Tool", "data": { "authors": ["ToriCalc"], "year": today.year, "title": "Dipole Moment Data Comparator", "retrieved_date": today.strftime('%Y-%m-%d'), "url": page_url }}
        library_citation = {"id": "chemicals_lib", "type": "Software", "source_type": "Calculation Engine", "data": { "authors": ["Caleb Bell, Yoel Rene Cortes-Pena, and Contributors"], "year": "2016-2024", "title": "Chemicals", "url": "https://github.com/CalebBell/chemicals"}}
        
        source_citations = {
            'CCCBDB': {"id": "nist_cccbdb", "type": "Database", "source_type": "Data Source", "data": {"authors": ["National Institute of Standards and Technology"], "year": 2015, "title": "Computational Chemistry Comparison and Benchmark Database SRD 101", "url": "http://cccbdb.nist.gov/"}},
            'MULLER': {"id": "muller_2012", "type": "Journal Article", "source_type": "Data Source", "data": {"authors": ["Muller, K.", "Mokrushina, L.", "Arlt, W."], "year": 2012, "title": "Second-Order Group Contribution Method for the Determination of the Dipole Moment", "journal": "Journal of Chemical & Engineering Data", "volume": "57", "issue": "4", "pages": "1231-1236"}},
            'POLING': {"id": "poling_2000", "type": "Book", "source_type": "Data Source", "data": {"authors": ["Poling, B. E.", "Prausnitz, J. M.", "O'Connell, J. P."], "year": 2000, "title": "The Properties of Gases and Liquids, 5th Edition", "journal": "McGraw-Hill"}},
            'PSI4_2022A': {"id": "psi4_2012", "type": "Software", "source_type": "Data Source", "data": {"authors": ["Turney, J. M., et al."], "year": 2012, "title": "Psi4: An Open-Source Ab Initio Electronic Structure Program", "url": "https://doi.org/10.1002/wcms.93"}}
        }
        final_citations = [tool_citation, library_citation]
        for method in available_methods:
            if method in source_citations: final_citations.append(source_citations[method])

        return jsonify({
            "chemical_name": metadata.common_name or chemical_name,
            "cas": cas_rn_str,
            "results": results,
            "citations": final_citations
        })
    except ValueError as e: return jsonify({"error": str(e)}), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": "An unexpected server error occurred."}), 500