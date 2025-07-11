// chapman_enskog_calculator.js: Handles form, unit conversion, AJAX, MathJax, copy buttons for Chapman-Enskog calculator
(function(){
  // --- Lennard-Jones Parameters (Appendix B) ---
  const LENNARD_JONES_PARAMS = {
    'Air': { sigma: 3.711, epsilon_k: 78.6 },
    'Ar': { sigma: 3.542, epsilon_k: 93.3 },
    'CO2': { sigma: 3.941, epsilon_k: 195.2 },
    'H2': { sigma: 2.827, epsilon_k: 59.7 },
    'H2O': { sigma: 2.641, epsilon_k: 809.1 },
    'He': { sigma: 2.551, epsilon_k: 10.22 },
    'N2': { sigma: 3.798, epsilon_k: 71.4 },
    'O2': { sigma: 3.467, epsilon_k: 106.7 },
    'Kr': { sigma: 3.655, epsilon_k: 178.9 },
    'Ne': { sigma: 2.820, epsilon_k: 32.8 },
    'Xe': { sigma: 4.047, epsilon_k: 231.0 },
    'AsH3': { sigma: 4.145, epsilon_k: 259.8 },
    'BCl3': { sigma: 5.127, epsilon_k: 337.7 },
    'BF3': { sigma: 4.198, epsilon_k: 186.3 },
    'B(OCH3)3': { sigma: 5.503, epsilon_k: 396.7 },
    'Br2': { sigma: 4.296, epsilon_k: 507.9 },
    'CCl4': { sigma: 5.947, epsilon_k: 322.7 },
    'CF4': { sigma: 4.662, epsilon_k: 134.0 },
    'CHCl3': { sigma: 5.389, epsilon_k: 340.2 },
    'CH2Cl2': { sigma: 4.898, epsilon_k: 356.3 },
    'CH3Br': { sigma: 4.118, epsilon_k: 449.2 },
    'CH3Cl': { sigma: 4.182, epsilon_k: 350 },
    'CH3OH': { sigma: 3.626, epsilon_k: 481.8 },
    'CH4': { sigma: 3.758, epsilon_k: 148.6 },
    'CO': { sigma: 3.690, epsilon_k: 91.7 },
    'COS': { sigma: 4.130, epsilon_k: 336.0 },
    'CS2': { sigma: 4.483, epsilon_k: 467 },
    'C2H2': { sigma: 4.033, epsilon_k: 231.8 },
    'C2H4': { sigma: 4.163, epsilon_k: 224.7 },
    'C2H6': { sigma: 4.443, epsilon_k: 215.7 },
    'C2H5Cl': { sigma: 4.898, epsilon_k: 300 },
    'C2H5OH': { sigma: 4.530, epsilon_k: 362.6 },
    'C2N2': { sigma: 4.361, epsilon_k: 348.6 },
    'CH3OCH3': { sigma: 4.307, epsilon_k: 395.0 },
    'CH2CHCH3': { sigma: 4.678, epsilon_k: 298.9 },
    'CH3CCH': { sigma: 4.761, epsilon_k: 251.8 },
    'C3H6': { sigma: 4.807, epsilon_k: 248.9 },
    'C3H8': { sigma: 5.118, epsilon_k: 237.1 },
    'nC3H7OH': { sigma: 4.549, epsilon_k: 576.7 },
    'CH3COCH3': { sigma: 4.600, epsilon_k: 560.2 },
    'CH3COOCH3': { sigma: 4.936, epsilon_k: 469.8 },
    'nC4H10': { sigma: 4.687, epsilon_k: 531.4 },
    'isoC4H10': { sigma: 5.278, epsilon_k: 330.1 },
    'C2H5OC2H5': { sigma: 5.678, epsilon_k: 313.8 },
    'CH3COOC2H5': { sigma: 5.205, epsilon_k: 521.3 },
    'nC5H12': { sigma: 5.784, epsilon_k: 341.1 },
    'C(CH3)4': { sigma: 6.464, epsilon_k: 193.4 },
    'C6H6': { sigma: 5.349, epsilon_k: 412.3 },
    'C6H12': { sigma: 6.182, epsilon_k: 297.1 },
    'nC6H14': { sigma: 5.949, epsilon_k: 399.3 },
    'Cl2': { sigma: 4.217, epsilon_k: 316.0 },
    'F2': { sigma: 3.357, epsilon_k: 112.6 },
    'HBr': { sigma: 3.353, epsilon_k: 449 },
    'HCN': { sigma: 3.630, epsilon_k: 569.1 },
    'HCl': { sigma: 3.339, epsilon_k: 344.7 },
    'HF': { sigma: 3.148, epsilon_k: 330 },
    'HI': { sigma: 4.211, epsilon_k: 288.7 },
    'H2O2': { sigma: 4.196, epsilon_k: 289.3 },
    'H2S': { sigma: 3.623, epsilon_k: 301.1 },
    'Hg': { sigma: 2.969, epsilon_k: 750 },
    'HgBr2': { sigma: 5.080, epsilon_k: 686.2 },
    'HgCl2': { sigma: 4.550, epsilon_k: 750 },
    'HgI2': { sigma: 5.625, epsilon_k: 695.6 },
    'I2': { sigma: 5.160, epsilon_k: 474.2 },
    'NH3': { sigma: 2.900, epsilon_k: 558.3 },
    'NO': { sigma: 3.492, epsilon_k: 116.7 },
    'NOCl': { sigma: 4.112, epsilon_k: 395.3 },
    'N2O': { sigma: 3.828, epsilon_k: 232.4 },
    'PH3': { sigma: 3.981, epsilon_k: 251.5 },
    'SF6': { sigma: 5.128, epsilon_k: 222.1 },
    'SO2': { sigma: 4.112, epsilon_k: 335.4 },
    'SiF4': { sigma: 4.880, epsilon_k: 171.9 },
    'SiH4': { sigma: 4.084, epsilon_k: 207.6 },
    'SnBr4': { sigma: 6.388, epsilon_k: 563.7 },
    'UF6': { sigma: 5.967, epsilon_k: 236.8 }
  };

  const API_ENDPOINT = '/calculators/diffusion/calculate_chapman_enskog';
  
  const toInternal = {
    T: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'K') return v;
      if(u === 'degC') return v + 273.15;
      if(u === 'degF') return (v - 32) * 5/9 + 273.15;
    },
    P: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'bar') return v;
      if(u === 'atm') return v * 1.01325;
      if(u === 'Pa') return v / 100000;
      if(u === 'kPa') return v / 100;
      if(u === 'psi') return v / 14.5038;
      if(u === 'mmHg') return v / 750.062;
    },
    M: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'g/mol') return v;
      if(u === 'kg/mol') return v * 1000;
    },
    sigma: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'A') return v;
      if(u === 'nm') return v * 10;
      if(u === 'm') return v * 1e10;
    },
    epsilon_k: (v, u) => { // u is always K for epsilon_k
        if (isNaN(v)) return NaN;
        return v;
    }
  };
  
  const fromInternal = {
    T: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'K') return v;
      if(u === 'degC') return v - 273.15;
      if(u === 'degF') return (v - 273.15) * 9/5 + 32;
    },
    P: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'bar') return v;
      if(u === 'atm') return v / 1.01325;
      if(u === 'Pa') return v * 100000;
      if(u === 'kPa') return v * 100;
      if(u === 'psi') return v * 14.5038;
      if(u === 'mmHg') return v * 750.062;
    },
    M: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'g/mol') return v;
      if(u === 'kg/mol') return v / 1000;
    },
    sigma: (v, u) => {
      if (isNaN(v)) return NaN;
      if(u === 'A') return v;
      if(u === 'nm') return v / 10;
      if(u === 'm') return v / 1e10;
    },
    epsilon_k: (v, u) => { // u is always K
        if (isNaN(v)) return NaN;
        return v;
    }
  };
  
  const paramRanges = {
    T: { min: 300, max: 2000, unit: 'K' }, // Adjusted min slightly above 0 for practical positive check
    P: { min: 0.00001, max: 10, unit: 'bar' }, // Adjusted min slightly above 0
    // No specific ranges for M, sigma, epsilon_k beyond being positive
    M_A: {min: 0.001}, M_B: {min: 0.001}, sigma_A: {min: 0.001}, sigma_B: {min: 0.001}, epsilon_k_A: {min: 0.001}, epsilon_k_B: {min: 0.001}
  };

  // Moved placeholderMap declaration earlier
  const placeholderMap = {
    T_ce: {K: 'e.g., 298.15', degC: 'e.g., 25', degF: 'e.g., 77'},
    P_ce: {bar: 'e.g., 1.01325', atm: 'e.g., 1', Pa: 'e.g., 101325', kPa: 'e.g., 101.325', psi: 'e.g., 14.7', mmHg: 'e.g., 760'},
    M_A_ce: {'g/mol': 'e.g., 28.01 (N₂)', 'kg/mol': 'e.g., 0.02801'},
    M_B_ce: {'g/mol': 'e.g., 32.00 (O₂)', 'kg/mol': 'e.g., 0.03200'},
    sigma_A_ce: {A: 'e.g., 3.798 (N₂)', nm: 'e.g., 0.3798', m: 'e.g., 3.798e-10'},
    sigma_B_ce: {A: 'e.g., 3.467 (O₂)', nm: 'e.g., 0.3467', m: 'e.g., 3.467e-10'},
    epsilon_k_A_ce: {K: 'e.g., 71.4 (N₂)'}, // Unit is fixed to K
    epsilon_k_B_ce: {K: 'e.g., 106.7 (O₂)'}  // Unit is fixed to K
  };

  // DOM Elements
  const form = document.getElementById('chapman-enskog-form');
  const resultPlaceholder = document.getElementById('result-placeholder-ce');
  const resultDisplayContainer = document.getElementById('result-display-ce');
  const resultValueDisplay = document.getElementById('result-value-ce');
  const resultUnitDisplay = document.getElementById('result-unit-ce');
  const resultHr = document.getElementById('result-hr-ce');
  const outputUnitSelect = document.getElementById('output_unit_ce');
  const latexDiv = document.getElementById('latex-steps-ce');
  const citationCeDisplay = document.getElementById('citation-ce-display');
  const applicabilityDiv = document.getElementById('applicability-notes-ce');
  const toolCitationCeDisplay = document.getElementById('tool-citation-ce-display');
  const copyLatexBtn = document.getElementById('copy-latex-btn-ce');
  const copyCitationBtn = document.getElementById('copy-citation-btn-ce');
  const copyToolCitationBtn = document.getElementById('copy-tool-citation-btn-ce');
  const warningEl = document.getElementById('range-warning-ce');

  let currentResult = { value: null, unit: 'cm^2/s' };
  let latexRaw = '';
  let paperCitationData = null;
  let toolCitationData = null;

  // Helper to populate molecule dropdowns
  function populateMoleculeDropdown(selectId) {
    const select = document.getElementById(selectId);
    if (!select) return;
    select.innerHTML = '<option value="">-- Select Molecule (Optional) --</option>';
    Object.entries(LENNARD_JONES_PARAMS).forEach(([key, val]) => {
      let label = key;
      if (key === 'N2') label = 'Nitrogen (N₂)'; else if (key === 'O2') label = 'Oxygen (O₂)';
      else if (key === 'H2') label = 'Hydrogen (H₂)'; else if (key === 'CO2') label = 'Carbon Dioxide (CO₂)';
      else if (key === 'H2O') label = 'Water (H₂O)'; else if (key === 'Ar') label = 'Argon (Ar)';
      else if (key === 'Air') label = 'Air';
      select.innerHTML += `<option value="${key}">${label}</option>`;
    });
  }
  
  function escapeHtml(unsafe) {
    return unsafe
         .replace(/&/g, "&amp;")
         .replace(/</g, "&lt;")
         .replace(/>/g, "&gt;")
         .replace(/"/g, "&quot;")
         .replace(/'/g, "&#039;");
  }

  form.addEventListener('submit', function(e) { e.preventDefault(); calculateAndDisplay(); });

  document.addEventListener('DOMContentLoaded', () => {
    if (resultPlaceholder) resultPlaceholder.style.display = 'block';
    if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
    if (resultHr) resultHr.style.display = 'none';
    if (latexDiv) latexDiv.innerHTML = "Step-by-step calculation details will appear here after calculation.";

    populateMoleculeDropdown('molecule_select_A_ce');
    populateMoleculeDropdown('molecule_select_B_ce');

    const molSelA = document.getElementById('molecule_select_A_ce');
    if (molSelA) {
      molSelA.addEventListener('change', function() {
        const key = molSelA.value;
        if (key && LENNARD_JONES_PARAMS[key]) {
          form.sigma_A_ce.value = LENNARD_JONES_PARAMS[key].sigma;
          form.epsilon_k_A_ce.value = LENNARD_JONES_PARAMS[key].epsilon_k;
          form.sigma_A_ce_unit.value = 'A';
          form.epsilon_k_A_ce_unit.value = 'K';
          calculateAndDisplay();
        }
      });
      form.sigma_A_ce.addEventListener('input', function() { molSelA.value = ''; calculateAndDisplay(); });
      form.epsilon_k_A_ce.addEventListener('input', function() { molSelA.value = ''; calculateAndDisplay(); });
    }

    const molSelB = document.getElementById('molecule_select_B_ce');
    if (molSelB) {
      molSelB.addEventListener('change', function() {
        const key = molSelB.value;
        if (key && LENNARD_JONES_PARAMS[key]) {
          form.sigma_B_ce.value = LENNARD_JONES_PARAMS[key].sigma;
          form.epsilon_k_B_ce.value = LENNARD_JONES_PARAMS[key].epsilon_k;
          form.sigma_B_ce_unit.value = 'A';
          form.epsilon_k_B_ce_unit.value = 'K';
          calculateAndDisplay();
        }
      });
      form.sigma_B_ce.addEventListener('input', function() { molSelB.value = ''; calculateAndDisplay(); });
      form.epsilon_k_B_ce.addEventListener('input', function() { molSelB.value = ''; calculateAndDisplay(); });
    }

    updateCitationDisplays();
    fetch(API_ENDPOINT, { method: 'GET' })
      .then(res => res.json())
      .then(json => {
        if (json.paper_citation_data) paperCitationData = json.paper_citation_data;
        if (json.tool_citation_data) toolCitationData = json.tool_citation_data;
        if (json.applicability_notes && applicabilityDiv) applicabilityDiv.textContent = json.applicability_notes;
        updateCitationDisplays();
      }).catch(err => console.error("Error fetching initial citation data:", err));

    // calculateAndDisplay(); // Avoid initial calculation on page load if fields are empty

    form.querySelectorAll('input[type="number"]').forEach(inp => inp.addEventListener('input', calculateAndDisplay));
    form.querySelectorAll('select').forEach(sel => {
        // The output unit select is handled by this generic listener for recalculation
        sel.addEventListener('change', calculateAndDisplay);
    });

    const paperFormatSel = document.getElementById('paper-citation-format-ce');
    const toolFormatSel = document.getElementById('tool-citation-format-ce');
    if (paperFormatSel) paperFormatSel.addEventListener('change', updateCitationDisplays);
    if (toolFormatSel) toolFormatSel.addEventListener('change', updateCitationDisplays);

    if (copyLatexBtn) copyLatexBtn.addEventListener('click', () => copyWithFeedback(copyLatexBtn, latexRaw));
    if (copyCitationBtn) copyCitationBtn.addEventListener('click', () => {
        if(citationCeDisplay) copyWithFeedback(copyCitationBtn, citationCeDisplay.innerText || citationCeDisplay.textContent);
    });
    if (copyToolCitationBtn) copyToolCitationBtn.addEventListener('click', () => {
        if(toolCitationCeDisplay) copyWithFeedback(copyToolCitationBtn, toolCitationCeDisplay.innerText || toolCitationCeDisplay.textContent);
    });

    form.addEventListener('reset', function() {
        if (resultPlaceholder) resultPlaceholder.style.display = 'block';
        if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
        if (resultValueDisplay) resultValueDisplay.textContent = '';
        if (resultUnitDisplay) resultUnitDisplay.textContent = '';
        if (resultHr) resultHr.style.display = 'none';
        if (latexDiv) latexDiv.innerHTML = 'Step-by-step calculation details will appear here.';
        if (warningEl) { warningEl.innerHTML = ''; warningEl.style.display = 'none'; }
        currentResult = { value: null, unit: 'cm^2/s' };
        latexRaw = '';
        const latexCollapseTarget = document.getElementById('latexStepsContent');
        if (latexCollapseTarget && latexCollapseTarget.classList.contains('show')) {
            new bootstrap.Collapse(latexCollapseTarget, { toggle: false }).hide();
        }
        ['collapseCitation', 'collapseToolCitation', 'collapseApplicability'].forEach(id => {
            const target = document.getElementById(id);
            if (target && target.classList.contains('show')) {
                new bootstrap.Collapse(target, { toggle: false }).hide();
            }
        });
        ['T_ce', 'P_ce', 'M_A_ce', 'M_B_ce', 'sigma_A_ce', 'sigma_B_ce', 'epsilon_k_A_ce', 'epsilon_k_B_ce'].forEach(input => {
            const inputEl = document.getElementById(input);
            const unitEl = document.getElementById(input + '_unit');
            if (unitEl) { // If there's a unit element, update placeholder based on its current value
                 updateRangePlaceholder(input, unitEl.value);
            } else if (inputEl && placeholderMap[input] && placeholderMap[input]['K']) { // Esp for epsilon_k fixed unit
                 inputEl.placeholder = placeholderMap[input]['K'];
            }
        });
    });

    // Initialize placeholders for all relevant inputs on load
    ['T_ce', 'P_ce', 'M_A_ce', 'M_B_ce', 'sigma_A_ce', 'sigma_B_ce', 'epsilon_k_A_ce', 'epsilon_k_B_ce'].forEach(input => {
        const unitEl = document.getElementById(input + '_unit');
        if (unitEl) {
            updateRangePlaceholder(input, unitEl.value);
        } else { // For epsilon_k which has fixed K unit and hidden select
            updateRangePlaceholder(input, 'K');
        }
    });
  });

  function getInputs() {
    return {
      T_ce: toInternal.T(Number(form.T_ce.value), form.T_ce_unit.value),
      P_ce: toInternal.P(Number(form.P_ce.value), form.P_ce_unit.value),
      M_A_ce: toInternal.M(Number(form.M_A_ce.value), form.M_A_ce_unit.value),
      M_B_ce: toInternal.M(Number(form.M_B_ce.value), form.M_B_ce_unit.value),
      sigma_A_ce: toInternal.sigma(Number(form.sigma_A_ce.value), form.sigma_A_ce_unit.value),
      sigma_B_ce: toInternal.sigma(Number(form.sigma_B_ce.value), form.sigma_B_ce_unit.value),
      epsilon_k_A_ce: toInternal.epsilon_k(Number(form.epsilon_k_A_ce.value), form.epsilon_k_A_ce_unit.value),
      epsilon_k_B_ce: toInternal.epsilon_k(Number(form.epsilon_k_B_ce.value), form.epsilon_k_B_ce_unit.value)
    };
  }

  function getOutputUnit() {
    return outputUnitSelect.value;
  }

  function formatCitation(data, format) {
    if (!data) return '';
    let textContent = '';
    if (data.authors && data.title) {
        const authors = data.authors || 'N/A'; const year = data.year || 'N/A';
        const title = data.title || 'N/A'; const edition = data.edition || '';
        const publisher = data.publisher || 'N/A'; const doi = data.doi || '';
        if (format === 'APA') textContent = `${authors} (${year}). *${title}*${edition ? ' (' + edition + ').' : '.'} ${publisher}.`;
        else if (format === 'ACS') textContent = `${authors}. *${title}*${edition ? ', ' + edition : ''}; ${publisher}, ${year}.`;
        else if (format === 'BibTeX') return `<pre>@book{ref_${year.replace(/[^a-zA-Z0-9]/g, "")}{${Math.random().toString(36).substring(2,8)}}_paper,
  author    = {${authors}},
  title     = {${title}},
  year      = {${year}},
  ${edition ? `edition   = {${edition}},` : ''}
  publisher = {${publisher}},
  ${doi ? `doi       = {${doi}}` : ''}
}</pre>`;
        else textContent = `${authors} (${year}). *${title}*. ${edition ? edition + ". " : ""}${publisher}.`;
    } else if (data.tool_author && data.tool_title) { 
        const author = data.tool_author || 'N/A'; const year = data.tool_year || 'N/A';
        const title = data.tool_title || 'N/A'; const retrieved_date = data.retrieved_date || new Date().toLocaleDateString();
        const url = data.url || window.location.href;
        if (format === 'APA') textContent = `${author}. (${year}). *${title}* [Web Application]. ToriCalc. Retrieved ${retrieved_date}, from ${url}`;
        else if (format === 'ACS') textContent = `${author}. ${title} [Online]; ToriCalc: ${url} (accessed ${retrieved_date}).`;
        else if (format === 'BibTeX') return `<pre>@misc{ToriCalc_${year.replace(/[^a-zA-Z0-9]/g, "")}_${Math.random().toString(36).substring(2,8)}_tool,
  author    = {${author}},
  title     = {${title}},
  year      = {${year}},
  howpublished = {Web Application},
  note      = {Retrieved ${retrieved_date} from ${url}},
  publisher = {ToriCalc}
}</pre>`;
        else textContent = `${author}. (${year}). *${title}* [Web Application]. ToriCalc. Retrieved ${retrieved_date}, from ${url}.`;
    } else return "Citation data not available or in unexpected format.";
    return textContent;
  }

  function updateCitationDisplays() {
    const paperFormat = document.getElementById('paper-citation-format-ce').value;
    const toolFormat = document.getElementById('tool-citation-format-ce').value;
    if (citationCeDisplay) {
        const formattedPaperCitation = formatCitation(paperCitationData, paperFormat);
        if (paperFormat === 'BibTeX') citationCeDisplay.innerHTML = formattedPaperCitation;
        else citationCeDisplay.textContent = formattedPaperCitation;
    }
    if (toolCitationCeDisplay) {
        const formattedToolCitation = formatCitation(toolCitationData, toolFormat);
         if (toolFormat === 'BibTeX') toolCitationCeDisplay.innerHTML = formattedToolCitation;
        else toolCitationCeDisplay.textContent = formattedToolCitation;
    }
    const doiSpan = document.getElementById('doi-link-ce');
    if (doiSpan) {
    if (paperCitationData && paperCitationData.doi) {
      doiSpan.innerHTML = `<a href="https://doi.org/${paperCitationData.doi}" target="_blank">https://doi.org/${paperCitationData.doi}</a>`;
        } else doiSpan.innerHTML = '';
    }
  }

  async function calculateAndDisplay() {
    const inputElementsIds = ['T_ce','P_ce','M_A_ce','M_B_ce','sigma_A_ce','sigma_B_ce','epsilon_k_A_ce','epsilon_k_B_ce'];
    let allRequiredFilledAndValid = true;
    let firstErrorMessage = "";

    const data = getInputs(); // This now returns converted values

    // Validate core inputs (Temperature and Pressure must be positive after conversion)
    if (isNaN(data.T_ce) || data.T_ce <= 0) {
        firstErrorMessage = 'Temperature must be a positive value (e.g., > 0 K after conversion).';
        allRequiredFilledAndValid = false;
    } else if (isNaN(data.P_ce) || data.P_ce <= 0) {
        firstErrorMessage = 'Pressure must be a positive value.';
        allRequiredFilledAndValid = false;
    }
    // Check other fields for being numbers and positive if applicable
    inputElementsIds.forEach(id => {
        if (id !== 'T_ce' && id !== 'P_ce') { // T_ce and P_ce already checked
            const value = data[id];
            const paramMin = paramRanges[id.replace('_ce','')]?.min;
            if (isNaN(value) || (paramMin !== undefined && value < paramMin) ) {
                if (allRequiredFilledAndValid) { // Capture first error for other fields
                    firstErrorMessage = `Field ${id.replace('_ce','')} must be a valid positive number.`;
                }
                allRequiredFilledAndValid = false;
            }
        }
    });
    
    // Check if form fields are empty - this is a basic check for presence, not validity which is done above
    const allRequiredDomFieldsHaveValue = inputElementsIds.every(id => form[id].value.trim() !== '');
    if(!allRequiredDomFieldsHaveValue && allRequiredFilledAndValid){
        firstErrorMessage = 'Please fill in all required fields.';
        allRequiredFilledAndValid = false;
    }

    if (!allRequiredFilledAndValid) {
      if (warningEl) { warningEl.textContent = firstErrorMessage; warningEl.style.display = 'block'; }
      if (resultPlaceholder) resultPlaceholder.style.display = 'block';
      if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
      if (resultHr) resultHr.style.display = 'none';
      if (latexDiv) latexDiv.innerHTML = "Error in input. Please check values.";
      return;
    }
    
    if (warningEl) warningEl.style.display = 'none';

    const nonBlockingWarnings = [];
    if (data.T_ce < paramRanges.T.min || data.T_ce > paramRanges.T.max) {
      nonBlockingWarnings.push(`Temperature ${form.T_ce.value}${form.T_ce_unit.value} (${data.T_ce.toFixed(2)} K) is outside typical validation range ${paramRanges.T.min}-${paramRanges.T.max} K.`);
    }
    if (data.P_ce > paramRanges.P.max || data.P_ce < paramRanges.P.min ) {
      nonBlockingWarnings.push(`Pressure ${form.P_ce.value}${form.P_ce_unit.value} (${data.P_ce.toFixed(3)} bar) is outside typical validation range ${paramRanges.P.min}-${paramRanges.P.max} bar.`);
    }
    if (nonBlockingWarnings.length && warningEl) {
      warningEl.innerHTML = nonBlockingWarnings.join('<br>');
      warningEl.style.display = 'block';
    } else if (warningEl) {
      warningEl.style.display = 'none';
    }

    try {
      const payload = { ...data, output_unit: getOutputUnit() };
      const resp = await fetch(API_ENDPOINT, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(payload)
      });
      const json = await resp.json();

      if (!resp.ok || json.error || json.detail) {
        const errorMsg = json.error || json.detail || (resp.statusText + " (HTTP " + resp.status + ")") || 'Unknown server error';
        if (warningEl) { warningEl.textContent = 'Error: ' + errorMsg; warningEl.style.display = 'block';}
        if (resultPlaceholder) resultPlaceholder.style.display = 'block';
        if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
        if (resultHr) resultHr.style.display = 'none';
        if (latexDiv) latexDiv.innerHTML = "Error during calculation from server.";
        return;
      }
      
      currentResult.value = Number(json.value_base);
      currentResult.unit = json.unit_base;
      latexRaw = json.latex_steps;
      updateResultDisplay(json.D_AB_formatted, json.output_unit_symbol);

      if (latexDiv) {
        latexDiv.innerHTML = ""; 
        const tempDiv = document.createElement('div'); tempDiv.innerHTML = latexRaw;
        latexDiv.textContent = tempDiv.textContent || tempDiv.innerText || "";
        if (window.MathJax && window.MathJax.typesetPromise) {
          window.MathJax.typesetPromise([latexDiv]).catch(err => {
            console.error('MathJax typesetting failed:', err);
            latexDiv.innerHTML = "Error rendering math. Raw LaTeX: <br><pre>" + escapeHtml(latexRaw) + "</pre>";
          });
        } else if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Queue) {
          window.MathJax.Hub.Queue(["Typeset", window.MathJax.Hub, latexDiv]);
        } else {
            latexDiv.innerHTML = "MathJax not loaded. Raw LaTeX: <br><pre>" + escapeHtml(latexRaw) + "</pre>";
        }
      }
      if (json.paper_citation_data) paperCitationData = json.paper_citation_data;
      if (json.tool_citation_data) toolCitationData = json.tool_citation_data;
      updateCitationDisplays();
      if (json.applicability_notes && applicabilityDiv) applicabilityDiv.textContent = json.applicability_notes;
      if (resultPlaceholder) resultPlaceholder.style.display = 'none';
      if (resultDisplayContainer) resultDisplayContainer.style.display = 'block';
      if (resultHr) resultHr.style.display = 'block';
    } catch (err) {
      console.error("Calculation or Display Error:", err);
      if (warningEl) { warningEl.textContent = 'Client-side Error: ' + err.message; warningEl.style.display = 'block'; }
      if (resultPlaceholder) resultPlaceholder.style.display = 'block';
      if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
      if (resultHr) resultHr.style.display = 'none';
      if (latexDiv) latexDiv.innerHTML = "An unexpected error occurred on the client.";
    }
  }

  // This function now primarily takes the pre-formatted values from backend
  // or falls back to client-side conversion if those aren't provided.
  function updateResultDisplay(formattedValue, formattedUnitSymbol) {
    if (formattedValue === null || formattedValue === undefined || formattedUnitSymbol === null || formattedUnitSymbol === undefined ) {
        // This case might occur if calculation fails before backend formatting
        // or if only output unit is changed and we need to re-convert client side
        if (currentResult.value !== null && currentResult.unit !== null) {
            const targetDisplayUnit = getOutputUnit();
            const factors = { 'cm^2/s': 1, 'm^2/s': 1e-4, 'mm^2/s': 100, 'ft^2/hr': 387.5008 };
            let baseValueInCm2s = currentResult.value;
            // Ensure baseValue is in cm^2/s before converting to targetDisplayUnit
            if (currentResult.unit !== 'cm^2/s' && factors[currentResult.unit]) {
                baseValueInCm2s = currentResult.value / factors[currentResult.unit];
            } else if (!factors[currentResult.unit]) {
                 console.warn(`Base unit ${currentResult.unit} for currentResult not recognized. Assuming cm^2/s.`);
            }

            let displayVal;
            let displayUnitSym = targetDisplayUnit;
            if (factors[targetDisplayUnit]) {
                displayVal = (baseValueInCm2s * factors[targetDisplayUnit]).toPrecision(6);
            } else {
                console.warn(`Unsupported output unit conversion: ${targetDisplayUnit}. Displaying in cm^2/s.`);
                displayVal = baseValueInCm2s.toPrecision(6);
                displayUnitSym = 'cm^2/s';
            }
            if (resultValueDisplay) resultValueDisplay.textContent = displayVal;
            if (resultUnitDisplay) resultUnitDisplay.textContent = displayUnitSym;
        } else {
            // No current result and no formatted values, show placeholder
            if (resultPlaceholder) resultPlaceholder.style.display = 'block';
            if (resultDisplayContainer) resultDisplayContainer.style.display = 'none';
            if (resultHr) resultHr.style.display = 'none';
            return;
        }
    } else {
        // Use backend-provided formatted values
        if (resultValueDisplay) resultValueDisplay.textContent = formattedValue;
        if (resultUnitDisplay) resultUnitDisplay.textContent = formattedUnitSymbol;
    }

    if (resultPlaceholder) resultPlaceholder.style.display = 'none';
    if (resultDisplayContainer) resultDisplayContainer.style.display = 'block';
    if (resultHr) resultHr.style.display = 'block';
  }

  function convertInputValueOnUnitChange(inputId, oldUnit, newUnit) {
    const inputEl = document.getElementById(inputId);
    if (!inputEl || !inputEl.value) return;
    let value = parseFloat(inputEl.value);
    if (isNaN(value)) return; // Don't attempt to convert if not a number
    
    let conversionType;
    if (inputId === 'T_ce') conversionType = 'T';
    else if (inputId === 'P_ce') conversionType = 'P';
    else if (inputId === 'M_A_ce' || inputId === 'M_B_ce') conversionType = 'M';
    else if (inputId === 'sigma_A_ce' || inputId === 'sigma_B_ce') conversionType = 'sigma';
    else if (inputId === 'epsilon_k_A_ce' || inputId === 'epsilon_k_B_ce') conversionType = 'epsilon_k';
    else return;

    const internalValue = toInternal[conversionType](value, oldUnit);
    const newValueInNewUnit = fromInternal[conversionType](internalValue, newUnit);
    
    if (isNaN(newValueInNewUnit)) { inputEl.value = ''; return; }

    if (Math.abs(newValueInNewUnit) > 0 && (Math.abs(newValueInNewUnit) < 0.001 || Math.abs(newValueInNewUnit) > 100000)) {
        inputEl.value = newValueInNewUnit.toExponential(4);
    } else if (Math.abs(newValueInNewUnit) >= 0) { 
        inputEl.value = parseFloat(newValueInNewUnit.toPrecision(6)); 
    } else {
        inputEl.value = '';
    }
  }
  
  function updateRangePlaceholder(inputId, unitValue) {
    const inputEl = document.getElementById(inputId);
    if (!inputEl) return;
    let conversionType;
    if (inputId === 'T_ce') conversionType = 'T';
    else if (inputId === 'P_ce') conversionType = 'P';
    else { // For other inputs, use the generic placeholder from placeholderMap
        if (placeholderMap[inputId] && placeholderMap[inputId][unitValue]) {
            inputEl.placeholder = placeholderMap[inputId][unitValue];
        } else if (placeholderMap[inputId] && placeholderMap[inputId]['K']) { // for epsilon_k
            inputEl.placeholder = placeholderMap[inputId]['K'];
        } else {
            inputEl.placeholder = 'Enter value';
        }
        return;
    }
    
    const range = paramRanges[conversionType];
    let placeholderText = '';
    if (range) {
        if (range.min !== undefined && range.max !== undefined) {
      const minConverted = fromInternal[conversionType](range.min, unitValue);
      const maxConverted = fromInternal[conversionType](range.max, unitValue);
            if(!isNaN(minConverted) && !isNaN(maxConverted)) placeholderText = `Range: ${minConverted.toFixed(1)} - ${maxConverted.toFixed(1)} ${unitValue}`;
        } else if (range.max !== undefined) {
      const maxConverted = fromInternal[conversionType](range.max, unitValue);
            if(!isNaN(maxConverted)) placeholderText = `Max: ${maxConverted.toFixed(1)} ${unitValue}`;
        } else if (range.min !== undefined) {
             const minConverted = fromInternal[conversionType](range.min, unitValue);
            if(!isNaN(minConverted)) placeholderText = `Min: ${minConverted.toFixed(1)} ${unitValue}`;
        }
    }

    if (placeholderText) {
        inputEl.placeholder = placeholderText;
    } else if (placeholderMap[inputId] && placeholderMap[inputId][unitValue]) {
         inputEl.placeholder = placeholderMap[inputId][unitValue];
    } else {
        inputEl.placeholder = 'Enter value';
    }
  }
  
  const currentUnits = {};
  ['T_ce', 'P_ce', 'M_A_ce', 'M_B_ce', 'sigma_A_ce', 'sigma_B_ce', 'epsilon_k_A_ce', 'epsilon_k_B_ce'].forEach(inputId => {
    const unitElId = inputId + '_unit';
    const unitEl = document.getElementById(unitElId);
    if (unitEl) {
        currentUnits[inputId] = unitEl.value;
        unitEl.addEventListener('change', () => {
      const oldUnit = currentUnits[inputId];
            const newUnit = unitEl.value;
            if (form[inputId].value) {
      convertInputValueOnUnitChange(inputId, oldUnit, newUnit);
            }
      currentUnits[inputId] = newUnit;
      updateRangePlaceholder(inputId, newUnit);
            // calculateAndDisplay(); // Already handled by generic select listener
        });
    }
  });

  function copyWithFeedback(btn, text) {
    if (!text || String(text).trim() === "") {
        const origHtml = btn.innerHTML;
        btn.innerHTML = '<i class="fas fa-times me-2"></i>Nothing to copy';
        btn.classList.add('btn-warning'); btn.classList.remove('btn-outline-secondary', 'btn-outline-primary', 'btn-success', 'btn-danger');
        setTimeout(() => {
            btn.innerHTML = origHtml; btn.classList.remove('btn-warning');
            if(btn.id === 'copy-latex-btn-ce') btn.classList.add('btn-outline-primary'); else btn.classList.add('btn-outline-secondary');
        }, 2000);
      return;
    }
    navigator.clipboard.writeText(String(text)).then(() => {
      const origHtml = btn.innerHTML;
      btn.innerHTML = '<i class="fas fa-check me-2"></i>Copied!';
      btn.classList.add('btn-success'); btn.classList.remove('btn-outline-secondary', 'btn-outline-primary', 'btn-danger', 'btn-warning');
      setTimeout(() => {
        btn.innerHTML = origHtml; btn.classList.remove('btn-success');
        if(btn.id === 'copy-latex-btn-ce') btn.classList.add('btn-outline-primary'); else btn.classList.add('btn-outline-secondary');
      }, 1500);
    }).catch(err => {
        console.error('Failed to copy text: ', err);
        const origHtml = btn.innerHTML;
        btn.innerHTML = '<i class="fas fa-times me-2"></i>Copy Failed';
        btn.classList.add('btn-danger'); btn.classList.remove('btn-outline-secondary', 'btn-outline-primary', 'btn-success', 'btn-warning');
        setTimeout(() => {
            btn.innerHTML = origHtml; btn.classList.remove('btn-danger');
            if(btn.id === 'copy-latex-btn-ce') btn.classList.add('btn-outline-primary'); else btn.classList.add('btn-outline-secondary');
        }, 2000);
    });
  }
})();
