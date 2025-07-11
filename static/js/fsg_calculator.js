// fsg_calculator.js: Handles form, unit conversion, AJAX, MathJax, copy buttons for FSG calculator
(function(){
  // Utility: unit conversion functions
  const toInternal = {
    T: (v, u) => {
      if(u === 'K') return v;
      if(u === 'C') return v + 273.15;
      if(u === 'F') return (v - 32) * 5/9 + 273.15;
    },
    P: (v, u) => {
      if(u === 'atm') return v;
      if(u === 'bar') return v / 1.01325;
      if(u === 'Pa') return v / 101325;
      if(u === 'kPa') return v / 101.325;
      if(u === 'psi') return v / 14.6959;
    },
    M: (v, u) => {
      if(u === 'g/mol') return v;
      if(u === 'kg/mol') return v * 1000;
    },
    sigma: (v, u) => {
      if(u === 'A') return v;
      if(u === 'nm') return v * 10;
      if(u === 'm') return v * 1e10;
    },
    epsilon_k: (v, u) => v // Only K for now
  };
  
  // Unit conversion both ways for auto-converting inputs
  const fromInternal = {
    T: (v, u) => {
      if(u === 'K') return v;
      if(u === 'C') return v - 273.15;
      if(u === 'F') return (v - 273.15) * 9/5 + 32;
    },
    P: (v, u) => {
      if(u === 'atm') return v;
      if(u === 'bar') return v * 1.01325;
      if(u === 'Pa') return v * 101325;
      if(u === 'kPa') return v * 101.325;
      if(u === 'psi') return v * 14.6959;
    },
    M: (v, u) => {
      if(u === 'g/mol') return v;
      if(u === 'kg/mol') return v / 1000;
    },
    sigma: (v, u) => {
      if(u === 'A') return v;
      if(u === 'nm') return v / 10;
      if(u === 'm') return v / 1e10;
    },
    epsilon_k: (v, u) => v // Only K for now
  };
  
  // Parameter ranges (in base units)
  const paramRanges = {
    T: { min: 300, max: 2000, unit: 'K' },
    P: { max: 10, unit: 'bar' },
    // Other params have wider ranges
  };

  // DOM
  const form = document.getElementById('fsg-form');
  // Intercept form submit to prevent clearing inputs
  form.addEventListener('submit', function(e) { e.preventDefault(); calculateAndDisplay(); });
  const outputSection = document.getElementById('output-section');
  const resultValue = document.getElementById('result-value');
  const outputUnitSelect = document.getElementById('output_unit');
  const latexDiv = document.getElementById('latex-steps');
  const citationDiv = document.getElementById('citation');
  const toolCitationDiv = document.getElementById('tool-citation');
  const applicabilityDiv = document.getElementById('applicability-notes');
  const copyLatexBtn = document.getElementById('copy-latex-btn');
  const copyCitationBtn = document.getElementById('copy-citation-btn');
  const copyToolCitationBtn = document.getElementById('copy-tool-citation-btn');
  const warningEl = document.getElementById('range-warning');

  let baseValue = null, baseUnit = 'cm^2/s', latexRaw = '';

  // Helper: get all input values (converted)
  function getInputs() {
    return {
      T: toInternal.T(Number(form.T.value), form.T_unit.value),
      P: toInternal.P(Number(form.P.value), form.P_unit.value),
      M_A: toInternal.M(Number(form.M_A.value), form.M_A_unit.value),
      M_B: toInternal.M(Number(form.M_B.value), form.M_B_unit.value),
      sigma_A: toInternal.sigma(Number(form.sigma_A.value), form.sigma_A_unit.value),
      sigma_B: toInternal.sigma(Number(form.sigma_B.value), form.sigma_B_unit.value),
      epsilon_k_A: toInternal.epsilon_k(Number(form.epsilon_k_A.value), form.epsilon_k_A_unit.value),
      epsilon_k_B: toInternal.epsilon_k(Number(form.epsilon_k_B.value), form.epsilon_k_B_unit.value)
    };
  }

  // Helper: get output unit
  function getOutputUnit() {
    return outputUnitSelect.value;
  }

  // AJAX calculation
  async function calculateAndDisplay() {
    const data = getInputs();
    let error = null;
    // Basic front-end validation
    for (const key in data) {
      if (isNaN(data[key]) || data[key] === null) {
        error = 'Please enter valid numeric values for all fields.';
        break;
      }
    }
    if (error) {
      resultValue.textContent = error;
      latexDiv.textContent = '';
      outputSection.style.display = 'block';
      return;
    }
    // Range warnings for T and P inputs
    const warnings = [];
    // Temperature
    const tRaw = Number(form.T.value);
    const tUnit = form.T_unit.value;
    const tMin = fromInternal.T(paramRanges.T.min, tUnit);
    const tMax = fromInternal.T(paramRanges.T.max, tUnit);
    if (tRaw < tMin || tRaw > tMax) {
      warnings.push(`Temperature ${tRaw}${tUnit} is outside valid range ${tMin.toFixed(1)}-${tMax.toFixed(1)}${tUnit}.`);
    }
    // Pressure
    const pRaw = Number(form.P.value);
    const pUnit = form.P_unit.value;
    const pMax = fromInternal.P(paramRanges.P.max, pUnit);
    if (pRaw > pMax) {
      warnings.push(`Pressure ${pRaw}${pUnit} exceeds valid maximum ${pMax.toFixed(1)}${pUnit}.`);
    }
    if (warnings.length) {
      warningEl.innerHTML = warnings.join('<br>');
      warningEl.style.display = 'block';
    } else {
      warningEl.style.display = 'none';
    }
    try {
      const resp = await fetch('/calculate_fsg', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ ...data, output_unit: getOutputUnit() })
      });
      const json = await resp.json();
      if (!resp.ok || json.error) {
        resultValue.textContent = 'Error: ' + (json.error || 'Unknown error');
        latexDiv.textContent = '';
        // Don't clear citation/applicability/tool citation
        outputSection.style.display = 'block';
        return;
      }
      baseValue = Number(json.value);
      baseUnit = json.unit;
      latexRaw = json.latex_steps;
      resultValue.textContent = baseValue.toPrecision(6) + ' ' + baseUnit;
      latexDiv.textContent = latexRaw;
      if (window.MathJax) {
        // Ensure MathJax is fully loaded before typesetting
        if (MathJax.typesetPromise) {
          MathJax.typesetPromise([latexDiv]).catch(err => {
            console.error('MathJax typesetting failed:', err);
          });
        } else if (MathJax.Hub && MathJax.Hub.Queue) {
          // Fallback for MathJax 2.x
          MathJax.Hub.Queue(["Typeset", MathJax.Hub, latexDiv]);
        }
      }
      citationDiv.textContent = json.citation;
      applicabilityDiv.textContent = json.applicability_notes;
      toolCitationDiv.textContent = json.tool_citation;
      outputSection.style.display = 'block';
    } catch (err) {
      resultValue.textContent = 'Error: ' + err.message;
      latexDiv.textContent = '';
      outputSection.style.display = 'block';
    }
  }

  // Helper: convert input value when unit changes
  function convertInputValueOnUnitChange(inputId, oldUnit, newUnit) {
    const inputEl = document.getElementById(inputId);
    if (!inputEl || !inputEl.value || isNaN(inputEl.value)) return;
    
    let value = parseFloat(inputEl.value);
    let conversionType;
    
    if (inputId === 'T') conversionType = 'T';
    else if (inputId === 'P') conversionType = 'P';
    else if (inputId === 'M_A' || inputId === 'M_B') conversionType = 'M';
    else if (inputId === 'sigma_A' || inputId === 'sigma_B') conversionType = 'sigma';
    else if (inputId === 'epsilon_k_A' || inputId === 'epsilon_k_B') conversionType = 'epsilon_k';
    else return;
    
    // Convert to internal first, then to new unit
    const internalValue = toInternal[conversionType](value, oldUnit);
    const newValue = fromInternal[conversionType](internalValue, newUnit);
    
    // Format based on size (avoid scientific notation for readability)
    inputEl.value = Math.abs(newValue) < 0.01 || Math.abs(newValue) > 10000 
      ? newValue.toExponential(4) 
      : newValue.toPrecision(6);
  }
  
  // Update range placeholders
  function updateRangePlaceholder(inputId, unitValue) {
    const inputEl = document.getElementById(inputId);
    if (!inputEl || !paramRanges[inputId]) return;
    
    const range = paramRanges[inputId];
    let conversionType;
    
    if (inputId === 'T') conversionType = 'T';
    else if (inputId === 'P') conversionType = 'P';
    else return; // Only T and P have specific ranges
    
    let placeholder = '';
    
    if (range.min != null && range.max != null) {
      const minConverted = fromInternal[conversionType](range.min, unitValue);
      const maxConverted = fromInternal[conversionType](range.max, unitValue);
      placeholder = `Range: ${minConverted.toFixed(1)} - ${maxConverted.toFixed(1)} ${unitValue}`;
    } else if (range.max) {
      const maxConverted = fromInternal[conversionType](range.max, unitValue);
      placeholder = `Range: <${maxConverted.toFixed(1)} ${unitValue}`;
    }
    
    if (placeholder) {
      inputEl.placeholder = placeholder;
    }
  }
  
  // Track current units for conversion
  const currentUnits = {};
  document.querySelectorAll('select[id$="_unit"]').forEach(select => {
    const inputId = select.id.replace('_unit', '');
    currentUnits[inputId] = select.value;
    
    select.addEventListener('change', () => {
      const oldUnit = currentUnits[inputId];
      const newUnit = select.value;
      convertInputValueOnUnitChange(inputId, oldUnit, newUnit);
      currentUnits[inputId] = newUnit;
      updateRangePlaceholder(inputId, newUnit);
      calculateAndDisplay();
    });
    
    // Initialize range tooltips
    updateRangePlaceholder(inputId, select.value);
  });
  
  // Auto-calc on input change
  form.querySelectorAll('input[type="number"]').forEach(input => {
    input.addEventListener('input', calculateAndDisplay);
  });

  // Output unit change triggers conversion
  outputUnitSelect.addEventListener('change', calculateAndDisplay);

  // Copy buttons with user feedback
  function copyWithFeedback(btn, text) {
    navigator.clipboard.writeText(text).then(() => {
      const orig = btn.textContent;
      btn.textContent = 'Copied!';
      setTimeout(() => { btn.textContent = orig; }, 1200);
    });
  }
  copyLatexBtn.addEventListener('click', () => {
    copyWithFeedback(copyLatexBtn, latexRaw);
  });
  copyCitationBtn.addEventListener('click', () => {
    copyWithFeedback(copyCitationBtn, citationDiv.textContent);
  });
  copyToolCitationBtn.addEventListener('click', () => {
    copyWithFeedback(copyToolCitationBtn, toolCitationDiv.textContent);
  });

  // Dynamic placeholders for units with ranges where applicable
  const placeholderMap = {
    T: {K: 'Range: 300-2000 K', C: 'Range: 26.85-1726.85 °C', F: 'Range: 80.33-3140.33 °F'},
    P: {atm: 'Range: <=10 atm', bar: 'Range: <=10.13 bar', Pa: 'Range: <=1013250 Pa', kPa: 'Range: <=1013.25 kPa', psi: 'Range: <=147 psi'},
    M_A: {'g/mol': 'e.g., 28', 'kg/mol': 'e.g., 0.028'},
    M_B: {'g/mol': 'e.g., 44', 'kg/mol': 'e.g., 0.044'},
    sigma_A: {A: 'e.g., 3.8', nm: 'e.g., 0.38', m: 'e.g., 3.8e-10'},
    sigma_B: {A: 'e.g., 3.9', nm: 'e.g., 0.39', m: 'e.g., 3.9e-10'},
    epsilon_k_A: {K: 'e.g., 71.4'},
    epsilon_k_B: {K: 'e.g., 195.2'}
  };

  [
    ['T', 'T_unit'],
    ['P', 'P_unit'],
    ['M_A', 'M_A_unit'],
    ['M_B', 'M_B_unit'],
    ['sigma_A', 'sigma_A_unit'],
    ['sigma_B', 'sigma_B_unit'],
    ['epsilon_k_A', 'epsilon_k_A_unit'],
    ['epsilon_k_B', 'epsilon_k_B_unit']
  ].forEach(([input, unit]) => {
    const inputEl = document.getElementById(input);
    const unitEl = document.getElementById(unit);
    if (inputEl && unitEl && placeholderMap[input]) {
      unitEl.addEventListener('change', () => {
        inputEl.placeholder = placeholderMap[input][unitEl.value] || '';
      });
      // Set initial placeholder
      inputEl.placeholder = placeholderMap[input][unitEl.value] || '';
    }
  });

  // Trigger calculation when Calculate button is clicked, without clearing inputs
  const calculateBtn = document.getElementById('calculate-btn');
  calculateBtn.addEventListener('click', calculateAndDisplay);
})();
