document.addEventListener('DOMContentLoaded', () => {
    let appState = {};
    const DOMElements = {
        findButton: document.getElementById('find-chemical-btn'),
        nameInput: document.getElementById('chemical-name-input'),
        finderStatus: document.getElementById('finder-status'),
        solverSections: document.getElementById('solver-sections-container'),
        citationsSection: document.getElementById('citations-master-section'),
        citationsList: document.getElementById('citations-list'),
        formatButtons: document.querySelectorAll('.format-btn')
    };

    /**
     * A robust function to render LaTeX. It queues the typesetting command
     * to run only after MathJax is fully initialized, preventing race conditions.
     * @param {HTMLElement} element The container element with new math to render.
     */
    function typesetMath(element) {
        if (window.MathJax) {
            window.MathJax.startup.promise.then(() => {
                window.MathJax.typesetPromise([element]).catch(console.error);
            });
        }
    }

    const citationManager = {
        citations: new Map(),
        add(id, cit) { if (id && cit) this.citations.set(id, cit); this.render(); },
        clear() { this.citations.clear(); },
        render(format) {
            const currentFormat = format || document.querySelector('.format-btn.active')?.dataset.format || 'default';
            DOMElements.citationsSection.style.display = 'block';
            const staticCits = [
                { id: 'toricalc', type: 'Tool', data: { authors: ['ToriCalc'], year: new Date().getFullYear(), title: 'Phase Change Super-Solver' } },
                { id: 'chemicals', type: 'Library', data: { authors: ['Bell, C., Cortes-Pena, Y.R., & Contributors'], year: '2016-2024', title: 'Chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)' } }
            ];
            let html = staticCits.map(cit => this.format(cit, currentFormat)).join('');
            this.citations.forEach(cit => html += this.format(cit, currentFormat));
            DOMElements.citationsList.innerHTML = html;
        },
        format(cit, format) {
            const d = cit.data;
            const authorStr = Array.isArray(d.authors) ? d.authors.join(' and ') : d.author || '';
            const authorApa = Array.isArray(d.authors) ? d.authors.join(', ') : d.author || '';
            const authorAcs = Array.isArray(d.authors) ? d.authors.join('; ') : d.author || '';
            let content = '';
            let plainText = '';

            switch (format) {
                case 'apa':
                    content = `${authorApa} (${d.year}). ${d.title}.`;
                    if (d.journal && d.pages) content += ` <i>${d.journal}</i>, <i>${d.volume}</i>, ${d.pages}.`;
                    break;
                case 'acs':
                    content = `${authorAcs} ${d.title}. <em>${d.journal}</em> <strong>${d.year}</strong>, <em>${d.volume}</em>, ${d.pages}.`;
                    break;
                case 'bibtex':
                    content = `<pre class="bibtex">@article{${cit.id},\n  author  = {${authorStr}},\n  title   = {${d.title}},\n  journal = {${d.journal}},\n  year    = {${d.year}},\n  volume  = {${d.volume}},\n  pages   = {${d.pages}}\n}</pre>`;
                    break;
                default:
                    let title = d.journal ? `"${d.title}"` : `<em>${d.title}</em>`;
                    content = `<strong>[${cit.id}]</strong> ${authorApa} (${d.year}). ${title}`;
                    if (d.journal && d.issue && d.pages) content += ` <em>${d.journal}</em>, ${d.volume}(${d.issue}), pp. ${d.pages}.`;
            }
            
            // Generate plain text version for copying, stripping HTML tags
            plainText = content.replace(/<[^>]+>/g, '');

            // Add the copy button to the HTML content
            return `<div class="citation-item d-flex justify-content-between align-items-start">
                        <div>${content}</div>
                        <button class="btn btn-outline-secondary btn-sm copy-citation-btn" data-clipboard-text="${escape(plainText)}">Copy</button>
                    </div>`;
        }
    };

    const listeners = {
        init() {
            DOMElements.findButton.addEventListener('click', actions.findChemical);
            DOMElements.nameInput.addEventListener('keyup', e => { if (e.key === 'Enter') actions.findChemical(); });
            
            DOMElements.formatButtons.forEach(btn => btn.addEventListener('click', (e) => {
                DOMElements.formatButtons.forEach(b => b.classList.remove('active'));
                e.target.classList.add('active');
                citationManager.render(e.target.dataset.format);
            }));

            // Event delegation for all dynamic elements
            DOMElements.solverSections.addEventListener('change', e => {
                if (e.target.matches('.std-method-select')) actions.handleStandardCalc(e.target.dataset.key);
                if (e.target.matches('.std-unit-select')) actions.updateUnits(e.target.dataset.key);
                if (e.target.id === 'hvap-method') actions.renderHvapInputs();
            });
            DOMElements.solverSections.addEventListener('click', e => {
                if (e.target.id === 'hvap-calculate-btn') actions.handleHvapCalculation();
                if (e.target.id === 'hvap-get-ref-btn') actions.handleGetRefPoint();
                if (e.target.id === 'hvap-copy-latex-btn') navigator.clipboard.writeText(e.target.dataset.latex);
            });
            
            // Event delegation for the new "Copy Citation" buttons
            DOMElements.citationsList.addEventListener('click', e => {
                if (e.target.matches('.copy-citation-btn')) {
                    const textToCopy = unescape(e.target.dataset.clipboardText);
                    navigator.clipboard.writeText(textToCopy);
                    e.target.textContent = 'Copied!';
                    setTimeout(() => { e.target.textContent = 'Copy'; }, 2000);
                }
            });
        }
    };

    const actions = {
        async findChemical() {
            const chemicalName = DOMElements.nameInput.value.trim();
            if (!chemicalName) return;
            DOMElements.finderStatus.innerHTML = '<div class="spinner-border spinner-border-sm" role="status"></div>';
            DOMElements.solverSections.innerHTML = '';
            DOMElements.citationsSection.style.display = 'none';
            citationManager.clear();
            try {
                const response = await fetch('/api/phase_change/find_methods', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ name: chemicalName }) });
                const data = await response.json();
                if (!response.ok) throw new Error(data.error);
                appState = data; appState.results = {};
                DOMElements.finderStatus.textContent = `Found: ${appState.chemical_name} (CAS: ${appState.cas})`;
                actions.renderAllSolvers();
            } catch (error) { DOMElements.finderStatus.textContent = `Error: ${error.message}`; }
        },
        renderAllSolvers() {
            DOMElements.solverSections.innerHTML = `<div id="tb-container"></div><div id="tm-container"></div><div id="hfus-container"></div><div id="hvap-container"></div>`;
            actions.renderStandardCard('tb', 'Boiling Point (T<sub>b</sub>)', appState.available_methods.tb, ['K', '°C', '°F']);
            actions.renderStandardCard('tm', 'Melting Point (T<sub>m</sub>)', appState.available_methods.tm, ['K', '°C', '°F']);
            actions.renderStandardCard('hfus', 'Heat of Fusion (H<sub>fus</sub>)', appState.available_methods.hfus, ['J/mol', 'kJ/mol']);
            actions.renderHvapCard('hvap', 'Heat of Vaporization (H<sub>vap</sub>)', appState.hvap_correlations);
        },
        renderStandardCard(key, title, methods, units) {
            if (!methods || methods.length === 0) return;
            const container = document.getElementById(`${key}-container`);
            const options = methods.map(m => `<option value="${m}">${m}</option>`).join('');
            const unitOptions = units.map(u => `<option value="${u}">${u}</option>`).join('');
            container.innerHTML = `<div class="card solver-card"><div class="card-header solver-card-header">${title}</div><div class="card-body solver-card-body"><div class="row align-items-center"><div class="col-md-5 mb-2 mb-md-0"><label class="form-label">Method:</label><select class="form-select std-method-select" data-key="${key}" id="${key}-method">${options}</select></div><div class="col-md-7"><label class="form-label">Result:</label><div class="d-flex align-items-center"><p class="fs-4 me-3 mb-0" id="${key}-result">-</p><select class="form-select form-select-sm w-auto std-unit-select" data-key="${key}" id="${key}-units">${unitOptions}</select></div></div></div><div class="glass-box" id="${key}-glass-box" style="display:none;"></div></div></div>`;
            actions.handleStandardCalc(key);
        },
        async handleStandardCalc(key) {
            const method = document.getElementById(`${key}-method`).value;
            const resultEl = document.getElementById(`${key}-result`);
            const glassBoxEl = document.getElementById(`${key}-glass-box`);
            glassBoxEl.style.display = 'none';
            resultEl.innerHTML = '<div class="spinner-border spinner-border-sm"></div>';
            try {
                const endpoint = (method === 'JOBACK') ? '/api/phase_change/joback_calculate' : '/api/phase_change/get_value';
                const payload = (method === 'JOBACK') ? { chemical_name: appState.chemical_name, property: key } : { cas: appState.cas, property: key, method: method };
                const response = await fetch(endpoint, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
                const data = await response.json();
                if (!response.ok) throw new Error(data.error);
                appState.results[key] = { value: data.value };
                actions.updateUnits(key);
                glassBoxEl.innerHTML = data.glass_box_html || `Source: ${method}`;
                glassBoxEl.style.display = 'block';
                if (data.citation) citationManager.add(data.citation.id, data.citation);
                typesetMath(glassBoxEl);
            } catch (error) { resultEl.textContent = 'Error'; glassBoxEl.innerHTML = `<p class="text-danger">${error.message}</p>`; glassBoxEl.style.display = 'block'; }
        },
        renderHvapCard(key, title, correlations) {
            if (!correlations || correlations.length === 0) return;
            const container = document.getElementById(`${key}-container`);
            const options = correlations.map(c => `<option value="${c.name}">${c.name}</option>`).join('');
            container.innerHTML = `<div class="card solver-card"><div class="card-header solver-card-header">${title}</div><div class="card-body solver-card-body"><div class="row"><div class="col-md-5 mb-3 mb-md-0"><label for="hvap-method" class="form-label">Correlation:</label><select class="form-select" id="hvap-method">${options}</select></div><div class="col-md-7"><label class="form-label">Required Inputs:</label><div class="hvap-inputs row" id="hvap-inputs-container"></div></div></div><div class="mt-3"><button class="btn btn-success" id="hvap-calculate-btn">Calculate H<sub>vap</sub></button></div><div id="hvap-result-container" class="mt-3"></div><div class="glass-box" id="hvap-glass-box" style="display:none;"></div></div></div>`;
            actions.renderHvapInputs();
        },
        renderHvapInputs() {
            const correlationName = document.getElementById('hvap-method').value;
            const container = document.getElementById('hvap-inputs-container');
            const correlation = appState.hvap_correlations.find(c => c.name === correlationName);
            if (!correlation) { container.innerHTML = ''; return; }
            let inputsHtml = correlation.required_inputs.map(inputKey => {
                const unit = { Tb: 'K', Tc: 'K', Pc_pa: 'Pa', T: 'K', T_ref: 'K', Hvap_ref: 'J/mol' }[inputKey] || '';
                return `<div class="col-6 mb-2"><label for="hvap-in-${inputKey}" class="form-label">${inputKey} ${unit ? `(${unit})` : ''}</label><input type="number" step="any" class="form-control form-control-sm" id="hvap-in-${inputKey}"></div>`;
            }).join('');
            if (correlation.name === 'Watson') {
                inputsHtml += `<div class="col-12 mt-2"><button class="btn btn-sm btn-outline-primary" id="hvap-get-ref-btn">Use T<sub>b</sub> as Reference Point</button></div>`;
            }
            container.innerHTML = inputsHtml;
            actions.autoPopulateHvapInputs();
        },
        autoPopulateHvapInputs() {
            const { prerequisite_data, results } = appState;
            const setValue = (id, value) => {
                const el = document.getElementById(id);
                if (el && value != null) el.value = parseFloat(value.toPrecision(6));
            };
            setValue('hvap-in-Tb', results.tb ? results.tb.value : prerequisite_data.tb);
            setValue('hvap-in-Tc', prerequisite_data.tc);
            setValue('hvap-in-Pc_pa', prerequisite_data.pc_pa);
            setValue('hvap-in-omega', prerequisite_data.omega);
            const tInput = document.getElementById('hvap-in-T');
            if (tInput && !tInput.value) setValue('hvap-in-T', results.tb ? results.tb.value : prerequisite_data.tb);
        },
        async handleGetRefPoint() {
            const tRefInput = document.getElementById('hvap-in-T_ref');
            const hvapRefInput = document.getElementById('hvap-in-Hvap_ref');
            if (!tRefInput || !hvapRefInput) return;
            const currentTb = appState.results.tb ? appState.results.tb.value : appState.prerequisite_data.tb;
            if (!currentTb) { alert("Boiling Point (Tb) must be available first."); return; }
            tRefInput.value = ''; hvapRefInput.value = '';
            tRefInput.placeholder = 'Calculating...'; hvapRefInput.placeholder = 'Calculating...';
            tRefInput.disabled = true; hvapRefInput.disabled = true;
            try {
                const payload = { correlation: 'Riedel', inputs: { Tb: currentTb, Tc: appState.prerequisite_data.tc, Pc_pa: appState.prerequisite_data.pc_pa } };
                const response = await fetch('/api/phase_change/calculate_hvap', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
                const data = await response.json();
                if (!response.ok) throw new Error(data.error);
                tRefInput.value = parseFloat(currentTb.toPrecision(6));
                hvapRefInput.value = parseFloat(data.value.toPrecision(6));
                tRefInput.readOnly = true; hvapRefInput.readOnly = true;
            } catch(error) {
                alert(`Could not fetch reference point: ${error.message}`);
            } finally { tRefInput.placeholder = ''; hvapRefInput.placeholder = ''; tRefInput.disabled = false; hvapRefInput.disabled = false; }
        },
        async handleHvapCalculation() {
          const correlation = document.getElementById('hvap-method').value;
          const resultContainer = document.getElementById('hvap-result-container');
          const glassBoxEl = document.getElementById('hvap-glass-box');
          const button = document.getElementById('hvap-calculate-btn');
          
          resultContainer.innerHTML = '';
          glassBoxEl.style.display = 'block';
          glassBoxEl.innerHTML = '<div class="spinner-border spinner-border-sm"></div> Calculating...';
          button.disabled = true;

          const inputs = {};
          document.querySelectorAll('#hvap-inputs-container input').forEach(el => {
              inputs[el.id.replace('hvap-in-', '')] = el.value;
          });

          try {
              const response = await fetch('/api/phase_change/calculate_hvap', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ correlation, inputs }) });
              const data = await response.json();
              if (!response.ok) throw new Error(data.error);

              // CORRECTED: Use the same HTML structure and classes as the other cards
              resultContainer.innerHTML = `
                  <div class="d-flex align-items-center justify-content-between">
                      <div class="d-flex align-items-center gap-3">
                          <p class="fs-4 mb-0" id="hvap-result">-</p>
                          <select class="form-select form-select-sm w-auto std-unit-select" id="hvap-units" data-key="hvap">
                              <option value="J/mol">J/mol</option>
                              <option value="kJ/mol">kJ/mol</option>
                              <option value="BTU/lb-mol">BTU/lb-mol</option>
                          </select>
                      </div>
                      <button class="btn btn-sm btn-secondary" id="hvap-copy-latex-btn" data-latex="${escape(data.latex_steps)}">Copy LaTeX</button>
                  </div>`;

              // CORRECTED: Save the result to the main state and call the unified update function
              appState.results['hvap'] = { value: data.value };
              actions.updateUnits('hvap');
              
              glassBoxEl.innerHTML = data.plain_steps;
              typesetMath(glassBoxEl);
              if (data.citation) citationManager.add(data.citation.id, data.citation);

          } catch (error) { 
              glassBoxEl.innerHTML = `<p class="text-danger">${error.message}</p>`; 
          }
          finally { button.disabled = false; }
      },
          updateUnits(key) {
            // This now handles all four properties, including Hvap
            if (!appState.results[key]) return;

            const resultEl = document.getElementById(`${key}-result`);
            const unit = document.getElementById(`${key}-units`).value;
            const baseValue = appState.results[key].value;
            if (baseValue === null || isNaN(baseValue)) { resultEl.textContent = '-'; return; }
            
            let displayValue = baseValue;
            if (['tb', 'tm'].includes(key)) {
                if (unit === '°C') displayValue = baseValue - 273.15;
                else if (unit === '°F') displayValue = (baseValue - 273.15) * 9/5 + 32;
            } else if (key === 'hfus' || key === 'hvap') { // Logic for Hfus and Hvap
                if (unit === 'kJ/mol') displayValue = baseValue / 1000;
                else if (unit === 'BTU/lb-mol') displayValue = baseValue * 0.0004303; // More precise factor
            }
            resultEl.textContent = parseFloat(displayValue.toPrecision(6));
        },
        updateHvapDisplay() {
            const resultEl = document.getElementById('hvap-final-result');
            const targetUnit = document.getElementById('hvap-unit-select').value;
            const baseValueJmol = parseFloat(resultEl.dataset.baseValueJmol);
            if (isNaN(baseValueJmol)) return;
            let displayValue = baseValueJmol;
            if (targetUnit === 'kJ/mol') displayValue = baseValueJmol / 1000;
            else if (targetUnit === 'BTU/lb-mol') displayValue = baseValueJmol * 0.4303;
            resultEl.textContent = displayValue.toFixed(2);
        }
    };
    
    listeners.init();
});