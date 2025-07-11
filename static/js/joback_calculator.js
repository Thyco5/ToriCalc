// static/js/joback_calculator.js
(function() {
    console.log("Joback Calculator JS: Script Start");

    // --- Configuration & Constants ---
    const API_ENDPOINT = '/calculate_joback_critical';
    const TB_ESTIMATE_API_ENDPOINT = '/critical-properties/joback/estimate_tb';
    const NAME_TO_SMILES_API_ENDPOINT = '/api/name_to_smiles';
    
    const popularMolecules = [
        { name: "Methane", smiles: "C", formula: "CH4" }, { name: "Ethane", smiles: "CC", formula: "C2H6" },
        { name: "Propane", smiles: "CCC", formula: "C3H8" }, { name: "n-Butane", smiles: "CCCC", formula: "C4H10" },
        { name: "Isobutane (2-Methylpropane)", smiles: "CC(C)C", formula: "C4H10" }, { name: "n-Pentane", smiles: "CCCCC", formula: "C5H12" },
        { name: "Isopentane (2-Methylbutane)", smiles: "CCC(C)C", formula: "C5H12" }, { name: "Neopentane (2,2-Dimethylpropane)", smiles: "CC(C)(C)C", formula: "C5H12" },
        { name: "n-Hexane", smiles: "CCCCCC", formula: "C6H14" }, { name: "n-Heptane", smiles: "CCCCCCC", formula: "C7H16" },
        { name: "n-Octane", smiles: "CCCCCCCC", formula: "C8H18" }, { name: "Ethene (Ethylene)", smiles: "C=C", formula: "C2H4" },
        { name: "Propene (Propylene)", smiles: "CC=C", formula: "C3H6" }, { name: "1-Butene", smiles: "CCC=C", formula: "C4H8" },
        { name: "2-Butene (cis/trans mix)", smiles: "CC=CC", formula: "C4H8" }, { name: "Isobutene (2-Methylpropene)", smiles: "CC(=C)C", formula: "C4H8" },
        { name: "1,3-Butadiene", smiles: "C=CC=C", formula: "C4H6" }, { name: "Ethyne (Acetylene)", smiles: "C#C", formula: "C2H2" },
        { name: "Propyne", smiles: "CC#C", formula: "C3H4" }, { name: "1-Butyne", smiles: "CCC#C", formula: "C4H6" },
        { name: "2-Butyne", smiles: "CC#CC", formula: "C4H6" }, { name: "Cyclopropane", smiles: "C1CC1", formula: "C3H6" },
        { name: "Cyclobutane", smiles: "C1CCC1", formula: "C4H8" }, { name: "Cyclopentane", smiles: "C1CCCC1", formula: "C5H10" },
        { name: "Cyclohexane", smiles: "C1CCCCC1", formula: "C6H12" }, { name: "Methylcyclohexane", smiles: "CC1CCCCC1", formula: "C7H14" },
        { name: "Benzene", smiles: "c1ccccc1", formula: "C6H6" }, { name: "Toluene (Methylbenzene)", smiles: "Cc1ccccc1", formula: "C7H8" },
        { name: "Ethylbenzene", smiles: "CCc1ccccc1", formula: "C8H10" }, { name: "o-Xylene", smiles: "Cc1ccccc1C", formula: "C8H10" },
        { name: "m-Xylene", smiles: "Cc1cccc(C)c1", formula: "C8H10" }, { name: "p-Xylene", smiles: "Cc1ccc(C)cc1", formula: "C8H10" },
        { name: "Styrene (Vinylbenzene)", smiles: "C=Cc1ccccc1", formula: "C8H8" }, { name: "Naphthalene", smiles: "c1ccc2ccccc2c1", formula: "C10H8" },
        { name: "Methanol", smiles: "CO", formula: "CH3OH" }, { name: "Ethanol", smiles: "CCO", formula: "C2H5OH" },
        { name: "1-Propanol (n-Propanol)", smiles: "CCCO", formula: "C3H7OH" }, { name: "2-Propanol (Isopropanol)", smiles: "CC(O)C", formula: "C3H7OH" },
        { name: "1-Butanol (n-Butanol)", smiles: "CCCCO", formula: "C4H9OH" }, { name: "Phenol", smiles: "Oc1ccccc1", formula: "C6H5OH" },
        { name: "Ethylene Glycol", smiles: "OCCO", formula: "C2H6O2" }, { name: "Dimethyl Ether", smiles: "COC", formula: "C2H6O" },
        { name: "Diethyl Ether", smiles: "CCOCC", formula: "C4H10O" }, { name: "Tetrahydrofuran (THF)", smiles: "C1CCOC1", formula: "C4H8O" },
        { name: "Anisole (Methoxybenzene)", smiles: "COc1ccccc1", formula: "C7H8O" }, { name: "Formaldehyde (Methanal)", smiles: "C=O", formula: "CH2O" },
        { name: "Acetaldehyde (Ethanal)", smiles: "CC=O", formula: "C2H4O" }, { name: "Propanal", smiles: "CCC=O", formula: "C3H6O" },
        { name: "Benzaldehyde", smiles: "O=Cc1ccccc1", formula: "C7H6O" }, { name: "Acetone (Propanone)", smiles: "CC(=O)C", formula: "C3H6O" },
        { name: "2-Butanone (Methyl Ethyl Ketone, MEK)", smiles: "CCC(=O)C", formula: "C4H8O" }, { name: "Cyclohexanone", smiles: "O=C1CCCCC1", formula: "C6H10O" },
        { name: "Acetophenone", smiles: "CC(=O)c1ccccc1", formula: "C8H8O" }, { name: "Formic Acid (Methanoic Acid)", smiles: "C(=O)O", formula: "CH2O2" },
        { name: "Acetic Acid (Ethanoic Acid)", smiles: "CC(=O)O", formula: "C2H4O2" }, { name: "Propionic Acid (Propanoic Acid)", smiles: "CCC(=O)O", formula: "C3H6O2" },
        { name: "Benzoic Acid", smiles: "O=C(O)c1ccccc1", formula: "C7H6O2" }, { name: "Methyl Formate", smiles: "COC=O", formula: "C2H4O2" },
        { name: "Methyl Acetate", smiles: "COC(=O)C", formula: "C3H6O2" }, { name: "Ethyl Acetate", smiles: "CCOC(=O)C", formula: "C4H8O2" },
        { name: "Methylamine", smiles: "CN", formula: "CH3NH2" }, { name: "Ethylamine", smiles: "CCN", formula: "C2H5NH2" },
        { name: "Dimethylamine", smiles: "CNC", formula: "(CH3)2NH" }, { name: "Trimethylamine", smiles: "CN(C)C", formula: "(CH3)3N" },
        { name: "Aniline (Phenylamine)", smiles: "Nc1ccccc1", formula: "C6H5NH2" }, { name: "Pyridine", smiles: "n1ccccc1", formula: "C5H5N" },
        { name: "Chloromethane (Methyl Chloride)", smiles: "CCl", formula: "CH3Cl" }, { name: "Dichloromethane (Methylene Chloride)", smiles: "ClCCl", formula: "CH2Cl2" },
        { name: "Chloroform (Trichloromethane)", smiles: "ClC(Cl)Cl", formula: "CHCl3" }, { name: "Carbon Tetrachloride", smiles: "ClC(Cl)(Cl)Cl", formula: "CCl4" },
        { name: "Vinyl Chloride (Chloroethene)", smiles: "C=CCl", formula: "C2H3Cl" }, { name: "Bromoethane (Ethyl Bromide)", smiles: "CCBr", formula: "C2H5Br" },
        { name: "Trifluoromethane (Fluoroform)", smiles: "FC(F)F", formula: "CHF3" }, { name: "Methanethiol (Methyl Mercaptan)", smiles: "CS", formula: "CH3SH" },
        { name: "Ethanethiol (Ethyl Mercaptan)", smiles: "CCS", formula: "CH3CH2SH" }, { name: "Dimethyl Sulfide", smiles: "CSC", formula: "C2H6S" },
        { name: "Thiophene", smiles: "c1ccsc1", formula: "C4H4S" }, { name: "Acetonitrile (Methyl Cyanide)", smiles: "CC#N", formula: "C2H3N" },
        { name: "Nitromethane", smiles: "C[N+](=O)[O-]", formula: "CH3NO2" }, { name: "Water", smiles: "O", formula: "H2O"},
        { name: "Ammonia", smiles: "N", formula: "NH3"}
    ];

    // --- Input Unit Conversion Utilities ---
    const toInternalUnits = {
        T: (v, u) => { if (u === 'K') return parseFloat(v); if (u === 'C' || u === 'degC') return parseFloat(v) + 273.15; if (u === 'F' || u === 'degF') return (parseFloat(v) - 32) * 5 / 9 + 273.15; console.warn(`toInternalUnits.T: Unknown unit ${u}`); return NaN; }
    };
    const fromInternalUnits = {
        T: (v, u) => { if (u === 'K') return v; if (u === 'C' || u === 'degC') return v - 273.15; if (u === 'F' || u === 'degF') return (v - 273.15) * 9 / 5 + 32; console.warn(`fromInternalUnits.T: Unknown unit ${u}`); return NaN; }
    };

    // --- Output Unit Conversion Utilities ---
    const convertTemperature = (valueK, toUnit) => { if (isNaN(valueK)) return "---"; if (toUnit === 'K') return valueK.toFixed(2); if (toUnit === 'C') return (valueK - 273.15).toFixed(2); if (toUnit === 'F') return ((valueK - 273.15) * 9/5 + 32).toFixed(2); return "---"; };
    const convertPressure = (valueBar, toUnit) => { if (isNaN(valueBar)) return "---"; if (toUnit === 'bar') return valueBar.toFixed(2); if (toUnit === 'atm') return (valueBar * 0.986923).toFixed(2); if (toUnit === 'Pa') return (valueBar * 100000).toExponential(3); if (toUnit === 'kPa') return (valueBar * 100).toFixed(2); if (toUnit === 'psi') return (valueBar * 14.5038).toFixed(2); return "---"; };
    const convertMolarVolume = (valueCm3Mol, toUnit) => { if (isNaN(valueCm3Mol)) return "---"; if (toUnit === 'cm3/mol') return valueCm3Mol.toFixed(1); if (toUnit === 'm3/mol') return (valueCm3Mol * 1e-6).toExponential(3); if (toUnit === 'L/mol') return (valueCm3Mol * 0.001).toFixed(3); return "---"; };

    // --- DOM Element References (Updated for new HTML structure) ---
    let form, tbInput, tbUnitSelect, smilesInput, 
        resultsPlaceholderEl, resultsContentEl, // Updated for new results structure
        tcEl, pcEl, vcEl,
        latexStepsEl, paperCitationDisplayEl, toolCitationDisplayEl, applicabilityEl,
        paperCitationFormatSelect, toolCitationFormatSelect,
        copyLatexBtn, copyPaperCitationBtn, copyToolCitationBtn,
        drawMoleculeBtn, ketcherModalElement, ketcherModal, ketcherFrame, useStructureBtn,
        popularMoleculeSelect, formulaInput, loadFormulaBtn,
        estimateTbBtn, jobackTbEstimationCardWrapper, jobackTbEstimationLatexEl, // jobackTbEstimationCard is now jobackTbEstimationCardWrapper
        outputUnitTcSelect, outputUnitPcSelect, outputUnitVcSelect,
        chemicalNameInput, getNameSmilesBtn, warningDisplayEl; // Added warningDisplayEl

    // --- State Variables ---
    let rawLatexSteps = '', paperCitationData = null, toolCitationData = null, currentTbUnit = 'K', debounceTimer;
    let baseTcInK = NaN, basePcInBar = NaN, baseVcInCm3Mol = NaN;

    // --- Helper Functions ---
    function formatCitation(data, format = 'Default') {
        if (!data) return "Citation data not available.";
        let authors = data.authors || data.author || "N.A."; let year = data.year || "n.d."; let title = data.title || "Untitled Work";
        if (format === 'APA') { return `${authors} (${year}). *${title}*.${data.edition ? ` (${data.edition}).` : ''} ${data.publisher ? `${data.publisher}.` : ''} ${data.journal ? `${data.journal}, *${data.volume}*(${data.issue}), ${data.pages}.` : ''} ${data.doi ? `https://doi.org/${data.doi}` : ''} ${data.retrieved_date ? `Retrieved ${data.retrieved_date}, from ${data.url}`: ''}`; }
        else if (format === 'ACS') { return `${authors}. ${title}${data.edition ? `, ${data.edition}` : ''}; ${data.journal ? `*${data.journal}* ` : ''}${data.publisher ? `${data.publisher}:` : ''} ${year}${data.volume ? `, *${data.volume}* (${data.issue})` : ''}${data.pages ? `, ${data.pages}` : ''}.${data.doi ? ` DOI: ${data.doi}` : ''} ${data.retrieved_date ? `Retrieved ${data.retrieved_date}, from ${data.url}`: ''}`; }
        else if (format === 'BibTeX') {
            let type = data.journal ? "@article" : "@book"; if (data.retrieved_date) type = "@misc";
            let bibKeyPart = authors.split(',')[0].replace(/[\s\W]+/g, '') || 'Anon'; let bibKey = `${bibKeyPart}${year}${title.substring(0,10).replace(/[\s\W]+/g, '')}`;
            let bibtex = `${type}{${bibKey},\n  author    = "${authors.replace(/"/g, "''")}",\n  title     = "{${title.replace(/"/g, "''")}}",\n`;
            if(data.publisher && type !== "@misc") bibtex += `  publisher = "{${data.publisher.replace(/"/g, "''")}}",\n`;
            if(data.journal && type === "@article") bibtex += `  journal   = "{${data.journal.replace(/"/g, "''")}}",\n`;
            if(data.volume && type === "@article") bibtex += `  volume    = "${data.volume}",\n`; if(data.issue && type === "@article") bibtex += `  issue     = "${data.issue}",\n`;
            if(data.pages && type === "@article") bibtex += `  pages     = "${data.pages.replace(/"/g, "''")}",\n`; bibtex += `  year      = "${year}"`;
            if(data.doi) bibtex += `,\n  doi       = "${data.doi.replace(/"/g, "''")}"`;
            if(data.url && type === "@misc") bibtex += `,\n  howpublished = "\\url{${data.url}}",\n  note = "Retrieved ${data.retrieved_date}"`;
            else if(data.url) bibtex += `,\n  url       = "${data.url}"`;
            bibtex += `\n}`; return bibtex;
        }
        if (data.retrieved_date) { return `${data.author}. (${data.year}). ${data.title}. Retrieved ${data.retrieved_date}, from ${data.url}`; }
        return `${data.authors} (${data.year}). "${data.title}". ${data.journal ? `${data.journal}, ${data.volume||''}(${data.issue||''}), ${data.pages}.` : `${data.edition ? `${data.edition}, ` : ''}${data.publisher}.`} ${data.doi ? `DOI: ${data.doi}`: ''}`;
    }

    function displayError(message) {
        if (warningDisplayEl) {
            warningDisplayEl.textContent = message;
            warningDisplayEl.style.display = 'block';
        }
        // Hide results content, show placeholder
        if (resultsContentEl) resultsContentEl.style.display = 'none';
        if (resultsPlaceholderEl) resultsPlaceholderEl.style.display = 'block';
        
        // Clear specific result fields if they were showing calculation errors before
        if (tcEl) tcEl.textContent = '---';
        if (pcEl) pcEl.textContent = '---';
        if (vcEl) vcEl.textContent = '---';
        if (latexStepsEl) latexStepsEl.innerHTML = `<p class="text-danger">Error: ${message}</p>`;
    }

    function clearError() {
        if (warningDisplayEl) {
            warningDisplayEl.textContent = '';
            warningDisplayEl.style.display = 'none';
        }
    }

    function displayConvertedResults() {
        if (tcEl && outputUnitTcSelect) { tcEl.textContent = isNaN(baseTcInK) ? "---" : convertTemperature(baseTcInK, outputUnitTcSelect.value); } else if (tcEl) { tcEl.textContent = "---"; }
        if (pcEl && outputUnitPcSelect) { pcEl.textContent = isNaN(basePcInBar) ? "---" : convertPressure(basePcInBar, outputUnitPcSelect.value); } else if (pcEl) { pcEl.textContent = "---"; }
        if (vcEl && outputUnitVcSelect) { vcEl.textContent = isNaN(baseVcInCm3Mol) ? "---" : convertMolarVolume(baseVcInCm3Mol, outputUnitVcSelect.value); } else if (vcEl) { vcEl.textContent = "---"; }
    }

    function updateDisplayFromData(apiData) {
        console.log("Joback: updateDisplayFromData received:", apiData);
        if (!resultsPlaceholderEl || !resultsContentEl || !tcEl || !pcEl || !vcEl || !latexStepsEl || !paperCitationDisplayEl || !toolCitationDisplayEl || !applicabilityEl) {
            console.error("Joback: Critical display elements missing in updateDisplayFromData."); 
            displayError("A critical display element is missing. Contact support.");
            return;
        }
        clearError(); // Clear any previous errors

        if (apiData.error) {
            displayError(apiData.error); // Use new displayError function
            rawLatexSteps = '';
            baseTcInK = NaN; basePcInBar = NaN; baseVcInCm3Mol = NaN;
        } else {
            if (resultsPlaceholderEl) resultsPlaceholderEl.style.display = 'none';
            if (resultsContentEl) resultsContentEl.style.display = 'block';

            baseTcInK = (apiData.Tc && apiData.Tc !== "Error") ? parseFloat(apiData.Tc) : NaN;
            basePcInBar = (apiData.Pc && apiData.Pc !== "Error") ? parseFloat(apiData.Pc) : NaN;
            baseVcInCm3Mol = (apiData.Vc && apiData.Vc !== "Error") ? parseFloat(apiData.Vc) : NaN;
            
            rawLatexSteps = apiData.latex_steps || '';
            if (window.MathJax && rawLatexSteps) {
                latexStepsEl.innerHTML = rawLatexSteps; 
                MathJax.typesetPromise([latexStepsEl]).catch(err => {
                    console.error('Joback: MathJax typesetting error:', err);
                    latexStepsEl.innerHTML = `<p class="text-danger">Error rendering MathJax. Raw LaTeX: <pre>${rawLatexSteps}</pre></p>`;
                });
            } else if (!rawLatexSteps) {
                latexStepsEl.textContent = "No calculation steps to display.";
            } else { 
                latexStepsEl.textContent = rawLatexSteps; // Fallback if MathJax isn't available
            }
        }
        displayConvertedResults(); 

        paperCitationData = apiData.paper_citation_data || paperCitationData;
        toolCitationData = apiData.tool_citation_data || toolCitationData;
        if (paperCitationDisplayEl && paperCitationData) { paperCitationDisplayEl.textContent = formatCitation(paperCitationData, paperCitationFormatSelect.value); }
        else if (paperCitationDisplayEl) { paperCitationDisplayEl.textContent = "Paper citation data unavailable."; }
        if (toolCitationDisplayEl && toolCitationData) { toolCitationDisplayEl.textContent = formatCitation(toolCitationData, toolCitationFormatSelect.value); }
        else if (toolCitationDisplayEl) { toolCitationDisplayEl.textContent = "Tool citation data unavailable."; }
        if (applicabilityEl && apiData.applicability_notes) { applicabilityEl.textContent = apiData.applicability_notes; }
        else if (applicabilityEl && !apiData.error) { applicabilityEl.textContent = "Applicability notes unavailable."; }
    }
    
    function initializePageUIState() {
        console.log("Joback: initializePageUIState called");
        clearError();
        if (resultsPlaceholderEl) resultsPlaceholderEl.style.display = 'block';
        if (resultsContentEl) resultsContentEl.style.display = 'none';
        
        baseTcInK = NaN; basePcInBar = NaN; baseVcInCm3Mol = NaN;
        displayConvertedResults();

        if (latexStepsEl) latexStepsEl.innerHTML = 'Step-by-step calculation details will appear here.';
        // Collapse main calculation steps if open
        const latexCollapseTarget = document.getElementById('jobackLatexStepsContent');
        if (latexCollapseTarget && latexCollapseTarget.classList.contains('show')) {
            const bsCollapse = new bootstrap.Collapse(latexCollapseTarget, { toggle: false });
            bsCollapse.hide();
        }

        if (jobackTbEstimationCardWrapper) jobackTbEstimationCardWrapper.style.display = 'none';
        if (jobackTbEstimationLatexEl) jobackTbEstimationLatexEl.innerHTML = '';
        
        // Collapse accordion items
        ['#jobackCollapseCitation', '#jobackCollapseToolCitation', '#jobackCollapseApplicability'].forEach(selector => {
            const target = document.querySelector(selector);
            if (target && target.classList.contains('show')) {
                const bsCollapse = new bootstrap.Collapse(target, { toggle: false });
                bsCollapse.hide();
            }
        });
        fetchInitialData(); // Fetch default citations and applicability notes
    }

    async function fetchInitialData() {
        console.log("Joback: fetchInitialData for critical props");
        try {
            const resp = await fetch(API_ENDPOINT, { method: 'GET' }); 
            if (resp.ok) {
                const data = await resp.json();
                paperCitationData = data.paper_citation_data || null; toolCitationData = data.tool_citation_data || null;
                if (paperCitationDisplayEl && paperCitationData) paperCitationDisplayEl.textContent = formatCitation(paperCitationData, paperCitationFormatSelect.value);
                else if (paperCitationDisplayEl) paperCitationDisplayEl.textContent = 'Paper citation details will load with results.';
                if (toolCitationDisplayEl && toolCitationData) toolCitationDisplayEl.textContent = formatCitation(toolCitationData, toolCitationFormatSelect.value);
                else if (toolCitationDisplayEl) toolCitationDisplayEl.textContent = 'Tool citation details will load with results.';
                if (applicabilityEl && data.applicability_notes) applicabilityEl.textContent = data.applicability_notes;
                else if(applicabilityEl) applicabilityEl.textContent = 'Applicability notes will load with results.';
            } else {
                const errorText = await resp.text(); 
                console.error("Joback: Failed to fetch initial critical props data:", resp.status, errorText);
                displayError("Could not load initial page data (citations/notes). Please try refreshing.");
            }
        } catch (error) {
            console.error("Joback: Error fetching initial critical props data:", error);
            displayError("Network error loading initial page data. Please try refreshing.");
        }
    }
    
    function triggerCalculation() {
        console.log("Joback: triggerCalculation (critical props) CALLED");
        clearError();
        if (!tbInput || !smilesInput || !tbUnitSelect) { displayError("Critical Error: Input fields are missing from the page."); return; }
        const tbValueStr = tbInput.value.trim(); const currentSelectedTbUnit = tbUnitSelect.value; const smiles = smilesInput.value.trim();

        if (!smiles && !tbValueStr) {
            // Neither SMILES nor Tb is entered; reset to placeholder state without error
            baseTcInK = NaN; basePcInBar = NaN; baseVcInCm3Mol = NaN; displayConvertedResults(); 
            if (latexStepsEl) latexStepsEl.innerHTML = 'Step-by-step calculation details will appear here.';
            if (jobackTbEstimationCardWrapper) jobackTbEstimationCardWrapper.style.display = 'none';
            if (resultsPlaceholderEl) resultsPlaceholderEl.style.display = 'block';
            if (resultsContentEl) resultsContentEl.style.display = 'none';
            return;
        }
        if (!tbValueStr && smiles) { displayError("Normal Boiling Point (T<sub>b</sub>) is required."); return; }
        if (tbValueStr && !smiles) { displayError("SMILES string is required."); return; }

        const tbRaw = parseFloat(tbValueStr);
        if (isNaN(tbRaw) || tbRaw <= 0) { displayError("T<sub>b</sub> must be a positive number."); return; }
        const tbInKelvin = toInternalUnits.T(tbRaw, currentSelectedTbUnit);
        if (isNaN(tbInKelvin)) { displayError("Invalid T<sub>b</sub> unit or value for conversion."); return; }

        if(resultsPlaceholderEl) resultsPlaceholderEl.style.display = 'none';
        if(resultsContentEl) resultsContentEl.style.display = 'block'; // Show results area, Tc etc will show calculating

        tcEl.textContent = 'Calculating...'; pcEl.textContent = 'Calculating...'; vcEl.textContent = 'Calculating...';
        if(latexStepsEl) latexStepsEl.innerHTML = '<p class="text-info">Calculating steps...</p>';

        fetch(API_ENDPOINT, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ Tb: tbInKelvin, smiles: smiles }) })
        .then(async resp => {
            const responseText = await resp.text(); // Always get text first for better error inspection
            console.log("Joback: Raw API Response:", responseText);
            if (!resp.ok) {
                console.error("Joback: API response not OK:", resp.status);
                let errorMsg = `Server error (${resp.status}).`;
                try { 
                    const errorData = JSON.parse(responseText); 
                    errorMsg = errorData.error || errorData.detail || errorMsg; // FastAPI uses detail
                } catch (e) { 
                    console.warn("Could not parse error JSON from server:", responseText); 
                }
                updateDisplayFromData({ error: errorMsg }); // updateDisplayFromData handles tcEl, etc.
                return;
            }
            try { 
                const data = JSON.parse(responseText); 
                updateDisplayFromData(data); 
            } catch (e) { 
                updateDisplayFromData({ error: "Malformed data from server. Check console." }); 
                console.error("JSON Parse Error:", e, "Raw:", responseText); 
            }
        })
        .catch(error => { 
            console.error("Fetch error for critical props:", error); 
            displayError("Network error. Check console and connection."); 
        });
    }
    
    function handleTbUnitChange() {
        if (!tbInput || !tbUnitSelect) return;
        const tbValueStr = tbInput.value.trim();
        if (tbValueStr && !isNaN(parseFloat(tbValueStr))) {
            const tbOldValue = parseFloat(tbValueStr); const tbInKelvin = toInternalUnits.T(tbOldValue, currentTbUnit); 
            const tbNewValue = fromInternalUnits.T(tbInKelvin, tbUnitSelect.value);
            if (!isNaN(tbNewValue)) { tbInput.value = (Math.abs(tbNewValue) < 0.01 && tbNewValue !== 0) || Math.abs(tbNewValue) > 100000 ? tbNewValue.toExponential(3) : tbNewValue.toFixed(2); }
        }
        currentTbUnit = tbUnitSelect.value; 
        // No automatic recalculation on Tb unit change alone if Tb value is present. Let user initiate.
        // triggerCalculationWithDebounce(); 
    }
    
    function setupCopyButton(btn, textProviderOrElement) {
        if (!btn) return;
        btn.addEventListener('click', () => {
            let textToCopy = '';
            if (typeof textProviderOrElement === 'function') {
                textToCopy = textProviderOrElement();
            } else if (textProviderOrElement && (textProviderOrElement.textContent || textProviderOrElement.innerText)) {
                textToCopy = textProviderOrElement.textContent || textProviderOrElement.innerText;
            }
            
            if (textToCopy && typeof textToCopy === 'string' && textToCopy.trim()) {
                navigator.clipboard.writeText(textToCopy.trim()).then(() => {
                    const originalIcon = btn.innerHTML; // Assumes button might have an icon
                    btn.disabled = true; 
                    btn.innerHTML = '<i class="fas fa-check me-1"></i>Copied!';
                    btn.classList.add('btn-success'); btn.classList.remove('btn-outline-secondary');
                    setTimeout(() => { 
                        btn.innerHTML = originalIcon;
                        btn.disabled = false; 
                        btn.classList.remove('btn-success'); btn.classList.add('btn-outline-secondary');
                    }, 1500);
                }).catch(err => { 
                    console.error('Copy failed:', err); 
                    const originalIcon = btn.innerHTML;
                    btn.innerHTML = '<i class="fas fa-times me-1"></i>Failed';
                    btn.classList.add('btn-danger'); btn.classList.remove('btn-outline-secondary');
                    setTimeout(() => { 
                        btn.innerHTML = originalIcon;
                        btn.classList.remove('btn-danger'); btn.classList.add('btn-outline-secondary');
                     }, 2000);
                });
            } else { 
                // alert("Nothing to copy."); 
                const originalIcon = btn.innerHTML;
                btn.innerHTML = '<i class="fas fa-exclamation-triangle me-1"></i>Empty';
                btn.classList.add('btn-warning'); btn.classList.remove('btn-outline-secondary');
                setTimeout(() => { 
                    btn.innerHTML = originalIcon;
                    btn.classList.remove('btn-warning'); btn.classList.add('btn-outline-secondary');
                 }, 2000);
            }
        });
    }

    function loadSmilesAndCalculate(smiles) {
        if (smilesInput && smiles) {
            smilesInput.value = smiles;
            if (jobackTbEstimationCardWrapper) jobackTbEstimationCardWrapper.style.display = 'none';
            if (jobackTbEstimationLatexEl) jobackTbEstimationLatexEl.innerHTML = '';
            clearTimeout(debounceTimer); 
            triggerCalculation(); // Calculate immediately after loading SMILES if Tb is present
        } else if (smilesInput) { displayError("Could not find SMILES for selected molecule/formula."); }
    }
    
    const triggerCalculationWithDebounce = () => { clearTimeout(debounceTimer); debounceTimer = setTimeout(triggerCalculation, 500); };

    document.addEventListener('DOMContentLoaded', () => {
        // Cache DOM elements
        form = document.getElementById('joback_form');
        tbInput = document.getElementById('Tb_joback');
        tbUnitSelect = document.getElementById('Tb_joback_unit');
        smilesInput = document.getElementById('smiles_input_joback');
        
        resultsPlaceholderEl = document.getElementById('joback_results_placeholder');
        resultsContentEl = document.getElementById('joback_results_content');
        warningDisplayEl = document.getElementById('joback-warning-display');

        tcEl = document.getElementById('result_tc_joback');
        pcEl = document.getElementById('result_pc_joback');
        vcEl = document.getElementById('result_vc_joback');
        latexStepsEl = document.getElementById('joback_latex_steps');
        paperCitationDisplayEl = document.getElementById('joback_citation_display');
        toolCitationDisplayEl = document.getElementById('joback_tool_citation_display');
        applicabilityEl = document.getElementById('joback_applicability_notes');
        paperCitationFormatSelect = document.getElementById('paper_citation_format_joback');
        toolCitationFormatSelect = document.getElementById('tool_citation_format_joback');
        copyLatexBtn = document.getElementById('copy_latex_btn_joback');
        copyPaperCitationBtn = document.getElementById('copy_citation_btn_joback');
        copyToolCitationBtn = document.getElementById('copy_tool_citation_btn_joback');
        drawMoleculeBtn = document.getElementById('drawMoleculeBtn');
        ketcherModalElement = document.getElementById('ketcherModal');
        if (ketcherModalElement) ketcherModal = new bootstrap.Modal(ketcherModalElement);
        ketcherFrame = document.getElementById('ketcherFrame');
        useStructureBtn = document.getElementById('useStructureBtn');
        popularMoleculeSelect = document.getElementById('popularMoleculeSelect');
        formulaInput = document.getElementById('formulaInput');
        loadFormulaBtn = document.getElementById('loadFormulaBtn');
        estimateTbBtn = document.getElementById('estimateTbBtnJoback');
        jobackTbEstimationCardWrapper = document.getElementById('jobackTbEstimationCardWrapper'); 
        jobackTbEstimationLatexEl = document.getElementById('joback_tb_estimation_latex_steps');
        outputUnitTcSelect = document.getElementById('outputUnit_tc_joback');
        outputUnitPcSelect = document.getElementById('outputUnit_pc_joback');
        outputUnitVcSelect = document.getElementById('outputUnit_vc_joback');
        chemicalNameInput = document.getElementById('chemicalNameInput');
        getNameSmilesBtn = document.getElementById('getNameSmilesBtn');

        if (!form) { console.error("Joback: #joback_form NOT FOUND! Aborting script."); return; }

        if (popularMoleculeSelect) {
            popularMolecules.forEach(mol => { const option = document.createElement('option'); option.value = mol.smiles; option.textContent = `${mol.name} (${mol.formula || 'N/A'})`; popularMoleculeSelect.appendChild(option); });
            popularMoleculeSelect.addEventListener('change', (e) => { if (e.target.value) { loadSmilesAndCalculate(e.target.value); chemicalNameInput.value = ''; formulaInput.value = ''; }});
        }

        if (loadFormulaBtn && formulaInput) {
            loadFormulaBtn.addEventListener('click', () => {
                const enteredFormula = formulaInput.value.trim().toUpperCase();
                if (!enteredFormula) { displayError("Please enter a chemical formula."); formulaInput.focus(); return; }
                clearError();
                const foundMolecule = popularMolecules.find(mol => (mol.formula || "").toUpperCase() === enteredFormula);
                if (foundMolecule) { loadSmilesAndCalculate(foundMolecule.smiles); formulaInput.value = ''; popularMoleculeSelect.value=''; chemicalNameInput.value = ''; }
                else { displayError(`Formula "${enteredFormula}" not found in popular list. Please provide SMILES or draw.`); }
            });
        }
        
        if (getNameSmilesBtn && chemicalNameInput && smilesInput) {
            getNameSmilesBtn.addEventListener('click', async () => {
                const chemicalName = chemicalNameInput.value.trim();
                if (!chemicalName) { displayError("Please enter a chemical name."); chemicalNameInput.focus(); return; }
                clearError();
                const originalButtonText = getNameSmilesBtn.innerHTML; // Store full HTML including icon
                getNameSmilesBtn.disabled = true; getNameSmilesBtn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Fetching...';
                try {
                    const response = await fetch(NAME_TO_SMILES_API_ENDPOINT, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ chemical_name: chemicalName }) });
                    const data = await response.json();
                    if (response.ok && data.smiles) { loadSmilesAndCalculate(data.smiles); chemicalNameInput.value = ''; popularMoleculeSelect.value=''; formulaInput.value = '';} 
                    else { displayError(`Could not get SMILES: ${data.error || 'Unknown server error.'}`); }
                } catch (error) { console.error("Name to SMILES fetch error:", error); displayError("Error connecting to name resolver service.");
                } finally { getNameSmilesBtn.disabled = false; getNameSmilesBtn.innerHTML = originalButtonText; }
            });
        }

        if (tbUnitSelect) { currentTbUnit = tbUnitSelect.value; tbUnitSelect.addEventListener('change', handleTbUnitChange); }
        if (tbInput) tbInput.addEventListener('input', triggerCalculationWithDebounce);
        if (smilesInput) smilesInput.addEventListener('input', () => { 
             // Clear Tb estimation card if user types SMILES manually
            if (jobackTbEstimationCardWrapper) jobackTbEstimationCardWrapper.style.display = 'none';
            if (jobackTbEstimationLatexEl) jobackTbEstimationLatexEl.innerHTML = '';
            triggerCalculationWithDebounce(); 
        });
        form.addEventListener('submit', (e) => { e.preventDefault(); clearTimeout(debounceTimer); triggerCalculation(); });
        form.addEventListener('reset', () => {
            setTimeout(() => { // Allow browser reset to clear fields first
                initializePageUIState();
                smilesInput.value = '';
                tbInput.value = '';
                tbUnitSelect.value = 'K'; currentTbUnit = 'K';
                formulaInput.value = '';
                chemicalNameInput.value = '';
                if(popularMoleculeSelect) popularMoleculeSelect.value = '';
                // Any other specific fields to clear for Joback
            }, 0);
        });

        if (paperCitationFormatSelect && paperCitationDisplayEl) paperCitationFormatSelect.addEventListener('change', () => paperCitationData && (paperCitationDisplayEl.textContent = formatCitation(paperCitationData, paperCitationFormatSelect.value)));
        if (toolCitationFormatSelect && toolCitationDisplayEl) toolCitationFormatSelect.addEventListener('change', () => toolCitationData && (toolCitationDisplayEl.textContent = formatCitation(toolCitationData, toolCitationFormatSelect.value)));
        
        setupCopyButton(copyLatexBtn, () => rawLatexSteps);
        setupCopyButton(copyPaperCitationBtn, paperCitationDisplayEl); // Pass element directly
        setupCopyButton(copyToolCitationBtn, toolCitationDisplayEl); // Pass element directly
        
        if (useStructureBtn && ketcherFrame && smilesInput && ketcherModal) {
            useStructureBtn.addEventListener('click', async () => {
                try {
                    const ketcherInstance = ketcherFrame.contentWindow.ketcher;
                    if (ketcherInstance && typeof ketcherInstance.getSmiles === 'function') {
                        const smilesValue = await ketcherInstance.getSmiles(); // Ketcher API might return null for empty canvas
                        if (smilesValue && smilesValue.trim() !== '') {
                             loadSmilesAndCalculate(smilesValue); 
                             ketcherModal.hide(); 
                             // Clear other molecule inputs
                             if(popularMoleculeSelect) popularMoleculeSelect.value = '';
                             if(formulaInput) formulaInput.value = '';
                             if(chemicalNameInput) chemicalNameInput.value = '';
                        }
                        else { alert("No structure drawn or structure is invalid."); }
                    } else { alert("Ketcher drawer API not available."); }
                } catch (error) { console.error("Ketcher SMILES error:", error); alert("Error using drawn structure."); }
            });
        }
        
        if (estimateTbBtn && smilesInput && tbInput && tbUnitSelect) {
            estimateTbBtn.addEventListener('click', async () => {
                const currentSmiles = smilesInput.value.trim();
                if (!currentSmiles) { displayError("Enter/draw SMILES to estimate Tb."); smilesInput.focus(); return; }
                clearError();
                const originalTbPlaceholder = tbInput.placeholder;
                const originalButtonHtml = estimateTbBtn.innerHTML; // Store full HTML including icon
                tbInput.value = ''; tbInput.placeholder = 'Estimating Tb...'; 
                estimateTbBtn.disabled = true; estimateTbBtn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Estimating...';
                
                if(jobackTbEstimationCardWrapper) jobackTbEstimationCardWrapper.style.display = 'none'; 
                if(jobackTbEstimationLatexEl) jobackTbEstimationLatexEl.innerHTML = '';
        
                try {
                    const response = await fetch(TB_ESTIMATE_API_ENDPOINT, { 
                        method: 'POST', headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ smiles: currentSmiles })
                    });
                    const data = await response.json();
                    if (response.ok && data.estimated_tb_k !== undefined && data.estimated_tb_k !== null) {
                        tbInput.value = data.estimated_tb_k.toFixed(2);
                        tbUnitSelect.value = 'K'; currentTbUnit = 'K'; 
                        if (data.tb_estimation_latex_steps && jobackTbEstimationLatexEl && jobackTbEstimationCardWrapper) {
                            jobackTbEstimationLatexEl.innerHTML = data.tb_estimation_latex_steps;
                            jobackTbEstimationCardWrapper.style.display = 'block'; 
                            if (window.MathJax && window.MathJax.typesetPromise) { MathJax.typesetPromise([jobackTbEstimationLatexEl]).catch(err => console.error('MathJax Tb steps error:', err)); }
                            else if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Queue) { window.MathJax.Hub.Queue(["Typeset",MathJax.Hub,jobackTbEstimationLatexEl]);}
                        }
                        clearTimeout(debounceTimer); triggerCalculation(); 
                    } else {
                        displayError(`Error estimating Tb: ${data.error || 'Estimation failed or not applicable.'}`);
                        tbInput.placeholder = originalTbPlaceholder;
                    }
                } catch (error) { 
                    console.error("Network/JSON parse error estimating Tb:", error);
                    displayError("Could not connect or parse server response for Tb estimation.");
                    tbInput.placeholder = originalTbPlaceholder;
                } finally {
                    estimateTbBtn.disabled = false; estimateTbBtn.innerHTML = originalButtonHtml;
                    if (tbInput.placeholder === 'Estimating Tb...') tbInput.placeholder = originalTbPlaceholder;
                }
            });
        } else { if(!estimateTbBtn) console.warn("Joback: Estimate Tb button not found."); }

        if (outputUnitTcSelect) outputUnitTcSelect.addEventListener('change', displayConvertedResults);
        if (outputUnitPcSelect) outputUnitPcSelect.addEventListener('change', displayConvertedResults);
        if (outputUnitVcSelect) outputUnitVcSelect.addEventListener('change', displayConvertedResults);
        
        initializePageUIState(); // Set initial state of the page (placeholders, hidden sections, etc.)
    });
})();