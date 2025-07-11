// static/js/constantinou_gani_calculator.js
(function() {
    "use strict"; // Helps catch common coding errors
    console.log("Constantinou-Gani Calculator JS: Script Loaded and Revamped");

    // --- API Endpoints ---
    const API_ENDPOINT_CG_CRITICAL = '/api/calculate_cg_critical';
    const NAME_TO_SMILES_API_ENDPOINT = '/api/name_to_smiles'; // Ensure this endpoint exists and works

    // --- Data ---
    const popularMoleculesCG = [ // Renamed to avoid potential conflicts
        { name: "Methane", smiles: "C", formula: "CH4" }, { name: "Ethane", smiles: "CC", formula: "C2H6" },
        { name: "Propane", smiles: "CCC", formula: "C3H8" }, { name: "Ethanol", smiles: "CCO", formula: "C2H5OH" },
        { name: "Acetic Acid", smiles: "CC(=O)O", formula: "C2H4O2" }, { name: "Benzene", smiles: "c1ccccc1", formula: "C6H6" },
        { name: "Acetone", smiles: "CC(=O)C", formula: "C3H6O" },
        { name: "Phenol", smiles: "Oc1ccccc1", formula: "C6H5OH" },
        { name: "2,3-Dimethylbutane", smiles: "CC(C)C(C)C", formula: "C6H14"},
        { name: "Benzaldehyde", smiles: "O=Cc1ccccc1", formula: "C7H6O"},
        { name: "Perfluoropropane", smiles: "FC(F)(F)C(F)(F)C(F)(F)F", formula: "C3F8"},
        { name: "Pyrrolidinone", smiles: "O=C1CCCN1", formula: "C4H7NO"}
        // Add more relevant examples if needed
    ];

    // --- Unit Conversion Functions ---
    const convertTemperature = (valueK, toUnit) => { if (isNaN(valueK)) return "---"; if (toUnit === 'K') return valueK.toFixed(2); if (toUnit === 'C') return (valueK - 273.15).toFixed(2); if (toUnit === 'F') return ((valueK - 273.15) * 9/5 + 32).toFixed(2); return "---"; };
    const convertPressure = (valueBar, toUnit) => { if (isNaN(valueBar)) return "---"; if (toUnit === 'bar') return valueBar.toFixed(2); if (toUnit === 'atm') return (valueBar * 0.986923).toFixed(2); if (toUnit === 'Pa') return (valueBar * 100000).toExponential(3); if (toUnit === 'kPa') return (valueBar * 100).toFixed(2); if (toUnit === 'psi') return (valueBar * 14.5038).toFixed(2); return "---"; };
    const convertMolarVolume = (valueCm3Mol, toUnit) => { if (isNaN(valueCm3Mol)) return "---"; if (toUnit === 'cm3/mol') return valueCm3Mol.toFixed(1); if (toUnit === 'm3/mol') return (valueCm3Mol * 1e-6).toExponential(3); if (toUnit === 'L/mol') return (valueCm3Mol * 0.001).toFixed(3); return "---"; };

    // --- DOM Element Variables ---
    let formCG, smilesInputCG, calculationOrderSelectCG,
        resultsPlaceholderCG, resultsContentCG, // For new results card structure
        tcElCG, pcElCG, vcElCG, omegaElCG, tbCalcElCG,
        latexStepsElCG, // latexCardCG removed, steps are inside results card now
        paperCitationDisplayElCG, toolCitationDisplayElCG, applicabilityElCG,
        paperCitationSelectCG, paperCitationFormatSelectCG, toolCitationFormatSelectCG,
        copyLatexBtnCG, copyPaperCitationBtnCG, copyToolCitationBtnCG,
        drawMoleculeBtnCG, ketcherModalElementCG, ketcherModalCG, ketcherFrameCG, useStructureBtnCG,
        popularMoleculeSelectCG, formulaInputCG, loadFormulaBtnCG,
        chemicalNameInputCG, getNameSmilesBtnCG,
        outputUnitTcSelectCG, outputUnitPcSelectCG, outputUnitVcSelectCG, outputUnitTbCalcSelectCG,
        warningDisplayElCG; // For the new error display area
        // infoCardCG (accordion container) might be useful if its visibility needs to be controlled, but not strictly for individual items.

    // --- State Variables ---
    let rawLatexStepsCG = '';
    let paperCitationsListCG = []; // Stores the list of paper citations
    let toolCitationDataCG = null; // Stores the single tool citation object
    let debounceTimerCG;
    // Base results stored in standard units (K, bar, cm³/mol)
    let baseTcInK_CG = NaN, basePcInBar_CG = NaN, baseVcInCm3Mol_CG = NaN, baseOmega_CG = NaN, baseTbCalcInK_CG = NaN; // Added Omega and Tb_calc

    // --- Citation Formatting ---
    function formatCitation(data, format = 'Default') {
        // (Using the improved version from the AI suggestion - seems more robust)
        if (!data) return "Citation data not available.";
        let authors = data.authors || data.author || "N.A.";
        let year = data.year || "n.d.";
        let title = data.title || "Untitled Work";
        // Basic escaping to prevent accidental HTML injection if title contains < or >
        title = title.replace(/</g, "<").replace(/>/g, ">");

        if (format === 'APA') { return `${authors} (${year}). *${title}*.${data.edition ? ` (${data.edition}).` : ''} ${data.journal ? ` ${data.journal}, *${data.volume}*(${data.issue}), ${data.pages}.` : (data.publisher ? ` ${data.publisher}.` : '')} ${data.doi ? ` https://doi.org/${data.doi}` : ''} ${data.retrieved_date ? ` Retrieved ${data.retrieved_date}, from ${data.url}`: ''}`; }
        else if (format === 'ACS') { return `${authors}. ${title}${data.edition ? `, ${data.edition}` : ''}; ${data.journal ? `*${data.journal}* ` : ''}${data.publisher ? `${data.publisher}: ` : ''}${year}${data.volume ? `; *${data.volume}* (${data.issue})` : ''}${data.pages ? `, pp ${data.pages}` : ''}.${data.doi ? ` DOI: ${data.doi}` : ''} ${data.retrieved_date ? `Retrieved ${data.retrieved_date}, from ${data.url}`: ''}`; }
        else if (format === 'BibTeX') {
            let type = data.journal ? "@article" : (data.publisher ? "@book" : "@misc");
            if (data.retrieved_date && type !== "@article" && type !== "@book") type = "@misc";
            let bibKeyPart = (authors.split(',')[0].split(' ').pop() || 'Anon').replace(/[\s\W_]+/g, '');
            let titlePart = title.split(' ')[0].replace(/[\s\W_]+/g, '');
            let bibKey = `${bibKeyPart}${year}${titlePart}`;
            // Escape double quotes within fields for BibTeX
            const escapeBib = (str) => str.replace(/"/g, "''");
            let bibtex = `${type}{${bibKey},\n  author    = "${escapeBib(authors)}",\n  title     = "{${escapeBib(title)}}",\n`;
            if(data.publisher && type === "@book") bibtex += `  publisher = "{${escapeBib(data.publisher)}}",\n`;
            if(data.journal && type === "@article") bibtex += `  journal   = "{${escapeBib(data.journal)}}",\n`;
            if(data.volume && type === "@article") bibtex += `  volume    = "${data.volume}",\n`;
            if(data.issue && type === "@article") bibtex += `  number    = "${data.issue}",\n`; // BibTeX uses 'number'
            if(data.pages && type === "@article") bibtex += `  pages     = "${data.pages.replace('--', '-')}",\n`; // Ensure single dash
            bibtex += `  year      = "${year}"`;
            if(data.doi) bibtex += `,\n  doi       = "${escapeBib(data.doi)}"`;
            if(data.url && type === "@misc") bibtex += `,\n  howpublished = "\\url{${data.url}}",\n  note      = "{Retrieved ${data.retrieved_date}}"`; // Added braces for note
            else if(data.url) bibtex += `,\n  url       = "${data.url}"`;
            bibtex += `\n}`; return bibtex;
        }
        // Default format (simplified)
        let citationString = `${authors} (${year}). "${title}".`;
        if (data.journal) citationString += ` *${data.journal}*`;
        if (data.volume) citationString += `, *${data.volume}*`;
        if (data.issue) citationString += `(${data.issue})`;
        if (data.pages) citationString += `, ${data.pages}`;
        if (!data.journal && data.publisher) citationString += `. ${data.publisher}`;
        if (data.doi) citationString += `. DOI: ${data.doi}`;
        if (data.retrieved_date) citationString += `. Retrieved ${data.retrieved_date}, from ${data.url}`;
        return citationString;
    }

    // --- Error Handling (New/Improved) ---
    function displayErrorCG(message) {
        if (warningDisplayElCG) {
            warningDisplayElCG.innerHTML = message; // Use innerHTML for <sub> etc.
            warningDisplayElCG.style.display = 'block';
        }
        if (resultsContentCG) resultsContentCG.style.display = 'none';
        if (resultsPlaceholderCG) resultsPlaceholderCG.style.display = 'block';
        
        if (tcElCG) tcElCG.textContent = '---'; if (pcElCG) pcElCG.textContent = '---';
        if (vcElCG) vcElCG.textContent = '---'; if (omegaElCG) omegaElCG.textContent = '---';
        if (tbCalcElCG) tbCalcElCG.textContent = '---';
        if (latexStepsElCG) latexStepsElCG.innerHTML = `<p class="text-danger">Error: ${message}</p>`;
    }

    function clearErrorCG() {
        if (warningDisplayElCG) {
            warningDisplayElCG.textContent = '';
            warningDisplayElCG.style.display = 'none';
        }
    }

    // --- Display Logic (adapted for new structure) ---
    function displayConvertedResultsCG() {
        if (tcElCG && outputUnitTcSelectCG) tcElCG.textContent = convertTemperature(baseTcInK_CG, outputUnitTcSelectCG.value); else if (tcElCG) tcElCG.textContent = "---";
        if (pcElCG && outputUnitPcSelectCG) pcElCG.textContent = convertPressure(basePcInBar_CG, outputUnitPcSelectCG.value); else if (pcElCG) pcElCG.textContent = "---";
        if (vcElCG && outputUnitVcSelectCG) vcElCG.textContent = convertMolarVolume(baseVcInCm3Mol_CG, outputUnitVcSelectCG.value); else if (vcElCG) vcElCG.textContent = "---";
        if (omegaElCG) omegaElCG.textContent = isNaN(baseOmega_CG) ? "---" : baseOmega_CG.toFixed(4);
        if (tbCalcElCG && outputUnitTbCalcSelectCG) tbCalcElCG.textContent = convertTemperature(baseTbCalcInK_CG, outputUnitTbCalcSelectCG.value); else if (tbCalcElCG) tbCalcElCG.textContent = "---";
    }

    function displaySelectedPaperCitation() {
        if (paperCitationsListCG && paperCitationsListCG.length > 0 && paperCitationSelectCG && paperCitationDisplayElCG && paperCitationFormatSelectCG) {
            const selectedIndex = parseInt(paperCitationSelectCG.value, 10);
            if (!isNaN(selectedIndex) && paperCitationsListCG[selectedIndex]) {
                const formattedCitation = formatCitation(paperCitationsListCG[selectedIndex], paperCitationFormatSelectCG.value);
                if (paperCitationFormatSelectCG.value === 'BibTeX') {
                    paperCitationDisplayElCG.innerHTML = formattedCitation.replace(/\n/g, '<br>');
                } else {
                    paperCitationDisplayElCG.textContent = formattedCitation;
                }
            } else {
                 paperCitationDisplayElCG.textContent = "Selected citation not found.";
            }
        } else if (paperCitationDisplayElCG) {
            paperCitationDisplayElCG.textContent = "Paper citation unavailable.";
        }
    }

    function updateDisplayFromDataCG(apiData) {
        console.log("CG: updateDisplayFromData received:", apiData);
        if (!resultsPlaceholderCG || !resultsContentCG || !tcElCG || !pcElCG || !vcElCG || !omegaElCG || !tbCalcElCG ||
            !latexStepsElCG || !paperCitationDisplayElCG || !toolCitationDisplayElCG || !applicabilityElCG || !paperCitationSelectCG) {
            console.error("CG: Critical display elements missing in updateDisplayFromDataCG.");
            displayErrorCG("A critical display element is missing. Contact support.");
            return;
        }
        clearErrorCG();

        if (apiData.error) {
            console.error("CG: API Error:", apiData.error);
            displayErrorCG(apiData.error);
            rawLatexStepsCG = ''; 
            baseTcInK_CG = NaN; basePcInBar_CG = NaN; baseVcInCm3Mol_CG = NaN; baseOmega_CG = NaN; baseTbCalcInK_CG = NaN;
            
            // Clear paper citation select if error and it was populated based on results
            if (paperCitationSelectCG.options.length > 0 && (!apiData.paper_citations || apiData.paper_citations.length === 0)) {
                 paperCitationSelectCG.innerHTML = '<option value="">N/A</option>';
                 paperCitationSelectCG.disabled = true;
            }
        } else {
            if (resultsPlaceholderCG) resultsPlaceholderCG.style.display = 'none';
            if (resultsContentCG) resultsContentCG.style.display = 'block';

            baseTcInK_CG = (apiData.Tc && apiData.Tc !== "Error") ? parseFloat(apiData.Tc) : NaN;
            basePcInBar_CG = (apiData.Pc && apiData.Pc !== "Error") ? parseFloat(apiData.Pc) : NaN;
            baseVcInCm3Mol_CG = (apiData.Vc && apiData.Vc !== "Error") ? parseFloat(apiData.Vc) : NaN;
            baseOmega_CG = (apiData.omega && apiData.omega !== "Error") ? parseFloat(apiData.omega) : NaN;
            baseTbCalcInK_CG = (apiData.Tb_calc && apiData.Tb_calc !== "Error") ? parseFloat(apiData.Tb_calc) : NaN;

            rawLatexStepsCG = apiData.latex_steps || '';
            if (window.MathJax && rawLatexStepsCG) {
                latexStepsElCG.innerHTML = rawLatexStepsCG;
                MathJax.typesetPromise([latexStepsElCG]).catch(err => {
                    console.error('CG: MathJax typesetting error:', err);
                    latexStepsElCG.innerHTML = `<p class="text-danger">Error rendering MathJax. Raw LaTeX: <pre>${rawLatexStepsCG}</pre></p>`;
                });
            } else {
                latexStepsElCG.textContent = rawLatexStepsCG || "No calculation steps to display.";
            }

            paperCitationsListCG = apiData.paper_citations || paperCitationsListCG;
            toolCitationDataCG = apiData.tool_citation_data || toolCitationDataCG;

            if (paperCitationsListCG && paperCitationsListCG.length > 0) {
                paperCitationSelectCG.disabled = false;
                if (paperCitationSelectCG.options.length !== paperCitationsListCG.length || paperCitationSelectCG.options[0]?.value === "") { // Repopulate if different or was N/A
                     paperCitationSelectCG.innerHTML = ''; 
                     paperCitationsListCG.forEach((citation, index) => {
                        const option = document.createElement('option');
                        option.value = index;
                        option.textContent = citation.short_ref || `Citation ${index + 1}${citation.year ? ` (${citation.year})` : ''}`;
                        paperCitationSelectCG.appendChild(option);
                    });
                }
            } else {
                paperCitationSelectCG.innerHTML = '<option value="">N/A</option>';
                paperCitationSelectCG.disabled = true;
            }
            displaySelectedPaperCitation();

            if (toolCitationDisplayElCG && toolCitationDataCG && toolCitationFormatSelectCG) {
                 const formattedToolCitation = formatCitation(toolCitationDataCG, toolCitationFormatSelectCG.value);
                 if (toolCitationFormatSelectCG.value === 'BibTeX') {
                    toolCitationDisplayElCG.innerHTML = formattedToolCitation.replace(/\n/g, '<br>');
                 } else {
                    toolCitationDisplayElCG.textContent = formattedToolCitation;
                 }
            } else if (toolCitationDisplayElCG) {
                 toolCitationDisplayElCG.textContent = "Tool citation unavailable.";
            }

            if (applicabilityElCG) {
                 applicabilityElCG.textContent = apiData.applicability_notes || "Applicability notes unavailable.";
            }
        }
        displayConvertedResultsCG();
    }
    
    function initializePageUIStateCG() {
        console.log("CG: initializePageUIStateCG called");
        clearErrorCG();
        if (resultsPlaceholderCG) resultsPlaceholderCG.style.display = 'block';
        if (resultsContentCG) resultsContentCG.style.display = 'none';

        baseTcInK_CG = NaN; basePcInBar_CG = NaN; baseVcInCm3Mol_CG = NaN; baseOmega_CG = NaN; baseTbCalcInK_CG = NaN;
        displayConvertedResultsCG();

        if (latexStepsElCG) latexStepsElCG.innerHTML = 'Step-by-step calculation details will appear here.';
        const latexCollapseTarget = document.getElementById('cgLatexStepsContent');
        if (latexCollapseTarget && latexCollapseTarget.classList.contains('show')) {
            const bsCollapse = new bootstrap.Collapse(latexCollapseTarget, { toggle: false });
            bsCollapse.hide();
        }

        ['#cgCollapseCitation', '#cgCollapseToolCitation', '#cgCollapseApplicability'].forEach(selector => {
            const target = document.querySelector(selector);
            if (target && target.classList.contains('show')) {
                const bsCollapse = new bootstrap.Collapse(target, { toggle: false });
                bsCollapse.hide();
            }
        });
        fetchInitialDataCG();
    }

    async function fetchInitialDataCG() {
        console.log("CG: fetchInitialData called");
        try {
            const resp = await fetch(API_ENDPOINT_CG_CRITICAL, { method: 'GET' });
            if (resp.ok) {
                const data = await resp.json();
                console.log("CG: Initial data fetched:", data);
                paperCitationsListCG = data.paper_citations || [];
                toolCitationDataCG = data.tool_citation_data;

                if (paperCitationSelectCG && paperCitationsListCG.length > 0) {
                    paperCitationSelectCG.disabled = false;
                    paperCitationSelectCG.innerHTML = ''; 
                    paperCitationsListCG.forEach((citation, index) => {
                        const option = document.createElement('option');
                        option.value = index;
                        option.textContent = citation.short_ref || `Citation ${index + 1}${citation.year ? ` (${citation.year})` : ''}`;
                        paperCitationSelectCG.appendChild(option);
                    });
                    displaySelectedPaperCitation();
                } else if (paperCitationSelectCG) {
                    paperCitationSelectCG.innerHTML = '<option value="">N/A</option>';
                    paperCitationSelectCG.disabled = true;
                    if(paperCitationDisplayElCG) paperCitationDisplayElCG.textContent = "Paper citation unavailable.";
                }
                
                if (toolCitationDisplayElCG && toolCitationDataCG && toolCitationFormatSelectCG) {
                    const formattedToolCitation = formatCitation(toolCitationDataCG, toolCitationFormatSelectCG.value);
                    if (toolCitationFormatSelectCG.value === 'BibTeX') toolCitationDisplayElCG.innerHTML = formattedToolCitation.replace(/\n/g, '<br>');
                    else toolCitationDisplayElCG.textContent = formattedToolCitation;
                } else if (toolCitationDisplayElCG) {
                    toolCitationDisplayElCG.textContent = "Tool citation details will load with results.";
                }
                if (applicabilityElCG) {
                    applicabilityElCG.textContent = data.applicability_notes || "Applicability notes will load with results.";
                }
            } else {
                const errorText = await resp.text();
                console.error("CG: Failed to fetch initial data:", resp.status, errorText);
                displayErrorCG("Could not load initial page data. Please try refreshing.");
            }
        } catch (error) {
            console.error("CG: Error fetching initial data:", error);
            displayErrorCG("Network error loading initial page data. Please try refreshing.");
        }
    }
    
    function triggerCalculationCG() {
        console.log("CG: triggerCalculationCG CALLED");
        clearErrorCG();
        if (!smilesInputCG || !calculationOrderSelectCG) {
            displayErrorCG("Critical Error: Input fields (SMILES, Order) are missing.");
            return;
        }
        const smiles = smilesInputCG.value.trim();
        const order = calculationOrderSelectCG.value;

        if (!smiles) { displayErrorCG("SMILES string is required."); return; }

        if (resultsPlaceholderCG) resultsPlaceholderCG.style.display = 'none';
        if (resultsContentCG) resultsContentCG.style.display = 'block';
        
        tcElCG.textContent = 'Calculating...'; pcElCG.textContent = 'Calculating...'; vcElCG.textContent = 'Calculating...';
        omegaElCG.textContent = 'Calculating...'; tbCalcElCG.textContent = 'Calculating...';
        if (latexStepsElCG) latexStepsElCG.innerHTML = '<p class="text-info">Calculating steps...</p>';

        fetch(API_ENDPOINT_CG_CRITICAL, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: smiles, order: parseInt(order, 10) })
        })
        .then(async resp => {
            const responseText = await resp.text();
            console.log("CG: Raw API Response:", responseText);
            if (!resp.ok) {
                console.error("CG: API response not OK:", resp.status);
                let errorMsg = `Server error (${resp.status}).`;
                try { const errorData = JSON.parse(responseText); errorMsg = errorData.error || errorData.detail || errorMsg;}
                catch (e) { console.warn("CG: Could not parse error JSON from server:", responseText); }
                updateDisplayFromDataCG({ error: errorMsg });
                return;
            }
            try { const data = JSON.parse(responseText); updateDisplayFromDataCG(data); }
            catch (e) { updateDisplayFromDataCG({ error: "Malformed data from server. Check console." }); console.error("CG: JSON Parse Error:", e, "Raw:", responseText); }
        })
        .catch(error => {
            console.error("Fetch error for CG critical props:", error);
            displayErrorCG("Network error. Check console and connection.");
        });
    }
    
    const triggerCalculationWithDebounceCG = () => { clearTimeout(debounceTimerCG); debounceTimerCG = setTimeout(triggerCalculationCG, 500); };

    function loadSmilesAndCalculateCG(smiles) {
        if (smilesInputCG && smiles) {
            smilesInputCG.value = smiles;
            clearTimeout(debounceTimerCG);
            triggerCalculationCG();
        } else if (smilesInputCG) {
            displayErrorCG("Could not find SMILES for selected molecule/formula.");
        }
    }

    function setupCopyButton(btn, textProviderOrElement) {
        if (!btn) return;
        btn.addEventListener('click', () => {
            let textToCopy = '';
            if (typeof textProviderOrElement === 'function') {
                textToCopy = textProviderOrElement();
            } else if (textProviderOrElement && typeof textProviderOrElement.textContent !== 'undefined') { // Check textContent first
                textToCopy = textProviderOrElement.textContent;
            } else if (textProviderOrElement && typeof textProviderOrElement.innerText !== 'undefined') { // Fallback to innerText
                textToCopy = textProviderOrElement.innerText;
            }
            
            if (textToCopy && typeof textToCopy === 'string' && textToCopy.trim()) {
                navigator.clipboard.writeText(textToCopy.trim()).then(() => {
                    const originalIcon = btn.innerHTML;
                    btn.disabled = true; 
                    btn.innerHTML = '<i class="fas fa-check me-1"></i>Copied!';
                    btn.classList.add('btn-success'); btn.classList.remove('btn-outline-secondary', 'btn-warning', 'btn-danger');
                    setTimeout(() => { 
                        btn.innerHTML = originalIcon;
                        btn.disabled = false; 
                        btn.classList.remove('btn-success'); btn.classList.add('btn-outline-secondary');
                    }, 1500);
                }).catch(err => { 
                    console.error('Copy failed:', err); 
                    const originalIcon = btn.innerHTML;
                    btn.innerHTML = '<i class="fas fa-times me-1"></i>Failed';
                    btn.classList.add('btn-danger'); btn.classList.remove('btn-outline-secondary', 'btn-warning', 'btn-success');
                    setTimeout(() => { 
                        btn.innerHTML = originalIcon;
                        btn.classList.remove('btn-danger'); btn.classList.add('btn-outline-secondary');
                     }, 2000);
                });
            } else { 
                const originalIcon = btn.innerHTML;
                btn.innerHTML = '<i class="fas fa-exclamation-triangle me-1"></i>Empty';
                btn.classList.add('btn-warning'); btn.classList.remove('btn-outline-secondary', 'btn-success', 'btn-danger');
                setTimeout(() => { 
                    btn.innerHTML = originalIcon;
                    btn.classList.remove('btn-warning'); btn.classList.add('btn-outline-secondary');
                 }, 2000);
            }
        });
    }

    document.addEventListener('DOMContentLoaded', () => {
        // Cache DOM elements
        formCG = document.getElementById('cg_form');
        smilesInputCG = document.getElementById('smiles_input_cg');
        calculationOrderSelectCG = document.getElementById('calculationOrderCG');
        resultsPlaceholderCG = document.getElementById('cg_results_placeholder');
        resultsContentCG = document.getElementById('cg_results_content');
        warningDisplayElCG = document.getElementById('cg-warning-display');
        tcElCG = document.getElementById('result_tc_cg');
        pcElCG = document.getElementById('result_pc_cg');
        vcElCG = document.getElementById('result_vc_cg');
        omegaElCG = document.getElementById('result_omega_cg');
        tbCalcElCG = document.getElementById('result_tb_calc_cg');
        latexStepsElCG = document.getElementById('cg_latex_steps');
        paperCitationDisplayElCG = document.getElementById('cg_citation_display');
        toolCitationDisplayElCG = document.getElementById('cg_tool_citation_display');
        applicabilityElCG = document.getElementById('cg_applicability_notes');
        paperCitationSelectCG = document.getElementById('paperCitationSelectCG');
        paperCitationFormatSelectCG = document.getElementById('paper_citation_format_cg');
        toolCitationFormatSelectCG = document.getElementById('tool_citation_format_cg');
        copyLatexBtnCG = document.getElementById('copy_latex_btn_cg');
        copyPaperCitationBtnCG = document.getElementById('copy_citation_btn_cg');
        copyToolCitationBtnCG = document.getElementById('copy_tool_citation_btn_cg');
        drawMoleculeBtnCG = document.getElementById('drawMoleculeBtnCG');
        ketcherModalElementCG = document.getElementById('ketcherModalCG');
        if (ketcherModalElementCG) ketcherModalCG = new bootstrap.Modal(ketcherModalElementCG);
        ketcherFrameCG = document.getElementById('ketcherFrameCG');
        useStructureBtnCG = document.getElementById('useStructureBtnCG');
        popularMoleculeSelectCG = document.getElementById('popularMoleculeSelectCG');
        formulaInputCG = document.getElementById('formulaInputCG');
        loadFormulaBtnCG = document.getElementById('loadFormulaBtnCG');
        chemicalNameInputCG = document.getElementById('chemicalNameInputCG');
        getNameSmilesBtnCG = document.getElementById('getNameSmilesBtnCG');
        outputUnitTcSelectCG = document.getElementById('outputUnit_tc_cg');
        outputUnitPcSelectCG = document.getElementById('outputUnit_pc_cg');
        outputUnitVcSelectCG = document.getElementById('outputUnit_vc_cg');
        outputUnitTbCalcSelectCG = document.getElementById('outputUnit_tb_calc_cg');

        if (!formCG) { console.error("CG: #cg_form NOT FOUND! Aborting script."); return; }

        if (popularMoleculeSelectCG) {
            popularMoleculesCG.forEach(mol => { const option = document.createElement('option'); option.value = mol.smiles; option.textContent = `${mol.name} (${mol.formula || 'N/A'})`; popularMoleculeSelectCG.appendChild(option); });
            popularMoleculeSelectCG.addEventListener('change', (e) => { if (e.target.value) { loadSmilesAndCalculateCG(e.target.value); chemicalNameInputCG.value = ''; formulaInputCG.value = ''; }});
        }

        if (loadFormulaBtnCG && formulaInputCG) {
            loadFormulaBtnCG.addEventListener('click', () => {
                const enteredFormula = formulaInputCG.value.trim().toUpperCase();
                if (!enteredFormula) { displayErrorCG("Please enter a chemical formula."); formulaInputCG.focus(); return; }
                clearErrorCG();
                const foundMolecule = popularMoleculesCG.find(mol => (mol.formula || "").toUpperCase() === enteredFormula);
                if (foundMolecule) { loadSmilesAndCalculateCG(foundMolecule.smiles); formulaInputCG.value = ''; popularMoleculeSelectCG.value=''; chemicalNameInputCG.value = ''; }
                else { displayErrorCG(`Formula "${enteredFormula}" not found in popular list. Please provide SMILES or draw.`); }
            });
        }
        
        if (getNameSmilesBtnCG && chemicalNameInputCG && smilesInputCG) {
            getNameSmilesBtnCG.addEventListener('click', async () => {
                const chemicalName = chemicalNameInputCG.value.trim();
                if (!chemicalName) { displayErrorCG("Please enter a chemical name."); chemicalNameInputCG.focus(); return; }
                clearErrorCG();
                const originalButtonText = getNameSmilesBtnCG.innerHTML;
                getNameSmilesBtnCG.disabled = true; getNameSmilesBtnCG.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Fetching...';
                try {
                    const response = await fetch(NAME_TO_SMILES_API_ENDPOINT, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ chemical_name: chemicalName }) });
                    const data = await response.json();
                    if (response.ok && data.smiles) { loadSmilesAndCalculateCG(data.smiles); chemicalNameInputCG.value = ''; popularMoleculeSelectCG.value=''; formulaInputCG.value = '';}
                    else { displayErrorCG(`Could not get SMILES: ${data.error || 'Unknown server error.'}`); }
                } catch (error) { console.error("Name to SMILES fetch error (CG):", error); displayErrorCG("Error connecting to name resolver service.");
                } finally { getNameSmilesBtnCG.disabled = false; getNameSmilesBtnCG.innerHTML = originalButtonText; }
            });
        }

        if (smilesInputCG) smilesInputCG.addEventListener('input', triggerCalculationWithDebounceCG);
        if (calculationOrderSelectCG) calculationOrderSelectCG.addEventListener('change', triggerCalculationWithDebounceCG);
        
        formCG.addEventListener('submit', (e) => { e.preventDefault(); clearTimeout(debounceTimerCG); triggerCalculationCG(); });
        formCG.addEventListener('reset', () => {
            setTimeout(() => {
                initializePageUIStateCG();
                smilesInputCG.value = '';
                if(calculationOrderSelectCG) calculationOrderSelectCG.value = '1'; // Default to 1st order
                formulaInputCG.value = '';
                chemicalNameInputCG.value = '';
                if(popularMoleculeSelectCG) popularMoleculeSelectCG.value = '';
            }, 0);
        });

        if (paperCitationSelectCG) paperCitationSelectCG.addEventListener('change', displaySelectedPaperCitation);
        if (paperCitationFormatSelectCG) paperCitationFormatSelectCG.addEventListener('change', displaySelectedPaperCitation);
        if (toolCitationFormatSelectCG && toolCitationDisplayElCG) {
            toolCitationFormatSelectCG.addEventListener('change', () => {
                if (toolCitationDataCG) {
                    const formatted = formatCitation(toolCitationDataCG, toolCitationFormatSelectCG.value);
                    if (toolCitationFormatSelectCG.value === 'BibTeX') toolCitationDisplayElCG.innerHTML = formatted.replace(/\n/g, '<br>');
                    else toolCitationDisplayElCG.textContent = formatted;
                }
            });
        }
        
        setupCopyButton(copyLatexBtnCG, () => rawLatexStepsCG);
        setupCopyButton(copyPaperCitationBtnCG, paperCitationDisplayElCG);
        setupCopyButton(copyToolCitationBtnCG, toolCitationDisplayElCG);
        
        if (useStructureBtnCG && ketcherFrameCG && smilesInputCG && ketcherModalCG) {
            useStructureBtnCG.addEventListener('click', async () => {
                try {
                    const ketcherInstance = ketcherFrameCG.contentWindow.ketcher;
                    if (ketcherInstance && typeof ketcherInstance.getSmiles === 'function') {
                        const smilesValue = await ketcherInstance.getSmiles();
                        if (smilesValue && smilesValue.trim() !== '') {
                             loadSmilesAndCalculateCG(smilesValue); 
                             ketcherModalCG.hide(); 
                             if(popularMoleculeSelectCG) popularMoleculeSelectCG.value = '';
                             if(formulaInputCG) formulaInputCG.value = '';
                             if(chemicalNameInputCG) chemicalNameInputCG.value = '';
                        }
                        else { alert("No structure drawn or structure is invalid."); }
                    } else { alert("Ketcher drawer API not available."); }
                } catch (error) { console.error("Ketcher SMILES error (CG):", error); alert("Error using drawn structure."); }
            });
        }

        if (outputUnitTcSelectCG) outputUnitTcSelectCG.addEventListener('change', displayConvertedResultsCG);
        if (outputUnitPcSelectCG) outputUnitPcSelectCG.addEventListener('change', displayConvertedResultsCG);
        if (outputUnitVcSelectCG) outputUnitVcSelectCG.addEventListener('change', displayConvertedResultsCG);
        if (outputUnitTbCalcSelectCG) outputUnitTbCalcSelectCG.addEventListener('change', displayConvertedResultsCG);
        
        initializePageUIStateCG();
    });
})();