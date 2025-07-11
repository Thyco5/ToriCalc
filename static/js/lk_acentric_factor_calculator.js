// static/js/lk_acentric_factor_calculator.js

document.addEventListener('DOMContentLoaded', function() {
    console.log("LK Acentric Factor Calculator JS loaded.");

    // --- DOM Element References ---
    const form = document.getElementById('lk-omega-form');
    const boilingPointInput = document.getElementById('boiling-point');
    const boilingPointUnitsSelect = document.getElementById('boiling-point-units');
    const criticalTempInput = document.getElementById('critical-temp');
    const criticalTempUnitsSelect = document.getElementById('critical-temp-units');
    const criticalPressureInput = document.getElementById('critical-pressure');
    const criticalPressureUnitsSelect = document.getElementById('critical-pressure-units');
    
    // --- Output and Error Elements ---
    const resultContainer = document.getElementById('result-container');
    const placeholderContainer = document.getElementById('placeholder-container');
    const omegaResultElem = document.getElementById('omega-result');
    const resultSourceElem = document.getElementById('result-source-text');
    const errorContainer = document.getElementById('lk-error-container');
    
    const stepsCard = document.getElementById('calculation-steps-card');
    const stepsContainer = document.getElementById('steps-container');
    
    const citationsContentElem = document.getElementById('citations-content');
    const copyLatexBtn = document.getElementById('copy-latex-btn');
    
    // --- CITATION FORMATTING LOGIC ---
    const toolCitationContainer = document.getElementById('tool-citation-container');
    const sourceCitationsContainer = document.getElementById('source-citations-container');
    const citationControlButtons = document.querySelectorAll('.citation-controls .format-btn');

    let currentCitations = []; // To store the raw citation data from the backend
    let rawLatexString = ''; // Variable to store the raw LaTeX from the API

    // --- Unit Conversion Functions ---
    // Converts a temperature from a given unit TO Kelvin
    function convertTemperatureToKelvin(value, unit) {
        if (unit === 'C') return value + 273.15;
        if (unit === 'F') return (value - 32) * 5/9 + 273.15;
        return value; // Already in Kelvin
    }

    // Converts a temperature FROM Kelvin TO a given unit
    function convertTemperatureFromKelvin(value, unit) {
        if (unit === 'C') return value - 273.15;
        if (unit === 'F') return (value - 273.15) * 9/5 + 32;
        return value; // Kelvin
    }

    // Converts a pressure from a given unit TO Pascals
    function convertPressureToPascals(value, unit) {
        if (unit === 'bar') return value * 1e5;
        if (unit === 'kPa') return value * 1000;
        if (unit === 'atm') return value * 101325;
        if (unit === 'psi') return value * 6894.76;
        return value; // Already in Pa
    }

    // Converts a pressure FROM Pascals TO a given unit
    function convertPressureFromPascals(value, unit) {
        if (unit === 'bar') return value / 1e5;
        if (unit === 'kPa') return value / 1000;
        if (unit === 'atm') return value / 101325;
        if (unit === 'psi') return value / 6894.76;
        return value; // Pascals
    }

    // --- Utility Functions ---
    function showError(message) {
        errorContainer.textContent = message;
        errorContainer.style.display = 'block';
    }

    function hideError() {
        errorContainer.style.display = 'none';
    }

    function copyToClipboard(text) {
        navigator.clipboard.writeText(text).then(() => {
            const originalText = copyLatexBtn.innerHTML;
            copyLatexBtn.innerHTML = '<i class="fas fa-check me-2"></i>Copied!';
            setTimeout(() => { copyLatexBtn.innerHTML = originalText; }, 2000);
        }, (err) => {
            alert('Failed to copy LaTeX.');
            console.error('Clipboard copy failed: ', err);
        });
    }

    // --- Dynamic Unit Conversion for Inputs ---
    function createUnitConverter(inputElem, unitSelectElem, toBaseConverter, fromBaseConverter) {
        let previousUnit = unitSelectElem.value;
        unitSelectElem.addEventListener('change', () => {
            const currentValue = parseFloat(inputElem.value);
            if (!isNaN(currentValue)) {
                const baseValue = toBaseConverter(currentValue, previousUnit);
                const newValue = fromBaseConverter(baseValue, unitSelectElem.value);
                inputElem.value = newValue.toFixed(4); // Use a reasonable number of decimal places
            }
            previousUnit = unitSelectElem.value;
        });
    }

    createUnitConverter(boilingPointInput, boilingPointUnitsSelect, convertTemperatureToKelvin, convertTemperatureFromKelvin);
    createUnitConverter(criticalTempInput, criticalTempUnitsSelect, convertTemperatureToKelvin, convertTemperatureFromKelvin);
    createUnitConverter(criticalPressureInput, criticalPressureUnitsSelect, convertPressureToPascals, convertPressureFromPascals);
    
    // --- API Call and UI Update ---
    async function handleFormSubmit(event) {
        event.preventDefault();
        hideError();
        resetOutputUI();

        const tbNum = parseFloat(boilingPointInput.value);
        const tcNum = parseFloat(criticalTempInput.value);
        const pcNum = parseFloat(criticalPressureInput.value);

        if (isNaN(tbNum) || isNaN(tcNum) || isNaN(pcNum)) {
            showError('Error: Please fill in all fields with valid numbers.');
            return;
        }

        const tbInKelvin = convertTemperatureToKelvin(tbNum, boilingPointUnitsSelect.value);
        const tcInKelvin = convertTemperatureToKelvin(tcNum, criticalTempUnitsSelect.value);
        const pcInBar = convertPressureFromPascals(convertPressureToPascals(pcNum, criticalPressureUnitsSelect.value), 'bar');

        const payload = { 
            Tb: tbInKelvin, 
            Tc: tcInKelvin, 
            Pc_bar: pcInBar 
        };

        try {
            const response = await fetch('/api/calculate_lk_omega', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload),
            });
            const data = await response.json();
            if (!response.ok) throw new Error(data.error || `HTTP error! status: ${response.status}`);
            updateUIWithResults(data);
        } catch (error) {
            showError(`Calculation Error: ${error.message}`);
            resetOutputUI();
        }
    }

    // --- UI Update Functions ---
    function updateUIWithResults(data) {
        // Result
        omegaResultElem.textContent = data.omega;
        resultSourceElem.textContent = `Source: ${data.source}`;
        resultContainer.style.display = 'block';
        placeholderContainer.style.display = 'none';

        // Store raw LaTeX and show copy button
        rawLatexString = data.latex_steps;
        copyLatexBtn.style.display = 'inline-block';

        // Display plain HTML steps immediately
        stepsContainer.innerHTML = data.plain_steps;
        stepsCard.style.display = 'block';

        // Define the function that will render MathJax
        const renderMath = () => {
            stepsContainer.innerHTML = data.latex_steps;
            if (window.MathJax) {
                MathJax.typesetPromise([stepsContainer]).catch(err => console.error('MathJax typeset failed:', err));
            }
        };

        // Check if MathJax is already ready, otherwise wait for our custom event
        if (window.MathJax && window.MathJax.startup.ready) {
            renderMath();
        } else {
            document.body.addEventListener('mathjax-ready', renderMath, { once: true });
        }

        // Citations
        currentCitations = data.citations;
        displayCitations('default');
    }

    function resetOutputUI() {
        resultContainer.style.display = 'none';
        placeholderContainer.style.display = 'block';
        stepsCard.style.display = 'none';
        copyLatexBtn.style.display = 'none';
        rawLatexString = '';
        stepsContainer.innerHTML = '';
        toolCitationContainer.innerHTML = '';
        sourceCitationsContainer.innerHTML = `<p class="text-muted small">Source citations will appear here after calculation.</p>`;
        
        // Reset citation buttons to default
        citationControlButtons.forEach(btn => btn.classList.remove('active'));
        document.querySelector('.format-btn[data-format="default"]').classList.add('active');
    }

    const formatters = {
        default: (cit) => {
            const d = cit.data;
            let formatted = `<strong>[${cit.source_type}]</strong><br>`;
            const authorStr = Array.isArray(d.authors) ? d.authors.join(', ') : d.author;
            formatted += `${authorStr} (${d.year}). `;
            formatted += `<em>${d.title}</em>.`;
            if (d.journal) formatted += ` ${d.journal}, <strong>${d.volume}</strong>(${d.issue}), pp. ${d.pages}.`;
            if (d.url) formatted += ` <a href="${d.url}" target="_blank" rel="noopener noreferrer">Source Link</a>`;
            if (d.retrieved_date) formatted += ` Retrieved ${d.retrieved_date}.`;
            return `<div class="citation-text mb-2">${formatted}</div>`;
        },
        apa: (cit) => {
            const d = cit.data;
            let formatted = '';
            const authorStr = Array.isArray(d.authors) ? d.authors.map(a => a.split(',').reverse().join(' ')).join(', ') : d.author;
            formatted = `${authorStr} (${d.year}).`;
            if (d.title) formatted += ` <em>${d.title}</em>.`;
            if (d.journal) {
                 formatted += ` ${d.journal}, <em>${d.volume}</em>(${d.issue}), ${d.pages}.`;
            } else if (d.type === 'Software') {
                 formatted += ` [${d.type}].`;
            }
            if (d.retrieved_date) {
                formatted += ` Retrieved ${d.retrieved_date}, from ${d.url}`;
            } else if (d.url) {
                formatted += ` ${d.url}`;
            }
            return `<p class="citation-text">${formatted}</p>`;
        },
        bibtex: (cit) => {
            const d = cit.data;
            let bibtex = '';
            const authorStr = Array.isArray(d.authors) ? d.authors.join(' and ') : `{${d.author}}`;

            if (cit.type === 'Journal Article') {
                bibtex = `@article{${cit.id},\n  author  = {${authorStr}},\n  title   = {${d.title}},\n  journal = {${d.journal}},\n  year    = {${d.year}},\n  volume  = {${d.volume}},\n  pages   = {${d.pages}}\n}`;
            } else { // Software or Web App
                bibtex = `@misc{${cit.id},\n  author       = {${authorStr}},\n  title        = {${d.title}},\n  year         = {${d.year}},\n  url          = {${d.url}}`;
                if(d.retrieved_date) bibtex += `,\n  note = {Accessed: ${d.retrieved_date}}`
                bibtex += `\n}`;
            }
            return `<pre class="citation-bibtex">${bibtex}</pre>`;
        }
    };

    function displayCitations(format = 'default') {
        toolCitationContainer.innerHTML = '';
        sourceCitationsContainer.innerHTML = '';

        currentCitations.forEach(cit => {
            const formattedHtml = formatters[format](cit);
            if (cit.source_type === 'Tool') {
                toolCitationContainer.innerHTML += formattedHtml;
            } else {
                sourceCitationsContainer.innerHTML += formattedHtml;
            }
        });
    }

    citationControlButtons.forEach(button => {
        button.addEventListener('click', () => {
            citationControlButtons.forEach(btn => btn.classList.remove('active'));
            button.classList.add('active');
            displayCitations(button.dataset.format);
        });
    });

    // --- Event Listeners ---
    form.addEventListener('submit', handleFormSubmit);
    form.addEventListener('reset', () => {
        hideError();
        resetOutputUI();
    });
    copyLatexBtn.addEventListener('click', () => {
        if (rawLatexString) copyToClipboard(rawLatexString);
    });
}); 