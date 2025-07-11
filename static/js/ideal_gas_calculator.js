document.addEventListener('DOMContentLoaded', function() {
    const findButton = document.getElementById('find-btn-main');
    const nameInput = document.getElementById('name-input-main');
    const resultsContainer = document.getElementById('results-container');
    let chartInstance = null;
    let currentData = null; 

    function typesetMath(element) {
        if (window.MathJax && window.MathJax.startup?.promise) {
            window.MathJax.startup.promise.then(() => {
                window.MathJax.typesetPromise([element]).catch(console.error);
            });
        }
    }
    
    const citationManager = {
        citations: new Map(),
        add(id, cit) { if(id && cit) this.citations.set(id, cit); this.render(); },
        clear() { this.citations.clear(); },
        render(format) {
            const displayEl = document.getElementById('citation-display');
            if(!displayEl) return;
            const currentFormat = format || document.querySelector('.format-btn.active')?.dataset.format || 'default';
            const staticCits = [
                {id: 'toricalc_cp', type: 'Tool', data: {authors: ['ToriCalc'], year: new Date().getFullYear(), title: 'Ideal Gas Properties Calculator'}},
                {id: 'chemicals_lib', type: 'Library', data: {authors: ['Bell, C., et al.'], year: '2016-2024', title: 'Chemicals Library'}}
            ];
            let html = staticCits.map(cit => this.format(cit, currentFormat)).join('');
            this.citations.forEach(cit => html += this.format(cit, currentFormat));
            displayEl.innerHTML = html;
        },
        format(cit, format) {
            const d = cit.data;
            const authorStr = Array.isArray(d.authors) ? d.authors.join(' and ') : d.author || '';
            const authorApa = Array.isArray(d.authors) ? d.authors.join(', ') : d.author || '';
            const authorAcs = Array.isArray(d.authors) ? d.authors.join('; ') : d.author || '';
            let content = '';
            
            switch (format) {
                case 'apa': content = `${authorApa} (${d.year}). <em>${d.title}</em>. ${d.journal || ''}`; break;
                case 'acs': content = `${authorAcs} ${d.title}. <em>${d.journal}</em>, <strong>${d.year}</strong>.`; break;
                case 'bibtex': content = `<pre class="bibtex">@misc{${cit.id},\n  author = {${authorStr}},\n  title = {${d.title}},\n  year = {${d.year}}\n}</pre>`; break;
                default:
                    let title = d.journal ? `"${d.title}"` : `<em>${d.title}</em>`;
                    content = `<strong>[${cit.id}]</strong> ${authorApa} (${d.year}). ${title}`;
            }
            const plainText = content.replace(/<[^>]+>/g, '');
            return `<div class="citation-item d-flex justify-content-between align-items-center"><div>${content}</div><button class="btn btn-outline-secondary btn-sm py-0 copy-citation-btn" data-clipboard-text="${escape(plainText)}">Copy</button></div>`;
        }
    };

    findButton.addEventListener('click', handleLookup);
    nameInput.addEventListener('keyup', e => { if (e.key === 'Enter') handleLookup(); });
    
    resultsContainer.addEventListener('click', e => {
        if(e.target.matches('#plot-btn')) plotChart();
        if(e.target.matches('.format-btn')) {
            document.querySelectorAll('.format-btn').forEach(btn => btn.classList.remove('active'));
            e.target.classList.add('active');
            citationManager.render(e.target.dataset.format);
        }
        if (e.target.matches('.copy-citation-btn')) {
            navigator.clipboard.writeText(unescape(e.target.dataset.clipboardText));
            e.target.textContent = 'Copied!'; setTimeout(() => { e.target.textContent = 'Copy'; }, 2000);
        }
        if (e.target.matches('.copy-latex-btn')) {
            navigator.clipboard.writeText(unescape(e.target.dataset.latex));
            e.target.textContent = 'Copied!'; setTimeout(() => { e.target.textContent = 'Copy LaTeX'; }, 2000);
        }
    });

    async function handleLookup() {
        const chemicalName = nameInput.value.trim();
        if (!chemicalName) return;
        resultsContainer.innerHTML = '<div class="d-flex align-items-center"><div class="spinner-border spinner-border-sm me-2"></div>Searching...</div>';
        findButton.disabled = true;
        citationManager.clear();
        try {
            const response = await fetch('/api/ideal_gas/get_comparison_data', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ chemical_name: chemicalName }) });
            currentData = await response.json();
            if (!response.ok) throw new Error(currentData.error);
            renderResults(currentData);
            plotChart();
            if(currentData.citations) {
                currentData.citations.forEach(c => citationManager.add(c.id, c));
            }
        } catch (err) {
            resultsContainer.innerHTML = `<div class="alert alert-danger p-2"><strong>Error:</strong> ${err.message}</div>`;
        } finally {
            findButton.disabled = false;
        }
    }

    function renderResults(data) {
        let perryHtml = data.perry_results.error ? `<p class="text-danger small">${data.perry_results.error}</p>` : `<div class="glass-box"><p class="small mb-2"><strong>Perry's Equation:</strong> $C_p = a + bT + c/T^2 + dT^2$</p><div class="table-responsive"><table class="table table-sm table-bordered mb-2"><thead><tr><th>a</th><th>b</th><th>c</th><th>d</th></tr></thead><tbody><tr><td>${data.perry_results.coeffs.a.toExponential(3)}</td><td>${data.perry_results.coeffs.b.toExponential(3)}</td><td>${data.perry_results.coeffs.c.toExponential(3)}</td><td>${data.perry_results.coeffs.d.toExponential(3)}</td></tr></tbody></table></div><button class="btn btn-sm btn-outline-secondary copy-latex-btn" data-latex="${escape(data.perry_results.latex_steps)}">Copy LaTeX</button></div><p class="small text-muted mt-2 mb-0">Valid from ${data.perry_results.T_low} K to ${data.perry_results.T_high} K. Note: Original units in cal/mol/K.</p>`;
        let jobackHtml = data.joback_results.error ? `<p class="text-danger small">${data.joback_results.error}</p>` : `<div class="glass-box" id="joback-glass-box">${data.joback_results.plain_steps}</div><button class="btn btn-sm btn-outline-secondary mt-2 copy-latex-btn" data-latex="${escape(data.joback_results.latex_steps)}">Copy LaTeX</button>`;

        resultsContainer.innerHTML = `
            <h4>Comparison for: <strong>${data.chemical_name}</strong></h4>
            <div class="row mt-3">
                <div class="col-lg-6 mb-3"><div class="card h-100"><div class="card-header"><h6>Database Method (Perry's)</h6></div><div class="card-body small">${perryHtml}</div></div></div>
                <div class="col-lg-6 mb-3"><div class="card h-100"><div class="card-header"><h6>Estimation (Joback "Glass Box")</h6></div><div class="card-body small">${jobackHtml}</div></div></div>
            </div>
            <div class="card mt-3"><div class="card-body"><h5>Interactive Comparison Plot</h5><div class="row g-2 align-items-center"><div class="col-auto"><label class="col-form-label-sm">T(K):</label></div><div class="col"><input type="number" class="form-control form-control-sm" id="tmin-input" value="300"></div><div class="col-auto">to</div><div class="col"><input type="number" class="form-control form-control-sm" id="tmax-input" value="1000"></div><div class="col-auto"><button class="btn btn-secondary btn-sm" id="plot-btn">Re-Plot</button></div></div><canvas id="cp-chart" class="mt-3"></canvas></div></div>
            <div class="mt-4"><hr><h5>Citations</h5><div class="citation-controls"><span class="me-2">Format:</span><div class="btn-group btn-group-sm"><button class="btn btn-outline-secondary format-btn active" data-format="default">Default</button><button class="btn btn-outline-secondary format-btn" data-format="apa">APA</button><button class="btn btn-outline-secondary format-btn" data-format="acs">ACS</button><button class="btn btn-outline-secondary format-btn" data-format="bibtex">BibTeX</button></div></div><div id="citation-display" class="mt-2 citation-display"></div></div>`;
        
        typesetMath(resultsContainer);
    }

    function plotChart() {
        if (!currentData) return;
        const tMin = parseFloat(document.getElementById('tmin-input').value);
        const tMax = parseFloat(document.getElementById('tmax-input').value);
        if (isNaN(tMin) || isNaN(tMax) || tMin >= tMax) { alert("Invalid temperature range."); return; }
        
        const labels = [], jobackData = [], perryData = [];
        const step = (tMax - tMin) / 50;
        for (let T = tMin; T <= tMax; T += step) {
            labels.push(T.toFixed(1));
            if (!currentData.joback_results.error) {
                const jc = currentData.joback_results.coeffs;
                jobackData.push(jc.A + jc.B*T + jc.C*T**2 + jc.D*T**3);
            }
            if (!currentData.perry_results.error) {
                const pc = currentData.perry_results.coeffs;
                perryData.push((pc.a + pc.b*T + pc.c/(T**2) + pc.d*T**2) * 4.184);
            }
        }
        
        const datasets = [];
        if (jobackData.length > 0) datasets.push({ label: 'Joback Estimate', data: jobackData, borderColor: 'rgba(255, 99, 132, 1)', tension: 0.1, fill: false });
        if (perryData.length > 0) datasets.push({ label: `Perry's Polynomial (Database)`, data: perryData, borderColor: 'rgba(54, 162, 235, 1)', tension: 0.1, fill: false });
        
        const ctx = document.getElementById('cp-chart').getContext('2d');
        if (chartInstance) chartInstance.destroy();
        chartInstance = new Chart(ctx, { type: 'line', data: { labels: labels, datasets: datasets }, options: { scales: { x: { title: { display: true, text: 'Temperature (K)' } }, y: { title: { display: true, text: 'Cp (J/mol·K)' } } } } });
    }
});