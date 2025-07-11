// static/js/joback_calculator.js
(function() {
  console.log("Joback Calculator JS Loaded"); // <<< DO YOU SEE THIS IN CONSOLE?

  // ... (consts, utils, DOM refs) ...
  const form = document.getElementById('joback_form');
  const tbInput = document.getElementById('Tb_joback');
  const tbUnitSelect = document.getElementById('Tb_joback_unit');
  const smilesInput = document.getElementById('smiles_input_joback');
  // ... other DOM refs ...

  if (!form) {
      console.error("Joback form #joback_form not found! Aborting script."); // <<< CHECK FOR THIS
      return;
  }
  // ... (helper functions like formatCitation, updateDisplayFromData, initializeOutputSection, getAndValidateInputs, etc.) ...

  async function triggerCalculation() {
      console.log("Joback: triggerCalculation called"); // <<< DO YOU SEE THIS?
      // ... (validation logic) ...
      // ... (fetch call) ...
  }

  // Event Listener for Tb Unit Change
  if (tbUnitSelect && tbInput) {
      tbUnitSelect.addEventListener('change', function() {
          console.log("Joback: Tb unit changed"); // <<< DO YOU SEE THIS?
          // ... (conversion logic) ...
          currentTbUnit = this.value;
          triggerCalculation();
      });
  }
  // ... (other listeners like format selects, copy buttons) ...

  // ---- CRITICAL LISTENERS FOR CALCULATION ----
  form.addEventListener('submit', function(e) {
      console.log("Joback: Form submitted (Calculate button clicked)"); // <<< DO YOU SEE THIS?
      e.preventDefault();
      triggerCalculation();
  });

  // Auto-calculate on input change
  if (tbInput) {
      tbInput.addEventListener('input', () => {
          console.log("Joback: Tb input changed"); // <<< DO YOU SEE THIS?
          triggerCalculation();
      });
  }
  if (smilesInput) {
      smilesInput.addEventListener('input', () => {
          console.log("Joback: SMILES input changed"); // <<< DO YOU SEE THIS?
          triggerCalculation();
      });
  }
  // ---- END CRITICAL LISTENERS ----

  document.addEventListener('DOMContentLoaded', () => {
      console.log("Joback: DOMContentLoaded"); // <<< DO YOU SEE THIS?
      initializeOutputSection();
      // Other initializations...
       if (tbUnitSelect) { // Store initial Tb unit for conversion logic
          currentTbUnit = tbUnitSelect.value;
      }
  });

})();