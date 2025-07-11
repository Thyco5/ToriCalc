document.addEventListener('DOMContentLoaded', function () {
    const calculatorsDropdownElement = document.querySelector('.calculators-dropdown');
    const mainCalculatorsLink = document.getElementById('calculatorsDropdown');
    
    // Items that have vertically expanding submenus
    const verticalExpandItems = document.querySelectorAll('.calculators-dropdown .vertical-expand-item');

    // 1. Allow main 'Calculators' link to navigate, Bootstrap JS handles dropdown toggle.
    if (mainCalculatorsLink) {
        // The event.preventDefault() was removed from here.
        // Bootstrap's JS will handle showing/hiding the main dropdown on click.
        // Our custom CSS will handle showing it on hover.
    }

    function closeAllVerticalExpansions(exceptThisItem = null) {
        verticalExpandItems.forEach(item => {
            if (item !== exceptThisItem) {
                item.classList.remove('open');
                const trigger = item.querySelector('.vertical-expand-trigger');
                if (trigger) {
                    trigger.classList.remove('active-section');
                    trigger.setAttribute('aria-expanded', 'false');
                }
            }
        });
    }

    verticalExpandItems.forEach(item => {
        const trigger = item.querySelector('.vertical-expand-trigger');
        // const submenu = item.querySelector('.dropdown-menu-vertical-expansion'); // Not directly used for show/hide, CSS handles it.

        if (trigger) {
            trigger.addEventListener('mouseenter', function () {
                closeAllVerticalExpansions(item); // Close others, keep this one if already open by click
                item.classList.add('open');
                trigger.classList.add('active-section');
                trigger.setAttribute('aria-expanded', 'true');
            });

            // Clicking the trigger will also toggle its own expansion panel
            trigger.addEventListener('click', function (event) {
                event.preventDefault();
                event.stopPropagation(); // Prevent Bootstrap from closing main dropdown
                
                const currentlyOpen = item.classList.contains('open');
                closeAllVerticalExpansions(); // Close all first for exclusive behavior on click
                
                if (!currentlyOpen) { // If it wasn't open, open it now
                    item.classList.add('open');
                    trigger.classList.add('active-section');
                    trigger.setAttribute('aria-expanded', 'true');
                } else {
                    // If it was open and we clicked it again, it's now closed by closeAllVerticalExpansions
                    // and its active states are removed.
                    trigger.setAttribute('aria-expanded', 'false'); // Ensure it's false if it was toggled off
                }
            });
        }
    });

    // Clear highlights and close panels when mouse leaves the entire Calculators dropdown area
    if (calculatorsDropdownElement) {
        calculatorsDropdownElement.addEventListener('mouseleave', function(e) {
            setTimeout(() => {
                if (!calculatorsDropdownElement.matches(':hover')) {
                    // Check if any child within the dropdown is focused (e.g. keyboard nav)
                    let isChildFocused = false;
                    if (document.activeElement) {
                        isChildFocused = calculatorsDropdownElement.contains(document.activeElement);
                    }
                    if (!isChildFocused) { // Only close if no child is focused
                         closeAllVerticalExpansions();
                    }
                }
            }, 200); 
        });
    }
}); 