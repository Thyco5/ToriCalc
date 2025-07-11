document.addEventListener('DOMContentLoaded', function () {
    const navItemDropdown = document.querySelector('.nav-item-dropdown');
    const calculatorsLink = navItemDropdown ? navItemDropdown.querySelector('.nav-links') : null;
    const dropdownMenu = navItemDropdown ? navItemDropdown.querySelector('.dropdown-menu') : null;
    const submenuContainers = document.querySelectorAll('.dropdown-submenu-container');

    const allDropdownItemsWithSubmenus = dropdownMenu ? Array.from(dropdownMenu.querySelectorAll('.dropdown-item, .dropdown-submenu-container')) : [];

    function closeAllCalculatorSubSections() {
        submenuContainers.forEach(container => {
            const submenu = container.querySelector('.submenu');
            container.classList.remove('open'); // For in-panel
            container.classList.remove('active-dropdown-item');

            if (submenu && !container.dataset.submenuType === 'in-panel') { // only for flyouts
                submenu.classList.remove('open-by-js');
                submenu.classList.add('hidden'); // Fallback if animations not fully set up
                submenu.style.display = 'none';
                submenu.style.left = '';
                submenu.style.right = '';
                submenu.style.top = '';
                submenu.classList.remove('flipped-left');
            }
            // For in-panel, the .open class removal handles hiding via max-height
        });
        // Remove active class from direct dropdown items that are not submenu containers
        dropdownMenu.querySelectorAll('.dropdown-item:not(.dropdown-submenu-container)').forEach(item => {
            item.classList.remove('active-dropdown-item');
        });
    }

    function positionFlyoutSubmenu(submenuElement, parentElement) {
        if (!submenuElement || !parentElement) return;

        submenuElement.classList.remove('hidden', 'flipped-left');
        submenuElement.classList.add('open-by-js'); // Use for animation
        // submenuElement.style.display = 'block'; // CSS will handle with .open-by-js
        
        // Default fly-out position (to the right)
        submenuElement.style.left = '100%';
        submenuElement.style.right = 'auto';
        submenuElement.style.top = '0';

        const rect = submenuElement.getBoundingClientRect();
        const parentRect = parentElement.getBoundingClientRect();
        const viewportWidth = window.innerWidth;
        const viewportHeight = window.innerHeight;

        if (rect.right > viewportWidth) {
            submenuElement.style.left = 'auto';
            submenuElement.style.right = '100%';
            submenuElement.classList.add('flipped-left');
        }

        if (rect.bottom > viewportHeight) {
            const newTop = parentRect.top - (rect.bottom - viewportHeight) - 5;
            submenuElement.style.top = Math.max(0, newTop - parentRect.top) + 'px';
        }

        const finalRect = submenuElement.getBoundingClientRect();
        if (finalRect.left < 0 && submenuElement.classList.contains('flipped-left')) {
            submenuElement.style.left = '100%';
            submenuElement.style.right = 'auto';
            submenuElement.classList.remove('flipped-left');
        }
    }

    if (calculatorsLink) {
        calculatorsLink.addEventListener('click', function(event) {
            event.preventDefault();
        });
    }

    if (navItemDropdown && dropdownMenu) {
        navItemDropdown.addEventListener('mouseenter', () => {
            // dropdownMenu.style.display = 'block'; // Replaced by class for animation
            dropdownMenu.classList.add('open-by-js');
        });

        navItemDropdown.addEventListener('mouseleave', () => {
            setTimeout(() => {
                let isHoveringSubmenuRelated = false;
                if (navItemDropdown.matches(':hover')) isHoveringSubmenuRelated = true;
                // Check if still hovering any part of the open dropdown or its submenus
                if (dropdownMenu.classList.contains('open-by-js') && dropdownMenu.matches(':hover')) {
                    isHoveringSubmenuRelated = true;
                }

                if (!isHoveringSubmenuRelated) {
                    // dropdownMenu.style.display = 'none'; // Replaced by class for animation
                    dropdownMenu.classList.remove('open-by-js');
                    closeAllCalculatorSubSections();
                }
            }, 150);
        });
    }

    allDropdownItemsWithSubmenus.forEach(item => {
        const trigger = item.classList.contains('dropdown-submenu-container') ? item.querySelector('.submenu-trigger') : item;
        const submenuContainer = item.classList.contains('dropdown-submenu-container') ? item : null;
        const submenu = submenuContainer ? submenuContainer.querySelector('.submenu') : null;

        trigger.addEventListener('mouseenter', () => {
            closeAllCalculatorSubSections(); // Close all others first

            if (submenuContainer && submenu) {
                submenuContainer.classList.add('active-dropdown-item');
                if (submenuContainer.dataset.submenuType === 'in-panel') {
                    submenuContainer.classList.add('open'); // Triggers CSS for max-height
                } else {
                    // It's a fly-out
                    positionFlyoutSubmenu(submenu, submenuContainer);
                }
            } else if (!submenuContainer) {
                // This is a direct dropdown item without a submenu
                item.classList.add('active-dropdown-item');
            }
        });

        // No specific mouseleave for individual items if covered by parent
        // Mouseleave for submenuContainer is tricky due to in-panel vs flyout
    });


    const menuIcon = document.querySelector('.menu-icon');
    const navMenu = document.querySelector('.nav-menu');

    if (menuIcon && navMenu) {
        menuIcon.addEventListener('click', () => {
            navMenu.classList.toggle('active');
            const isExpanded = navMenu.classList.contains('active');
            menuIcon.setAttribute('aria-expanded', isExpanded.toString());
            // Add CSS for .nav-menu.active like:
            // .nav-menu.active {
            //   display: flex; flex-direction: column; position: absolute;
            //   top: 70px; left: 0; right: 0; background: #fff; z-index: 999;
            // }
        });
    }

    const currentPath = window.location.pathname;
    const navLinks = document.querySelectorAll('.nav-links, .dropdown-item, .submenu-item');

    navLinks.forEach(link => {
        if (link.getAttribute('href') === currentPath && currentPath !== "#") {
            link.classList.add('active');
            const mainDropdown = link.closest('.nav-item-dropdown');
            if (mainDropdown) {
                mainDropdown.querySelector('.nav-links').classList.add('active');
            }
            const submenuContainer = link.closest('.dropdown-submenu-container');
            if (submenuContainer) {
                 const submenuTrigger = submenuContainer.querySelector('.submenu-trigger');
                 if(submenuTrigger) submenuTrigger.classList.add('active');
                 if (submenuContainer.dataset.submenuType === 'in-panel' && link.closest('.in-panel-submenu')){
                    submenuContainer.classList.add('open'); // Auto-open in-panel if child is active
                 }
            }
        }
    });

    // ARIA attributes 
    if (navItemDropdown) {
        const mainCalculatorsLink = navItemDropdown.querySelector('.nav-links');
        if (mainCalculatorsLink) {
            mainCalculatorsLink.setAttribute('aria-haspopup', 'true');
            mainCalculatorsLink.setAttribute('aria-expanded', 'false');
        }
        if (dropdownMenu) {
            dropdownMenu.setAttribute('aria-label', 'Calculators Menu');
            // Update aria-expanded based on .open-by-js class
            // Using MutationObserver would be more robust for class changes
            const observer = new MutationObserver(() => {
                mainCalculatorsLink?.setAttribute('aria-expanded', dropdownMenu.classList.contains('open-by-js').toString());
            });
            observer.observe(dropdownMenu, { attributes: true, attributeFilter: ['class'] });
        }
    }
    submenuContainers.forEach(container => {
        const trigger = container.querySelector('.submenu-trigger');
        const submenu = container.querySelector('.submenu');
        if (trigger) {
            trigger.setAttribute('aria-haspopup', 'true');
            trigger.setAttribute('aria-expanded', 'false');
        }
        if (submenu) {
            submenu.setAttribute('aria-label', trigger ? trigger.textContent.replace(/<[^>]*>/g, '').trim() + ' Submenu' : 'Submenu');
            const observer = new MutationObserver(() => {
                const isOpen = container.classList.contains('open') || submenu.classList.contains('open-by-js');
                trigger?.setAttribute('aria-expanded', isOpen.toString());
            });
            if (container.dataset.submenuType === 'in-panel') {
                observer.observe(container, { attributes: true, attributeFilter: ['class'] });
            } else {
                observer.observe(submenu, { attributes: true, attributeFilter: ['class'] });
            }
        }
    });
}); 