// config-dirty.js — Shiny message handler for sidebar dirty-state badge
//
// Receives: { sidebar_id: "config_association", dirty_count: 3 }
// Toggles a .config-dirty-badge span inside the sidebar title element.
// The sidebar title selector targets the bslib sidebar panel's header.

Shiny.addCustomMessageHandler("config_dirty_update", function(msg) {
    var sid  = msg.sidebar_id;
    var n    = msg.dirty_count;

    // The bslib sidebar toggle button has aria-controls = "{id}"
    // The title text span is inside the sidebar panel itself
    // Target: .bslib-sidebar-layout > .sidebar > .sidebar-title
    // Fall back to searching by sidebar id attribute
    var container = document.querySelector(
        '[data-bslib-sidebar-id="' + sid + '"], #' + sid
    );

    // Also find the sidebar toggle button (for badge on the button)
    var toggleBtn = document.querySelector('[aria-controls="' + sid + '"]');

    if (!container && !toggleBtn) return;

    // Manage badge on the toggle button
    if (toggleBtn) {
        var badge = toggleBtn.querySelector('.config-dirty-badge');
        if (n > 0) {
            if (!badge) {
                badge = document.createElement('span');
                badge.className = 'config-dirty-badge';
                toggleBtn.appendChild(badge);
            }
            badge.textContent = 'modified';
        } else {
            if (badge) badge.remove();
        }
    }

    // Also update the sidebar panel title if accessible
    if (container) {
        var titleEl = container.querySelector('.config-sidebar-title');
        if (titleEl) {
            var panelBadge = titleEl.querySelector('.config-dirty-badge');
            if (n > 0) {
                if (!panelBadge) {
                    panelBadge = document.createElement('span');
                    panelBadge.className = 'config-dirty-badge ms-2';
                    titleEl.appendChild(panelBadge);
                }
                panelBadge.textContent = 'modified';
            } else {
                if (panelBadge) panelBadge.remove();
            }
        }
    }
});
