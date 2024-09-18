# vim: set expandtab ts=4 sw=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2022 Regents of the University of California. All rights reserved.
# The ChimeraX application is provided pursuant to the ChimeraX license
# agreement, which covers academic and commercial uses. For more details, see
# <https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html>
#
# This particular file is part of the ChimeraX library. You can also
# redistribute and/or modify it under the terms of the GNU Lesser General
# Public License version 2.1 as published by the Free Software Foundation.
# For more details, see
# <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>
#
# THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
# EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. ADDITIONAL LIABILITY
# LIMITATIONS ARE DESCRIBED IN THE GNU LESSER GENERAL PUBLIC LICENSE
# VERSION 2.1
#
# This notice must be embedded in or attached to all copies, including partial
# copies, of the software or any revisions or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.toolshed import BundleAPI

class _MapEraserAPI(BundleAPI):

    @staticmethod
    def start_tool(session, tool_name):
        from .eraser import map_eraser_panel, MapEraser
        p = map_eraser_panel(session)
        # Bind mouse button when panel shown.
        mm = session.ui.mouse_modes
        mm.bind_mouse_mode('right', [], MapEraser(session))
        return p

    @staticmethod
    def initialize(session, bundle_info):
        """Register map eraser mouse mode."""
        if session.ui.is_gui:
            from . import eraser
            eraser.register_mousemode(session)

    @staticmethod
    def finish(session, bundle_info):
        # TODO: remove mouse mode
        pass

    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        from . import eraser
        eraser.register_volume_erase_command(logger)

bundle_api = _MapEraserAPI()
