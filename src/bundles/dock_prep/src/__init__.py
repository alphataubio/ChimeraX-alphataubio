# vim: set expandtab shiftwidth=4 softtabstop=4:

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

# --- Public API ---
from .cmd import dock_prep_caller

# All modules involved in the DockPrep pipeline provide these variables/functions
from .cmd import dock_prep_arg_info
from .prep import prep as run_for_dock_prep, handle_memorization, MEMORIZE_USE, MEMORIZE_SAVE, MEMORIZE_NONE

# --- Toolshed / Session-Init Functions ---
from chimerax.core.toolshed import BundleAPI
from chimerax.save_command import SaverInfo
from chimerax.open_command import OpenerInfo
from chimerax.core.commands import ModelsArg
from .dms import save_dms
from .sph import save_sph, open_sph

class DockPrepAPI(BundleAPI):

    @staticmethod
    def register_command(command_name, logger):
        from .cmd import register_command
        register_command(logger)

    @staticmethod
    def start_tool(session, tool_name):
        from .tool import DockPrepTool
        DockPrepTool(session)

    @staticmethod
    def run_provider(session, name, mgr):
        if mgr == session.open_command:
            class SphOpener(OpenerInfo):
                def open(self, session, data, file_name, **kw):
                    return open_sph(session, file_name, **kw)
            return SphOpener()

        if mgr == session.save_command:
            class DmsSaver(SaverInfo):
                def save(self, session, path, *, models=None):
                    save_dms(session, path, models)

                @property
                def save_args(self):
                    return {'models': ModelsArg}

            class SphSaver(SaverInfo):
                def save(self, session, path, *, models=None):
                    save_sph(session, path, models)

                @property
                def save_args(self):
                    return {'models': ModelsArg}

            if name == "DMS":
                return DmsSaver()
            elif name == "SPH":
                return SphSaver()
            else:
                raise ValueError(f"Unknown save format: {name}")

        raise ValueError(f"Unknown provider manager: {mgr}")

bundle_api = DockPrepAPI()


