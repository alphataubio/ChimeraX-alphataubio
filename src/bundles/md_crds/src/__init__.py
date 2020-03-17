# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.toolshed import BundleAPI

class _MDCrdsBundleAPI(BundleAPI):
    
    from chimerax.atomic import StructureArg

    @staticmethod
    def open_file(session, path, format_name, structure_model=None, replace=True):
        if structure_model is None:
            from chimerax.core.errors import UserError
            raise UserError("Must specify a structure model to read the coordinates into")
        from .read_coords import read_coords
        num_coords = read_coords(session, path, structure_model, format_name, replace=replace)
        if replace:
            return [], "Replaced existing frames of %s with  %d new frames" % (structure_model,
                num_coords)
        return [], "Added %d frames to %s" % (num_coords, model)

    @staticmethod
    def run_provider(session, name, mgr, **kw):
        print("MD run_provider, name:", repr(name))
        from chimerax.open import OpenerInfo
        class MDOpenerInfo(OpenerInfo):
            print("at definition time, name is", name)
            def open(self, session, data, file_name, *, structure_model=None, md_type=name, replace=True, **kw):
                if structure_model is None:
                    from chimerax.core.errors import UserError
                    raise UserError("Must specify a structure model to read the coordinates into")
                from .read_coords import read_coords
                print("at call time, md_type is", md_type)
                num_coords = read_coords(session, data, structure_model, md_type, replace=replace)
                if replace:
                    return [], "Replaced existing frames of %s with  %d new frames" % (structure_model,
                        num_coords)
                return [], "Added %d frames to %s" % (num_coords, structure_model)
            
            @property
            def open_args(self):
                from chimerax.atomic import StructureArg
                from chimerax.core.commands import BoolArg
                return {
                    'structure_model': StructureArg,
                    'replace': BoolArg
                }
        return MDOpenerInfo()

    @staticmethod
    def save_file(session, path, format_name, models=None):
        from chimerax import atomic
        if models is None:
            models = atomic.all_structures(session)
        else:
            models = [m for m in models if isinstance(m, atomic.Structure)]
        if len(models) == 0:
            from chimerax.core.errors import UserError
            raise UserError("Must specify models to write DCD coordinates")
        # Check that all models have same number of atoms.
        nas = set(m.num_atoms for m in models)
        if len(nas) != 1:
            from chimerax.core.errors import UserError
            raise UserError("Models have different number of atoms")
        from .write_coords import write_coords
        write_coords(session, path, format_name, models)

bundle_api = _MDCrdsBundleAPI()
