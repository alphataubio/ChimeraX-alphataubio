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

from chimerax.core.errors import UserError

default_criteria = "dchp"
def swap_aa(session, residues, res_type, *, angle_slop=None, bfactor=None, criteria=default_criteria,
    density=None, dist_slop=None, hbond_allowance=None, ignore_other_models=False, lib=None, log=True,
    preserve=None, relax=True, retain=False, score_method="num", overlap_cutoff=None):
    ''' Command to swap amino acid side chains '''

    residues = _check_residues(residues)

    if type(criteria) == str:
        for c in criteria:
            if c not in "dchp":
                raise UserError("Unknown criteria: '%s'" % c)
    elif preserve is not None:
        raise UserError("'preserve' not compatible with Nth-most-probable criteria")

    if lib is None:
        lib = _get_lib(session)

    if log:
        session.logger.info("Using %s library" % lib)

    from . import swap_res
    swap_res.swap_aa(session, residues, res_type, bfactor=bfactor, clash_hbond_allowance=hbond_allowance,
        clash_score_method=score_method, clash_overlap_cutoff=overlap_cutoff,
        criteria=criteria, density=density, hbond_angle_slop=angle_slop,
        hbond_dist_slop=dist_slop, ignore_other_models=ignore_other_models, lib=lib, log=log,
        preserve=preserve, hbond_relax=relax, retain=retain)

from chimerax.core.state import StateManager
class _RotamerStateManager(StateManager):
    def __init__(self, session, base_residue, rotamers):
        self.init_state_manager(session, "residue rotamers")
        self.session = session
        self.base_residue = base_residue
        self.rotamers = list(rotamers) # don't want auto-shrinking of a Collection
        self.group = session.models.add_group(rotamers, name="%s rotamers"
            % base_residue.string(omit_structure=True), parent=base_residue.structure)
        from chimerax.atomic import get_triggers
        self.handler = get_triggers().add_handler('changes', self._changes_cb)
        from chimerax.core.triggerset import TriggerSet
        self.triggers = TriggerSet()
        self.triggers.add_trigger('fewer rotamers') # but not zero
        self.triggers.add_trigger('self destroyed')

    def destroy(self):
        self.handler.remove()
        if self.group.id is not None:
            self.session.models.close([self.group])
        self.group = self.base_residue = self.rotamers = self.session = None
        super().destroy()

    def reset_state(self, session):
        self.triggers.activate_trigger('self destroyed', self)
        self.destroy()

    def _changes_cb(self, trigger_name, changes):
        if changes.num_deleted_residues() == 0:
            return
        remaining = [rot for rot in self.rotamers if not rot.deleted]
        if self.base_residue.deleted:
            self.triggers.activate_trigger('self destroyed', self)
            self.destroy()
            return
        remaining = [rot for rot in self.rotamers if not rot.deleted]
        if len(remaining) < len(self.rotamers):
            if remaining:
                self.rotamers = remaining
                self.triggers.activate_trigger('fewer rotamers', self)
            else:
                self.triggers.activate_trigger('self destroyed', self)
                self.destroy()

def rotamers(session, residues, res_type, *, lib=None):
    ''' Command to display possible side-chain rotamers '''

    residues = _check_residues(residues)

    if lib is None:
        lib = _get_lib(session)

    ret_val = []
    from . import swap_res
    from chimerax.atomic import AtomicStructures
    from chimerax.core.objects import Objects
    for r in residues:
        if res_type == "same":
            r_type = r.name
        else:
            r_type = res_type.upper()
        rotamers = swap_res.get_rotamers(session, r, res_type=r_type, lib=lib)
        mgr = _RotamerStateManager(session, r, rotamers)
        if session.ui.is_gui:
            from .tool import RotamerDialog
            RotamerDialog(session,
                "%s Side-Chain Rotamers" % r, mgr, res_type, session.rotamers.library(lib))
        ret_val.append(mgr)
        rot_structs = AtomicStructures(rotamers)
        from chimerax.std_commands.color import color
        color(session, Objects(atoms=rot_structs.atoms), color="byelement")
    return ret_val

def _check_residues(residues):
    residues = [r for r in residues if r.polymer_type == r.PT_AMINO]
    if not residues:
        raise UserError("No amino acid residues specified for swapping")
    return residues

def _get_lib(session):
    available_libs = session.rotamers.library_names()
    for lib_name in available_libs:
        if "Dunbrack" in lib_name:
            lib = lib_name
            break
    else:
        if available_libs:
            lib = list(available_libs)[0]
        else:
            raise UserError("No rotamer libraries installed!")
    return lib

def register_command(command_name, logger):
    from chimerax.core.commands import CmdDesc, register, StringArg, BoolArg, NonNegativeIntArg, Or
    from chimerax.core.commands import NonNegativeFloatArg, DynamicEnum, ListOf, FloatArg, EnumOf
    from chimerax.atomic import ResiduesArg
    from chimerax.map import MapArg
    desc = CmdDesc(
        required = [('residues', ResiduesArg), ('res_type', StringArg)],
        keyword = [
            ('angle_slop', FloatArg),
            ('bfactor', FloatArg),
            ('criteria', Or(ListOf(NonNegativeIntArg), StringArg)),
            ('density', MapArg),
            ('dist_slop', FloatArg),
            ('hbond_allowance', FloatArg),
            ('ignore_other_models', BoolArg),
            ('lib', DynamicEnum(logger.session.rotamers.library_names)),
            ('log', BoolArg),
            ('preserve', NonNegativeFloatArg),
            ('relax', BoolArg),
            ('retain', BoolArg),
            ('score_method', EnumOf(('sum', 'num'))),
            ('overlap_cutoff', FloatArg),
        ],
        synopsis = 'Swap amino acid side chain(s)'
    )
    register("swapaa", desc, swap_aa, logger=logger)

    desc = CmdDesc(
        required = [('residues', ResiduesArg), ('res_type', StringArg)],
        keyword = [
            ('lib', DynamicEnum(logger.session.rotamers.library_names)),
        ],
        synopsis = 'Show possible side-chain rotamers'
    )
    register("swapaa interactive", desc, rotamers, logger=logger)