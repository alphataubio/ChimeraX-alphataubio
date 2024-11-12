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

import traceback

from chimerax.atomic import AtomicStructure, Element

def determine_element_from_mass(mass, *, consider_hydrogens=True):
    H = Element.get_element('H')
    nearest = None
    for high in range(1, Element.NUM_SUPPORTED_ELEMENTS+1):
        if Element.get_element(high).mass > mass:
            break
    else:
        high = Element.NUM_SUPPORTED_ELEMENTS

    if high == 1:
        return H

    if consider_hydrogens:
        max_hyds = 6
    else:
        max_hyds = 0

    for num_hyds in range(max_hyds+1):
        adj_mass = mass - num_hyds * H.mass
        low_mass = Element.get_element(high-1).mass
        while low_mass > adj_mass and high > 1:
            high -= 1
            low_mass = Element.get_element(high-1).mass
        high_mass = Element.get_element(high).mass
        low_diff = abs(adj_mass - low_mass)
        high_diff = abs(adj_mass - high_mass)
        if low_diff < high_diff:
            diff = low_diff
            element = high-1
        else:
            diff = high_diff
            element = high
        if nearest is None or diff < nearest[1]:
            nearest = (element, diff)
    return Element.get_element(nearest[0])

def read_data(session, path, file_name, *, auto_style=True, coords=None, **kw):
    from chimerax.core.errors import UserError, CancelOperation
    import os
    if coords is None:
        if session.ui.is_gui and not session.in_script:
            from Qt.QtWidgets import QFileDialog
            # Don't use a native dialog so that the caption is actually shown;
            # otherwise the dialog is totally mystifying
            coords, types = QFileDialog.getOpenFileName(caption="Specify DUMP file for DATA",
                directory=os.path.dirname(path), options=QFileDialog.DontUseNativeDialog)
            if not coords:
                raise CancelOperation("No coordinates file specified for DATA")
            session.logger.info("Coordinates file: %s" % coords)
        else:
            raise UserError("'coords' keyword with coordinate-file argument must be supplied")
    from chimerax.data_formats import NoFormatError
    try:
        data_fmt = session.data_formats.open_format_from_file_name(coords)
    except NoFormatError as e:
        raise UserError("Cannot determine format of coordinates file '%s' from suffix" % coords)

    structure = AtomicStructure(session, name=os.path.basename(file_name), auto_style=auto_style)

    try:
        from chimerax.atomic.struct_edit import add_atom, add_bond
        from numpy import array, float64

        data = open(path, "r")
        data.readline()
        data.readline()

        # READ NUMBER OF ATOMS AND BONDS

        line = data.readline()
        while line != '\n':
          tokens = line.split()

          if tokens[1] == 'atoms':
            atoms = int(tokens[0])
          elif tokens[1] == 'bonds':
            bonds = int(tokens[0])

          line = data.readline()

        session.logger.info( f"{atoms} atoms {bonds} bonds")

        # SKIP UNTIL MASSES SECTION

        line = data.readline()
        while not line.startswith("Masses"): line = data.readline()
        line = data.readline() # SKIP BLANK LINE

        # PARSE MASSES

        masses = {}

        tokens = data.readline().split()
        while tokens and tokens[0].isdigit():
          masses[int(tokens[0])] = float(tokens[1])
          tokens = data.readline().split()

        # print( masses )

        # SKIP UNTIL ATOMS SECTION

        line = data.readline()
        while not line.startswith("Atoms"): line = data.readline()
        line = data.readline() # SKIP BLANK LINE

        # PARSE ATOMS

        atoms_list = []
        atoms = {}

        tokens = data.readline().split()
        while tokens:
          tag = int(tokens[0])
          mol = int(tokens[1])
          type = int(tokens[2])
          xyz = array([float(tokens[4]),float(tokens[5]),float(tokens[6])], dtype=float64)
          residue = structure.find_residue(" ", mol)
          if residue is None: residue = structure.new_residue(str(mol), " ", mol)
          element = determine_element_from_mass(masses[type])
          atoms_list.append([tag, element, residue, xyz])
          tokens = data.readline().split()

        atoms_list.sort(key=lambda atom:atom[0])

        for atom in atoms_list:
          atoms[atom[0]] = add_atom(str(atom[0]), atom[1], atom[2], atom[3], serial_number=atom[0])

        # SKIP UNTIL BONDS SECTION

        line = data.readline()
        while not line.startswith("Bonds"): line = data.readline()
        line = data.readline() # SKIP BLANK LINE

        # PARSE BONDS

        tokens = data.readline().split()
        while tokens:
          tag1 = int(tokens[2])
          tag2 = int(tokens[3])
          # FIXME: handle tag1 and/or tag2 not found
          add_bond(atoms[tag1], atoms[tag2])
          tokens = data.readline().split()

        data.close()

    except Exception as e:
        print(traceback.format_exc())
        raise UserError("Problem reading/processing DATA file '%s': %s" % (path, e))

    from .read_dump import read_dump
    read_dump(session, coords, structure, data_fmt.nicknames[0], replace=True, **kw)
    return [structure], ""
