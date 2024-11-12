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

from chimerax.core.errors import UserError, LimitationError

def read_dump(session, file_name, model, format_name, *, replace=True, start=1, step=1, end=None):
    from numpy import array, float64

    dump = open(file_name, "r")
    dump.readline()
    timestep = int(dump.readline().split()[0])
    dump.readline()
    num_atoms = int(dump.readline().split()[0])
    for j in range(5): dump.readline()

    coords_list = []
    done = False
    i = 0

    while not done:

      coords_list.append([])

      for j in range(num_atoms):
        # FIXME: handle dump format other than id type mol x y z
        tokens = dump.readline().split()
        tag = int(tokens[0])
        type = int(tokens[1])
        mol = int(tokens[2])
        x,y,z = float(tokens[3]),float(tokens[4]),float(tokens[5])
        coords_list[i].append([tag,x,y,z])

      coords_list[i].sort(key=lambda atom:atom[0])
      i += 1
      if dump.readline():
        for j in range(8): dump.readline()
      else:
        done = True

    coords = array(coords_list, dtype=float64)[:,:,1:]
    dump.close()

    if model.num_atoms != num_atoms:
        raise UserError("Specified structure has %d atoms"
            " whereas the coordinates are for %d atoms" % (model.num_atoms, num_atoms))
    start, step, end = process_limit_args(session, start, step, end, len(coords))
    if start > 0 or step > 1 or end < len(coords):
        coords = coords[start:end:step]

    model.add_coordsets(coords, replace=replace)
    return len(coords)

def process_limit_args(session, start, step, end, num_coords):
    if end is None:
        end = num_coords
    if end < start:
        raise UserError("'start' (%d) must be less than or equal to 'end' (%d)" % (start, end))
    start -= 1
    if start > 0 or step > 1 or end < num_coords:
        session.logger.info("start: %d, step: %d, end %d" % (start+1, step, end))
    return start, step, end
