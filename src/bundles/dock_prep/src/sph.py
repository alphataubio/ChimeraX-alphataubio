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

from chimerax.core.models import Model, Surface
from chimerax.atomic import AtomicStructure, Atom
from chimerax.core.colors import random_colors
from chimerax.core.errors import UserError
import numpy as np

from chimerax.atomic import Atom
from chimerax.core.errors import UserError
from sklearn.cluster import DBSCAN

from chimerax.atomic import Structure, Residue, Element
from chimerax.core.errors import UserError

def open_sph(session, file_name, **kw):
    """
    Load a DOCK 6 SPH file into ChimeraX.
    """
    spheres = []
    
    with open(file_name, 'r') as f:
        lines = f.readlines()

    if not lines:
        raise UserError(f"SPH file {file_name} is empty.")

    if lines[0].strip().startswith("cluster"):
        # Format 1: Clustered Spheres
        cluster_id = None
        for line in lines:
            parts = line.split()
            if "cluster" in line:
                cluster_id = int(parts[1])
            elif len(parts) >= 5:
                sphere_id = int(parts[0])
                x, y, z, radius = map(float, parts[1:5])
                spheres.append((sphere_id, x, y, z, radius, cluster_id))
    else:
        # Format 2: Simple Spheres
        for line in lines[1:]:  # Skip header
            parts = line.split()
            if len(parts) >= 4:
                cluster_id = int(parts[0])
                x, y, z, radius = map(float, parts[1:5])
                spheres.append((cluster_id, x, y, z, radius, cluster_id))

    if not spheres:
        raise UserError(f"No spheres found in {file_name}")

    # **Create a new structure model for the spheres**
    model = Structure(session, name=f"DOCK 6 Spheres ({len(set(c[5] for c in spheres))} clusters)")

    # **Track clusters as separate residues**
    cluster_residues = {}

    for sphere_id, x, y, z, radius, cluster_id in spheres:
        if cluster_id not in cluster_residues:
            cluster_residues[cluster_id] = model.new_residue("CLUSTER", "A", cluster_id)

        residue = cluster_residues[cluster_id]
        atom = model.new_atom(f"S{sphere_id}", Element.get_element(1))  # Use Hydrogen as placeholder
        atom.coord = (x, y, z)
        atom.radius = radius
        residue.add_atom(atom)

    session.logger.info(f"Loaded {len(spheres)} spheres from {file_name}")

    # **FIX: DO NOT call `session.models.add([model])` here**
    # Let ChimeraX handle adding it automatically

    # **Return model correctly to avoid duplication**
    return [model], f"Loaded {len(spheres)} spheres from {file_name}"













import numpy as np
from sklearn.cluster import DBSCAN
from chimerax.core.models import Surface
#from chimerax.surface import calculate_surface
#from chimerax.geometry import sphere_fit

def save_sph(session, filename, receptor_model, probe_radius=1.4, dbscan_eps=1.4, dbscan_min_pts=5):
    """
    Generate DOCK-compatible spheres using ChimeraX SES and DBSCAN clustering.

    Parameters:
    - session: ChimeraX session object
    - filename: Path to save the .sph file
    - receptor_model: ChimeraX molecular model
    - probe_radius: Probe radius for SES calculation (default: 1.4 Å)
    - dbscan_eps: DBSCAN clustering radius (default: 1.4 Å)
    - dbscan_min_pts: Minimum points per DBSCAN cluster (default: 5)
    """

    logger = session.logger  # Use ChimeraX session logger

    if not receptor_model:
        logger.error("Invalid receptor model: None provided.")
        return

    # Ensure the input is an atomic model
    if not hasattr(receptor_model, 'atoms'):
        logger.error("Provided receptor model is not an atomic structure.")
        return

    # Generate Solvent-Excluded Surface (SES)
    ses_surface = calculate_surface(session, receptor_model, method="ses", probe_radius=probe_radius)

    if ses_surface is None or not hasattr(ses_surface, 'vertices'):
        logger.error("Failed to generate SES. No surface model found.")
        return

    # Extract SES vertex coordinates
    vertices = np.array([v.xyz for v in ses_surface.vertices])

    if len(vertices) == 0:
        logger.error("No SES vertices found. Check receptor model and SES parameters.")
        return

    # Run DBSCAN clustering
    dbscan = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_pts)
    labels = dbscan.fit_predict(vertices)

    # Identify the largest clusters
    unique_labels, counts = np.unique(labels[labels >= 0], return_counts=True)
    largest_clusters = unique_labels[np.argsort(-counts)]  # Sort clusters by size

    spheres = []

    # Fit spheres to clusters
    for cluster_id in largest_clusters:
        cluster_points = vertices[labels == cluster_id]
        center, radius = sphere_fit(cluster_points)

        # Apply DOCK sphgen-like filtering
        if 1.4 <= radius <= 2.5:  # Typical DOCK sphgen radius range
            spheres.append((center, radius))

    # Save spheres in DOCK-compatible .sph format
    with open(filename, "w") as f:
        f.write("DOCK spheres generated using ChimeraX SES + DBSCAN\n")
        f.write(f"{len(spheres)}\n")
        for i, (center, radius) in enumerate(spheres):
            f.write(f"{i + 1} {center[0]:.3f} {center[1]:.3f} {center[2]:.3f} {radius:.3f}\n")

    logger.info(f"Saved {len(spheres)} spheres to {filename}.")
