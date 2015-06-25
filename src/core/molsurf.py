# vi: set expandtab shiftwidth=4 softtabstop=4:
"""
molsurf -- Compute molecular surfaces
=====================================
"""

from .generic3d import Generic3DModel
from . import cli

class MolecularSurface(Generic3DModel):

    def __init__(self, enclose_atoms, show_atoms, probe_radius, grid_spacing,
                 name, color, visible_patches, sharp_boundaries):
        
        Generic3DModel.__init__(self, name)

        self.atoms = enclose_atoms
        self.show_atoms = show_atoms	# Atoms for surface patch to show
        self.probe_radius = probe_radius
        self.grid_spacing = grid_spacing
        self.name = name
        self.color = color
        self.visible_patches = visible_patches
        self.sharp_boundaries = sharp_boundaries

        self._vertex_to_atom = None
        self._max_radius = None
        self._atom_colors = None

    def new_parameters(self, show_atoms, probe_radius, grid_spacing,
                       visible_patches, sharp_boundaries,
                       color = None, transparency = None):

        shape_change = (probe_radius != self.probe_radius or
                        grid_spacing != self.grid_spacing or
                        sharp_boundaries != self.sharp_boundaries)
        shown_changed = (show_atoms.hash() != self.show_atoms.hash() or
                         visible_patches != self.visible_patches)

        self.show_atoms = show_atoms	# Atoms for surface patch to show
        self.probe_radius = probe_radius
        self.grid_spacing = grid_spacing
        self.visible_patches = visible_patches
        self.sharp_boundaries = sharp_boundaries

        if shape_change:
            self.vertices = None
            self.normals = None
            self.triangles = None
            self.preserve_colors()
            self.vertex_colors = None
            self._vertex_to_atom = None
            self._max_radius = None
        elif shown_changed:
            self.triangle_mask = self.calc_triangle_mask()
                
    @property
    def atom_count(self):
        return len(self.atoms)

    def calculate_surface_geometry(self):
        if not self.vertices is None:
            return              # Geometry already computed
        atoms = self.atoms
        xyz = atoms.coords
        r = atoms.radii
        self._max_radius = r.max()
        from .surface import ses_surface_geometry
        va, na, ta = ses_surface_geometry(xyz, r, self.probe_radius, self.grid_spacing)
        if self.sharp_boundaries:
            from .surface import sharp_edge_patches
            va, na, ta, v2a = sharp_edge_patches(va, na, ta, self.vertex_to_atom_map(va), xyz)
            self._vertex_to_atom = v2a
#        for i in range(3):
#            vsa, nsa, tsa, v2a = sharp_edge_patches(vsa, nsa, tsa, v2a, xyz)
        self.vertices = va
        self.normals = na
        self.triangles = ta
        self.triangle_mask = self.calc_triangle_mask()
        self.restore_colors()

    def calc_triangle_mask(self):
        tmask = self.patch_display_mask(self.show_atoms)
        if self.visible_patches is None:
            return tmask

        from . import surface
        if self.sharp_boundaries:
            # With sharp boundaries triangles are not connected.
            vmap = surface.unique_vertex_map(self.vertices)
            tri = vmap[self.triangles]
        else:
            tri = self.triangles
        m = surface.largest_blobs_triangle_mask(self.vertices, tri, tmask,
                                                blob_count = self.visible_patches,
                                                rank_metric = 'area rank')
        return m

    def vertex_to_atom_map(self, vertices = None):
        v2a = self._vertex_to_atom
        if v2a is None:
            xyz1 = self.vertices if vertices is None else vertices
            xyz2 = self.atoms.coords
            max_dist = 1.1 * (self.probe_radius + self._max_radius)
            from . import geometry
            i1, i2, nearest1 = geometry.find_closest_points(xyz1, xyz2, max_dist)
            from numpy import empty, int32
            v2a = empty((len(xyz1),), int32)
            v2a[i1] = nearest1
            self._vertex_to_atom = v2a
        return v2a

    def patch_display_mask(self, patch_atoms):
        surf_atoms = self.atoms
        if len(patch_atoms) == len(surf_atoms):
            return None
        shown_atoms = surf_atoms.mask(patch_atoms)
        return self.atom_triangle_mask(shown_atoms)

    def atom_triangle_mask(self, atom_mask):
        v2a = self.vertex_to_atom_map()
        shown_vertices = atom_mask[v2a]
        t = self.triangles
        from numpy import logical_and, empty, bool
        shown_triangles = empty((len(t),), bool)
        logical_and(shown_vertices[t[:,0]], shown_vertices[t[:,1]], shown_triangles)
        logical_and(shown_triangles, shown_vertices[t[:,2]], shown_triangles)
        return shown_triangles

    def preserve_colors(self):
        vc = self.vertex_colors
        if vc is None:
            return

        # Preserve surface colors per atom.
        v2a = self.vertex_to_atom_map()
        from numpy import empty, uint8
        acolor = empty((len(self.atoms),4), uint8)
        acolor[:] = self.color
        acolor[v2a] = vc
        self._atom_colors = acolor

    def restore_colors(self):
        acolor = self._atom_colors
        if acolor is None:
            return

        v2a = self.vertex_to_atom_map()
        self.vertex_colors = acolor[v2a,:]
        self._atom_colors = None

    def first_intercept(self, mxyz1, mxyz2, exclude = None):
        # Pick atom associated with surface patch
        from .graphics import Drawing
        p = Drawing.first_intercept(self, mxyz1, mxyz2, exclude)
        if p is None:
            return None
        t = p.triangle_number
        v = self.triangles[t,0]
        v2a = self.vertex_to_atom_map()
        a = v2a[v]
        atom = self.atoms[a]
        from .structure import PickedAtom
        pa = PickedAtom(atom, p.distance)
        return pa

    def update_selection(self):
        asel = self.atoms.selected
        tmask = self.atom_triangle_mask(asel)
        self.selected_triangles_mask = tmask

# -------------------------------------------------------------------------------------
#
def surface_command(session, atoms = None, enclose = None, include = None,
                    probe_radius = 1.4, grid_spacing = 0.5,
                    color = None, transparency = None, visible_patches = None,
                    sharp_boundaries = True,
                    nthread = None, replace = True, hide = False, close = False):
    '''
    Compute and display solvent excluded molecular surfaces.
    '''
    atoms = check_atoms(atoms, session) # Warn if no atoms specifed

    if close:
        close_surfaces(atoms, session.models)
        return []

    if replace:
        all_surfs = dict((s.atoms.hash(), s) for s in session.models.list(type = MolecularSurface))
    else:
        all_surfs = {}
    
    surfs = []
    new_surfs = []
    if enclose is None:
        atoms, all_small = remove_solvent_ligands_ions(atoms, include)
        for m, chain_id, show_atoms in atoms.by_chain:
            if all_small:
                enclose_atoms = atoms.filter(atoms.chain_ids == chain_id)
            else:
                matoms = m.atoms
                chain_atoms = matoms.filter(matoms.chain_ids == chain_id)
                enclose_atoms = remove_solvent_ligands_ions(chain_atoms, include)[0]
            s = all_surfs.get(enclose_atoms.hash())
            if s is None:
                name = '%s_%s SES surface' % (m.name, chain_id)
                rgba = surface_rgba(color, transparency, chain_id)
                s = MolecularSurface(enclose_atoms, show_atoms, probe_radius, grid_spacing,
                                     name, rgba, visible_patches, sharp_boundaries)
                new_surfs.append((s,m))
            else:
                s.new_parameters(show_atoms, probe_radius, grid_spacing,
                                 visible_patches, sharp_boundaries)
                update_color(s, color, transparency)
            surfs.append(s)
    else:
        enclose_atoms, eall_small = remove_solvent_ligands_ions(enclose, include)
        show_atoms = enclose_atoms if atoms is None else atoms.intersect(enclose_atoms)
        s = all_surfs.get(enclose_atoms.hash())
        if s is None:
            mols = enclose.unique_structures
            parent = mols[0] if len(mols) == 1 else session.models.drawing
            name = 'Surface %s' % enclose.spec
            rgba = surface_rgba(color, transparency)
            s = MolecularSurface(enclose_atoms, show_atoms, probe_radius, grid_spacing,
                                 name, rgba, visible_patches, sharp_boundaries)
            new_surfs.append((s,parent))
        else:
            s.new_parameters(show_atoms, probe_radius, grid_spacing,
                             visible_patches, sharp_boundaries)
            update_color(s, color, transparency)
        surfs.append(s)

    # Close overlapping surfaces.
    if replace:
        other_surfs = set(all_surfs.values()) - set(surfs)
        from .molecule import concatenate
        surf_atoms = concatenate([s.atoms for s in surfs])
        osurfs = surfaces_overlapping_atoms(other_surfs, surf_atoms)
        if osurfs:
            session.models.close(osurfs)

    # Compute surfaces using multiple threads
    args = [(s,) for s in surfs]
    args.sort(key = lambda s: s[0].atom_count, reverse = True)      # Largest first for load balancing
    from . import threadq
    threadq.apply_to_list(lambda s: s.calculate_surface_geometry(), args, nthread)

    # Add new surfaces to open models list.
    for s, parent in new_surfs:
        session.models.add([s], parent = parent)

    # Make sure replaced surfaces are displayed.
    for s in surfs:
        s.display = True

    # Hide some surfaces
    if hide:
        hide_surfaces(atoms, session.models)

    return surfs

def register_surface_command():
    from .structure import AtomsArg
    from . import cli, color
    _surface_desc = cli.CmdDesc(
        optional = [('atoms', AtomsArg)],
        keyword = [('enclose', AtomsArg),
                   ('include', AtomsArg),
                   ('probe_radius', cli.FloatArg),
                   ('grid_spacing', cli.FloatArg),
                   ('color', color.ColorArg),
                   ('transparency', cli.FloatArg),
                   ('visible_patches', cli.IntArg),
                   ('sharp_boundaries', cli.BoolArg),
                   ('nthread', cli.IntArg),
                   ('replace', cli.BoolArg),
                   ('hide', cli.NoArg),
                   ('close', cli.NoArg)],
        synopsis = 'create molecular surface')
    cli.register('surface', _surface_desc, surface_command)

def check_atoms(atoms, session):
    if atoms is None:
        from .structure import all_atoms
        atoms = all_atoms(session)
        if len(atoms) == 0:
            raise cli.AnnotationError('No atomic models open.')
        atoms.spec = 'all atoms'
    elif len(atoms) == 0:
        raise cli.AnnotationError('No atoms specified by %s' % (atoms.spec,))
    return atoms

def remove_solvent_ligands_ions(atoms, keep = None):
    '''Remove solvent, ligands and ions unless that removes all atoms
    in which case don't remove any.'''
    # TODO: Properly identify solvent, ligands and ions.
    # Currently simply remove every atom is does not belong to a chain.
    fatoms = atoms.filter(atoms.in_chains)
    if keep:
        fatoms = fatoms.merge(atoms.intersect(keep))
    all_small = (len(fatoms) == 0)
    if all_small:
        return atoms, all_small
    return fatoms, all_small

def surface_rgba(color, transparency, chain_id = None):
    if color is None:
        if chain_id is None:
            from numpy import array, uint8
            rgba8 = array((180,180,180,255), uint8)
        else:
            from . import structure
            rgba8 = structure.chain_rgba8(chain_id)
    else:
        rgba8 = color.uint8x4()
    if not transparency is None:
        opacity = int(255*(100.0-transparency)/100.0)
        rgba8[3] = opacity
    return rgba8

def update_color(surf, color, transparency):
    if color is None:
        if not transparency is None:
            opacity = int(255*(100.0-transparency)/100.0)
            vcolors = surf.vertex_colors
            if vcolors is None:
                rgba = surf.color
                rgba[3] = opacity
                surf.color = rgba
            else:
                vcolors[:,3] = opacity
                surf.vertex_colors = vcolors
    else:
        rgba = color.uint8x4()
        if not transparency is None:
            opacity = int(255*(100.0-transparency)/100.0)
            rgba[3] = opacity
        surf.color = rgba
        surf.vertex_colors = None
        print ('set surf color', rgba)

def surfaces_overlapping_atoms(surfs, atoms):
    si = atoms.intersects_each([s.atoms for s in surfs])
    osurfs = [s for s,i in zip(surfs,si) if i]
    return osurfs

def show_surface(name, va, na, ta, color = (180,180,180,255)):
    surf = MolecularSurface(name)
    surf.geometry = va, ta
    surf.normals = na
    surf.color = color
    return surf

def molecule_surface(mol, probe_radius = 1.4, grid_spacing = 0.5):
    a = mol.atoms
    from numpy import array
    a = a.filter(array(a.residues.names) != 'HOH')     # exclude waters
    xyz, r = a.coords, a.radii
    from .surface import ses_surface_geometry
    va, na, ta = ses_surface_geometry(xyz, r, probe_radius, grid_spacing)
    from numpy import array, uint8
    color = array((180,180,180,255), uint8)
    surf = show_surface(mol.name + ' surface', va, na, ta, color, mol.position)
    return surf

def surfaces_with_atoms(atoms, models):
    surfs = []
    for m in list(atoms.unique_structures) + [models.drawing]:
        for s in m.child_drawings():
            if isinstance(s, MolecularSurface):
                if len(atoms.intersect(s.atoms)) > 0:
                    surfs.append(s)
    return surfs

def hide_surfaces(atoms, models):
    for s in surfaces_with_atoms(atoms, models):
        s.display = False

def close_surfaces(atoms, models):
    surfs = surfaces_with_atoms(atoms, models)
    if surfs:
        models.close(surfs)

def sasa_command(session, atoms = None, probe_radius = 1.4):
    '''
    Compute solvent accessible surface area.
    Only the specified atoms are considered.
    '''
    atoms = check_atoms(atoms, session)
    r = atoms.radii
    r += probe_radius
    from . import surface
    areas = surface.spheres_surface_area(atoms.scene_coords, r)
    area = areas.sum()
    msg = 'Solvent accessible area for %s = %.5g' % (atoms.spec, area)
    log = session.logger
    log.info(msg)
    log.status(msg)

def register_sasa_command():
    from .structure import AtomsArg
    _sasa_desc = cli.CmdDesc(
        optional = [('atoms', AtomsArg)],
        keyword = [('probe_radius', cli.FloatArg),],
        synopsis = 'compute solvent accessible surface area')
    cli.register('sasa', _sasa_desc, sasa_command)

def buriedarea_command(session, atoms1, with_atoms2 = None, probe_radius = 1.4):
    '''
    Compute solvent accessible surface area.
    Only the specified atoms are considered.
    '''
    if with_atoms2 is None:
        raise cli.AnnotationError('Require "with" keyword: buriedarea #1 with #2')
    atoms2 = with_atoms2

    ni = len(atoms1.intersect(atoms2))
    if ni > 0:
        raise cli.AnnotationError('Two sets of atoms must be disjoint, got %d atoms in %s and %s'
                                  % (ni, atoms1.spec, atoms2.spec))

    ba, a1a, a2a, a12a = buried_area(atoms1, atoms2, probe_radius)

    # Report result
    msg = 'Buried area between %s and %s = %.5g' % (atoms1.spec, atoms2.spec, ba)
    log = session.logger
    log.status(msg)
    msg += ('\n  area %s = %.5g, area %s = %.5g, area both = %.5g'
            % (atoms1.spec, a1a, atoms2.spec, a2a, a12a))
    log.info(msg)

def buried_area(a1, a2, probe_radius):
    from .surface import spheres_surface_area
    xyz1, r1 = atom_spheres(a1, probe_radius)
    a1a = spheres_surface_area(xyz1, r1).sum()
    xyz2, r2 = atom_spheres(a2, probe_radius)
    a2a = spheres_surface_area(xyz2, r2).sum()
    from numpy import concatenate
    xyz12, r12 = concatenate((xyz1,xyz2)), concatenate((r1,r2))
    a12a = spheres_surface_area(xyz12, r12).sum()
    ba = 0.5 * (a1a + a2a - a12a)
    return ba, a1a, a2a, a12a

def atom_spheres(atoms, probe_radius = 1.4):
    xyz = atoms.coords
    r = atoms.radii.copy()
    r += probe_radius
    return xyz, r

def register_buriedarea_command():
    from .structure import AtomsArg
    _buriedarea_desc = cli.CmdDesc(
        required = [('atoms1', AtomsArg)],
        keyword = [('with_atoms2', AtomsArg),
                   ('probe_radius', cli.FloatArg),],
        synopsis = 'compute buried area')
    cli.register('buriedarea', _buriedarea_desc, buriedarea_command)

def scolor(session, atoms = None, color = None, byatom = False):
    surfs = session.models.list(type = MolecularSurface)
    if not atoms is None:
        surfs = [s for s in surfs if s.atoms.intersects(atoms)]

    for s in surfs:
        nv = len(s.vertices)
        vcolors = s.vertex_colors
        if vcolors is None:
            from numpy import empty, uint8
            vcolors = empty((nv,4), uint8)
            vcolors[:] = s.color
        if atoms is None:
            v = slice(nv)
        else:
            ai = s.atoms.mask(atoms)
            v2a = s.vertex_to_atom_map()
            v = ai[v2a]
        if byatom:
            v2a = s.vertex_to_atom_map()
            vcolors[v] = s.atoms.colors[v2a[v],:]
        elif not color is None:
            vcolors[v] = color.uint8x4()
        s.vertex_colors = vcolors

def register_scolor_command():
    from .structure import AtomsArg
    from . import cli, color
    _scolor_desc = cli.CmdDesc(
        optional = [('atoms', cli.Or(AtomsArg, cli.EmptyArg)),
                    ('color', color.ColorArg),],
        keyword = [('byatom', cli.NoArg)],
        synopsis = 'color surfaces')
    cli.register('scolor', _scolor_desc, scolor)
