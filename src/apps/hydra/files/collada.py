def read_collada_surfaces(path, session, color = (178,178,178,255)):

    from os.path import basename
    from ..models import Model
    s = Model(basename(path))

    from collada import Collada
    c = Collada(path)
    from ..geometry.place import Place
    splist = surface_pieces_from_nodes(c.scene.nodes, s, color, Place(), {})
    set_positions_and_colors(splist)

    ai = c.assetInfo
    if ai:
        s.collada_unit_name = ai.unitname
        s.collada_contributors = ai.contributors

    return s

def surface_pieces_from_nodes(nodes, surf, color, place, ginst):

    # TODO: Copy collada hierarchy instead of flattening it.
    #       Code was originally written when only 2-level hierarchy was supported.
    splist = []
    from collada.scene import GeometryNode, Node
    from ..geometry.place import Place
    for n in nodes:
        if isinstance(n, GeometryNode):
            s2m = dict((m.symbol,m.target) for m in n.materials)
            g = n.geometry
            if g.id in ginst:
                add_geometry_instance(g.primitives, ginst[g.id], place, color, s2m)
            else:
                ginst[g.id] = spieces = geometry_surface_pieces(g.primitives, place, color, s2m, surf)
                splist.extend(spieces)
        elif isinstance(n, Node):
            pl = place * Place(n.matrix[:3,:])
            splist.extend(surface_pieces_from_nodes(n.children, surf, color, pl, ginst))
    return splist

def geometry_surface_pieces(primitives, place, color, s2m, surf):

    from collada import polylist, triangleset

    splist = []
    for p in primitives:
        if isinstance(p, polylist.Polylist):
            p = p.triangleset()
        if not isinstance(p, triangleset.TriangleSet):
            continue        # Skip line sets.

        t = p.vertex_index            # N by 3 array of vertex indices for triangles
        t = t.copy()                  # array from pycollada is not contiguous.
        v = p.vertex                  # M by 3 array of floats for vertex positions
        ni = p.normal_index           # N by 3 array of normal indices for triangles
        n = p.normal		      # M by 3 array of floats for vertex normals

        # Collada allows different normals on the same vertex in different triangles,
        # but Hydra only allows one normal per vertex.
        from numpy import empty
        vn = empty(v.shape, n.dtype)
        vn[t.ravel(),:] = n[ni.ravel(),:]

        sp = surf.new_drawing()
        sp.geometry = v, t
        sp.normals = vn
        c = material_color(s2m.get(p.material), color)
        sp.color_list = [c]
        sp.position_list = [place]

        splist.append(sp)

    return splist

def add_geometry_instance(primitives, spieces, place, color, s2m):

    for p,sp in zip(primitives, spieces):
        c = material_color(s2m.get(p.material), color)
        sp.position_list.append(place)
        sp.color_list.append(c)

def material_color(material, color):

    if material is None:
        return color
    e = material.effect
    if e is None:
        return color
    c = e.diffuse
    if c is None:
        return color
    c8bit = tuple(int(255*r) for r in c)
    return c8bit

# For drawings with multiple instances make colors a numpy uint8 array.
def set_positions_and_colors(drawings):
    from ..geometry.place import Places
    for d in drawings:
        clist = d.color_list
        from numpy import array, uint8
        d.colors = array(clist, uint8).reshape((len(clist),4))
        d.positions = Places(d.position_list)
