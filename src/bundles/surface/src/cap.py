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

def update_clip_caps(view):
    from chimerax.core.core_settings import settings
    if not settings.clipping_surface_caps:
        return
    cp = view.clip_planes
    planes = cp.planes()
    cpos = view.camera.position
    for p in planes:
        p.update_direction(cpos)
    update = (cp.changed or (view.shape_changed and planes))
    # TODO: Update caps only on specific drawings whose shape changed.
    if update:
        drawings = view.drawing.all_drawings()
        show_surface_clip_caps(planes, drawings, offset = settings.clipping_cap_offset)
        cp.changed = False
        view.redraw_needed = True

def show_surface_clip_caps(planes, drawings, offset = 0.01):
    for p in planes:
        for d in drawings:
            # Clip only drawings that have "clip_cap" attribute true.
            if (hasattr(d, 'clip_cap') and d.clip_cap and
                d.triangles is not None and not hasattr(d, 'clip_cap_owner')):
                varray, narray, tarray = compute_cap(d, p, offset)
                set_cap_drawing_geometry(d, p.name, varray, narray, tarray)

    # Remove caps for clip planes that are gone.
    plane_names = set(p.name for p in planes)
    for cap in drawings:
        if hasattr(cap, 'clip_cap_owner') and cap.clip_plane_name not in plane_names:
            d = cap.clip_cap_owner
            del d._clip_cap_drawings[cap.clip_plane_name]
            from chimerax.core.models import Model
            if isinstance(cap, Model):
                cap.session.models.remove([cap])
            else:
                cap.parent.remove_drawing(cap)

def remove_clip_caps(drawings):
    for cap in drawings:
        if hasattr(cap, 'clip_cap_owner'):
            d = cap.clip_cap_owner
            del d._clip_cap_drawings[cap.clip_plane_name]
            from chimerax.core.models import Model
            if isinstance(cap, Model):
                cap.session.models.remove([cap])
            else:
                cap.parent.remove_drawing(cap)

def compute_cap(drawing, plane, offset):
    # Undisplay cap for drawing with no geometry shown.
    d = drawing
    if (not d.display or
        (d.triangle_mask is not None and
         d.triangle_mask.sum() < len(d.triangle_mask))):
        return None, None, None

    # Handle surfaces with duplicate vertices, such as molecular
    # surfaces with sharp edges between atoms.
    if hasattr(d, 'joined_triangles'):
        t = d.joined_triangles
    else:
        t = d.triangles

    # Compute cap geometry.
    # TODO: Cap instances
    if len(d.positions) > 1:
        varray, tarray, pnormal = compute_instances_cap(d, t, plane, offset)
    else:
        dp = d.scene_position.inverse()
        pnormal = dp.apply_without_translation(plane.normal)
        from chimerax.core.geometry import inner_product
        poffset = inner_product(pnormal, dp*plane.plane_point) + offset + getattr(d, 'clip_offset', 0)
        from . import compute_cap
        varray, tarray = compute_cap(pnormal, poffset, d.vertices, t)

    if tarray is None or len(tarray) == 0:
        return None, None, None
    narray = varray.copy()
    narray[:] = pnormal

    return varray, narray, tarray

def compute_instances_cap(drawing, triangles, plane, offset):
    d = drawing
    doffset = offset + getattr(d, 'clip_offset', 0)
    geom = []
    # TODO: Handle two hierarchy levels of instancing.
    pp = drawing.parent.scene_position.inverse()
    parent_ppoint = pp*plane.plane_point
    parent_pnormal = pp.apply_without_translation(plane.normal)

    # TODO: Optimize by testing if plane intercepts bounding sphere.
    b = d.bounds(positions = False)
    if b is None:
        return None, None, None
        
    dpos = d.positions.masked(d.display_positions)
    ipos = box_positions_intersecting_plane(dpos, b, parent_ppoint, parent_pnormal)
    if len(ipos) == 0:
        return None, None, None
    for pos in ipos:
        pinv = pos.inverse()
        pnormal = pinv.apply_without_translation(parent_pnormal)
        from chimerax.core.geometry import inner_product
        poffset = inner_product(pnormal, pinv*parent_ppoint) + doffset
        from . import compute_cap
        ivarray, itarray = compute_cap(pnormal, poffset, d.vertices, triangles)
        pos.move(ivarray)
        geom.append((ivarray, itarray))
    varray, tarray = concatenate_geometry(geom)
    return varray, tarray, parent_pnormal

def box_positions_intersecting_plane(positions, b, origin, normal):
    c, r = b.center(), b.radius()
    pc = positions * c
    pc -= origin
    from numpy import dot, abs
    dist = abs(dot(pc,normal))
    bint = (dist <= r)
    ipos = positions.masked(bint)
    return ipos

def concatenate_geometry(geom):
    from numpy import concatenate
    varray = concatenate(tuple(v for v,t in geom))
    tarray = concatenate(tuple(t for v,t in geom))
    voffset = ts = 0
    for v,t in geom:
        nt = len(t)
        tarray[ts:ts+nt,:] += voffset
        ts += nt
        voffset += len(v)
    return varray, tarray

def set_cap_drawing_geometry(drawing, plane_name, varray, narray, tarray):
    d = drawing
    # Set cap drawing geometry.
    if not hasattr(d, '_clip_cap_drawings'):
        d._clip_cap_drawings = {}
    mcap = d._clip_cap_drawings.get(plane_name, None)     # Find cap drawing
    if mcap and mcap.was_deleted:
        mcap = None

    if mcap:
        cm = mcap
    elif varray is None:
        return
    else:
        cap_name = 'cap ' + plane_name
        if len(d.positions) == 1:
            cm = new_cap(d, cap_name)
        else:
            cm = new_cap(d.parent, cap_name)
            cm.pickable = False	  # Don't want pick of one cap to pick all instance caps.
        cm.clip_plane_name = plane_name
        cm.clip_cap_owner = d
        d._clip_cap_drawings[plane_name] = cm
        cm.color = d.color

    cm.set_geometry(varray, narray, tarray)

def new_cap(drawing, cap_name):
    from chimerax.core.models import Model, Surface
    if isinstance(drawing, Model):
        # Make cap a model when capping a model so color can be set by command.
        c = Surface(cap_name, drawing.session)
        c.SESSION_SAVE = False
        drawing.add([c])
    else:
        c = drawing.new_drawing(cap_name)
    return c