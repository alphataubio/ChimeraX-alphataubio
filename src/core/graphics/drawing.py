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

'''
Drawing
=======
'''


class Drawing:
    '''
    A Drawing represents a tree of objects each consisting of a set of
    triangles in 3 dimensional space.  Drawings are used to draw molecules,
    density maps, geometric shapes and other models.  A Drawing has a name,
    a unique id number which is a positive integer, it can be displayed
    or hidden, has a placement in space, or multiple copies can be placed
    in space, and a drawing can be selected.  The coordinates, colors,
    normal vectors and other geometric and display properties are managed
    by the Drawing objects.

    A drawing can have child drawings.  The purpose of child drawings is
    for convenience in adding, removing, displaying and selecting parts
    of a scene. Child drawings are created by the new_drawing() method.

    Multiple copies of a drawing be drawn with specified positions and
    colors. Copy positions can be specified by a shift and scale factor but
    no rotation, useful for copies of spheres.  Each copy can be displayed
    or hidden, selected or unselected.

    The basic data defining the triangles is an N by 3 array of vertices
    (float32 numpy array) and an array that defines triangles as 3 integer
    index values into the vertex array which define the 3 corners of
    a triangle.  The triangle array is of shape T by 3 and is a numpy
    int32 array.  The filled triangles or a mesh consisting of just the
    triangle edges can be shown.  The vertices can be individually colored
    with linear interpolation of colors across triangles, or all triangles
    can be given the same color, or a 2-dimensional texture can be used to
    color the triangles with texture coordinates assigned to the vertices.
    Transparency values can be assigned to the vertices. Individual
    triangles or triangle edges in mesh display style can be hidden.
    An N by 3 float array gives normal vectors, one normal per vertex,
    for lighting calculations.

    Rendering of drawings is done with OpenGL.
    '''
    def __init__(self, name):

        self._redraw_needed = None
        self._attribute_changes = set()         # Attribute names of changed data for updating buffers

        self.name = name
        "Name of this drawing."

        from ..geometry import Places
        # Copies of drawing are placed at these positions:
        self._positions = Places()
        from numpy import array, uint8
        # Colors for each position, N by 4 uint8 numpy array:
        self._colors = array(((178, 178, 178, 255),), uint8)
        self._opaque_color_count = 1
        # bool numpy array, show only some positions:
        self._displayed_positions = None
        self._any_displayed_positions = True
        # bool numpy array, selected positions:
        self._selected_positions = None
        self._selected_triangles_mask = None  # bool numpy array
        self._child_drawings = []

        self._cached_geometry_bounds = None
        self._cached_position_bounds = None

        # Geometry and colors
        self._vertices = None		# N x 3 float32 numpy array
        self._triangles = None		# N x 3 int32 numpy array
        self._normals = None		# N x 3 float32 numpy array

        self._vertex_colors = None
        self._opaque_vertex_color_count = 0
        self.auto_recolor_vertices = None	# Function to call when geometry changes

        self._triangle_mask = None
        '''
        A mask that allows hiding some triangles, a numpy array of
        length M (# of triangles) of type bool. This is used for
        showing just a patch of a surface.
        '''
        self.auto_remask_triangles = None	# Function to call when geometry changes

        self._edge_mask = None
        '''
        A mask that allows hiding some edges, a numpy array of length M
        (# of triangles) of type uint8, where bits 0, 1, and 2 are whether
        to display each edge of the triangle.  This is used for
        square mesh density map display.
        '''

        self.display_style = self.Solid
        '''
        Display style can be Drawing.Solid, Drawing.Mesh or Drawing.Dot.
        Only one style can be used for a single Drawing instance.
        '''

        self.texture = None
        '''
        Texture to use in coloring the surface, a graphics.Texture object.
        Only 2-dimensional textures are supported.  Can be None.
        '''

        self.multitexture = None
        '''
        List of N textures to use, each applying to 1/N of the triangles.
        This is used for volumetric rendering by texturing a stack of rectangles.
        Only 2-dimensional textures are supported.  Can be None.
        '''
        self.multitexture_reverse_order = False
        """Whether to draw multitextured geometry in reverse order for handling transparency.
        Used by grayscale rendering for depth ordering."""

        self.texture_coordinates = None
        """Texture coordinates, an N by 2 numpy array of float32 values
        in range 0-1"""

        # 3d texture that modulates colors.
        self.ambient_texture = None
        '''
        A 3-dimensional texture that modulates the brightness of surface
        vertex colors.  Used for fast rendering of ambient occlusion
        lighting.
        '''

        # Drawing to texture coordinates.
        self.ambient_texture_transform = None
        """Transformation mapping vertex coordinates to ambient_texture
        coordinates, a geometry.Place object."""

        self.opaque_texture = False
        "Whether the texture for surface coloring is opaque or transparent."

        self.use_lighting = True
        """Whether to use lighting when rendering.  If false then a flat
        unshaded color will be shown."""

        # OpenGL drawing
        self._draw_shape = None
        self._draw_selection = None
        self._shader_opt = None                 # Cached shader options
        self._vertex_buffers = []               # Buffers used by both drawing and selection
        self._opengl_context = None		# For deleting buffers, to make context current

        self.was_deleted = False
        "Indicates whether this Drawing has been deleted."

    def redraw_needed(self, **kw):
        """Function called when the drawing has been changed to indicate
        that the graphics needs to be redrawn."""
        rn = self._redraw_needed
        if rn is not None:
            rn(self, **kw)

    @property
    def vertices(self):
        '''
        Vertices of the rendered geometry, a numpy N by 3 array of float32 values.
        Read-only. Set using set_geometry() method.
        '''
        return self._vertices

    @property
    def normals(self):
        '''
        Normal vectors of the rendered geometry, a numpy N by 3 array of float32 values.
        Read-only. Set using set_geometry() method.
        '''
        return self._normals

    @property
    def triangles(self):
        '''
        Vertex indices for the corners of each triangle making up the
        rendered geometry, a numpy M by 3 array of int32 values.
        Read-only.  Set using set_geometry() method.
        '''
        return self._triangles
    
    def _get_shape_changed(self):
        rn = self._redraw_needed
        return rn.shape_changed if rn else False
    def _set_shape_changed(self, changed):
        self.redraw_needed(shape_changed = True)
    shape_changed = property(_get_shape_changed, _set_shape_changed)
    '''Did this drawing or any drawing in the same tree change shape since the last redraw.'''
    
    def __setattr__(self, key, value):
        if key in self._effects_shader:
            self._shader_opt = None       # Cause shader update
            self.redraw_needed()
        if key in self._effects_buffers:
            self._attribute_changes.add(key)
            sc = key in ('_vertices', '_triangles', '_triangle_mask')
            if sc:
                self._cached_geometry_bounds = None
                self._cached_position_bounds = None
            else:
                sc = key in ('_displayed_positions', '_positions')
                if sc:
                    self._cached_position_bounds = None
            self.redraw_needed(shape_changed=sc)

        super(Drawing, self).__setattr__(key, value)

    # Display styles
    Solid = 'solid'
    "Display style showing filled triangles."
    Mesh = 'mesh'
    "Display style showing only edges of triangles."
    Dot = 'dot'
    "Display style showing only dots at triangle vertices."

    def child_drawings(self):
        '''Return the list of surface pieces.'''
        return self._child_drawings

    def all_drawings(self, displayed_only = False):
        '''Return all drawings including self and children at all levels.'''
        if displayed_only and not self.display:
            return []
        dlist = [self]
        for d in self.child_drawings():
            dlist.extend(d.all_drawings(displayed_only))
        return dlist

    def new_drawing(self, name=None):
        '''Create a new empty child drawing.'''
        d = Drawing(name)
        self.add_drawing(d)
        return d

    def add_drawing(self, d):
        '''Add a child drawing.'''
        d.set_redraw_callback(self._redraw_needed)
        cd = self._child_drawings
        cd.append(d)
        if hasattr(d, 'parent') and d.parent is not None:
            # Reparent drawing.
            d.parent.remove_drawing(d, delete=False)
        d.parent = self
        if self.display:
            self.redraw_needed(shape_changed=True)

    def remove_drawing(self, d, delete=True):
        '''Remove a specified child drawing.'''
        self._child_drawings.remove(d)
        del d.parent
        if delete:
            d.delete()
        self.redraw_needed(shape_changed=True, selection_changed=True)

    def remove_drawings(self, drawings, delete=True):
        '''Remove specified child drawings.'''
        dset = set(drawings)
        self._child_drawings = [d for d in self._child_drawings
                                if d not in dset]
        for d in drawings:
            del d.parent
        if delete:
            for d in drawings:
                d.delete()
        self.redraw_needed(shape_changed=True, selection_changed=True)

    def remove_all_drawings(self, delete=True):
        '''Remove all child drawings.'''
        cd = self.child_drawings()
        if cd:
            self.remove_drawings(cd, delete)

    @property
    def drawing_lineage(self):
        '''Return a sequence of drawings from the root down to the current drawing.'''
        if hasattr(self, 'parent'):
            return self.parent.drawing_lineage + [self]
        else:
            return [self]

    def set_redraw_callback(self, redraw_needed):
        self._redraw_needed = redraw_needed
        for d in self.child_drawings():
            d.set_redraw_callback(redraw_needed)

    def get_display(self):
        return self._any_displayed_positions and len(self._positions) > 0

    def set_display(self, display):
        dp = self._displayed_positions
        if dp is None:
            from numpy import empty, bool
            dp = empty((len(self._positions),), bool)
        dp[:] = display
        self._displayed_positions = dp		# Need this to trigger buffer update
        self._any_displayed_positions = display
        self.redraw_needed(shape_changed=True)

    display = property(get_display, set_display)
    '''Whether or not the surface is drawn.'''

    def get_display_positions(self):
        return self._displayed_positions

    def set_display_positions(self, position_mask):
        from numpy import array_equal
        dp = self._displayed_positions
        if ((position_mask is None and dp is None) or
            (position_mask is not None and dp is not None and
             position_mask is not dp and array_equal(position_mask, dp))):
            return
        self._displayed_positions = position_mask
        self._any_displayed_positions = (position_mask.sum() > 0)
        self.redraw_needed(shape_changed=True)

    display_positions = property(get_display_positions, set_display_positions)
    '''Mask specifying which copies are displayed. Can be None meaning all positions displayed'''

    @property
    def num_displayed_positions(self):
        dp = self.display_positions
        ndp = len(self.positions) if dp is None else dp.sum()
        return ndp

    @property
    def parents_displayed(self):
        for d in self.drawing_lineage[:-1]:
            if not d.display:
                return False
        return True

    def get_selected(self):
        sp = self._selected_positions
        tmask = self._selected_triangles_mask
        return (((sp is not None) and sp.sum() > 0) or
                ((tmask is not None) and tmask.sum() > 0))

    def set_selected(self, sel):
        if sel:
            sp = self._selected_positions
            if sp is None:
                from numpy import ones, bool
                self._selected_positions = ones(len(self.positions), bool)
            else:
                sp[:] = True
                self._selected_positions = sp # Need to set to track changes
        else:
            self._selected_positions = None
            self._selected_triangles_mask = None
        self.redraw_needed(selection_changed=True)

    selected = property(get_selected, set_selected)
    '''Whether or not the drawing is selected.
    Does not include or effect children.'''

    def get_selected_positions(self):
        return self._selected_positions

    def set_selected_positions(self, spos):
        self._selected_positions = spos
        self.redraw_needed(selection_changed=True)

    selected_positions = property(get_selected_positions,
                                  set_selected_positions)
    '''Mask specifying which drawing positions are selected.
    Does not include or effect children.'''

    def get_selected_triangles_mask(self):
        return self._selected_triangles_mask

    def set_selected_triangles_mask(self, tmask):
        self._selected_triangles_mask = tmask
        self.redraw_needed(selection_changed=True)

    selected_triangles_mask = property(get_selected_triangles_mask,
                                       set_selected_triangles_mask)
    '''Mask specifying which triangles are selected.'''

    def any_part_selected(self):
        '''Is any part of this Drawing or its children selected.'''
        if self.selected:
            return True
        for d in self.child_drawings():
            if d.any_part_selected():
                return True
        return False

    def clear_selection(self, include_children=True):
        '''Unselect this drawing and child drawings in if include_children is True.'''
        self.selected = False
        self.redraw_needed(selection_changed=True)
        if include_children:
            for d in self.child_drawings():
                d.clear_selection()

    def get_position(self):
        return self._positions[0]

    def set_position(self, pos):
        from ..geometry import Places
        self._positions = Places([pos])
        if (not self._displayed_positions is None
            and len(self._displayed_positions) != 1):
            self._displayed_positions = None
        self.redraw_needed(shape_changed=True)

    position = property(get_position, set_position)
    '''Position and orientation of the surface in space.'''

    def _get_scene_position(self):
        from ..geometry import product
        return product([d.position for d in self.drawing_lineage])
    def _set_scene_position(self, pos):
        self.positions = self.positions * (self.scene_position.inverse() * pos)
    scene_position = property(_get_scene_position, _set_scene_position)
    '''Position in scene coordinates.'''

    def get_scene_positions(self, displayed_only = False):
        p = self.get_positions(displayed_only)
        for d in reversed(self.drawing_lineage[:-1]):
            p = d.get_positions(displayed_only) * p
        return p
        
    def get_positions(self, displayed_only=False):
        if displayed_only:
            return self._positions.masked(self.display_positions)
        return self._positions

    def set_positions(self, positions):
        from ..geometry import Places
        if positions and not isinstance(positions, Places):
            raise ValueError('Got %s instead of Places' % str(type(positions)))
        self._positions = positions
        np = len(positions)
        if self._displayed_positions is not None and len(self._displayed_positions) != np:
            self._displayed_positions = None
        if self._selected_positions is not None and len(self._selected_positions) != np:
            self._selected_positions = None
        if len(self._colors) != np:
            from numpy import empty, uint8
            c = empty((np, 4), uint8)
            c[:,:] = self._colors[0,:] if len(self._colors) > 0 else 255
            self._colors = c
        self.redraw_needed(shape_changed=True)

    positions = property(get_positions, set_positions)
    '''
    Copies of the surface piece are placed using a 3 by 4 matrix with
    the first 3 columns giving a linear transformation, and the last
    column specifying a shift.
    '''

    def number_of_positions(self, displayed_only=False):
        '''Number of positions the Drawing is placed at.'''
        if displayed_only and not self.display:
            return 0
        dp = self.display_positions
        np = len(self.positions) if dp is None else dp.sum()
        return np

    def get_color(self):
        return self._colors[0]

    def set_color(self, rgba):
        from numpy import empty, uint8
        np = len(self._positions)
        c = empty((np, 4), uint8)
        c[:, :] = rgba
        self._colors = c
        self._opaque_color_count = (np if rgba[3] == 255 else 0)
        self.redraw_needed()

    color = property(get_color, set_color)
    '''Single color of drawing used when per-vertex coloring is not
    specified, 0-255 red, green, blue, alpha values.'''

    def get_colors(self, displayed_only=False):
        if displayed_only:
            dp = self.display_positions
            return self._colors if dp is None else self._colors[dp]
        return self._colors

    def set_colors(self, rgba):
        from numpy import ndarray, array, uint8
        c = rgba if isinstance(rgba, ndarray) else array(rgba, uint8)
        self._colors = c
        opc = opaque_count(c)
        tchange = (opc != self._opaque_color_count)
        self._opaque_color_count = opc
        self.redraw_needed(transparency_changed = tchange)

    colors = property(get_colors, set_colors)
    '''Color for each position used when per-vertex coloring is not
    specified.'''

    def get_vertex_colors(self, create = False, copy = False):
        vc = self._vertex_colors
        if vc is None:
            if create:
                nv = len(self.vertices)
                from numpy import empty, uint8
                vc = empty((nv,4), uint8)
                vc[:] = self.color
        elif copy:
            vc = vc.copy()
        return vc
    def set_vertex_colors(self, vcolors):
        self._vertex_colors = vcolors
        opvc = opaque_count(vcolors)
        tchange = (opvc != self._opaque_vertex_color_count)
        self._opaque_vertex_color_count = opvc
        self.auto_recolor_vertices = None
        self.redraw_needed(transparency_changed = tchange)
    vertex_colors = property(get_vertex_colors, set_vertex_colors)
    '''
    R, G, B, A color and transparency for each vertex, a numpy N by
    4 array of uint8 values, can be None in which case a single color
    (attribute color) is used for the object.
    '''

    def set_transparency(self, alpha):
        '''
        Set transparency to alpha (0-255). Applies to per-vertex colors if
        currently showing per-vertex colors otherwise single color.
        Does not effect child drawings.
        '''
        vcolors = self.vertex_colors
        if vcolors is None:
            c = self.colors
            c[:, 3] = alpha
            self.colors = c
        else:
            vcolors[:, 3] = alpha
            self.vertex_colors = vcolors
        self.redraw_needed(transparency_changed = True)

    def _transparency(self):
        if self.texture is not None or self.multitexture:
            any_opaque = self.opaque_texture
            any_transparent = not self.opaque_texture
        else:
            vc = self.vertex_colors
            if vc is None:
                oc = self._opaque_color_count
                any_opaque = (oc > 0)
                any_transparent = (oc < len(self._colors))
            else:
                oc = self._opaque_vertex_color_count
                any_opaque = (oc > 0)
                any_transparent = (oc < len(vc))
        return any_opaque, any_transparent

    def showing_transparent(self):
        '''Are any transparent objects being displayed. Includes all
        children.'''
        if self.display:
            if not self.empty_drawing():
                any_opaque, any_transp = self._transparency()
                if any_transp:
                    return True
            for d in self.child_drawings():
                if d.showing_transparent():
                    return True
        return False

    def set_geometry(self, vertices, normals, triangles):
        '''Set vertices, normals and triangles defining the shape to be drawn.'''
        self._vertices = vertices
        self._normals = normals
        self._triangles = triangles
        self._vertex_colors = None
        self._edge_mask = None
        self._triangle_mask = None
        self.redraw_needed(shape_changed=True)

        arv = self.auto_recolor_vertices
        if arv:
            arv()

        art = self.auto_remask_triangles
        if art:
            art()

    def empty_drawing(self):
        '''Does this drawing have no geometry? Does not consider child
        drawings.'''
        v,t = self.vertices, self.triangles
        return v is None or t is None or len(t) == 0

    def number_of_triangles(self, displayed_only=False):
        '''Return the number of triangles including all child drawings
        and all positions.'''
        np = self.number_of_positions(displayed_only)
        if np == 0:
            return 0
        t = self.triangles
        tc = 0 if t is None else np * len(t)
        for d in self.child_drawings():
            tc += np * d.number_of_triangles(displayed_only)
        return tc

    OPAQUE_DRAW_PASS = 'opaque'
    "Draw pass to render only opaque drawings."
    TRANSPARENT_DRAW_PASS = 'transparent'
    "Draw pass to render only transparent drawings."
    TRANSPARENT_DEPTH_DRAW_PASS = 'transparent depth'
    "Draw pass to render only the depth of transparent drawings."
    SELECTION_DRAW_PASS = 'selection'
    "Draw pass to render only the selected parts of drawings."

    def draw(self, renderer, place, draw_pass, selected_only=False):
        '''Draw this drawing and children using the given draw pass.'''

        if not self.display:
            return

        if not self.empty_drawing():
            self.draw_self(renderer, place, draw_pass, selected_only)

        if self.child_drawings():
            for p in self.positions:
                pp = place if p.is_identity() else place * p
                self._draw_children(renderer, pp, draw_pass, selected_only)

    def draw_self(self, renderer, place, draw_pass, selected_only=False):
        '''Draw this drawing without children using the given draw pass.'''

        if selected_only and not self.selected:
            return

        if (len(self.positions) == 1 and
                self.positions.shift_and_scale_array() is None):
            p = self.position
            pp = place if p.is_identity() else place * p
        else:
            pp = place
        renderer.set_model_matrix(pp)

        if draw_pass == self.OPAQUE_DRAW_PASS:
            any_opaque, any_transp = self._transparency()
            if any_opaque:
                self._draw_geometry(renderer, selected_only, opaque_only = any_transp)
        elif draw_pass in (self.TRANSPARENT_DRAW_PASS,
                           self.TRANSPARENT_DEPTH_DRAW_PASS):
            any_opaque, any_transp = self._transparency()
            if any_transp:
                self._draw_geometry(renderer, selected_only, transparent_only = any_opaque)
        elif draw_pass == self.SELECTION_DRAW_PASS:
            self._draw_geometry(renderer, selected_only)

    def _draw_children(self, renderer, place, draw_pass, selected_only=False):
        dlist = self.child_drawings()
        for d in dlist:
            d.draw(renderer, place, draw_pass, selected_only)

    def _draw_geometry(self, renderer, selected_only=False,
                       transparent_only=False, opaque_only=False):
        ''' Draw the geometry.'''

        if self.vertices is None:
            return

        self._opengl_context = renderer.opengl_context
        
        if len(self._vertex_buffers) == 0:
            self._create_vertex_buffers()

        # Update opengl buffers to reflect drawing changes
        self._update_buffers()

        ds = self._draw_selection if selected_only else self._draw_shape
        ds.activate_bindings()

        sopt = self._shader_options(transparent_only, opaque_only)
        r = renderer
        shader = r.shader(sopt)
        r.use_shader(shader)

        # Set color
        if self.vertex_colors is None and len(self._colors) == 1:
            r.set_single_color([c / 255.0 for c in self._colors[0]])

        t = self.texture
        if t is not None:
            t.bind_texture()

        at = self.ambient_texture
        if at is not None:
            at.bind_texture()
            r.set_ambient_texture_transform(self.ambient_texture_transform)

        # Draw triangles
        mtex = self.multitexture
        if mtex:
            ds.draw_multitexture(self.display_style, mtex, self.multitexture_reverse_order)
        else:
            ds.draw(self.display_style)

        if t is not None:
            t.unbind_texture()

    def _shader_options(self, transparent_only = False, opaque_only = False):
        sopt = self._shader_opt
        if sopt is None:
            sopt = 0
            from .opengl import Render
            if self.use_lighting:
                sopt |= Render.SHADER_LIGHTING
            if (self.vertex_colors is not None) or len(self._colors) > 1:
                sopt |= Render.SHADER_VERTEX_COLORS
            t = self.texture
            if t or self.multitexture:
                if t and t.is_cubemap:
                    sopt |= Render.SHADER_TEXTURE_CUBEMAP
                else:
                    sopt |= Render.SHADER_TEXTURE_2D
            if self.ambient_texture is not None:
                sopt |= Render.SHADER_TEXTURE_3D_AMBIENT
            if self.positions.shift_and_scale_array() is not None:
                sopt |= Render.SHADER_SHIFT_AND_SCALE
            elif len(self.positions) > 1:
                sopt |= Render.SHADER_INSTANCING
            self._shader_opt = sopt
        if transparent_only:
            from .opengl import Render
            sopt |= Render.SHADER_TRANSPARENT_ONLY
        if opaque_only:
            from .opengl import Render
            sopt |= Render.SHADER_OPAQUE_ONLY
        return sopt

    _effects_shader = set(
        ('use_lighting', '_vertex_colors', '_colors', 'texture',
         'ambient_texture', '_positions'))

    # Update the contents of vertex, element and instance buffers if associated
    #  arrays have changed.
    def _update_buffers(self):

        changes = self._attribute_changes
        if len(changes) == 0:
            return

        ds, dss = self._draw_shape, self._draw_selection

        # Update drawing and selection triangle buffers
        if ('_triangles' in changes or
            'display_style' in changes or
            '_triangle_mask' in changes or
            (self.display_style == self.Mesh and '_edge_mask' in changes) or
            '_selected_triangles_mask' in changes):
            ta = self.triangles
            style = self.display_style
            em = self._edge_mask
            tm = self._triangle_mask
            tmsel = self.selected_displayed_triangles_mask
            ds.update_element_buffer(ta, style, tm, em)
            dss.update_element_buffer(ta, style, tmsel, em)

        # Update instancing buffers
        p = self.positions
        if ('_colors' in changes or
            '_vertex_colors' in changes or
            '_positions' in changes or
            '_displayed_positions' in changes or
            '_selected_positions' in changes):
            c = self.colors if self._vertex_colors is None else None
            pm = self._position_mask()
            pmsel = self._position_mask(True)
            ds.update_instance_buffers(p, c, pm)
            dss.update_instance_buffers(p, c, pmsel)

        # Update buffers shared by drawing and selection
        for b in self._vertex_buffers:
            aname = b.buffer_attribute_name
            if aname in changes:
                data = getattr(self, aname)
                ds.update_vertex_buffer(b, data)
                dss.update_vertex_buffer(b, data)

        changes.clear()

    def _position_mask(self, selected_only=False):
        dp = self._displayed_positions        # bool array
        if selected_only:
            sp = self._selected_positions
            if sp is not None:
                import numpy
                dp = sp if dp is None else numpy.logical_and(dp, sp)
        return dp

    def bounds(self, positions=True):
        '''
        The bounds of drawing and displayed children and displayed positions
        in the parent model coordinate system.
        '''
        from ..geometry import bounds
        dbounds = [d.bounds() for d in self.child_drawings() if d.display and not hasattr(d, 'skip_bounds')]
        if dbounds:
            cb = bounds.union_bounds(dbounds)
            if positions:
                cb = bounds.copies_bounding_box(cb, self.get_positions(displayed_only=True))
        else:
            cb = None
        if self.empty_drawing():
            return cb
        if not positions:
            b = self._geometry_bounds()
            if cb:
                b = bounds.union_bounds((b, cb))
            return b
        if self._cached_position_bounds is not None:
            pb = self._cached_position_bounds
        else:
            b = self._geometry_bounds()
            pb = bounds.copies_bounding_box(
                    b, self.get_positions(displayed_only=True))
            self._cached_position_bounds = pb
        if cb is None:
            return pb
        b = bounds.union_bounds((pb, cb))
        return b

    def _geometry_bounds(self):
        '''
        Return the bounds of the drawing not including positions nor children
        in this drawing's coordinate system.
        '''
        cb = self._cached_geometry_bounds
        if cb is not None:
            return cb

        va = self.vertices
        if va is None:
            return None
        tmask = self._triangle_mask
        if not tmask is None:
            import numpy
            vshown = numpy.unique(self.triangles[tmask,:])
            va = va[vshown,:]
        if len(va) == 0:
            return None
        xyz_min = va.min(axis=0)
        xyz_max = va.max(axis=0)
        from ..geometry.bounds import Bounds
        b = Bounds(xyz_min, xyz_max)
        self._cached_geometry_bounds = b
        return b

    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        '''
        Find the first intercept of a line segment with the displayed part
        of this drawing and its children.  The end points are in the parent
        drawing coordinates and do not take account of this Drawings positions.
        If the exclude option is given it is a function that takes a drawing
        and returns true if this drawing and its children should be excluded.
        Return a Pick object for the intercepted item.
        The Pick object has a distance attribute giving the fraction (0-1)
        along the segment where the intersection occurs.
        For no intersection None is returned.  This routine is used for
        selecting objects, for identifying objects during mouse-over, and
        to determine the front-most point in the center of view to be used
        as the interactive center of rotation.
        '''
        if not self.display:
            return None
        if exclude is not None and exclude(self):
            return None

        pclosest = None
        if not self.empty_drawing():
            p = self._first_intercept_excluding_children(mxyz1, mxyz2)
            if p and (pclosest is None or p.distance < pclosest.distance):
                pclosest = p
        p = self.first_intercept_children(self.child_drawings(), mxyz1, mxyz2, exclude)
        if p and (pclosest is None or p.distance < pclosest.distance):
            pclosest = p
        return pclosest

    def first_intercept_children(self, child_drawings, mxyz1, mxyz2, exclude=None):
        '''
        Like first_intercept() but check for intercepts with just the specified children.
        '''
        if len(child_drawings) == 0:
            return None
        pclosest = None
        pos = [p.inverse() * (mxyz1, mxyz2) for p in self.positions]
        for d in child_drawings:
            if d.display and (exclude is None or not exclude(d)):
                for cxyz1, cxyz2 in pos:
                    p = d.first_intercept(cxyz1, cxyz2, exclude)
                    if p and (pclosest is None or p.distance < pclosest.distance):
                        pclosest = p
        return pclosest

    def _first_intercept_excluding_children(self, mxyz1, mxyz2):
        if self.empty_drawing():
            return None
        va = self.vertices
        ta = self.masked_triangles
        if ta.shape[1] != 3:
            # TODO: Intercept only for triangles, not lines or points.
            return None
        p = None
        from ..geometry import closest_triangle_intercept
        if self.positions.is_identity():
            fmin, tmin = closest_triangle_intercept(va, ta, mxyz1, mxyz2)
            if fmin is not None:
                p = PickedTriangle(fmin, tmin, 0, self)
        else:
            pos_nums = self.bounds_intercept_copies(self._geometry_bounds(), mxyz1, mxyz2)
            for i in pos_nums:
                cxyz1, cxyz2 = self.positions[i].inverse() * (mxyz1, mxyz2)
                fmin, tmin = closest_triangle_intercept(va, ta, cxyz1, cxyz2)
                if fmin is not None and (p is None or fmin < p.distance):
                    p = PickedTriangle(fmin, tmin, i, self)
        return p

    def bounds_intercept_copies(self, bounds, mxyz1, mxyz2):
        '''
        Return indices of positions where line segment intercepts displayed bounds.
        This is to optimize picking so that positions where no intercept occurs do not
        need to be checked to see what is picked.
        '''
        # Only check objects with bounding box close to line. 
        b = bounds
        if b is None:
            return []
        c, r = b.center(), b.radius()
        pc = self.positions * c
        from ..geometry import segment_intercepts_spheres
        bi = segment_intercepts_spheres(pc, r, mxyz1, mxyz2)
        dp = self._displayed_positions
        if dp is not None:
            from numpy import logical_and
            logical_and(bi, dp, bi)
        pos_nums = bi.nonzero()[0]
        return pos_nums

    def planes_pick(self, planes, exclude=None):
        '''
        Find the displayed drawing instances bounded by the specified planes
        for this drawing and its children.  Each plane is a 4-vector v with
        points in the pick region v0*x + v1*y + v2*z + v3 >= 0 using coordinate
        system of the parent drawing.  If a drawing has instances then only
        the center of each instance is considered and the whole instance is
        picked if the center is within the planes.  If a drawing has only one
        instance (self.positions has length 1) then the pick lists the individual
        triangles which have at least one vertex within all of the planes.
        If exclude is not None then it is a function called with a Drawing argument
        that returns true if this drawing and its children should be excluded
        from the pick.
        Return a list of Pick objects for the contained items.
        This routine is used for selecting objects in a frustum.
        '''
        if not self.display:
            return []
        if exclude is not None and exclude(self):
            return []

        picks = []
        if not self.empty_drawing():
            from ..geometry import points_within_planes
            if len(self.positions) > 1:
                # Use center of instances.
                b = self._geometry_bounds()
                if b:
                    pc = self.positions * b.center()
                    pmask = points_within_planes(pc, planes)
                    if pmask.sum() > 0:
                        dp = self.display_positions
                        if dp is not None:
                            # Pick displayed positions only
                            from numpy import logical_and
                            logical_and(pmask, dp, pmask)
                        picks.append(PickedInstance(pmask, self))
            else:
                # For non-instances pick using all vertices.
                from ..geometry import transform_planes
                pplanes = transform_planes(self.position, planes)
                vmask = points_within_planes(self.vertices, pplanes)
                if vmask.sum() > 0:
                    t = self.triangles
                    from numpy import logical_or, logical_and
                    tmask = logical_or(vmask[t[:,0]], vmask[t[:,1]])
                    logical_or(tmask, vmask[t[:,2]], tmask)
                    tm = self._triangle_mask
                    if tm is not None:
                        logical_and(tmask, tm, tmask)
                    if tmask.sum() > 0:
                        picks.append(PickedTriangles(tmask, self))
        from ..geometry import transform_planes
        for d in self.child_drawings():
            for p in self.positions:
                pplanes = transform_planes(p, planes)
                picks.extend(d.planes_pick(pplanes, exclude))
        return picks

    def __del__(self):
        if not self.was_deleted:
            # Release opengl resources.
            self.delete()
            
    def delete(self):
        '''
        Delete drawing and all child drawings.
        '''
        c = self._opengl_context
        if c:
            c.make_current()	# Make OpenGL context current for deleting OpenGL resources.
        self._delete_geometry()
        self.remove_all_drawings()

    def _delete_geometry(self):
        '''Release all the arrays and graphics memory associated with
        the surface piece.'''
        self._positions = None
        self._colors = None
        self._displayed_positions = None
        self.auto_recolor_vertices = None
        self.auto_remask_triangles = None
        self._vertices = None
        self._triangles = None
        self._normals = None
        self._edge_mask = None
        self._triangle_mask = None
        if self.texture:
            self.texture.delete_texture()
            self.texture = None
        if self.multitexture:
            for t in self.multitexture:
                t.delete_texture()
            self.multitexture = None
        self.texture_coordinates = None

        for b in self._vertex_buffers:
            b.delete_buffer()
        self._vertex_buffers = []

        for ds in (self._draw_shape, self._draw_selection):
            if ds:
                ds.delete()
        self._draw_shape = None
        self._draw_selection = None

        self._opengl_context = None
        
        self.was_deleted = True

    def _create_vertex_buffers(self):
        from . import opengl
        vbufs = (
            ('_vertices', opengl.VERTEX_BUFFER),
            ('_normals', opengl.NORMAL_BUFFER),
            ('_vertex_colors', opengl.VERTEX_COLOR_BUFFER),
            ('texture_coordinates', opengl.TEXTURE_COORDS_BUFFER),
        )

        self._vertex_buffers = vb = []
        for a, v in vbufs:
            b = opengl.Buffer(v)
            b.buffer_attribute_name = a
            vb.append(b)

        self._draw_shape = _DrawShape(vb)
        self._draw_selection = _DrawShape(vb)

    _effects_buffers = set(
        ('_vertices', '_normals', '_vertex_colors', 'texture_coordinates',
         '_triangles', 'display_style', '_displayed_positions', '_colors', '_positions',
         '_edge_mask', '_triangle_mask', '_selected_triangles_mask', '_selected_positions'))

    EDGE0_DISPLAY_MASK = 1
    ALL_EDGES_DISPLAY_MASK = 7
    '''Edge mask for displaying all three triangle edges (bits 0, 1, 2).'''

    def get_triangle_mask(self):
        return self._triangle_mask

    def set_triangle_mask(self, tmask):
        self._triangle_mask = tmask
        self.redraw_needed(shape_changed=True)
        self.auto_remask_triangles = None
        
    triangle_mask = property(get_triangle_mask, set_triangle_mask)
    '''
    The triangle mask is a 1-dimensional bool numpy array of
    length equal to the number of triangles used to control display
    of individual triangles.
    '''

    @property
    def selected_displayed_triangles_mask(self):
        '''Mask of selected and displayed triangles.'''
        tm = self._triangle_mask
        tmsel = self._selected_triangles_mask
        if tm is not None:
            # Combine selected and displayed triangle masks
            if tmsel is not None:
                from numpy import logical_and
                tmsel = logical_and(tmsel, tm)
            else:
                tmsel = tm
        return tmsel

    def get_edge_mask(self):
        return self._edge_mask

    def set_edge_mask(self, emask):
        self._edge_mask = emask
        self.redraw_needed(shape_changed=True)

    edge_mask = property(get_edge_mask, set_edge_mask)
    '''
    The edge mask is a 1-dimensional uint8 numpy array of
    length equal to the number of triangles.  The lowest 3 bits are used
    to control display of the 3 triangle edges in mesh mode.
    '''

    @property
    def masked_triangles(self):
        ta = self.triangles
        if ta is None:
            return None
        tm = self._triangle_mask
        return ta if tm is None else ta[tm,:]

    @property
    def num_masked_triangles(self):
        ta = self.triangles
        if ta is None:
            return 0
        tm = self._triangle_mask
        return len(ta) if tm is None else tm.sum()

    def x3d_needs(self, x3d_scene):
        if not self.display:
            return
        dlist = self.child_drawings()
        for d in dlist:
            d.x3d_needs(x3d_scene)
        if self.empty_drawing():
            return
        any_opaque, any_transp = self._transparency()
        from .. import x3d
        # x3d_scene.need(x3d.Components.Core, 2)  # Prototyping
        x3d_scene.need(x3d.Components.Grouping, 1)  # Group, Transform
        if any_transp and self.vertex_colors:
            x3d_scene.need(x3d.Components.Rendering, 4)  # ColorRGBA
        else:
            x3d_scene.need(x3d.Components.Rendering, 3)  # IndexedTriangleSet
        # x3d_scene.need(x3d.Components.Rendering, 5)  # ClipPlane
        x3d_scene.need(x3d.Components.Shape, 1)  # Appearance, Material, Shape
        # x3d_scene.need(x3d.Components.Shape, 2)  # LineProperties
        # x3d_scene.need(x3d.Components.Shape, 3)  # FillProperties
        # x3d_scene.need(x3d.Components.Geometry3D, 1)  # Cylinder, Sphere
        # x3d_scene.need(x3d.Components.Geometry3D, 4)  # Extrusion
        if self.texture or self.multitexture:
            x3d_scene.need(x3d.Components.Texturing, 1)  # PixelTexture

    def write_x3d(self, stream, x3d_scene, indent, place):
        if not self.display:
            return
        dlist = self.child_drawings()
        if dlist:
            for p in self.positions:
                pp = place if p.is_identity() else place * p
                for d in dlist:
                    d.write_x3d(stream, x3d_scene, indent, pp)
        if not self.empty_drawing():
            self.custom_x3d(stream, x3d_scene, indent, place)

    def reuse_unlit_appearance(self, stream, x3d_scene, indent, color, line_width, line_type):
        tab = ' ' * indent
        use, name = self.def_or_use((tuple(color), line_width, line_type), 'aup')
        if use == 'USE':
            print("%s<Appearance USE='%s'/>" % (tab, name), file=stream)
            return

        from graphics.linetype import LineType
        print("%s<Appearance DEF='%s'>" % (tab, name), file=stream)
        if line_width != 1 or line_type != LineType.Solid:
            print("%s <LineProperties" % tab, end='', file=stream)
            if line_width != 1:
                    print(" linewidthScaleFactor='%g'" % line_width, end='', file=stream)
            if line_type != LineType.Solid:
                    print(" linetype='%d'" % line_type.value, end='', file=stream)
            print("/>", file=stream)
        if color is not None:
            color.write_x3d(stream, indent + 1, False)
        print("%s</Appearance>" % tab, file=stream)

    def reuse_appearance(self, stream, x3d_scene, indent, color):
        if color is None:
            return
        tab = ' ' * indent
        use, name = x3d_scene.def_or_use(tuple(color), 'ap')
        if use == 'USE':
            print("%s<Appearance USE='%s'/>" % (tab, name), file=stream)
            return

        print("%s<Appearance DEF='%s'>" % (tab, name), file=stream)
        print("%s <Material ambientIntensity='1' diffuseColor='%g %g %g' specularColor='0.85 0.85 0.85' shininess='0.234375' transparency='%g'/>" % (tab, color[0] / 255, color[1] / 255, color[2] / 255, 1 - color[3] / 255), file=stream)
        print('%s</Appearance>' % tab, file=stream)

    def reuse_its(self, stream, x3d_scene, indent, def_use_tag, indices, colors, normals, any_transp):
        tab = ' ' * indent
        if def_use_tag is None:
            def_use = ''
        else:
            use, name = x3d_scene.def_or_use(def_use_tag, 'its')
            if use == 'USE':
                print("%s<IndexedTriangleSet USE='%s'/>" % (tab, name), file=stream)
                return
            def_use = "%s='%s' " % (use, name)

        def bulk_write(values, chunk_size, stream):
            # break up list of values into chucks to try to avoid bug on Mac OS X writing huge strings

            sep = ''
            for i in range(0, len(values), chunk_size):
                print('%s%s' % (sep, ' '.join(values[i:i + chunk_size])), end='', file=stream)
                sep = ' '

        indices = ['%g' % i for i in indices]
        print('%s<IndexedTriangleSet %sindex="%s"' % (tab, def_use, ' '.join(indices)), end='', file=stream)

        print(' solid="false"', end='', file=stream)
        if colors is None:
            print(' colorPerVertex="false"', end='', file=stream)
        if normals is None:
            print(' normalPerVertex="false"', end='', file=stream)
        print('>', file=stream)
        vertices = ['%g' % x for x in self.vertices.flatten()]
        print('%s <Coordinate point="' % tab, end='', file=stream)
        bulk_write(vertices, 3 * 1024, stream)
        print('"/>', file=stream)
        if normals is not None:
            normals = ['%g' % x for x in normals.flatten()]
            print('%s <Normal vector="' % tab, end='', file=stream)
            bulk_write(normals, 3 * 1024, stream)
            print('"/>', file=stream)
        if colors is not None:
            if any_transp:
                colors = ['%g' % x for x in (colors / 255).flatten()]
                print('%s <ColorRGBA color="' % tab, end='', file=stream)
                bulk_write(colors, 4 * 1024, stream)
                print('"/>', file=stream)
            else:
                colors = ['%g' % x for x in (colors[:, 0:3] / 255).flatten()]
                print('%s <Color color="' % tab, end='', file=stream)
                bulk_write(colors, 3 * 1024, stream)
                print('"/>', file=stream)
        print('%s</IndexedTriangleSet>' % tab, file=stream)

    def custom_x3d(self, stream, x3d_scene, indent, place):
        """Override this function for custom X3D

        This is a generic version and assumes that positions are orthogonal.
        """
        any_opaque, any_transp = self._transparency()
        # cases:
        #  1 position, 1 color
        #  multpiple positions, multiple colors
        #  (multpiple positions, 1 color)
        #  1 postion, per-vertex coloring
        #  multiple postions, per-vertex coloring
        has_ssa = self.positions.shift_and_scale_array() is not None
        tab = ' ' * indent
        def_use_tag = self  # always allow reuse
        print('%s<Group>' % tab, file=stream)
        colors = self.vertex_colors
        normals = self.normals
        indices = self.masked_triangles.flatten()
        from math import radians
        for p, c in zip(self.positions, self.colors):
            p = place if p.is_identity() else place * p
            if has_ssa:
                s = (p.matrix[0][0], p.matrix[1][1], p.matrix[2][2])
                t = p.translation()
                print('%s<Transform scale="%g %g %g" translation="%g %g %g">' % (tab, s[0], s[1], s[2], t[0], t[1], t[2]), file=stream)
            else:
                r = p.rotation_axis_and_angle()
                t = p.translation()
                print('%s<Transform rotation="%g %g %g %g" translation="%g %g %g">' % (tab, r[0][0], r[0][1], r[0][2], radians(r[1]), t[0], t[1], t[2]), file=stream)
            print('%s <Shape>' % tab, file=stream)
            self.reuse_appearance(stream, x3d_scene, indent + 2, c)
            self.reuse_its(stream, x3d_scene, indent + 2, def_use_tag, indices,
                           colors, normals, any_transp)
            print('%s </Shape>' % tab, file=stream)
            print('%s</Transform>' % tab, file=stream)
        print('%s</Group>' % tab, file=stream)

def opaque_count(rgba):
    if rgba is None:
        return 0
    from . import _graphics
    return _graphics.count_value(rgba[:,3], 255)

def draw_drawings(renderer, cvinv, drawings, opaque_only = False):
    '''
    Render opaque and transparent draw passes for a given set of drawings,
    and given camera view (inverse camera transform).
    '''
    r = renderer
    r.set_view_matrix(cvinv)
    from ..geometry import Place
    p = Place()
    _draw_multiple(drawings, r, p, Drawing.OPAQUE_DRAW_PASS)
    if not opaque_only and _any_transparent_drawings(drawings):
        r.draw_transparent(
            lambda: _draw_multiple(drawings, r, p, Drawing.TRANSPARENT_DEPTH_DRAW_PASS),
            lambda: _draw_multiple(drawings, r, p, Drawing.TRANSPARENT_DRAW_PASS))

def draw_opaque(renderer, drawings):
    from ..geometry import identity
    _draw_multiple(drawings, renderer, identity(), Drawing.OPAQUE_DRAW_PASS)

def draw_transparent(renderer, drawings):
    if _any_transparent_drawings(drawings):
        from ..geometry import identity
        p = identity()
        r = renderer
        r.draw_transparent(
            lambda: _draw_multiple(drawings, r, p, Drawing.TRANSPARENT_DEPTH_DRAW_PASS),
            lambda: _draw_multiple(drawings, r, p, Drawing.TRANSPARENT_DRAW_PASS))


def _draw_multiple(drawings, renderer, place, draw_pass):
    selected_only = (draw_pass == Drawing.SELECTION_DRAW_PASS)
    for d in drawings:
        d.draw(renderer, place, draw_pass, selected_only)


def _any_transparent_drawings(drawings):
    for d in drawings:
        if d.showing_transparent():
            return True
    return False


def draw_depth(renderer, cvinv, drawings, opaque_only = True):
    '''Render only the depth buffer (not colors).'''
    r = renderer
    r.disable_shader_capabilities(r.SHADER_LIGHTING |
                                  r.SHADER_TEXTURE_2D)
    draw_drawings(r, cvinv, drawings, opaque_only)
    r.disable_shader_capabilities(0)


def draw_overlays(drawings, renderer):
    '''Render drawings using an identity projection matrix with no
    depth test.'''
    r = renderer
    r.disable_shader_capabilities(r.SHADER_STEREO_360 |	# Avoid geometry shift
                                  r.SHADER_DEPTH_CUE |
                                  r.SHADER_SHADOWS |
                                  r.SHADER_MULTISHADOW |
                                  r.SHADER_CLIP_PLANES)
    r.set_projection_matrix(((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0),
                             (0, 0, 0, 1)))

    from ..geometry import place
    p0 = place.identity()
    r.set_view_matrix(p0)
    r.set_model_matrix(p0)
    r.enable_depth_test(False)
    r.enable_blending(False)
    _draw_multiple(drawings, r, p0, Drawing.OPAQUE_DRAW_PASS)
    r.enable_blending(True)
    _draw_multiple(drawings, r, p0, Drawing.TRANSPARENT_DRAW_PASS)
    r.enable_blending(False)
    r.enable_depth_test(True)
    r.disable_shader_capabilities(0)


def draw_selection_outline(renderer, cvinv, drawings):
    '''Draw the outlines of selected parts of the specified drawings.'''
    r = renderer
    r.set_view_matrix(cvinv)
    r.start_rendering_outline()
    from ..geometry import Place
    p = Place()
    _draw_multiple(drawings, r, p, Drawing.SELECTION_DRAW_PASS)
    r.finish_rendering_outline()


def draw_xor_rectangle(renderer, x1, y1, x2, y2, color, drawing = None):
    '''Draw rectangle outline on front buffer using xor mode.'''

    if drawing is None:
        d = Drawing('selection drag box')
        from numpy import array, int32, uint8
        t = array(((0,1,2), (0,2,3)), int32)
        d.set_geometry(None, None, t)
        d.edge_mask = array((3, 6), uint8)
        d.display_style = d.Mesh
        d.use_lighting = False
    else:
        d = drawing

    r = renderer
    s = r.pixel_scale()
    from numpy import array, float32
    v = array(((s*x1, s*y1, 0), (s*x2, s*y1, 0),
               (s*x2, s*y2, 0), (s*x1, s*y2, 0)),
              float32)
    d.set_geometry(v, None, d.triangles)
    d.color = color

    from ..geometry import identity
    p0 = identity()
    r.set_view_matrix(p0)

    r.draw_front_buffer(True)
    r.enable_depth_test(False)
    r.enable_xor(True)
    rdc = r.disable_capabilities
    r.disable_capabilities = rdc | r.SHADER_CLIP_PLANES

    w, h = r.render_size()
    from .camera import ortho
    r.set_projection_matrix(ortho(0, w, 0, h, -1, 1))

    d.draw(r, p0, d.OPAQUE_DRAW_PASS)

    r.disable_capabilities = rdc
    r.enable_xor(False)
    r.enable_depth_test(True)
    r.draw_front_buffer(False)
    r.flush()

    return d

def _element_type(display_style):
    from .opengl import Buffer
    if display_style == Drawing.Solid:
        t = Buffer.triangles
    elif display_style == Drawing.Mesh:
        t = Buffer.lines
    elif display_style == Drawing.Dot:
        t = Buffer.points
    return t


class _DrawShape:

    def __init__(self, vertex_buffers):

        # Arrays derived from positions, colors and geometry
        self.instance_shift_and_scale = None   # N by 4 array, (x, y, z, scale)
        self.instance_matrices = None	    # matrices for displayed instances
        self.instance_colors = None
        self.elements = None                # Triangles after mask applied
        self.masked_edges = None
        self._edge_mask = None
        self._tri_mask = None

        # Vertex buffer data
        self.vertices = None
        self.normals = None
        self.vertex_colors = None
        self.texture_coordinates = None

        # OpenGL rendering
        self.bindings = None    	      # Shader variable bindings in an opengl vertex array object
        self._buffers_need_update = set()     # Buffers that need data copied to opengl buffer object
        self.vertex_buffers = vertex_buffers
        self.element_buffer = None
        self.instance_buffers = []

    def __del__(self):
        self.delete()
            
    def delete(self):

        self.masked_edges = None
        self.instance_shift_and_scale = None
        self.instance_matrices = None
        self.instance_colors = None
        if self.element_buffer:
            self.element_buffer.delete_buffer()
            self.element_buffer = None
            for b in self.instance_buffers:
                b.delete_buffer()
            self.instance_buffers = []

        if self.bindings:
            self.bindings.delete_bindings()
            self.bindings = None

    def draw(self, display_style):

        eb = self.element_buffer
        if eb is None:
            return
        etype = _element_type(display_style)
        ni = self.instance_count()
        if ni > 0:
            eb.draw_elements(etype, ni)

    def draw_multitexture(self, display_style, textures, reverse_order):

        eb = self.element_buffer
        if eb is None:
            return
        etype = _element_type(display_style)
        ni = self.instance_count()
        if ni > 0:
            nt = len(textures)
            ne = eb.size() // nt
            torder = range(nt-1,-1,-1) if reverse_order else range(nt)
            for ti in torder:
                textures[ti].bind_texture()
                eb.draw_elements(etype, ni, count=ne, offset=ti*ne)
            # TODO: I put unbind outside loop for better efficiency.  But if textures use
            # multiple texture units, this will only unbind the last one.
            textures[nt-1].unbind_texture()

    def update_vertex_buffer(self, b, data):

        setattr(self, b.buffer_attribute_name, data)
        self.buffer_needs_update(b)

    def create_element_buffer(self):

        from . import opengl
        eb = opengl.Buffer(opengl.ELEMENT_BUFFER)
        eb.buffer_attribute_name = 'elements'
        return eb

    def update_element_buffer(self, triangles, style, triangle_mask, edge_mask):

        self.elements = e = self.masked_elements(triangles, style, triangle_mask, edge_mask)

        eb = self.element_buffer
        if len(e) > 0 and eb is None:
            self.element_buffer = eb = self.create_element_buffer()
        elif len(e) == 0 and eb:
            eb.delete_buffer()
            self.element_buffer = eb = None

        if eb:
            self.buffer_needs_update(eb)

    def masked_elements(self, triangles, style, tmask, edge_mask):

        ta = triangles
        if ta is None:
            return None
        if style == Drawing.Solid:
            if tmask is not None:
                ta = ta[tmask, :]
        elif style == Drawing.Mesh:
            from ._graphics import masked_edges
            if ta.shape[1] == 2:
                pass    # Triangles array already contains edges.
            elif edge_mask is None:
                kw = {} if tmask is None else {'triangle_mask': tmask}
                ta = masked_edges(ta, **kw)
            else:
                # TODO: Need to reset masked_edges if edge_mask changed.
                me = self.masked_edges
                if (me is None or edge_mask is not self._edge_mask or
                    tmask is not self._tri_mask):
                    kw = {}
                    if edge_mask is not None:
                        kw['edge_mask'] = edge_mask
                    if tmask is not None:
                        kw['triangle_mask'] = tmask
                    self.masked_edges = me = masked_edges(ta, **kw)
                    self._edge_mask, self._tri_mask = edge_mask, tmask
                ta = me
        elif style == Drawing.Dot:
            from numpy import unique
            ta = unique(ta)
        return ta

    def create_instance_buffers(self):

        from . import opengl
        ibufs = (
            ('instance_shift_and_scale', opengl.INSTANCE_SHIFT_AND_SCALE_BUFFER),
            ('instance_matrices', opengl.INSTANCE_MATRIX_BUFFER),
            ('instance_colors', opengl.INSTANCE_COLOR_BUFFER),
        )
        ib = []
        for a, v in ibufs:
            b = opengl.Buffer(v)
            b.buffer_attribute_name = a
            ib.append(b)
        return ib

    def update_instance_buffers(self, positions, colors, position_mask):

        self.update_instance_arrays(positions, colors, position_mask)

        ib = self.instance_buffers
        if len(ib) == 0:
            self.instance_buffers = ib = self.create_instance_buffers()

        for b in ib:
            self.buffer_needs_update(b)

    def update_instance_arrays(self, positions, colors, position_mask):
        sas = positions.shift_and_scale_array()
        np = len(positions)
        im = positions.opengl_matrices() if sas is None and np > 1 else None
        ic = colors if np > 1 or sas is not None else None
        if ic is not None and len(ic) != np:
            # If instance colors array is not same length as positions, resize colors.
            import numpy
            ic = numpy.resize(ic, (np,4))

        pm = position_mask
        if pm is not None:
            im = im[pm, :, :] if im is not None else None
            ic = ic[pm, :] if ic is not None else None
            sas = sas[pm, :] if sas is not None else None

        self.instance_matrices = im
        self.instance_shift_and_scale = sas
        self.instance_colors = ic

    def instance_count(self):
        im = self.instance_matrices
        isas = self.instance_colors
        if im is not None:
            ninst = len(im)
        elif isas is not None:
            ninst = len(isas)
        else:
            ninst = 1
        return ninst

    def buffer_needs_update(self, b):
        self._buffers_need_update.add(b)
        b.up_to_date = False

    def update_buffers(self):
        bu = self._buffers_need_update
        if len(bu) == 0:
            return

        # TODO: This hack makes vertex colors override instance colors.
        bufs = sorted(bu, key = lambda b: b.buffer_attribute_name)

        bi = self.bindings
        for b in bufs:
            if not b.up_to_date:
                data = getattr(self, b.buffer_attribute_name)
                b.update_buffer_data(data)
                b.up_to_date = True
            bi.bind_shader_variable(b)
        bu.clear()

    def activate_bindings(self):
        bi = self.bindings
        if bi is None:
            from . import opengl
            self.bindings = bi = opengl.Bindings()

        bi.activate()
        self.update_buffers()

class Pick:
    '''
    A picked object returned by first_intercept() method of the Drawing class.
    '''
    def __init__(self, distance = None):
        self.distance = distance

    def description(self):
        '''Text description of the picked object.'''
        return None

    def select(self, mode = 'add'):
        '''
        Cause this picked object to be selected ('add' mode), unselected ('subtract' mode)
        or toggle selected ('toggle' mode).
        '''
        pass

class PickedTriangle(Pick):
    '''
    A picked triangle of a drawing.
    '''
    def __init__(self, distance, triangle_number, copy_number, drawing):
        Pick.__init__(self, distance)
        tm = drawing.triangle_mask
        # Convert to from displayed triangle number to all triangles number.
        tnum = triangle_number if tm is None else tm.nonzero()[0][triangle_number]
        self.triangle_number = tnum
        self._copy = copy_number
        self._drawing = drawing

    def description(self):
        d = self.drawing()
        fields = [self.id_string()]
        if d.name is not None:
            fields.append(d.name)
        if len(d.positions) > 1:
            fields.append('copy %d' % self._copy)
        fields.append('triangle %d of %d' % (self.triangle_number, len(d.triangles)))
        desc = ' '.join(fields)
        return desc

    def id_string(self):
        '''
        A text identifier that can be used in commands to specified the
        picked Model. The id number is not a standard attribute
        of Drawing, only of Model which is a subclass of Drawing,
        and is a tuple of integers.
        '''
        d = self.drawing()
        while True:
            if hasattr(d, 'id') and d.id is not None:
                s = '#' + '.'.join(('%d' % id) for id in d.id)
                return s
            if hasattr(d, 'parent') and not d.parent is None:
                d = d.parent
            else:
                break
        return '?'

    def drawing(self):
        return self._drawing

    def select(self, mode = 'add'):
        d = self.drawing()
        pmask = d.selected_positions
        if pmask is None:
            from numpy import zeros, bool
            pmask = zeros((len(d.positions),), bool)
        c = self._copy
        if mode == 'add':
            s = 1
        elif mode == 'subtract':
            s = 0
        elif mode == 'toggle':
            s = not pmask[c]
        pmask[c] = s
        d.selected_positions = pmask

    def is_transparent(self):
        d = self.drawing()
        vc = d.vertex_colors
        if vc is None:
            return d.color[3] < 255
        t = self.triangle_number
        for v in d.triangles[t]:
            if vc[v,3] < 255:
                return True
        return False

class PickedTriangles(Pick):
    '''
    Picked triangles of a drawing.
    '''
    def __init__(self, tmask, drawing):
        Pick.__init__(self)
        self._triangles_mask = tmask
        self._drawing = drawing

    def description(self):
        desc = self._drawing.name
        tm = self._triangles_mask
        nt = tm.sum()
        if nt < len(tm):
            desc += ', %d of %d triangles' % (nt, len(tm))
        return desc

    def id_string(self):
        d = self.drawing()
        return d.id_string() if hasattr(d, 'id_string') else '?'

    def drawing(self):
        return self._drawing

    def select(self, mode = 'add'):
        d = self.drawing()
        if mode == 'add':
            s = True
        elif mode == 'subtract':
            s = False
        elif mode == 'toggle':
            s = (not d.selected)
        d.selected = s


class PickedInstance(Pick):
    '''
    A picked triangle of a drawing.
    '''
    def __init__(self, pmask, drawing):
        Pick.__init__(self)
        self._positions_mask = pmask
        self._drawing = drawing

    def description(self):
        desc = self._drawing.name
        pm = self._posititions_mask
        np = pm.sum()
        if np < len(pm):
            desc += ', %d of %d instances' % (np, len(pm))
        return desc

    def id_string(self):
        d = self.drawing()
        return d.id_string() if hasattr(d, 'id_string') else '?'

    def drawing(self):
        return self._drawing

    def select(self, mode = 'add'):
        d = self.drawing()
        pm = self._positions_mask
        pmask = d.selected_positions
        if pmask is None and mode != 'subtract':
            pmask = pm.copy()
            pmask[:] = d.selected
        if mode == 'add':
            pmask[pm] = 1
        elif mode == 'subtract':
            if pmask is not None:
                pmask[pm] = 0
        elif mode == 'toggle':
            from numpy import logical_xor
            logical_xor(pmask, pm, pmask)
        d.selected_positions = pmask


def rgba_drawing(drawing, rgba, pos=(-1, -1), size=(2, 2)):
    '''
    Make a drawing that is a single rectangle with a texture to show an
    RGBA image on it.
    '''
    from . import opengl
    t = opengl.Texture(rgba)
    d = _texture_drawing(t, pos, size, drawing)
    return d

def position_rgba_drawing(drawing, pos, size):
    '''
    Use specified position and size for rgba drawing, values in fractional window size.
    '''
    x,y = pos
    sx,sy = size
    from numpy import array, float32
    v = array(((x, y, 0),
               (x + sx, y, 0),
               (x + sx, y + sy, 0),
               (x, y + sy, 0)), float32)
    drawing.set_geometry(v, drawing.normals, drawing.triangles)

def _texture_drawing(texture, pos=(-1, -1), size=(2, 2), drawing=None):
    '''
    Make a drawing that is a single rectangle colored with a texture.
    '''
    d = drawing if drawing else Drawing('rgba')
    x, y = pos
    sx, sy = size
    from numpy import array, float32, int32
    va = array(((x, y, 0),
                (x + sx, y, 0),
                (x + sx, y + sy, 0),
                (x, y + sy, 0)), float32)
    ta = array(((0, 1, 2), (0, 2, 3)), int32)
    tc = array(((0, 0), (1, 0), (1, 1), (0, 1)), float32)
    d.set_geometry(va, None, ta)
    d.color = (255, 255, 255, 255)         # Modulates texture values
    d.use_lighting = False
    d.texture_coordinates = tc
    if d.texture is not None:
        d.texture.delete_texture()
    d.texture = texture
    return d

def resize_rgba_drawing(drawing, pos = (-1,-1), size = (2,2)):
    x, y = pos
    sx, sy = size
    from numpy import array, float32
    varray = array(((x, y, 0),
                    (x + sx, y, 0),
                    (x + sx, y + sy, 0),
                    (x, y + sy, 0)), float32)
    drawing.set_geometry(varray, drawing.normals, drawing.triangles)
    
def _draw_texture(texture, renderer):
    d = _texture_drawing(texture)
    d.opaque_texture = True
    draw_overlays([d], renderer)
    
def qimage_to_numpy(qi):
    shape = (qi.height(), qi.width(), 4)
    buf = qi.bits().asstring(qi.byteCount())
    from numpy import uint8, frombuffer
    bgra = frombuffer(buf, uint8).reshape(shape)
    rgba = bgra.copy()
    rgba[:,:,0] = bgra[:,:,2]
    rgba[:,:,2] = bgra[:,:,0]
    rgba = rgba[::-1,:,:]	# flip
    return rgba
