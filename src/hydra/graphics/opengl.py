'''
OpenGL classes
==============

All calls to OpenGL are made through this module.  Currently all OpenGL is done with PyOpenGL.

The Render class manages shader, view matrices, and lighting.  The Buffer class handles object
geometry (vertices, normals, triangles) and colors and texture coordinates.  The Bindings class
defines the connections between Buffers and shader program variables.  The Texture class manages 2D
texture storage.
'''

from OpenGL import GL

class Render:
    '''
    Manage shaders, viewing matrices and lighting parameters to render a scene.
    '''
    def __init__(self):
                
        self.shader_programs = {}
        self.current_shader_program = None

        self.default_capabilities = set((self.SHADER_LIGHTING,self.SHADER_VERTEX_COLORS))
        self.override_capabilities = {}

        self.current_viewport = None
        self.current_projection_matrix = None   # Used when switching shaders
        self.current_model_view_matrix = None   # Used when switching shaders
        self.current_model_matrix = None        # Used for optimizing model view matrix updates
        self.current_view_matrix = None         # Maps scene to camera coordinates

        self.lighting = Lighting()
        self.material = Material()              # Currently there is only a global material

        self.framebuffer_stack = [default_framebuffer()]
        self.mask_framebuffer = None
        self.outline_framebuffer = None
        self.shadow_map_framebuffer = None
        self._silhouette_framebuffer = None

        # Texture warp parameters
        self.warp_center = (0.5, 0.5)
        self.radial_warp_coefficients = (1,0,0,0)
        self.chromatic_warp_coefficients = (1,0,1,0)

        # 3D ambient texture transform from model coordinates to texture coordinates.
        self.ambient_texture_transform = None

        # Maps camera coordinates to shadow map texture coordinates.
        self.shadow_transform = None

        self.single_color = (1,1,1,1)

    def render_size(self):
        fb = self.current_framebuffer()
        x,y,w,h = fb.viewport
        return (w,h)

    # use_shader() option names
    SHADER_LIGHTING = 'USE_LIGHTING'
    SHADER_DEPTH_CUE = 'USE_DEPTH_CUE'
    SHADER_TEXTURE_2D = 'USE_TEXTURE_2D'
    SHADER_TEXTURE_3D_AMBIENT = 'USE_TEXTURE_3D_AMBIENT'
    SHADER_SHADOWS = 'USE_SHADOWS'
    SHADER_RADIAL_WARP = 'USE_RADIAL_WARP'
    SHADER_SHIFT_AND_SCALE = 'USE_INSTANCING_SS'
    SHADER_INSTANCING = 'USE_INSTANCING_44'
    SHADER_TEXTURE_MASK = 'USE_TEXTURE_MASK'
    SHADER_DEPTH_OUTLINE = 'USE_DEPTH_OUTLINE'
    SHADER_VERTEX_COLORS = 'USE_VERTEX_COLORS'

    def set_override_capabilities(self, ocap):
        self.override_capabilities = ocap

    def shader(self, options):
        '''
        Return a shader that supports the specified capabilities.
        The capabilities are specified as keyword options with boolean values.
        The available option names are given by the values of SHADER_LIGHTING,
        SHADER_DEPTH_CUE, SHADER_TEXTURE_2D, SHADER_TEXTURE_3D_AMBIENT,
        SHADER_SHADOWS, SHADER_RADIAL_WARP, SHADER_SHIFT_AND_SCALE,
        SHADER_INSTANCING, SHADER_TEXTURE_MASK, SHADER_DEPTH_OUTLINE, SHADER_VERTEX_COLORS
        '''
        capabilities = self.default_capabilities.copy()
        ocap = self.override_capabilities
        if ocap:
            options = options.copy()
            options.update(ocap)
        for opt,onoff in options.items():
            if onoff:
                capabilities.add(opt)
            else:
                capabilities.discard(opt)

        p = self.opengl_shader(capabilities)
        return p

    def use_shader(self, shader):
        '''
        Set the current shader.
        '''
        if shader == self.current_shader_program:
            return

#        print ('changed shader',
#               self.current_shader_program.capabilities if self.current_shader_program else None,
#               shader.capabilities)
        self.current_shader_program = shader
        c = shader.capabilities
        GL.glUseProgram(shader.program_id)
        if self.SHADER_LIGHTING in c:
            self.set_shader_lighting_parameters()
        if self.SHADER_DEPTH_CUE in c:
            self.set_depth_cue_parameters()
        self.set_projection_matrix()
        self.set_model_matrix()
        if self.SHADER_TEXTURE_2D in c or self.SHADER_TEXTURE_MASK in c or self.SHADER_DEPTH_OUTLINE in c:
            GL.glUniform1i(shader.uniform_id("tex2d"), 0)    # Texture unit 0.
        if self.SHADER_TEXTURE_3D_AMBIENT in c:
            GL.glUniform1i(shader.uniform_id("tex3d"), 0)    # Texture unit 0.
        if self.SHADER_SHADOWS in c:
            GL.glUniform1i(shader.uniform_id("shadow_map"), self.shadow_texture_unit)
            if not self.shadow_transform is None:
                GL.glUniformMatrix4fv(shader.uniform_id("shadow_transform"), 1, False, self.shadow_transform)
        if self.SHADER_RADIAL_WARP in c:
            self.set_radial_warp_parameters()
        if not self.SHADER_VERTEX_COLORS in c:
            self.set_single_color()

    shadow_texture_unit = 1

    def push_framebuffer(self, fb):
        self.framebuffer_stack.append(fb)
        fb.activate()
        self.set_viewport(*fb.viewport)

    def current_framebuffer(self):
        return self.framebuffer_stack[-1]

    def pop_framebuffer(self):
        s = self.framebuffer_stack
        pfb = s.pop()
        if len(s) == 0:
            raise RuntimeError('No framebuffer left on stack.')
        fb = s[-1]
        fb.activate()
        self.set_viewport(*fb.viewport)
        return pfb

    def rendering_to_screen(self):
        return len(self.framebuffer_stack) == 0

    def opengl_shader(self, capabilities = (SHADER_LIGHTING,)):
        'Private.  OpenGL shader program id.'

        ckey = tuple(sorted(capabilities))
        p = self.shader_programs.get(ckey)
        if not p is None:
            return p

        p = Shader(capabilities)
        self.shader_programs[ckey] = p

        return p
        
    def set_projection_matrix(self, pm = None):
        '''
        Set the shader to use the given 4x4 OpenGL projection matrix.
        If no matrix is specified use the last specified one.
        '''
        if pm is None:
            if self.current_projection_matrix is None:
                return
            pm = self.current_projection_matrix
        else:
            self.current_projection_matrix = pm
        p = self.current_shader_program
        if not p is None:
            GL.glUniformMatrix4fv(p.uniform_id('projection_matrix'), 1, False, pm)

    def set_view_matrix(self, vm):
        '''Set the camera view matrix, mapping scene to camera coordinates.'''
        self.current_view_matrix = vm
        self.current_model_matrix = None
        self.current_model_view_matrix = None

    def set_model_matrix(self, model_matrix = None):
        '''
        Set the shader model view using the given model matrix and previously set view matrix.
        If no matrix is specified, the shader gets the last used one model view matrix.
        '''

        if model_matrix is None:
            mv4 = self.current_model_view_matrix
            if mv4 is None:
                return
        else:
            # TODO: optimize check of same model matrix.
            cmm = self.current_model_matrix
            if cmm:
                if (model_matrix.is_identity() and cmm.is_identity()) or model_matrix.same(cmm):
                    return
            self.current_model_matrix = model_matrix
            # TODO: optimize matrix multiply.  Rendering bottleneck with 200 models open.
            mv4 = (self.current_view_matrix * model_matrix).opengl_matrix()
            self.current_model_view_matrix = mv4

        p = self.current_shader_program
        if not p is None:
            var_id = p.uniform_id('model_view_matrix')
            # Note: Getting about 5000 glUniformMatrix4fv() calls per second on 2013 Mac hardware.
            # This can be a rendering bottleneck for large numbers of models or instances.
            GL.glUniformMatrix4fv(var_id, 1, False, mv4)
            if not self.lighting.move_lights_with_camera:
                self.set_shader_lighting_parameters()

    def set_shader_lighting_parameters(self):
        'Private. Sets shader lighting variables using the lighting parameters object given in the contructor.'

        p = self.current_shader_program.program_id
        lp = self.lighting
        mp = self.material

        move = None if lp.move_lights_with_camera else self.current_view_matrix

        # Key light
        key_light_dir = GL.glGetUniformLocation(p, b"key_light_direction")
        kld = tuple(move.apply_without_translation(lp.key_light_direction)) if move else lp.key_light_direction
        GL.glUniform3f(key_light_dir, *kld)
        key_diffuse = GL.glGetUniformLocation(p, b"key_light_diffuse_color")
        ds = mp.diffuse_reflectivity
        kdc = tuple(ds*c for c in lp.key_light_color)
        GL.glUniform3f(key_diffuse, *kdc)

        # Key light specular
        key_specular = GL.glGetUniformLocation(p, b"key_light_specular_color")
        ss = mp.specular_reflectivity
        ksc = tuple(ss*c for c in lp.key_light_color)
        GL.glUniform3f(key_specular, *ksc)
        key_shininess = GL.glGetUniformLocation(p, b"key_light_specular_exponent")
        GL.glUniform1f(key_shininess, mp.specular_exponent)

        # Fill light
        fill_light_dir = GL.glGetUniformLocation(p, b"fill_light_direction")
        fld = tuple(move.apply_without_translation(lp.fill_light_direction)) if move else lp.fill_light_direction 
        GL.glUniform3f(fill_light_dir, *fld)
        fill_diffuse = GL.glGetUniformLocation(p, b"fill_light_diffuse_color")
        fdc = tuple(ds*c for c in lp.fill_light_color)
        GL.glUniform3f(fill_diffuse, *fdc)

        # Ambient light
        ambient = GL.glGetUniformLocation(p, b"ambient_color")
        ams = mp.ambient_reflectivity
        ac = tuple(ams*c for c in lp.ambient_light_color)
        GL.glUniform3f(ambient, *ac)

    def set_depth_cue_parameters(self):
        'Private. Sets shader depth variables using the lighting parameters object given in the contructor.'

        p = self.current_shader_program.program_id
        lp = self.lighting

        dc_distance = GL.glGetUniformLocation(p, b"depth_cue_distance")
        GL.glUniform1f(dc_distance, lp.depth_cue_distance)
        dc_darkest = GL.glGetUniformLocation(p, b"depth_cue_darkest")
        GL.glUniform1f(dc_darkest, lp.depth_cue_darkest)

    def set_single_color(self, color = None):
        '''
        Set the OpenGL shader color for shader single color mode.
        '''
        if not color is None:
            self.single_color = color
        p = self.current_shader_program.program_id
        c = GL.glGetUniformLocation(p, b"color")
        GL.glUniform4fv(c, 1, self.single_color)

    def set_radial_warp_parameters(self):
        p = self.current_shader_program.program_id
        wcenter = GL.glGetUniformLocation(p, b"warp_center")
        GL.glUniform2fv(wcenter, 1, self.warp_center)
        rcoef = GL.glGetUniformLocation(p, b"radial_warp")
        GL.glUniform4fv(rcoef, 1, self.radial_warp_coefficients)
        ccoef = GL.glGetUniformLocation(p, b"chromatic_warp")
        GL.glUniform4fv(ccoef, 1, self.chromatic_warp_coefficients)

    def set_ambient_texture_transform(self, tf):
        # Transform from model coordinates to ambient texture coordinates.
        p = self.current_shader_program
        m = tf.opengl_matrix()
        GL.glUniformMatrix4fv(p.uniform_id("ambient_tex3d_transform"), 1, False, m)

    def set_shadow_transform(self, tf):
        # Transform from camera coordinates to shadow map texture coordinates.
        self.shadow_transform = m = tf.opengl_matrix()
        p = self.current_shader_program
        if self.SHADER_SHADOWS in p.capabilities:
            GL.glUniformMatrix4fv(p.uniform_id("shadow_transform"), 1, False, m)

    def opengl_version(self):
        'String description of the OpenGL version for the current context.'
        return GL.glGetString(GL.GL_VERSION).decode('utf-8')

    def support_stereo(self):
        'Return if sequential stereo is supported.'
        return GL.glGetBoolean(GL.GL_STEREO)

    def initialize_opengl(self, width, height):
        'Create an initial vertex array object.'

        # OpenGL 3.2 core profile requires a bound vertex array object
        # for drawing, or binding shader attributes to VBOs.  Mac 10.8
        # gives an error if no VAO is bound when glCompileProgram() called.
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)

        fb = default_framebuffer()
        fb.width, fb.height = width, height
        self.set_viewport(0,0,width,height)

    def set_viewport(self, x, y, w, h):
        'Set the OpenGL viewport.'
        if (x,y,w,h) != self.current_viewport:
            GL.glViewport(x, y, w, h)
            self.current_viewport = (x,y,w,h)
        fb = self.current_framebuffer()
        fb.viewport = (x,y,w,h)

    def full_viewport(self):
        fb = self.current_framebuffer()
        self.set_viewport(0,0,fb.width,fb.height)

    def set_background_color(self, rgba):
        'Set the OpenGL clear color.'
        r,g,b,a = rgba
        GL.glClearColor(r, g, b, a)

    def draw_background(self):
        'Draw the background color and clear the depth buffer.'
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

    def enable_depth_test(self, enable):
        'Enable OpenGL depth testing.'
        if enable:
            GL.glEnable(GL.GL_DEPTH_TEST)
        else:
            GL.glDisable(GL.GL_DEPTH_TEST)

    def enable_blending(self, enable):
        'Enable OpenGL alpha blending.'
        if enable:
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        else:
            GL.glDisable(GL.GL_BLEND)

    def draw_transparent(self, draw_depth, draw):
        '''
        Render using single-layer transparency. This is a two-pass drawing.
        In the first pass is only sets the depth buffer, but not colors, and in
        the second path it draws the colors for pixels at or in front of the
        recorded depths.  The draw_depth and draw routines, taking no arguments
        perform the actual drawing, and are invoked by this routine after setting
        the appropriate OpenGL color and depth drawing modes.
        '''
        # Single layer transparency
        GL.glColorMask(GL.GL_FALSE, GL.GL_FALSE, GL.GL_FALSE, GL.GL_FALSE)
        draw_depth()
        GL.glColorMask(GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE)
        GL.glDepthFunc(GL.GL_LEQUAL)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        draw()
        GL.glDepthFunc(GL.GL_LESS)

    IMAGE_FORMAT_RGBA32 = 'rgba32'
    IMAGE_FORMAT_RGBA8 = 'rgba8'
    IMAGE_FORMAT_RGB32 = 'rgb32'

    def frame_buffer_image(self, w, h, format = IMAGE_FORMAT_RGBA8):
        '''
        Return the current frame buffer image as a numpy array of size (h,w) for 32-bit
        formats or (h,w,4) for 8-bit formats where w and h are the framebuffer width and height.
        Array index 0,0 is at the bottom left corner of the OpenGL viewport for RGB32 format
        and at the upper left corner for the other formats.  For 32-bit formats the array values
        are uint32 and contain 8-bit red, green, and blue values is the low 24 bits for RGB32,
        and 8-bit red, green, blue and alpha for RGBA32.  The RGBA8 format has uint8 values.
        '''

        if format == self.IMAGE_FORMAT_RGBA32:
            from numpy import empty, uint32
            rgba = empty((h,w),uint32)
            GL.glReadPixels(0,0,w,h,GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, rgba)
            return rgba
        elif format == self.IMAGE_FORMAT_RGB32:
            rgba = self.frame_buffer_image(w, h, self.IMAGE_FORMAT_RGBA32)
            rgba >>= 8
            rgb = rgba[::-1,:].copy()
            return rgb
        elif format == self.IMAGE_FORMAT_RGBA8:
            rgba = self.frame_buffer_image(w, h, self.IMAGE_FORMAT_RGBA32)
            from numpy import little_endian, uint8
            if little_endian:
                rgba.byteswap(True) # in place
                rgba8 = rgba.view(uint8).reshape((h,w,4))
                return rgba8

    def set_stereo_buffer(self, eye_num):
        '''Set the draw and read buffers for the left eye (0) or right eye (1).'''
        self.full_viewport()
        if not self.rendering_to_screen():
            return
        b = GL.GL_BACK_LEFT if eye_num == 0 else GL.GL_BACK_RIGHT
        GL.glDrawBuffer(b)
        GL.glReadBuffer(b)

    def set_mono_buffer(self):
        '''Set the draw and read buffers for mono rendering.'''
        self.full_viewport()
        if not self.rendering_to_screen():
            return
        b = GL.GL_BACK
        GL.glDrawBuffer(b)
        GL.glReadBuffer(b)

    def start_rendering_shadowmap(self, center, radius, camera_view, size = 1024, depth_bias = 0.005):

        # Set projection matrix to be orthographic and span all models.
        from ..geometry.place import translation, scale, orthonormal_frame
        pm = scale((1/radius,1/radius,-1/radius))*translation((0,0,radius))       # orthographic projection along z
        self.set_projection_matrix(pm.opengl_matrix())

        # Make a framebuffer for depth texture rendering
        fb = self.shadow_map_framebuffer
        if fb is None:
            dt = Texture()
            dt.initialize_depth((size,size))
            fb = Framebuffer(depth_texture = dt)
            if not fb.valid():
                return           # Requested size exceeds framebuffer limits
            self.shadow_map_framebuffer = fb

        # Compute the view matrix looking along the light direction.
        kl = self.lighting.key_light_direction
        ld = camera_view.apply_without_translation(kl) # Light direction in scene coords.
        lv = translation(center - radius*ld) * orthonormal_frame(-ld)   # Light view frame
        lvinv = lv.inverse()	# Scene to light view coordinates

        # Set the transform mapping camera to depth texture coordinates.
        ntf = translation((0.5,0.5,0.5-depth_bias))*scale(0.5)    # (-1,1) normalized device coords to (0,1) texture coords.
        stf = ntf * pm * lvinv                       # Scene to shadowmap coordinates

        # Make sure depth texture is not bound from previous drawing so that it is not
        # used for rendering shadows while the depth texture is being written.
        # TODO: The depth rendering should not render colors or shadows.
        dt = fb.depth_texture
        dt.unbind_texture(self.shadow_texture_unit)

        # Draw the models recording depth in light direction, i.e. calculate the shadow map.
        self.push_framebuffer(fb)
        self.draw_background()             # Clear depth buffer

        return lvinv, stf

    def finish_rendering_shadowmap(self):

        self.pop_framebuffer()

        # Bind the depth texture for computing shadows.
        # TODO: Use different texture unit to avoid conflict with texture coloring.
        fb = self.shadow_map_framebuffer
        return fb.depth_texture

    def start_rendering_outline(self):

        fb = self.current_framebuffer()
        mfb = self.make_mask_framebuffer()
        self.push_framebuffer(mfb)
        self.set_background_color((0,0,0,0))
        self.draw_background()
        # Use flat single color rendering.
        self.set_override_capabilities({self.SHADER_VERTEX_COLORS:False,
                                        self.SHADER_LIGHTING:False,
                                        self.SHADER_TEXTURE_2D:False,
                                        self.SHADER_TEXTURE_3D_AMBIENT:False,
                                        self.SHADER_SHADOWS:False,
                                        })
        self.set_depth_range(0,0.99999)      # Depth test GL_LEQUAL results in z-fighting
        self.copy_from_framebuffer(fb, color = False)      # Copy depth to outline framebuffer

    def finish_rendering_outline(self):

        self.pop_framebuffer()
        self.set_override_capabilities({})
        self.set_depth_range(0,1)
        t = self.mask_framebuffer.color_texture
        self.draw_texture_mask_outline(t)

    def make_mask_framebuffer(self):
        size = self.render_size()
        mfb = self.mask_framebuffer
        w,h = size
        if mfb and mfb.width == w and mfb.height == h:
            return mfb
        if mfb:
            mfb.delete()
        t = Texture()
        t.initialize_8_bit(size)
        self.mask_framebuffer = mfb = Framebuffer(color_texture = t)
        return mfb

    def make_outline_framebuffer(self, size):
        ofb = self.outline_framebuffer
        w,h = size
        if ofb and ofb.width == w and ofb.height == h:
            return ofb
        if ofb:
            ofb.delete()
        t = Texture()
        t.initialize_8_bit(size)
        self.outline_framebuffer = ofb = Framebuffer(color_texture = t, depth = False)
        return ofb

    def draw_texture_mask_outline(self, texture, color = (0,1,0,1)):

        # Draw to a new texture 4 shifted copies of texture and erase an unshifted copy to produce an outline.
        ofb = self.make_outline_framebuffer(texture.size)
        self.push_framebuffer(ofb)
        self.set_background_color((0,0,0,0))
        self.draw_background()

        # Render region with texture red > 0.
        # Texture map a full-screen quad to blend texture with frame buffer.
        tc = Texture_Window(self, self.SHADER_TEXTURE_MASK)
        texture.bind_texture()

        # Draw 4 shifted copies of mask
        w,h = texture.size
        dx, dy = 1.0/w, 1.0/h
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE_MINUS_SRC_ALPHA)
        self.set_texture_mask_color((1,1,1,1))
        for xs,ys in ((-dx,-dy), (dx,-dy), (dx,dy), (-dx,dy)):
            tc.draw(xshift = xs, yshift = ys)

        # Erase unshifted copy of mask
        GL.glBlendFunc(GL.GL_ONE_MINUS_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        tc.draw()

        # Now outline is in texture of the outline framebuffer
        outline = ofb.color_texture
        self.pop_framebuffer()

        # Draw outline on original framebuffer
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        outline.bind_texture()
        self.set_texture_mask_color(color)
        tc.draw()
        
    def set_texture_mask_color(self, color):

        p = self.current_shader_program.program_id
        mc = GL.glGetUniformLocation(p, b"color")
        GL.glUniform4fv(mc, 1, color)

    def set_depth_range(self, min, max):
#        GL.glDepthFunc(GL.GL_LEQUAL)   # Get z-fighting with screen depth copied to framebuffer object on Mac/Nvidia
        GL.glDepthRange(min, max)

    def start_silhouette_drawing(self):
        fb = self.silhouette_framebuffer(self.render_size())
        self.push_framebuffer(fb)

    def finish_silhouette_drawing(self, perspective_near_far_ratio):
        fb = self.pop_framebuffer()
        self.copy_from_framebuffer(fb, depth = False)
        self.draw_depth_outline(fb.depth_texture, perspective_near_far_ratio = perspective_near_far_ratio)

    def silhouette_framebuffer(self, size = None):
        sfb = self._silhouette_framebuffer
        if size is None:
            return sfb
        if sfb and (size[0] != sfb.width or size[1] != sfb.height):
            sfb.delete()
            sfb = None
        if sfb is None:
            dt = Texture()
            dt.initialize_depth(size, depth_compare_mode = False)
            self._silhouette_framebuffer = sfb = Framebuffer(depth_texture = dt)
        return sfb

    def draw_depth_outline(self, depth_texture, color = (0,0,0,1), depth_jump = 0.01, perspective_near_far_ratio = 1):

        # Render pixels with depth less than neighbor pixel by at least depth_jump
        # Texture map a full-screen quad to blend depth jump pixels with frame buffer.
        tc = Texture_Window(self, self.SHADER_DEPTH_OUTLINE)
        depth_texture.bind_texture()

        # Draw 4 shifted copies of mask
        w,h = depth_texture.size
        dx, dy = 1.0/w, 1.0/h
        self.enable_depth_test(False)
        self.enable_blending(True)
        GL.glDepthMask(False)   # Disable depth write
        self.set_depth_outline_color(color)
        for xs,ys in ((-dx,-dy), (dx,-dy), (dx,dy), (-dx,dy)):
            self.set_depth_outline_shift_and_jump(xs, ys, depth_jump, perspective_near_far_ratio)
            tc.draw()
        GL.glDepthMask(True)
        self.enable_blending(False)
        self.enable_depth_test(True)

    def set_depth_outline_color(self, color):

        p = self.current_shader_program.program_id
        mc = GL.glGetUniformLocation(p, b"color")
        GL.glUniform4fv(mc, 1, color)

    def set_depth_outline_shift_and_jump(self, xs, ys, depth_jump, perspective_near_far_ratio):

        p = self.current_shader_program.program_id
        mc = GL.glGetUniformLocation(p, b"depth_shift_and_jump")
        GL.glUniform4fv(mc, 1, (xs,ys,depth_jump,perspective_near_far_ratio))

    def copy_from_framebuffer(self, framebuffer, color = True, depth = True):
        # Copy current framebuffer contents to another framebuffer.
        # This leaves read and draw framebuffers set to the current framebuffer.
        cfb = self.current_framebuffer()
        GL.glBindFramebuffer(GL.GL_READ_FRAMEBUFFER, framebuffer.fbo)
        GL.glBindFramebuffer(GL.GL_DRAW_FRAMEBUFFER, cfb.fbo)
        what = GL.GL_COLOR_BUFFER_BIT if color else 0
        if depth:
            what |= GL.GL_DEPTH_BUFFER_BIT
        w,h = framebuffer.width, framebuffer.height
        GL.glBlitFramebuffer(0,0,w,h, 0,0,w,h, what, GL.GL_NEAREST)
        # Restore read buffer 
        GL.glBindFramebuffer(GL.GL_READ_FRAMEBUFFER, cfb.fbo)

    def finish_rendering(self):
        GL.glFinish()

class Framebuffer:

    def __init__(self, width = None, height = None,
                 color = True, color_texture = None,
                 depth = True, depth_texture = None):

        if not width is None and not height is None:
            w,h = width,height
        elif not color_texture is None:
            w,h = color_texture.size
        elif not depth_texture is None:
            w,h = depth_texture.size
        else:
            w = h = None
        
        self.width = w
        self.height = h
        self.viewport = None if w is None else (0,0,w,h)

        self.color_texture = color_texture
        self.color_rb = None if not color or color_texture or w is None else self.color_renderbuffer(w,h)
        self.depth_texture = depth_texture
        self.depth_rb = None if not depth or depth_texture or w is None else self.depth_renderbuffer(w,h)

        if w is None:
            fbo = 0
        elif not self.valid_size(w,h):
            fbo = None
        else:
            fbo = self.create_fbo(color_texture or self.color_rb,
                                  depth_texture or self.depth_rb)
        self.fbo = fbo

    def __del__(self):

        self.delete()

    def valid_size(self, width, height):

        max_rb_size = GL.glGetInteger(GL.GL_MAX_RENDERBUFFER_SIZE)
        max_tex_size = GL.glGetInteger(GL.GL_MAX_TEXTURE_SIZE)
        max_size = min(max_rb_size, max_tex_size)
        return width <= max_size and height <= max_size
            
    def color_renderbuffer(self, width, height):

        color_rb = GL.glGenRenderbuffers(1)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, color_rb)
        GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_RGB8, width, height)
        return color_rb

        return color_rb, depth_rb, fbo

    def depth_renderbuffer(self, width, height):

        depth_rb = GL.glGenRenderbuffers(1)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, depth_rb)
        GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT24, width, height)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, 0)
        return depth_rb

    def create_fbo(self, color_buf, depth_buf):

        fbo = GL.glGenFramebuffers(1)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo)

        if isinstance(color_buf, Texture):
            level = 0
            GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                      GL.GL_TEXTURE_2D, color_buf.id, level)
        else:
            GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                         GL.GL_RENDERBUFFER, color_buf)

        if isinstance(depth_buf, Texture):
            level = 0
            GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                      GL.GL_TEXTURE_2D, depth_buf.id, level)
        elif not depth_buf is None:
            GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                         GL.GL_RENDERBUFFER, depth_buf)

        status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
        if status != GL.GL_FRAMEBUFFER_COMPLETE:
            self.delete()
            # TODO: Need to rebind previous framebuffer.
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
            return None

        return fbo

    def valid(self):
        return not self.fbo is None

    def activate(self):
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo)

    def delete(self):
        if self.fbo is None:
            return

        if not self.color_rb is None:
            GL.glDeleteRenderbuffers(1, (self.color_rb,))
        if not self.depth_rb is None:
            GL.glDeleteRenderbuffers(1, (self.depth_rb,))
        GL.glDeleteFramebuffers(1, (self.fbo,))
        self.color_rb = self.depth_rb = self.fbo = None

_default_framebuffer = None
def default_framebuffer():
    global _default_framebuffer
    if _default_framebuffer is None:
        _default_framebuffer = fb = Framebuffer(color = False, depth = False)
    return _default_framebuffer

class Lighting:
    '''
    Lighting parameters specifying colors and directions of two lights:
    a key (main) light, and a fill light, as well as ambient light color.
    Directions are unit vectors in camera coordinates (x right, y up, z opposite camera view).
    Colors are R,G,B float values in the range 0-1.
    '''

    def __init__(self):

        from numpy import array, float32
        self.key_light_direction = array((.577,-.577,-.577), float32)    # Should have unit length
        '''Direction key light shines in.'''

        self.key_light_color = (1,1,1)
        '''Key light color.'''

        self.fill_light_direction = array((-.2,-.2,-.959), float32)        # Should have unit length
        '''Direction fill light shines in.'''

        self.fill_light_color = (.5,.5,.5)
        '''Fill light color.'''

        self.ambient_light_color = (1,1,1)
        '''Ambient light color.'''

        self.depth_cue_distance = 15.0  # Distance where dimming begins (Angstroms)
        self.depth_cue_darkest = 0.2    # Smallest dimming factor

        self.move_lights_with_camera = True
        '''Whether lights are attached to camera, or fixed in the scene.'''

class Material:
    '''
    Surface properties that control the reflection of light.
    '''
    def __init__(self):
        
        self.ambient_reflectivity = 0.3
        '''Fraction of ambient light reflected.  Ambient light comes from all directions
        and the amount reflected does not depend on the surface orientation of view direction.'''

        self.diffuse_reflectivity = 0.8
        '''Fraction of direction light reflected diffusely, that is
        depending on light angle to surface but not viewing direction.'''

        self.specular_reflectivity = 0.3
        '''Fraction of directional key light reflected specularly, that is
        depending how close reflected light direction is to the viewing direction.'''

        self.specular_exponent = 30
        '''Controls the spread of specular light. The specular exponent is a single float value
        used as an exponent e with reflected intensity scaled by cosine(a) ** e where a is the angle
        between the reflected light and the view direction. A typical value for e is 30.'''
        
class Bindings:
    '''
    Use an OpenGL vertex array object to save buffer bindings.
    The bindings are for a specific shader program since they use the shader variable ids.
    '''
    def __init__(self, shader = None):
        self.shader = shader
        self.vao_id = GL.glGenVertexArrays(1)
        self.bound_attr_ids = {}        # Maps buffer to list of ids

    def __del__(self):
        'Delete the OpenGL vertex array object.'
        GL.glDeleteVertexArrays(1, (self.vao_id,))

    def activate(self):
        'Activate the bindings by binding the OpenGL vertex array object.'
        GL.glBindVertexArray(self.vao_id)

    def bind_shader_variable(self, buffer):
        '''
        Bind the shader variable associated with buffer to the buffer data.
        This enables the OpenGL attribute array, and enables instancing if
        the buffer is instance data.
        '''
        buf_id = buffer.opengl_buffer
        btype = buffer.buffer_type
        if buf_id is None:
            # Unbind already bound variable
            for a in self.bound_attr_ids.get(buffer,[]):
                GL.glDisableVertexAttribArray(a)
            self.bound_attr_ids[buffer] = []
            if btype == GL.GL_ELEMENT_ARRAY_BUFFER:
                GL.glBindBuffer(btype, 0)
            return

        vname = buffer.shader_variable_name
        if vname is None:
            if btype == GL.GL_ELEMENT_ARRAY_BUFFER:
                # Element array buffer binding is saved in VAO.
                GL.glBindBuffer(btype, buf_id)
            return

        shader = self.shader
        if not buffer.shader_has_required_capabilities(shader):
            return

        attr_id = shader.attribute_id(vname)
        if attr_id == -1:
            raise NameError('Failed to find shader attribute %s\n in shader with capabilites = %s'
                            % (vname, str(shader.capabilities)))
        nattr = buffer.attribute_count()
        ncomp = buffer.component_count()
        from numpy import float32, uint8
        gtype = {float32:GL.GL_FLOAT,
                 uint8:GL.GL_UNSIGNED_BYTE}[buffer.value_type]
        normalize = GL.GL_TRUE if buffer.normalize else GL.GL_FALSE

        GL.glBindBuffer(btype, buf_id)
        if nattr == 1:
            GL.glVertexAttribPointer(attr_id, ncomp, gtype, normalize, 0, None)
            GL.glEnableVertexAttribArray(attr_id)
            glVertexAttribDivisor(attr_id, 1 if buffer.instance_buffer else 0)
            self.bound_attr_ids[buffer] = [attr_id]
        else:
            # Matrices use multiple vector attributes
            esize = buffer.array_element_bytes()
            abytes = ncomp * esize
            stride = nattr * abytes
            import ctypes
            for a in range(nattr):
                # Pointer arg must be void_p, not an integer.
                p = ctypes.c_void_p(a*abytes)
                GL.glVertexAttribPointer(attr_id+a, ncomp, gtype, normalize, stride, p)
                GL.glEnableVertexAttribArray(attr_id+a)
                glVertexAttribDivisor(attr_id+a, 1 if buffer.instance_buffer else 0)
            self.bound_attr_ids[buffer] = [attr_id+a for a in range(nattr)]
        GL.glBindBuffer(btype, 0)

        return attr_id

from numpy import uint8, uint32, float32

class Buffer_Type:
    def __init__(self, shader_variable_name,
                 buffer_type = GL.GL_ARRAY_BUFFER, value_type = float32,
                 normalize = False, instance_buffer = False, requires_capabilities = ()):
        self.shader_variable_name = shader_variable_name
        self.buffer_type = buffer_type
        self.value_type = value_type
        self.normalize = normalize
        self.instance_buffer = instance_buffer
        self.requires_capabilities = requires_capabilities      # Requires at least one of these.

# Buffer types with associated shader variable names
VERTEX_BUFFER = Buffer_Type('position')
NORMAL_BUFFER = Buffer_Type('normal', requires_capabilities = (Render.SHADER_LIGHTING,))
VERTEX_COLOR_BUFFER = Buffer_Type('vcolor', value_type = uint8, normalize = True,
                                  requires_capabilities = (Render.SHADER_VERTEX_COLORS,))
INSTANCE_SHIFT_AND_SCALE_BUFFER = Buffer_Type('instanceShiftAndScale', instance_buffer = True)
INSTANCE_MATRIX_BUFFER = Buffer_Type('instancePlacement', instance_buffer = True)
INSTANCE_COLOR_BUFFER = Buffer_Type('vcolor', instance_buffer = True, value_type = uint8, normalize = True,
                                    requires_capabilities = (Render.SHADER_VERTEX_COLORS,))
TEXTURE_COORDS_2D_BUFFER = Buffer_Type('tex_coord_2d', requires_capabilities = (Render.SHADER_TEXTURE_2D,
                                                                                Render.SHADER_TEXTURE_MASK,
                                                                                Render.SHADER_DEPTH_OUTLINE))
ELEMENT_BUFFER = Buffer_Type(None, buffer_type = GL.GL_ELEMENT_ARRAY_BUFFER, value_type = uint32)

class Buffer:
    '''
    Create an OpenGL buffer of vertex data such as vertex positions, normals, or colors,
    or per-instance data (e.g. color per sphere) or an element buffer for specifying which
    primitives (triangles, lines, points) to draw.  Vertex data buffers can be attached to
    a specific shader variable.
    '''
    def __init__(self, buffer_type):

        t = buffer_type
        self.shader_variable_name = t.shader_variable_name
        self.opengl_buffer = None
        self.buffered_array = None  # numpy array for vbo
        self.buffered_data = None   # data need not be numpy array
        self.value_type = t.value_type
        self.buffer_type = t.buffer_type
        self.normalize = t.normalize
        self.instance_buffer = t.instance_buffer
        self.requires_capabilities = t.requires_capabilities

    def __del__(self):
        self.delete_buffer()

    def delete_buffer(self):
        'Delete the OpenGL buffer object.'

        if self.opengl_buffer is None:
            return
        GL.glDeleteBuffers(1, [self.opengl_buffer])
        self.opengl_buffer = None
        self.buffered_array = None
        self.buffered_data = None

    def attribute_count(self):
        'Private.'
        # matrix attributes use multiple attribute ids
        barray = self.buffered_array
        if barray is None:
            return 0
        bshape = barray.shape
        nattr = 1 if len(bshape) == 2 else bshape[1]
        return nattr

    def component_count(self):
        'Private.'
        barray = self.buffered_array
        if barray is None:
            return 0
        ncomp = barray.shape[-1]
        return ncomp
  
    def array_element_bytes(self):
        'Private.'
        barray = self.buffered_array
        return 0 if barray is None else barray.itemsize

    def update_buffer_data(self, data):
        '''
        Update the buffer with data supplied by a numpy array and bind it to the
        associated shader variable.
        '''
        bdata = self.buffered_data
        if data is bdata:
            return False

        replace_buffer = (data is None or bdata is None or data.shape != bdata.shape)
        if replace_buffer:
            self.delete_buffer()

        if not data is None:
            b = GL.glGenBuffers(1) if replace_buffer else self.opengl_buffer
            btype = self.buffer_type
            GL.glBindBuffer(btype, b)
            d = data if data.dtype == self.value_type else data.astype(self.value_type)
            size = d.size * d.itemsize        # Bytes
            if replace_buffer:
                GL.glBufferData(btype, size, d, GL.GL_STATIC_DRAW)
            else:
                # TODO:PyOpenGL-20130502 has glBufferSubData() has python 3 bug, long undefined
                #   So use size None so size is computed from array.
                GL.glBufferSubData(btype, 0, None, d)
            GL.glBindBuffer(btype, 0)
            self.opengl_buffer = b
            self.buffered_array = d
            self.buffered_data = data

        return True

    # Element types for Buffer draw_elements()
    triangles = GL.GL_TRIANGLES
    lines = GL.GL_LINES
    points = GL.GL_POINTS

    def draw_elements(self, element_type = GL.GL_TRIANGLES, ninst = None):
        '''
        Draw primitives using this buffer as the element buffer.
        All the required buffers are assumed to be already bound using a
        vertex array object.
        '''
        # Don't bind element buffer since it is bound by VAO.
        ne = self.buffered_array.size
        if ninst is None:
            GL.glDrawElements(element_type, ne, GL.GL_UNSIGNED_INT, None)
        else:
            glDrawElementsInstanced(element_type, ne, GL.GL_UNSIGNED_INT, None, ninst)

    def shader_has_required_capabilities(self, shader):
        if not self.requires_capabilities:
            return True
        # Require at least one capability of list
        for cap in self.requires_capabilities:
            if cap in shader.capabilities:
                return True
        return False

def glDrawElementsInstanced(mode, count, etype, indices, ninst):
    'Private. Handle old or defective OpenGL instanced drawing.'
    if bool(GL.glDrawElementsInstanced):
        # OpenGL 3.1 required for this call.
        GL.glDrawElementsInstanced(mode, count, etype, indices, ninst)
    else:
        from OpenGL.GL.ARB.draw_instanced import glDrawElementsInstancedARB
        if not bool(glDrawElementsInstancedARB):
            # Mac 10.6 does not list draw_instanced as an extension using OpenGL 3.2
            from .pyopengl_draw_instanced import glDrawElementsInstancedARB
            glDrawElementsInstancedARB(mode, count, etype, indices, ninst)

def glVertexAttribDivisor(attr_id, d):
    'Private. Handle old or defective OpenGL attribute divisor.'
    if bool(GL.glVertexAttribDivisor):
        GL.glVertexAttribDivisor(attr_id, d)  # OpenGL 3.3
    else:
        from OpenGL.GL.ARB.instanced_arrays import glVertexAttribDivisorARB
        glVertexAttribDivisorARB(attr_id, d)

class Shader:
    '''Private. OpenGL shader program with specified capabilities.'''

    def __init__(self, capabilities):

        self.capabilities = capabilities
        self.program_id = self.compile_shader(capabilities)
        self.uniform_ids = {}
        self.attribute_ids = {}
        
    def uniform_id(self, name):
        uids = self.uniform_ids
        if name in uids:
            uid = uids[name]
        else:
            p = self.program_id
            uids[name] = uid = GL.glGetUniformLocation(p, name.encode('utf-8'))
        return uid

    def attribute_id(self, name):
        aids = self.attribute_ids
        if name in aids:
            aid = aids[name]
        else:
            p = self.program_id
            aids[name] = aid = GL.glGetAttribLocation(p, name.encode('utf-8'))
#            print('attrib id for %s is %d, shader %d, cap %s' % (name, aid, p, self.capabilities))
        return aid

    def compile_shader(self, capabilities):

        from os.path import dirname, join
        d = dirname(__file__)
        f = open(join(d,'vertexShader.txt'), 'r')
        vshader = insert_define_macros(f.read(), capabilities)
        f.close()

        f = open(join(d,'fragmentShader.txt'), 'r')
        fshader = insert_define_macros(f.read(), capabilities)
        f.close()

        from OpenGL.GL import shaders
        vs = shaders.compileShader(vshader, GL.GL_VERTEX_SHADER)
        fs = shaders.compileShader(fshader, GL.GL_FRAGMENT_SHADER)

        prog_id = shaders.compileProgram(vs, fs)

        # msg = (('Compiled shader %d,\n'
        #        ' capbilities %s,\n'
        #        ' vertex shader compile info log\n'
        #        ' %s\n'
        #        ' fragment shader compile info log\n'
        #        ' %s\n'
        #        ' program link info log\n'
        #        ' %s')
        #        % (prog_id, capabilities,
        #           GL.glGetShaderInfoLog(vs), GL.glGetShaderInfoLog(fs), GL.glGetProgramInfoLog(prog_id)))
        # print(msg)

        return prog_id

# Add #define lines after #version line of shader
def insert_define_macros(shader, capabilities):
    'Private. Puts "#define" statements in shader program templates to specify shader capabilities.'
    defs = '\n'.join('#define %s 1' % c for c in capabilities)
    v = shader.find('#version')
    eol = shader[v:].find('\n')+1
    s = shader[:eol] + defs + '\n' + shader[eol:]
    return s

class Texture:
    '''
    Create an OpenGL 1d, 2d, or 3d texture from a numpy array.  For a N dimensional texture
    the data array can be N or N+1 dimensional.  For example, for 2d shape (h,w,c) or (h,w)
    where w and h are the texture width and height and c is the number of color components.
    If the data array is 2-dimensional, the values must be 32-bit RGBA8.  If the data
    array is 3 dimensional the texture format is GL_RED, GL_RG, GL_RGB, or GL_RGBA depending
    on whether c is 1, 2, 3 or 4 and only value types of uint8 or float32 are allowed and
    texture of type GL_UNSIGNED_BYTE or GL_FLOAT is created.  Clamp to edge mode and
    nearest interpolation is set.  The c = 2 mode uses the second component as alpha and
    the first componet for red,green,blue.
    '''
    def __init__(self, data = None, dimension = 2):

        self.id = None
        self.dimension = dimension
        self.gl_target = (GL.GL_TEXTURE_1D, GL.GL_TEXTURE_2D, GL.GL_TEXTURE_3D)[dimension-1]

        if not data is None:
            size = tuple(data.shape[dimension-1::-1])
            format, iformat, tdtype, ncomp = self.texture_format(data)
            self.initialize_texture(size, format, iformat, tdtype, ncomp, data)

    def initialize_rgba(self, size):

        format = GL.GL_BGRA
        iformat = GL.GL_RGBA8
        tdtype = GL.GL_UNSIGNED_BYTE
        ncomp = 4
        self.initialize_texture(size, format, iformat, tdtype, ncomp)

    def initialize_8_bit(self, size):

        format = GL.GL_RED
        # TODO: PyOpenGL-20130502 does not have GL_R8.
        GL_R8 = 0x8229
        iformat = GL_R8
        tdtype = GL.GL_UNSIGNED_BYTE
        ncomp = 1
        self.initialize_texture(size, format, iformat, tdtype, ncomp)

    def initialize_depth(self, size, depth_compare_mode = True):

        format = GL.GL_DEPTH_COMPONENT
        iformat = GL.GL_DEPTH_COMPONENT24
        tdtype = GL.GL_FLOAT
        ncomp = 1
        self.initialize_texture(size, format, iformat, tdtype, ncomp,
                                depth_compare_mode = depth_compare_mode)

    def initialize_texture(self, size, format, iformat, tdtype, ncomp, data = None, depth_compare_mode = False):

        from OpenGL import GL
        self.id = t = GL.glGenTextures(1)
        self.size = size
        gl_target = self.gl_target
        GL.glBindTexture(gl_target, t)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        if data is None:
            import ctypes
            data = ctypes.c_void_p(0)
        dim = self.dimension
        if dim == 1:
            GL.glTexImage1D(gl_target, 0, iformat, size[0], 0, format, tdtype, data)
        elif dim == 2:
            GL.glTexImage2D(gl_target, 0, iformat, size[0], size[1], 0, format, tdtype, data)
        elif dim == 3:
            GL.glTexImage3D(gl_target, 0, iformat, size[0], size[1], size[2], 0, format, tdtype, data)

        GL.glTexParameterfv(gl_target, GL.GL_TEXTURE_BORDER_COLOR, (0,0,0,0))
        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_BORDER)
        if dim >= 2:
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_BORDER)
        if dim >= 3:
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_BORDER)

#        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
#        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
#        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_EDGE)

#        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
#        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(gl_target, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)

        if depth_compare_mode:
            # For GLSL sampler2dShadow objects to compare depth to r texture coord.
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_COMPARE_MODE, GL.GL_COMPARE_REF_TO_TEXTURE)
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_COMPARE_FUNC, GL.GL_LEQUAL);

        if ncomp == 1 or ncomp == 2:
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_SWIZZLE_G, GL.GL_RED)
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_SWIZZLE_B, GL.GL_RED)
        if ncomp == 2:
            GL.glTexParameteri(gl_target, GL.GL_TEXTURE_SWIZZLE_A, GL.GL_GREEN)
        GL.glBindTexture(gl_target, 0)

    def __del__(self):
        'Delete the OpenGL texture.'
        GL.glDeleteTextures((self.id,))

    def bind_texture(self, tex_unit = None):
        'Bind the OpenGL texture.'
        if tex_unit is None:
            GL.glBindTexture(self.gl_target, self.id)
        else:
            GL.glActiveTexture(GL.GL_TEXTURE0 + tex_unit)
            GL.glBindTexture(self.gl_target, self.id)
            GL.glActiveTexture(GL.GL_TEXTURE0)

    def unbind_texture(self, tex_unit = None):
        'Unbind the OpenGL texture.'
        if tex_unit is None:
            GL.glBindTexture(self.gl_target, 0)
        else:
            GL.glActiveTexture(GL.GL_TEXTURE0 + tex_unit)
            GL.glBindTexture(self.gl_target, 0)
            GL.glActiveTexture(GL.GL_TEXTURE0)

    def reload_texture(self, data):
        '''
        Replace the texture values in texture with OpenGL id using numpy array data.
        The data is interpreted the same as for the Texture constructor data argument.
        '''

        dim = self.dimension
        size = data.shape[dim-1::-1]
        format, iformat, tdtype, ncomp = self.texture_format(data)
        from OpenGL import GL
        gl_target = self.gl_target
        GL.glBindTexture(gl_target, self.id)
        if dim == 1:
            GL.glTexSubImage2D(gl_target, 0, 0, 0, size[0], format, tdtype, data)
        elif dim == 2:
            GL.glTexSubImage2D(gl_target, 0, 0, 0, size[0], size[1], format, tdtype, data)
        elif dim == 3:
            GL.glTexSubImage3D(gl_target, 0, 0, 0, size[0], size[1], size[2], format, tdtype, data)
        GL.glBindTexture(gl_target, 0)

    def texture_format(self, data):
        '''
        Return the OpenGL texture format, internal format, and texture value type
        that will be used by the glTexImageNd() function when creating a texture from
        a numpy array of colors.
        '''
        dim = self.dimension
        from OpenGL import GL
        if dim == 2 and len(data.shape) == dim and data.itemsize == 4:
            format = GL.GL_RGBA
            iformat = GL.GL_RGBA8
            tdtype = GL.GL_UNSIGNED_BYTE
            ncomp = 4
            return format, iformat, tdtype, ncomp

        ncomp = data.shape[dim] if len(data.shape) > dim else 1
        # TODO: Report pyopengl bug, GL_RG missing
        GL.GL_RG = 0x8227
        # luminance texture formats are not in opengl 3.
        format = {1:GL.GL_RED, 2:GL.GL_RG,
                  3:GL.GL_RGB, 4:GL.GL_RGBA}[ncomp]
        iformat = {1:GL.GL_RED, 2:GL.GL_RG,
                   3:GL.GL_RGB8, 4:GL.GL_RGBA8}[ncomp]
        dtype = data.dtype
        from numpy import uint8, float32
        if dtype == uint8:
            tdtype = GL.GL_UNSIGNED_BYTE
        elif dtype == float32:
            tdtype = GL.GL_FLOAT
        else:
            raise TypeError('Texture value type %s not supported' % str(dtype))
        return format, iformat, tdtype, ncomp

class Texture_Window:
    '''Draw a texture on a full window rectangle.'''
    def __init__(self, render, shader_mode):

        # Must have vao bound before compiling shader.
        self.vao = vao = Bindings(None)
        vao.activate()

        p = render.shader({shader_mode:True})
        render.use_shader(p)
        vao.shader = p

        self.vertex_buf = vb = Buffer(VERTEX_BUFFER)
        from numpy import array, float32, int32
        vb.update_buffer_data(array(((-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)), float32))
        vao.bind_shader_variable(vb)
        self.tex_coord_buf = tcb = Buffer(TEXTURE_COORDS_2D_BUFFER)
        tcb.update_buffer_data(array(((0,0,0),(1,0,0),(1,1,0),(0,1,0)), float32))
        vao.bind_shader_variable(tcb)
        self.element_buf = eb = Buffer(ELEMENT_BUFFER)
        eb.update_buffer_data(array(((0,1,2),(0,2,3)), int32))
        vao.bind_shader_variable(eb)    # Binds element buffer for rendering

    def __del__(self):
        self.vao = None
        for b in (self.vertex_buf, self.tex_coord_buf, self.element_buf):
            b.delete_buffer()

    def draw(self, xshift = 0, yshift = 0):
        xs, ys = xshift, yshift
        tcb = self.tex_coord_buf
        from numpy import array, float32
        tcb.update_buffer_data(array(((xs,ys,0),(1+xs,ys,0),(1+xs,1+ys,0),(xs,1+ys,0)), float32))
        GL.glDepthMask(False)   # Don't overwrite depth buffer
        eb = self.element_buf
        eb.draw_elements(eb.triangles)
        GL.glDepthMask(True)
