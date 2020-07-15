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

import ctypes
class HoloPlayCore:
    '''Access C library for LookingGlass render parameters.'''

    def __init__(self, library_path = None):
        self._c_funcs = self._open_c_library(library_path)

        self._device_params = {}		# Cached device parameters
        self._show_quilt = False
        self._string_length = 1000
        self._string = ctypes.create_string_buffer(self._string_length)
        self._string_p = ctypes.c_char_p(ctypes.addressof(self._string))

    def _open_c_library(self, library_path):
        if library_path is None:
            from sys import platform
            libname = {'darwin': 'libHoloPlayCore.dylib',
                       'win32': 'HoloPlayCore.dll',
                       'linux': 'libHoloPlayCore.so'}.get(platform)
            from os.path import join, dirname
            library_path = join(dirname(__file__), 'lib', libname)
        c_funcs = CFunctions(library_path)
        return c_funcs
        
    # license type argument to hpc_InitializeApp()
    hpc_LICENSE_NONCOMMERCIAL = 0
    hpc_LICENSE_COMMERCIAL = 1

    # Return error codes for hpc_InitializeApp()
    hpc_CLIERR_NOERROR = 0
    hpc_CLIERR_NOSERVICE = 1
    hpc_CLIERR_VERSIONERR = 2
    hpc_CLIERR_SERIALIZEERR = 3
    hpc_CLIERR_DESERIALIZEERR = 4
    hpc_CLIERR_MSGTOOBIG = 5
    hpc_CLIERR_SENDTIMEOUT = 6
    hpc_CLIERR_RECVTIMEOUT = 7
    hpc_CLIERR_PIPEERROR = 8
    hpc_CLIERR_APPNOTINITIALIZED = 9

    def hpc_InitializeApp(self, app_name, license_type):
        f = self._c_function('hpc_InitializeApp', args = (ctypes.c_char_p, ctypes.c_int), ret = ctypes.c_int)
        err_code = f(app_name.encode('utf-8'), license_type)
        return err_code

    def error_code_message(self, err_code):
        error_strings = {
            self.hpc_CLIERR_NOERROR: "HoloPlay no error",
            self.hpc_CLIERR_NOSERVICE: "HoloPlay Service not running",
            self.hpc_CLIERR_VERSIONERR: "Incompatible version of HoloPlay Service",
            self.hpc_CLIERR_SERIALIZEERR: "Client message could not be serialized",
            self.hpc_CLIERR_DESERIALIZEERR: "Client message could not be deserialized",
            self.hpc_CLIERR_MSGTOOBIG: "HoloPlay message too big",
            self.hpc_CLIERR_SENDTIMEOUT: "Interprocess pipe send timeout",
            self.hpc_CLIERR_RECVTIMEOUT: "Interprocess pipe receive timeout",
            self.hpc_CLIERR_PIPEERROR: "Interprocess pipe broken",
            self.hpc_CLIERR_APPNOTINITIALIZED: "HoloPlay app not initialized",
        }
        if err_code in error_strings:
            msg = error_strings.get(err_code)
        else:
            msg = 'HoloPlay unknown error %d' % err_code
        return msg
    
    def hpc_CloseApp(self):
        f = self._c_function('hpc_CloseApp')
        f()
        
    def hpc_GetHoloPlayCoreVersion(self):
        return self._string_property('hpc_GetHoloPlayCoreVersion')
        
    def hpc_GetHoloPlayServiceVersion(self):
        return self._string_property('hpc_GetHoloPlayCoreVersion')

    def hpc_GetNumDevices(self):
        f = self._c_function('hpc_GetNumDevices', ret = ctypes.c_int)
        return f()

    def hpc_GetDeviceHDMIName(self, device):
        # For 8.9" display on Mac got "LKG03cWgSCevy"
        return self._device_string('hpc_GetDeviceHDMIName', device)
    def hpc_GetDeviceType(self, device):
        # For 8.9" display on Mac got LKG03cWgSCevy "standard"
        return self._device_string('hpc_GetDeviceType', device)

    def hpc_GetDevicePropertyWinX(self, device):
        return self._device_int('hpc_GetDevicePropertyWinX', device)
    def hpc_GetDevicePropertyWinY(self, device):
        return self._device_int('hpc_GetDevicePropertyWinY', device)
    def hpc_GetDevicePropertyScreenW(self, device):
        return self._device_int('hpc_GetDevicePropertyScreenW', device)
    def hpc_GetDevicePropertyScreenH(self, device):
        return self._device_int('hpc_GetDevicePropertyScreenH', device)
    def hpc_GetDevicePropertyDisplayAspect(self, device):
        return self._device_float('hpc_GetDevicePropertyDisplayAspect', device)

    def hpc_GetDevicePropertyPitch(self, device):
        return self._device_float('hpc_GetDevicePropertyPitch', device)
    def hpc_GetDevicePropertyTilt(self, device):
        return self._device_float('hpc_GetDevicePropertyTilt', device)
    def hpc_GetDevicePropertyCenter(self, device):
        return self._device_float('hpc_GetDevicePropertyCenter', device)
    def hpc_GetDevicePropertySubp(self, device):
        return self._device_float('hpc_GetDevicePropertySubp', device)
    def hpc_GetDeviceViewCone(self, dev):
        return self.hpc_GetDevicePropertyFloat(dev, "/calibration/viewCone/value")
    def hpc_GetDevicePropertyFloat(self, device, prop_name):
        f = self._c_function('hpc_GetDevicePropertyFloat',
                            args = (ctypes.c_int, ctypes.c_char_p), ret = ctypes.c_float)
        return f(device, prop_name.encode('utf-8'))
    def hpc_GetDevicePropertyFringe(self, device):
        return self._device_float('hpc_GetDevicePropertyFringe', device)
    def hpc_GetDevicePropertyRi(self, device):
        return self._device_int('hpc_GetDevicePropertyRi', device)
    def hpc_GetDevicePropertyBi(self, device):
        return self._device_int('hpc_GetDevicePropertyBi', device)
    def hpc_GetDevicePropertyInvView(self, device):
        return self._device_int('hpc_GetDevicePropertyInvView', device)

    def info(self):

        ndev = self.hpc_GetNumDevices()
        lines = [
            "HoloPlay Core version %s" % self.hpc_GetHoloPlayCoreVersion(),
            "HoloPlay Service version %s" % self.hpc_GetHoloPlayServiceVersion(),
            "%d devices connected" % ndev,
        ]

        for dev in range(ndev):
            lines.extend([
                "Device information for display %d" % dev,
                "\tDevice name: %s" % self.hpc_GetDeviceHDMIName(dev),
                "\tDevice type: %s" % self.hpc_GetDeviceType(dev),

                "\nWindow parameters for display %s:" % dev,
                "\tPosition: (%d, %d)" % (self.hpc_GetDevicePropertyWinX(dev),
                                          self.hpc_GetDevicePropertyWinY(dev)),
                "\tSize: (%d, %d)" % (self.hpc_GetDevicePropertyScreenW(dev),
                                      self.hpc_GetDevicePropertyScreenH(dev)),
                "\tAspect ratio: %g" % self.hpc_GetDevicePropertyDisplayAspect(dev),

                "\nShader uniforms for display %d" % dev,
                "\tPitch: %g" % self.hpc_GetDevicePropertyPitch(dev),
                "\tTilt: %g" % self.hpc_GetDevicePropertyTilt(dev),
                "\tCenter: %g" % self.hpc_GetDevicePropertyCenter(dev),
                "\tSubpixel width: %g" % self.hpc_GetDevicePropertySubp(dev),
                "\tView cone: %g" % self.hpc_GetDeviceViewCone(dev),
                "\tFringe: %g" % self.hpc_GetDevicePropertyFringe(dev),
                "\tRI: %d" % self.hpc_GetDevicePropertyRi(dev),
                "\tBI: %d" % self.hpc_GetDevicePropertyBi(dev),
                "\tinvView: %d" % self.hpc_GetDevicePropertyInvView(dev),
            ])

        info = '\n'.join(lines)
        return info

    def set_shader_uniforms(self, shader, device, quilt_size):
        p = self.device_parameters(device)
        s = shader
        # Values in comments reported for 8.9" display by HoloPlayCore 1.1.1
        # (July 2020) on Mac, displaysize 2560 x 1600.
        s.set_float('pitch', p['pitch'])		# 354.686
        s.set_float('tilt', p['tilt'])			# -0.113296
        s.set_float('center', p['center'])		# -0.101902
        s.set_integer('invView', p['invView'])		# 1
        s.set_float('subp', p['subp'])			# 0.000130208
        s.set_float('displayAspect', p['aspect'])	# 1.6
        s.set_integer('ri', p['ri'])			# 0
        s.set_integer('bi', p['bi'])			# 2

        qs_width, qs_height, qs_rows, qs_columns = quilt_size	# 4096, 4096, 9, 5
        qs_total_views = qs_rows * qs_columns			# 45
        qs_view_width = qs_width // qs_columns			# 819		
        qs_view_height = qs_height // qs_rows			# 455
        s.set_vector3('tile', (qs_columns, qs_rows, qs_total_views))	# 5, 9, 45
        s.set_vector2('viewPortion', (qs_view_width * qs_columns / qs_width,
                                      qs_view_height * qs_rows / qs_height))	# 0.999755, 0.999755
        s.set_float('quiltAspect', p['aspect'])		# 1.6
        s.set_integer('overscan', 0)			# 0
        s.set_integer('quiltInvert', 0)			# 0
        s.set_integer('debug', (1 if self._show_quilt else 0))

    def device_parameters(self, device):
        d = device
        params = self._device_params.get(d)
        if params:
            return params
        params = {
            'pitch': self.hpc_GetDevicePropertyPitch(d),
            'tilt': self.hpc_GetDevicePropertyTilt(d),
            'center': self.hpc_GetDevicePropertyCenter(d),
            'invView': self.hpc_GetDevicePropertyInvView(d),
            'subp': self.hpc_GetDevicePropertySubp(d),
            'aspect': self.hpc_GetDevicePropertyDisplayAspect(d),
            'ri': self.hpc_GetDevicePropertyRi(d),
            'bi': self.hpc_GetDevicePropertyBi(d),
            
        }
        self._device_params[d] = params
        return params
        
    def _string_property(self, func_name):
        f = self._c_function(func_name, args = (ctypes.c_char_p, ctypes.c_int))
        f(self._string_p, self._string_length)
        bytes = self._string.value
        s = bytes.decode('ascii')
        return s

    def _device_string(self, func_name, device):
        f = self._c_function(func_name, args = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int))
        f(device, self._string_p, self._string_length)
        bytes = self._string.value
        s = bytes.decode('ascii')
        return s

    def _device_int(self, func_name, device):
        f = self._c_function(func_name, args = (ctypes.c_int,), ret = ctypes.c_int)
        return f(device)

    def _device_float(self, func_name, device):
        f = self._c_function(func_name, args = (ctypes.c_int,), ret = ctypes.c_float)
        return f(device)

    def _c_function(self, func_name, args = None, ret = None):
        return self._c_funcs.c_function(func_name, args=args, ret=ret)


# -----------------------------------------------------------------------------
#
class CFunctions:
    '''Access C functions from a shared library and create Python properties using these functions.'''

    def __init__(self, library_path):
        self._c_lib = ctypes.PyDLL(library_path)
        
    # -----------------------------------------------------------------------------
    #
    def c_function(self, func_name, args = None, ret = None):
        '''Look up a C function and set its argument types if they have not been set.'''
        f = getattr(self._c_lib, func_name)
        if args is not None and f.argtypes is None:
            f.argtypes = args
        if ret is not None:
            f.restype = ret
        return f