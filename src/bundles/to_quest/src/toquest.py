# vim: set expandtab ts=4 sw=4:

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

# -----------------------------------------------------------------------------
# User interface for sending scenes to Quest VR headset for viewing with
# LookSee app
#
from chimerax.core.tools import ToolInstance

# ------------------------------------------------------------------------------
#
class ToQuest(ToolInstance):
    SESSION_ENDURING = True

    # help = 'help:user/tools/markerplacement.html'

    def __init__(self, session, tool_name):
        ToolInstance.__init__(self, session, tool_name)

        self._settings = _ToQuestSettings(session, "send_to_quest")

        from chimerax.ui import MainToolWindow
        tw = MainToolWindow(self, close_destroys=True)
        self.tool_window = tw
        parent = tw.ui_area
        
        from chimerax.ui.widgets import vertical_layout
        layout = vertical_layout(parent, margins = (5,0,0,0))
        
        # Triangles settings
        tf = self._create_triangles_gui(parent)
        layout.addWidget(tf)
        
        # Send to Quest, Options, and Help buttons
        bf = self._create_action_buttons(parent)
        layout.addWidget(bf)

        # Options panel
        options = self._create_options_gui(parent)
        layout.addWidget(options)

        layout.addStretch(1)    # Extra space at end
        
        tw.manage(placement="side")

        vt = self.session.main_view.triggers
        h = vt.add_handler('shape changed', self._geometry_changed)
        self._geometry_change_handler = h

    # ---------------------------------------------------------------------------
    #
    def delete(self):
        vt = self.session.main_view.triggers
        vt.remove_handler(self._geometry_change_handler)
        self._geometry_change_handler = None
        ToolInstance.delete(self)
        
    # ---------------------------------------------------------------------------
    #
    @classmethod
    def get_singleton(self, session, create=True):
        from chimerax.core import tools
        return tools.get_singleton(session, ToQuest, 'Send to Quest', create=create)
        
    # ---------------------------------------------------------------------------
    #
    def _create_triangles_gui(self, parent):
        from Qt.QtWidgets import QFrame
        f = QFrame(parent)
        from chimerax.ui.widgets import vertical_layout, EntriesRow
        layout = vertical_layout(f, margins = (5,0,0,0))
        max_tri = 700000 if self._settings.limit_for_room_view else 900000
        tc = EntriesRow(f, '#', 'scene triangles', '    ', True, 'Maximum', max_tri)
        self._triangle_count = tcount = tc.labels[0]
        self._use_max_triangles, self._max_triangles = um, mt = tc.values
        mt.pixel_width = 70

        abt = EntriesRow(f, 'Simplify atoms', 100,
                         ('4', lambda: self._set_atom_triangles(4)),
                         ('10', lambda: self._set_atom_triangles(10)),
                         ('30', lambda: self._set_atom_triangles(30)),
                         ('100', lambda: self._set_atom_triangles(100)),
                         ('300', lambda: self._set_atom_triangles(300)),
                         '   bonds', 60,
                         ('12', lambda: self._set_bond_triangles(12)),
                         ('20', lambda: self._set_bond_triangles(20)),
                         ('60', lambda: self._set_bond_triangles(60)))
        self._atom_triangles, self._bond_triangles = at,bt = abt.values
        at.pixel_width, bt.pixel_width = 40,20
        at.return_pressed.connect(lambda *unused: self._set_atom_triangles())
        bt.return_pressed.connect(lambda *unused: self._set_bond_triangles())
        abt.frame.setStyleSheet('QPushButton {padding: 5px}') # Smaller buttons
        
        rt = EntriesRow(f, '          ribbon sides', 12,
                        ('4', lambda: self._set_ribbon_sides(4)),
                        ('8', lambda: self._set_ribbon_sides(8)),
                        ('12', lambda: self._set_ribbon_sides(12)),
                        '   divisions', 20,
                        ('2', lambda: self._set_ribbon_divisions(2)),
                        ('5', lambda: self._set_ribbon_divisions(5)),
                        ('10', lambda: self._set_ribbon_divisions(10)),
                        ('20', lambda: self._set_ribbon_divisions(20)))
        self._ribbon_sides, self._ribbon_divisions = rs,rd = rt.values
        rs.pixel_width = rd.pixel_width = 20
        rt.frame.setStyleSheet('QPushButton {padding: 5px}') # Smaller buttons
        rs.return_pressed.connect(lambda *unused: self._set_ribbon_sides())
        rd.return_pressed.connect(lambda *unused: self._set_ribbon_divisions())
        
        self._report_triangle_count()
        
        return f
    
    # ---------------------------------------------------------------------------
    #
    def _set_atom_triangles(self, tri = None):
        if tri is None:
            tri = self._atom_triangles.value
        from chimerax.core.commands import run
        run(self.session, f'graphics quality atomTriangles {tri}')
        self._atom_triangles.value = tri
        
    # ---------------------------------------------------------------------------
    #
    def _set_bond_triangles(self, tri = None):
        if tri is None:
            tri = self._bond_triangles.value
        from chimerax.core.commands import run
        run(self.session, f'graphics quality bondTriangles {tri}')
        self._bond_triangles.value = tri
        
    # ---------------------------------------------------------------------------
    #
    def _set_ribbon_sides(self, sides = None):
        if sides is None:
            sides = self._ribbon_sides.value
        from chimerax.core.commands import run
        run(self.session, f'graphics quality ribbonSides {sides}')
        self._ribbon_sides.value = sides
        
    # ---------------------------------------------------------------------------
    #
    def _set_ribbon_divisions(self, div = None):
        if div is None:
            div = self._ribbon_divisions.value
        from chimerax.core.commands import run
        run(self.session, f'graphics quality ribbonDivisions {div}')
        self._ribbon_divisions.value = div
    
    # ---------------------------------------------------------------------------
    #
    def _report_quality(self):
        from chimerax.atomic import Structure
        structs = [m for m in self.session.models.list(type = Structure) if m.visible]
        atri = max([len(s._atoms_drawing.triangles) for s in structs if s._atoms_drawing],
                   default = None)
        btri = max([len(s._bonds_drawing.triangles) for s in structs if s._bonds_drawing],
                   default = None)

        STYLE_ROUND = 1
        rside = max([s.ribbon_xs_mgr.params[STYLE_ROUND]['sides'] for s in structs],
                    default = None)
        from chimerax.atomic import structure_graphics_updater
        gu = structure_graphics_updater(self.session)
        lod = gu.level_of_detail
        rdiv = max([lod.ribbon_divisions(s.num_residues) for s in structs],
                   default = None)
        if atri:
            self._atom_triangles.value = atri
        if btri:
            self._bond_triangles.value = 2*btri		# btri is for half-bond
        if rside:
            self._ribbon_sides.value = rside
        if rdiv:
            self._ribbon_divisions.value = rdiv
    
    # ---------------------------------------------------------------------------
    #
    def _report_triangle_count(self):
        tcount = self._scene_triangles()
        tc = self._triangle_count
        tc.setText('%d' % tcount)
        tlimit = self._triangle_limit
        color = '' if tlimit is None or tcount <= tlimit else 'QLabel { color : red; }'
        tc.setStyleSheet(color);
        self._report_quality()
        return tcount
    
    # ---------------------------------------------------------------------------
    #
    def _scene_triangles(self):
        self.session.update_loop.update_graphics_now()  # Need to draw to get current ribbon
        models = [m for m in self.session.models.list() if len(m.id) == 1 and m.display]
        lines = []
        from chimerax.std_commands.graphics import _drawing_triangles
        tri = _drawing_triangles(models, lines)
        return tri

    # ---------------------------------------------------------------------------
    #
    def _geometry_changed(self, *unused):
        self._report_triangle_count()

    # ---------------------------------------------------------------------------
    #
    @property
    def _triangle_limit(self):
        if self._use_max_triangles.enabled:
            try:
                return self._max_triangles.value
            except:
                return None
        return None

    # ---------------------------------------------------------------------------
    #
    def _too_many_triangles(self, warn = True):
        tlimit = self._triangle_limit
        tcount = self._report_triangle_count()
        too_many = tlimit is not None and tcount > tlimit
        if warn and too_many:
            self.session.logger.error('Too many triangles to send to Quest, '
                                      f'{tcount} > {tlimit}')
        return too_many
    
    # ---------------------------------------------------------------------------
    #
    def _create_action_buttons(self, parent):
        from chimerax.ui.widgets import EntriesRow
        r = EntriesRow(parent,
                       ('Send to Quest', self._send_to_quest),
                       'filename', '',
                       ('Options', self._show_or_hide_options),
                       ('Help', self._help))
        self._filename = fn = r.values[0]
        fn.value = '    scene'
        return r.frame
    
    # ---------------------------------------------------------------------------
    #
    @property
    def _scene_filename(self):
        fn = self._filename.value.strip()
        if not fn:
            fn = 'scene.glb'
        elif not fn.endswith('.glb'):
            fn += '.glb'
        return fn
    
    # ---------------------------------------------------------------------------
    #
    def _send_to_quest(self):
        if self._too_many_triangles():
            return

        if not self._have_adb_path():
            self._need_adb()
            return
        
        # Save current scene.
        from os.path import expanduser, sep
        path = expanduser(f'~/Desktop/{self._scene_filename}').replace(sep, '/')
        from chimerax.core.commands import run
        run(self.session, f'save {path}')

        # Transfer scene file to Quest using adb
        adb = self._adb_path.value
        scenes_dir = f'/sdcard/Android/data/com.UCSF.LookSee/files'

        # all output is on stderr, but Windows needs all standard I/O to
        # be redirected if one is, so stdout is a pipe too
        args = [adb, "push", path, scenes_dir]
        cmd = ' '.join(args)
        self.session.logger.info(f'Running command: {cmd}')
        from subprocess import Popen, PIPE, DEVNULL
        p = Popen(args, stdin=DEVNULL, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out, err = out.decode('utf-8'), err.decode('utf-8')
        exit_code = p.returncode

        if exit_code != 0:
            lines = []
            if 'no devices' in out or 'no devices' in err:
                lines.extend([
                    'No VR headset was found by adb. '
                    'Try connecting a USB cable between the computer and headset '
                    'or using the adb connect <headset-ip-address> command from a shell.',
                     ''])
            if out:
                lines.extend(['stdout:', out])
            if err:
                lines.extend(['stderr:', err])
            self.session.logger.error('\n'.join(lines))
        else:
            self._settings.adb_executable_path = adb
            
    # ---------------------------------------------------------------------------
    #
    def _show_or_hide_options(self):
        self._options_panel.toggle_panel_display()
        
    # ---------------------------------------------------------------------------
    #
    def _create_options_gui(self, parent):
        from chimerax.ui.widgets import CollapsiblePanel
        self._options_panel = p = CollapsiblePanel(parent, title = None)
        f = p.content_area

        from chimerax.ui.widgets import EntriesRow

        # Path to adb executable
        ac = EntriesRow(f, 'adb executable', '', ('Browse', self._choose_adb_path))
        self._adb_path = adb = ac.values[0]
        adb.pixel_width = 350
        adb.value = self._settings.adb_executable_path

        # Whether to limit triangles for pass-through or black background
        # Black background can handle about 25% more triangles.
        rv = EntriesRow(f, False, 'Reduce maximum triangles for room view')
        self._room_view_max = rvm = rv.values[0]
        def set_max_tri(room_view):
            self._max_triangles.value = 700000 if room_view else 900000
            self._settings.limit_for_room_view = room_view
        rvm.changed.connect(set_max_tri)
        rvm.value = self._settings.limit_for_room_view

        return p
        
    # ---------------------------------------------------------------------------
    #
    def _choose_adb_path(self):
        parent = self.tool_window.ui_area
        from Qt.QtWidgets import QFileDialog
        path, ftype  = QFileDialog.getOpenFileName(parent, caption = f'adb executable')
        if path:
            self._adb_path.value = path
            self._settings.adb_executable_path = path

    # ---------------------------------------------------------------------------
    #
    def _have_adb_path(self):
        adb = self._adb_path.value.strip()
        from os import path
        return adb and path.exists(adb)

    # ---------------------------------------------------------------------------
    #
    def _need_adb(self):
        msg = (
            '<h2>ADB program needed to send files to Quest</h2>'
            '<p>'
            'To transfer files to the Quest you need the "adb" program installed '
            'on your computer and you need to give the path to that program by '
            'clicking the Options button in the Send to Quest panel '
            'and filling in the adb executable entry field. '
            'The adb command-line program is part of Android Studio SDK Platform-Tools '
            'obtained here: '
            '</p>'
            '<a href="https://developer.android.com/studio/releases/platform-tools">https://developer.android.com/studio/releases/platform-tools</a>'
            )
        self.session.logger.error(msg, is_html = True)
    
    # ---------------------------------------------------------------------------
    #
    def _help(self):
      from os.path import dirname, join
      help_url = 'file://' + join(dirname(__file__), 'sendtoquest.html')
      from chimerax.help_viewer import show_url
      show_url(self.session, help_url)

# -----------------------------------------------------------------------------
#
from chimerax.core.settings import Settings
class _ToQuestSettings(Settings):
    AUTO_SAVE = {
        'adb_executable_path': '',
        'limit_for_room_view': False,
    }
        
def to_quest_panel(session, create = True):
  return ToQuest.get_singleton(session, create)
