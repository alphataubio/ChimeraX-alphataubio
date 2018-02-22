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

"""OptionsPanel and CategorizedOptionsPanel are for interfaces that want to present a set of
options that aren't remembered in any way between different sessions.  The corresponding
Settings classes offer buttons to Save/Reset/Restore the option values to/from a Settings
instance (found in chimerax.core.settings).

The Categorized classes organize the presented options into categories, which the user can
switch between.
"""

from PyQt5.QtWidgets import QWidget, QFormLayout, QTabWidget, QBoxLayout, QGridLayout, \
    QPushButton, QCheckBox, QScrollArea

class OptionsPanel(QScrollArea):
    """OptionsPanel is a container for single-use (not savable) Options"""

    def __init__(self, parent=None, *, sorting=True, **kw):
        """sorting:
            False; options shown in order added
            True: options sorted alphabetically by name
            func: options sorted based on the provided key function
        """
        QScrollArea.__init__(self, parent, **kw)
        self._scrolled = QWidget()
        self.setWidget(self._scrolled)
        self._sorting = sorting
        self._options = []
        self._form = QFormLayout()
        self._form.setSizeConstraint(self._form.SetMinAndMaxSize)
        self._form.setVerticalSpacing(1)
        from PyQt5.QtCore import Qt
        self._form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self._scrolled.setLayout(self._form)

    def add_option(self, option):
        if self._sorting is None:
            insert_row = len(self._options)
        else:
            if self._sorting is True:
                test = lambda o1, o2: o1.name < o2.name
            else:
                test = lambda o1, o2: self._sorting(o1) < self._sorting(o2)
            for insert_row in range(len(self._options)):
                if test(option, self._options[insert_row]):
                    break
            else:
                insert_row = len(self._options)
        self._form.insertRow(insert_row, option.name, option.widget)
        self._options.insert(insert_row, option)
        if option.balloon:
            self._form.itemAt(insert_row,
                QFormLayout.LabelRole).widget().setToolTip(option.balloon)

    def options(self):
        return self._options

    def sizeHint(self):
        from PyQt5.QtCore import QSize
        form_size = self._form.minimumSize()
        return QSize(min(form_size.width(), 800), min(form_size.height(), 800))

class CategorizedOptionsPanel(QTabWidget):
    """CategorizedOptionsPanel is a container for single-use (not savable) Options sorted by category"""

    def __init__(self, parent=None, *, category_sorting=True, option_sorting=True, **kw):
        """sorting:
            False: categories/options shown in order added
            True: categories/options sorted alphabetically by name
            func: categories/options sorted based on the provided key function
        """
        QTabWidget.__init__(self, parent, **kw)
        self._category_sorting = category_sorting
        self._option_sorting = option_sorting
        self._category_to_panel = {}

    def add_option(self, category, option):
        """Add option to given category"""
        try:
            panel = self._category_to_panel[category]
        except KeyError:
            panel = OptionsPanel(sorting=self._option_sorting)
            self.add_tab(category, panel)
        panel.add_option(option)

    def add_tab(self, category, panel):
        """Only used if a tab needs to present a custom interface.

           The panel needs to offer a .options() method that returns its options."""
        self._category_to_panel[category] = panel
        if len(self._category_to_panel) == 1 or self._category_sorting == False:
            self.addTab(panel, category)
        else:
            cats = list(self._category_to_panel.keys()) + [category]
            if self._category_sorting is True:
                cats.sort()
            else:
                cats.sort(key=lambda cat: self._category_sorting(cat))
            self.insertTab(cats.index(category), panel, category)

    def categories(self):
        return self._category_to_panel.keys()

    def current_category(self):
        return self.tabText(self.currentIndex())

    def options(self, category):
        return self._category_to_panel[category].options()

class SettingsPanelBase(QWidget):
    def __init__(self, settings, owner_description, parent, option_sorting, multicategory,
            category_sorting=None, **kw):
        QWidget.__init__(self, parent, **kw)
        self.settings = settings
        self.multicategory = multicategory
        if multicategory:
            self.options_panel = CategorizedOptionsPanel(option_sorting=option_sorting,
                    category_sorting=category_sorting)
        else:
            self.options_panel = OptionsPanel(sorting=option_sorting)

        layout = QBoxLayout(QBoxLayout.TopToBottom)
        layout.setSpacing(5)
        layout.addWidget(self.options_panel, 1)
        layout.setContentsMargins(0,0,0,0)

        button_container = QWidget()
        bc_layout = QGridLayout()
        bc_layout.setContentsMargins(0, 0, 0, 0)
        bc_layout.setVerticalSpacing(5)
        if multicategory:
            self.current_check = QCheckBox("Buttons below apply to current section only")
            self.current_check.setToolTip("If checked, buttons only affect current section")
            self.current_check.setChecked(True)
            from .. import shrink_font
            shrink_font(self.current_check)
            from PyQt5.QtCore import Qt
            bc_layout.addWidget(self.current_check, 0, 0, 1, 3, Qt.AlignRight)
        save_button = QPushButton("Save")
        save_button.clicked.connect(self._save)
        save_button.setToolTip("Save as startup defaults")
        bc_layout.addWidget(save_button, 1, 0)
        reset_button = QPushButton("Reset")
        reset_button.clicked.connect(self._reset)
        reset_button.setToolTip("Reset to initial-installation defaults")
        bc_layout.addWidget(reset_button, 1, 1)
        restore_button = QPushButton("Restore")
        restore_button.clicked.connect(self._restore)
        restore_button.setToolTip("Restore from saved defaults")
        bc_layout.addWidget(restore_button, 1, 2)

        button_container.setLayout(bc_layout)
        layout.addWidget(button_container, 0)

        self.setLayout(layout)

    def _get_actionable_options(self):
        if self.multicategory:
            if self.current_check.isChecked():
                options = self.options_panel.options(self.options_panel.current_category())
            else:
                options = []
                for cat in self.options_panel.categories():
                    options.extend(self.options_panel.options(cat))
        else:
            options = self.options_panel.options()
        return options

    def _reset(self):
        from ...configfile import Value
        for opt in self._get_actionable_options():
            setting = opt.attr_name
            default_val = self.settings.PROPERTY_INFO[setting]
            if isinstance(default_val, Value):
                default_val = default_val.default
            # '==' on numpy objects doesn't return a boolean
            import numpy
            if not numpy.array_equal(opt.get(), default_val):
                opt.set(default_val)
                opt.make_callback()

    def _restore(self):
        for opt in self._get_actionable_options():
            setting = opt.attr_name
            restore_val = self.settings.saved_value(setting)
            # '==' on numpy objects doesn't return a boolean
            import numpy
            if not numpy.array_equal(opt.get(), default_val):
                opt.set(restore_val)
                opt.make_callback()

    def _save(self):
        save_settings = []
        for opt in self._get_actionable_options():
            setting = opt.attr_name
            # need to ensure "current value" is up to date before saving...
            setattr(self.settings, setting, opt.get())
            save_settings.append(setting)
        # We don't simply use settings.save() when all options are being saved
        # since there may be settings that aren't presented in the GUI
        self.settings.save(settings=save_settings)

class SettingsPanel(SettingsPanelBase):
    """SettingsPanel is a container for remember-able Options that work in conjunction with a
       Settings instance (found in chimerax.core.settings).

       The callback function for each option will have to call option.set_attribute(settings)
       in order to update the "current" value of the attribute in settings.
    """

    def __init__(self, settings, owner_description, parent=None, *, sorting=True, **kw):
        """'settings' is a Settings instance.  The remaining arguments are the same as
            for OptionsPanel
        """
        SettingsPanelBase.__init__(self, settings, owner_description, parent, sorting, multicategory=False, **kw)

    def add_option(self, option):
        self.options_panel.add_option(option)

class CategorizedSettingsPanel(SettingsPanelBase):
    """CategorizedSettingsPanel is a container for remember-able Options that work in conjunction
       with a Settings instance (found in chimerax.core.settings) and that are presented in
       categories.
    """

    def __init__(self, settings, owner_description, parent=None, *, category_sorting=True, option_sorting=True, **kw):
        """'settings' is a Settings instance.  The remaining arguments are the same as
            for CategorizedOptionsPanel
        """
        SettingsPanelBase.__init__(self, settings, owner_description, parent, option_sorting, multicategory=True,
            category_sorting=category_sorting, **kw)

    def add_option(self, category, option):
        self.options_panel.add_option(category, option)

    def add_tab(self, category, panel):
        """Same as CategorizedOptionsPanel.add_tab(...)"""
        self.options_panel.add_tab(category, panel)
