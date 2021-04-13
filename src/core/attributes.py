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

"""
Attributes: support for custom attibute saving in sessions
==========================================================

This module allows classes that can have arbitrary new attributes defined to preserve those attributes
in session files.  Such attributes can have a default value, and if provided then instances of the class
where the attribute value has not been explicitly set will return the default value instead of raising
AttributeError.
"""

# Custom attrs:
# At class definition, need to call a method of the manager to add in the registration machinery.

class RegistrationConflict(ValueError):
    pass

def type_attrs(t):
    """Return known attribute names for the class/type 't'"""
    from types import GetSetDescriptorType
    attrs = [name for name in dir(t)
        if name[0] != '_' and type(getattr(t, name)) in [property, GetSetDescriptorType]]
    attrs.extend(t._attr_registration.reg_attr_info.keys())
    attrs.sort()
    return attrs

from .state import State, StateManager

# something that can't be a default value, yet can be saved in sessions...
class _NoDefault(State):

    def take_snapshot(self, session, flags):
        return {'version':1}

    @staticmethod
    def restore_snapshot(session, data):
        return NO_DEFAULT
NO_DEFAULT = _NoDefault()

# methods that the manager will insert into the managed classes (inheritance won't work)
def _getattr_(self, attr_name, look_in_class=None):
    if look_in_class is None:
        look_in_class = self.__class__
    try:
        return look_in_class._attr_registration.get_attr(attr_name)
    except AttributeError:
        for base in look_in_class.__bases__:
            if hasattr(base, "__getattr__"):
                return base.__getattr__(self, attr_name, look_in_class=base)
        else:
            raise

@classmethod
def register_attr(cls, session, attr_name, registerer, *, default_value=NO_DEFAULT, attr_type=None,
        can_return_none=False):
    cls._attr_registration.register(session, attr_name, registerer, default_value, (attr_type, can_return_none))

# used within the class to hold the registration info
class AttrRegistration:
    def __init__(self, class_):
        self.reg_attr_info = {}
        self.class_ = class_
        self._session_attrs = {}
        self._ses_attr_counts = {}

    def get_attr(self, attr_name):
        try:
            registrant, default_value, type_info = self.reg_attr_info[attr_name]
        except KeyError:
            if attr_name in dir(self.class_):
                raise AttributeError("Execution of '%s' object's '%s' property raised AttributeError" % (self.class_.__name__, attr_name)) from None
            else:
                raise AttributeError("'%s' object has no attribute '%s'" % (self.class_.__name__, attr_name)) from None
        if default_value == NO_DEFAULT:
            raise AttributeError("'%s' object has no attribute '%s'" % (self.class_.__name__, attr_name))
        return default_value

    def register(self, session, attr_name, registrant, default_value, type_info):
        if attr_name in self.reg_attr_info:
            prev_registrant, prev_default, prev_type_info = self.reg_attr_info[attr_name]
            if prev_default == default_value and prev_type_info == type_info:
                return
            raise RegistrationConflict("Registration of attr '%s' with %s by %s conflicts with previous"
                " registration by %s" % (attr_name, self.class_.__name__, registrant, prev_registrant))
        self.reg_attr_info[attr_name] = (registrant, default_value, type_info)
        session_attrs = self._session_attrs.setdefault(session, set())
        if attr_name not in session_attrs:
            session_attrs.add(attr_name)

    # session functions; called from manager, not directly from session-saving mechanism,
    # so API varies from that for State class
    def reset_state(self, session):
        # don't nuke existing attributes, since tools that registered attributes with default
        # values may still depend on those attributes working
        return

    def take_snapshot(self, session, flags):
        ses_reg_attr_info = {}
        attr_names = self._session_attrs.get(session, [])
        for attr_name in attr_names:
            ses_reg_attr_info[attr_name] = self.reg_attr_info[attr_name]
        data = {'reg_attr_info': ses_reg_attr_info, 'version': 2}
        return data

    def restore_session_data(self, session, data):
        if data.get('version', 1) == 1:
            reg_attr_info = {}
            for attr_name, reg_info in data['reg_attr_info'].items():
                registrant, default_value, attr_type = reg_info
                reg_attr_info[attr_name] = (registrant, default_value, (attr_type, False))
        else:
            reg_attr_info = data['reg_attr_info']
        for attr_name, reg_info in reg_attr_info.items():
            self.register(session, attr_name, *reg_info)

@property
def has_custom_attrs(self):
    for attr_name in self.__class__._attr_registration.reg_attr_info.keys():
        if hasattr(self, attr_name):
            return True
    return False

@property
def custom_attrs(self):
    custom_attrs = []
    for attr_name, attr_info in self.__class__._attr_registration.reg_attr_info.items():
        if hasattr(self, attr_name):
            registrant, default_value, attr_type = attr_info
            val = getattr(self, attr_name)
            if default_value == NO_DEFAULT or val != default_value:
                custom_attrs.append((attr_name, val))
    return custom_attrs

def set_custom_attrs(self, ses_data):
    for attr_name, val in ses_data.get('custom attrs', []):
        setattr(self, attr_name, val)

# since attribute registration occurs at class-definition time, the manager may not
# exist yet, so provide for that
_mgr = None
_pending_classes = {}
def register_class(reg_class, instances_func, builtin_attr_info={}):
    """'reg_class' is the class that wants to be able to register custom attributes.
    'instances_func' is a function (taking 'session' as its only argument) that returns
        all existing Python instances of the class.  If an instance exists only in the
        C++ layer, the function should not create a Python instance.
    'builtin_attr_info' should only be provided by classes that want to automatically
        show up in some tools, such as Render By Attribute, and provides information about
        the class's builtin (i.e non-custom) attributes.  The dictionary value should
        be a mapping from attribute name to a list of possible types that attribute can
        return (including None).
    """
    if hasattr(reg_class, '_attr_registration'):
        return
    reg_class._attr_registration = AttrRegistration(reg_class)
    reg_class.__getattr__ = _getattr_
    reg_class.register_attr = register_attr
    reg_class.has_custom_attrs = has_custom_attrs
    reg_class.custom_attrs = custom_attrs
    reg_class.set_custom_attrs = set_custom_attrs
    if _mgr:
        _mgr.class_info[reg_class] = (instances_func, builtin_attr_info)
    else:
        _pending_classes[reg_class] = (instances_func, builtin_attr_info)

MANAGER_NAME = "attribute registration"
class RegAttrManager(StateManager):

    def __init__(self, session):
        self.class_info = {}
        self.class_info.update(_pending_classes)
        _pending_classes.clear()
        global _mgr
        _mgr = self
        self.init_state_manager(session, MANAGER_NAME)

    def attributes_returning(self, class_obj, return_types, *, none_okay=False):
        """Return list of attribute names for class 'class_obj' whose return types
           are the same or a subset of 'return_types'.  If 'none_okay' is True,
           then it is okay for the attribute to also possibly be None.
        """
        try:
            instances_func, builtin_attr_info = self.class_info[class_obj]
        except KeyError:
            raise ValueError("Class '%s' has not registered attribute information" % class_obj.__name__)

        # builtin properties
        matching_attr_names = []
        for attr_name, types in builtin_attr_info.items():
            matched_one = False
            for t in types:
                if t in return_types:
                    matched_one = True
                else:
                    if not (t is None and none_okay):
                        break
            else:
                if matched_one:
                    matching_attr_names.append(attr_name)

        # registered attributes
        for attr_name, attr_info in class_obj._attr_registration.reg_attr_info.items():
            registrant, default_value, type_info = attr_info
            attr_type, can_return_none = type_info
            if attr_type not in return_types:
                # this will also catch attr_type of None
                continue
            if can_return_none and not none_okay:
                continue
            matching_attr_names.append(attr_name)
        return matching_attr_names

    def has_attribute(self, class_obj, attr_name):
        """Does the class have a builtin or registered attribute with the given name?"""
        try:
            instances_func, builtin_attr_info = self.class_info[class_obj]
        except KeyError:
            raise ValueError("Class '%s' has not registered attribute information" % class_obj.__name__)

        return attr_name in builtin_attr_info or attr_name in class_obj._attr_registration.reg_attr_info

    @property
    def registered_classes(self):
        return self.class_info.keys()

    # session functions; there is one manager per session, and is only in charge of
    # remembering registrations from its session (instances save their own attrs)
    def reset_state(self, session):
        for reg_class in self.registered_classes:
            reg_class._attr_registration.reset_state(session)

    def take_snapshot(self, session, flags):
        # force save of registration instance session info

        # gather bundle information for the classes
        bundles = {}
        for rc in self.registered_classes:
            bundle = session.toolshed.find_bundle_for_class(rc)
            if not bundle:
                session.logger.warning("Attribute-registration manager could not find bundle for class '%s'"
                    % rc.__name__)
                continue
            bundles[rc] = bundle
        return {
            'version': 2,
            'registrations': {(bundles[rc].name, rc.__name__):
                rc._attr_registration.take_snapshot(session, flags) for rc in bundles},
            # force all Python instances of registered classes to be restored (instance itself
            # is responsible for actually restoring the attributes), so that they get saved/restored in
            # sessions even if there are no Python-layer references to them
            'instances': [[inst for  inst in inst_func(session)
                if inst.has_custom_attrs and getattr(inst, 'session', None) == session]
                for inst_func, builtin_info in self.class_info.values()]
        }

    @staticmethod
    def restore_snapshot(session, data):
        _mgr._restore_session_data(session, data)
        return _mgr

    def _restore_session_data(self, session, data):
        # Version 1 is from when the atomic bundle handled this, and lack the 'instances' key,
        # which was handled by a separate manager.  It is otherwise compatible with version 2.
        for class_info, registration in data['registrations'].items():
            if data['version'] == 1:
                bundle_name = "ChimeraX-Atomic"
                class_name = class_info
            else:
                bundle_name, class_name = class_info
            bundle = session.toolshed.find_bundle(bundle_name, session.logger)
            if not bundle:
                session.logger.warning("Could not find bundle '%s' needed to restore attribute-registration"
                    " information for class '%s'" % (bundle_name, class_name))
                continue
            try:
                rc = bundle.get_class(class_name, session.logger)
                if not rc:
                    raise ValueError("class not found")
            except:
                session.logger.warning("Could not find class '%s' in bundle '%s', needed to restore"
                    " attribute-registration information" % (class_name, bundle_name))
                continue
            rc._attr_registration.restore_session_data(session, registration)
        # don't have to do anything with 'instances'; having them restored was the whole point