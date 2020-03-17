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

class NoOpenerError(ValueError):
    pass

from chimerax.core.toolshed import ProviderManager
class OpenManager(ProviderManager):
    """Manager for open command"""

    def __init__(self, session):
        self.session = session
        self._openers = {}
        self._fetchers = {}
        self._compression_info = {}
        from chimerax.core.triggerset import TriggerSet
        self.triggers = TriggerSet()
        self.triggers.add_trigger("open command changed")

    def add_provider(self, bundle_info, name, *, type="open", want_path=False, check_path=True,
            batch=False, compression_suffixes=None, format_name=None, is_default=True, **kw):
        logger = self.session.logger

        bundle_name = _readable_bundle_name(bundle_info)
        is_default = bool_cvt(is_default, name, bundle_name, "is_default")
        want_path = bool_cvt(want_path, name, bundle_name, "want_path")
        check_path = bool_cvt(check_path, name, bundle_name, "check_path")
        if batch or not check_path:
            want_path = True
        type_description = "Open-command" if type == "open" else type.capitalize()
        if kw:
            logger.warning("%s provider '%s' supplied unknown keywords in provider description: %s"
                % (type_description, name, repr(kw)))
        if type == "open":
            try:
                data_format = self.session.data_formats[name]
            except KeyError:
                logger.warning("Open-command provider in bundle %s specified unknown data format '%s';"
                    " skipping" % (bundle_name, name))
                return
            if data_format in self._openers:
                logger.warning("Replacing opener for '%s' from %s bundle with that from %s bundle"
                    % (data_format.name, _readable_bundle_name(self._openers[data_format][0]), bundle_name))
            self._openers[data_format] = (bundle_info, name, want_path, check_path, batch)
        elif type == "fetch":
            if format_name is None:
                raise ValueError("Database fetch '%s' in bundle %s failed to specify file format name"
                    % (name, bundle_name))
            try:
                data_format = self.session.data_formats[format_name]
            except KeyError:
                raise ValueError("Database-fetch provider '%s' in bundle %s"
                    " specified unknown data format '%s'" % (name, bundle_name, format_name))
            if name in self._fetchers and format_name in self._fetchers[name]:
                logger.warning("Replacing fetcher for '%s' and format %s from %s bundle with that from %s"
                    " bundle" % (name, format_name,
                    _readable_bundle_name(self._fetchers[name][format_name][0]), bundle_name))
            self._fetchers.setdefault(name, {})[format_name] = (bundle_info, is_default)
            if is_default and len([fmt for fmt, info in self._fetchers[name].items() if info[1]]) > 1:
                logger.warning("Multiple default formats declared for database fetch '%s'" % name)
        elif type == "compression":
            suffixes = process_suffixes(compression_suffixes, "compression", name, logger, bundle_name)
            for suffix in suffixes:
                if suffix in self._compression_info:
                    prev_bi, prev_name = self._compression_info[suffix]
                    logger.warning("Duplicate decompression provider registered for file suffix '%s'."
                        " ('%s' from bundle %s and '%s' from bundle %s)" % (prev_name,
                        _readable_bundle_name(prev_bi), name, bundle_name))
                self._compression_info[suffix] = (bundle_info, name)
        else:
            logger.warning("Unknown provider type '%s' with name '%s' from bundle %s"
                % (type, name, bundle_name))

    def database_info(self, database_name):
        try:
            return self._fetchers[database_name]
        except KeyError:
            raise NoOpenerError("No such database '%s'" % database_name)

    @property
    def database_names(self):
        return list(self._fetchers.keys())

    def end_providers(self):
        self.triggers.activate_trigger("open command changed", self)
        if True:
            return
        in_common = []
        old_only = []
        from chimerax.core.io import _file_formats
        for old in _file_formats.keys():
            try:
                in_common.append(self.session.data_formats[old])
            except KeyError:
                print("old only:", old)
                old_only.append(old)
        for df in self.session.data_formats.formats:
            print("new:", df.name)
        print("%d data formats in common, %d new only, %d old only" %(len(in_common),
            len(self.session.data_formats) - len(in_common), len(old_only)))
        from random import choice
        print("Port format", choice(list(old_only_formats)+['http']))

    def fetch_args(self, database_name, *, format_name=None):
        try:
            db_formats = self._fetchers[database_name]
        except KeyError:
            raise NoOpenerError("No such database '%s'" % database_name)
        if format_name:
            try:
                bundle_info, is_default = db_formats[format_name]
            except KeyError:
                raise NoOpenerError("Format '%s' not supported for database '%s'.  Supported formats are: %s"
                    % (format_name, database_name, ", ".join(dbf for dbf in db_formats)))
        else:
            for format_name, info in db_formats.items():
                bundle_info, is_default = info
                if is_default:
                    break
            else:
                raise NoOpenerError("No default format for database '%s'.  Possible formats are: %s"
                    % (database_name, ", ".join(dbf for dbf in db_formats)))
        args = self.open_args(self.session.data_formats[format_name])
        args.update(bundle_info.run_provider(self.session, database_name, self, type="fetch").fetch_args)
        return args

    def open_data(self, path, **kw):
        from .cmd import provider_open
        return provider_open(self.session, [path], return_status=True, _add_to_file_history=False, **kw)

    def open_file(self, path, encoding=None):
        """Open possibly compressed file for reading.  If encoding is 'None', open as binary"""
        for suffix, comp_info in self._compression_info.items():
            if path.endswith(suffix):
                bundle_info, name = comp_info
                return bundle_info.decompress(session, name, path, encoding)
        return open(path, ('r' if encoding else 'rb'), encoding=encoding)

    def open_args(self, data_format):
        try:
            bundle_info, name, *args = self._openers[data_format]
        except KeyError:
            raise NoOpenerError("No opener registered for format '%s'" % data_format.name)
        return bundle_info.run_provider(self.session, name, self, type="open").open_args

    def open_info(self, data_format):
        try:
            bi, name, *args = self._openers[data_format]
            return (bi.run_provider(self.session, name, self, type="open"), name) + tuple(args)
        except KeyError:
            raise NoOpenerError("No opener registered for format '%s'" % data_format.name)

    def remove_compression_suffix(self, file_name):
        for suffix in self._compression_info.keys():
            if file_name.endswith(suffix):
                file_name = file_name[:-len(suffix)]
                break
        return file_name

def bool_cvt(val, name, bundle_name, var_name):
    if not isinstance(val, bool):
        try:
            val = eval(val.capitalize())
        except (ValueError, NameError):
            logger.warning("Database-fetch provider '%s' in bundle %s specified '%s' value"
                " (%s) that was neither 'true' nor 'false'" % (name, bundle_name, var_name, val))
    return val

def _readable_bundle_name(bundle_info):
    name = bundle_info.name
    if name.lower().startswith("chimerax"):
        return name[9:]
    return name

def process_suffixes(suffix_string, suffix_type, name, logger, bundle_name, none_okay=False):
    if not suffix_string:
        if not none_okay:
            logger.warning("%s provider for '%s' from bundle %s specified no %s suffixes"
                % (suffix_type.capitalize(), name, bundle_name, suffix_type))
        return []
    suffixes = []
    for suffix in suffix_string.split(','):
        if not suffix:
            logger.warning("Empty %s suffix found in suffix list of provider for '%s'"
                " from bundle %s" % (suffix_type, name, bundle_name))
            continue
        if suffix[0] != '.':
            logger.warning("%s suffixes must start with '.'.  '%s' found in suffix list of provider for"
                " '%s' in bundle %s does not." % (suffix_type, suffix, name, bundle_name))
        suffixes.append(suffix)
    return suffixes
