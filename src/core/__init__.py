# vi: set expandtab shiftwidth=4 softtabstop=4:
"""
chimera.core: collection of base chimera functionality
======================================================

"""
__copyright__ = (
    "Copyright \u00A9 2015 by the Regents of the University of California."
    "  All Rights Reserved."
    "  Free for non-commercial use."
    "  See http://www.cgl.ucsf.edu/chimera/ for license details."
)
_class_cache = {}


def get_class(class_name):
    """Return chimera.core class for instance in a session

    Parameters
    ----------
    class_name : str
        Class name
    """
    try:
        return _class_cache[class_name]
    except KeyError:
        pass
    if class_name == 'AtomicStructure':
        from . import atomic
        cls = atomic.AtomicStructure
    elif class_name == 'Generic3DModel':
        from . import generic3d
        cls = generic3d.Generic3DModel
    elif class_name == 'MolecularSurface':
        from . import molsurf
        cls = molsurf.MolecularSurface
    elif class_name == 'STLModel':
        from . import stl
        cls = stl.STLModel
    elif class_name == 'Job':
        from . import tasks
        cls = tasks.Job
    elif class_name == '_Input':
        from .ui import nogui
        cls = nogui._Input
    else:
        return None
    _class_cache[class_name] = cls
    return cls

def profile(func):
    def wrapper(*args, **kw):
        import cProfile, pstats, sys
        prof = cProfile.Profile()
        v = prof.runcall(func, *args, **kw)
        print(func.__name__, file=sys.__stderr__)
        p = pstats.Stats(prof, stream=sys.__stderr__)
        p.strip_dirs().sort_stats("cumulative", "time").print_callers(40)
        return v
    return wrapper
