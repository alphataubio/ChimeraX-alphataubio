# vim: set expandtab shiftwidth=4 softtabstop=4:
from .models import Model

from . import io
CATEGORY = io.GENERIC3D


class Generic3DModel(Model):
    """Commom base class for generic 3D data"""

    def take_snapshot(self, session, flags):
        from .state import CORE_STATE_VERSION
        from .graphics.gsession import DrawingState
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'drawing state': DrawingState(self).take_snapshot(session, flags),
            'version': CORE_STATE_VERSION,
        }
        return data

    @classmethod
    def restore_snapshot(cls, session, data):
        m = cls('name', session)
        m.set_state_from_snapshot(session, data['model state'])
        from .graphics.gsession import DrawingState
        DrawingState(m).set_state_from_snapshot(session, data['drawing state'])
        return m

    def reset_state(self, session):
        pass
