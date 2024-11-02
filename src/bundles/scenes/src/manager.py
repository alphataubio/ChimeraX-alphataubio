# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2022 Regents of the University of California. All rights reserved.
# The ChimeraX application is provided pursuant to the ChimeraX license
# agreement, which covers academic and commercial uses. For more details, see
# <https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html>
#
# This particular file is part of the ChimeraX library. You can also
# redistribute and/or modify it under the terms of the GNU Lesser General
# Public License version 2.1 as published by the Free Software Foundation.
# For more details, see
# <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>
#
# THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
# EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. ADDITIONAL LIABILITY
# LIMITATIONS ARE DESCRIBED IN THE GNU LESSER GENERAL PUBLIC LICENSE
# VERSION 2.1
#
# This notice must be embedded in or attached to all copies, including partial
# copies, of the software or any revisions or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.state import StateManager
from chimerax.graphics.gsession import ViewState
from .scene import Scene, SceneColors, SceneVisibility
from chimerax.core.triggerset import TriggerSet
from chimerax.core.models import REMOVE_MODELS
from chimerax.std_commands.view import _interpolate_views, _model_motion_centers
from .triggers import activate_trigger, ADDED, DELETED, EDITED


class SceneManager(StateManager):
    """
    Manager for scenes in ChimeraX.

    This class manages the creation, deletion, saving, restoring, and interpolation of scenes. It also handles the
    removal of models from scenes and provides methods to reset the state and take or restore snapshots.

    Attributes:
        version (int): The version of the SceneManager.
        ADDED (str): Trigger name for added scenes.
        DELETED (str): Trigger name for deleted scenes.
        scenes (dict): A dictionary mapping scene names to Scene objects.
        session: The current session.
    """

    version = 0

    def __init__(self, session):
        """
        Initialize the SceneManager.

        Args:
            session: The current session.
        """
        self.scenes: [Scene] = []
        self.session = session
        session.triggers.add_handler(REMOVE_MODELS, self._remove_models_cb)

    def scene_exists(self, scene_name: str) -> bool:
        """
        Check if a scene exists by name and return True if it does, False otherwise.
        """
        return scene_name in [scene.get_name() for scene in self.scenes]

    def delete_scene(self, scene_name):
        """
        Delete scene by name.
        """
        if self.scene_exists(scene_name):
            self.scenes = [scene for scene in self.scenes if scene.get_name() != scene_name]
            activate_trigger(DELETED, scene_name)
        else:
            self.session.logger.warning(f"Scene {scene_name} does not exist.")

    def edit_scene(self, scene_name):
        if self.scene_exists(scene_name):
            self.get_scene(scene_name).init_form_session()
            activate_trigger(EDITED, scene_name)
        else:
            self.session.logger.warning(f"Scene {scene_name} does not exist.")

    def clear(self):
        """
        Delete all scenes.
        """
        for scene in self.scenes:
            self.delete_scene(scene.get_name())

    def save_scene(self, scene_name):
        """
        Save the current state as a scene.
        """
        if self.scene_exists(scene_name):
            self.session.logger.warning(f"Scene {scene_name} already exists.")
            return
        self.scenes.append(Scene(self.session, scene_name))
        activate_trigger(ADDED, scene_name)
        return

    def restore_scene(self, scene_name):
        """
        Restore a scene by name.
        """
        if self.scene_exists(scene_name):
            self.get_scene(scene_name).restore_scene()
        return

    def interpolate_scenes(self, scene_name1, scene_name2, fraction):
        """
        Interpolate between two scenes.

        Args:
            scene_name1 (str): The name of the first scene.
            scene_name2 (str): The name of the second scene.
            fraction (float): The interpolation fraction (0.0 to 1.0).
        """
        if self.scene_exists(scene_name1) and self.scene_exists(scene_name2):
            scene1 = self.get_scene(scene_name1)
            scene2 = self.get_scene(scene_name2)

            # Check if scenes are compatible
            if not Scene.interpolatable(scene1, scene2):
                self.session.logger.warning(f"Cannot interpolate between scenes {scene_name1} and {scene_name2} because"
                                         f"they are incompatible.")
                return

            # Interpolate main view data
            ViewState.interpolate(
                self.session.view,
                scene1.main_view_data,
                scene2.main_view_data,
                fraction
            )

            # Use NamedViews to interpolate camera, clip planes, and model positions. See _InterpolateViews
            view1 = scene1.named_view
            view2 = scene2.named_view
            centers = _model_motion_centers(view1.positions, view2.positions)
            _interpolate_views(view1, view2, fraction, self.session.main_view, centers)

            # Interpolate scene color and visibility data
            SceneColors.interpolate(self.session, scene1.get_colors(), scene2.get_colors(), fraction)
            SceneVisibility.interpolate(self.session, scene1.get_visibility(), scene2.get_visibility(), fraction)
        return

    def _remove_models_cb(self, trig_name, models):
        """
        Callback for removing models from scenes.

        Args:
            trig_name (str): The name of the trigger.
            models: The models to remove.
        """
        for scene in self.scenes:
            scene.models_removed(models)

    # session methods
    def reset_state(self, session):
        """
        Reset the state of the SceneManager by removing all the scenes.
        """
        self.clear()

    @staticmethod
    def restore_snapshot(session, data):
        if data['version'] != SceneManager.version:
            raise ValueError("scenes restore_snapshot: unknown version in data: %d" % data['version'])
        mgr = session.scenes
        mgr._ses_restore(data)
        return mgr

    def take_snapshot(self, session, flags):
        # viewer_info is "session independent"
        return {
            'version': self.version,
            'scenes': [scene.take_snapshot(session, flags) for scene in self.scenes]
        }

    def _ses_restore(self, data):
        """
        Restore the SceneManager scenes attribute from session data.

        Args:
            data (dict): The session data.
        """
        self.clear()
        for scene_snapshot in data['scenes']:
            scene = Scene.restore_snapshot(self.session, scene_snapshot)
            self.scenes.append(scene)

    def get_scene(self, scene_name: str) -> Scene | None:
        """
        Get a scene by name. If the scene does not exist, return None.
        """
        for scene in self.scenes:
            if scene.get_name() == scene_name:
                return scene
        return None

    def get_scene_names(self) -> {str, Scene}:
        return [scene.get_name() for scene in self.scenes]
