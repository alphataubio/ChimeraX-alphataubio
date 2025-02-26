# === UCSF ChimeraX Copyright ===
# Copyright 2025 Regents of the University of California. All rights reserved.
# The ChimeraX application is provided pursuant to the ChimeraX license
# agreement, which covers academic and commercial uses. For more details, see
# <https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html>
#
# You can also
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

from chimerax.core.commands import register, CmdDesc, StringArg, FloatArg


def register_command(command_name, logger):
    if command_name == "scenes save":
        func = save_scene
        desc = save_scene_desc
    elif command_name == "scenes delete":
        func = delete_scene
        desc = delete_scene_desc
    elif command_name == "scenes edit":
        func = edit_scene
        desc = edit_scene_desc
    elif command_name == "scenes restore":
        func = restore_scene
        desc = restore_scene_desc
    elif command_name == "scenes list":
        func = list_scenes
        desc = list_scenes_desc
    else:
        raise ValueError("trying to register unknown command: %s" % command_name)
    register(command_name, desc, func)


def save_scene(session, scene_name):
    """Save the current scene as 'scene_name'."""
    session.scenes.save_scene(scene_name)


save_scene_desc = CmdDesc(
    required=[("scene_name", StringArg)],
    synopsis="Save the current scene as 'scene_name'."
)


def delete_scene(session, scene_name: str):
    """Delete the scene 'scene_name'."""
    session.scenes.delete_scene(scene_name)


delete_scene_desc = CmdDesc(
    required=[("scene_name", StringArg)],
    synopsis="Delete scene 'scene_name'."
)


def edit_scene(session, scene_name: str):
    """Edit the scene 'scene_name'."""
    session.scenes.edit_scene(scene_name)


edit_scene_desc = CmdDesc(
    required=[("scene_name", StringArg)],
    synopsis="Edit scene 'scene_name'."
)


def restore_scene(session, scene_name):
    """Restore the scene named 'scene_name'."""
    session.scenes.restore_scene(scene_name)


restore_scene_desc = CmdDesc(
    required=[("scene_name", StringArg)],
    synopsis="Restore the scene named 'scene_name'."
)


def list_scenes(session):
    """List all saved scenes."""
    for scene_name in session.scenes.get_scene_names():
        print(scene_name)


list_scenes_desc = CmdDesc(
    synopsis="List all saved scenes."
)