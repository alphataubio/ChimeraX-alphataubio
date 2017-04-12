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

def wait(session, frames=None):
    '''Wait before proceeding to the next command. Used in movie recording scripts.

    Parameters
    ----------
    frames : integer
       Wait until this many frames have been rendered before executing the next
       command in a command script.
    '''
    if frames is None:
        from . import motion
        while motion.motion_in_progress(session):
            draw_frame(session)
    else:
        for f in range(frames):
            draw_frame(session)

def draw_frame(session, limit_frame_rate = True):
    '''Draw a graphics frame.'''
    if limit_frame_rate:
        from time import time, sleep
        t0 = time()

    v = session.main_view
    v.redraw_needed = True  # Trigger frame rendered callbacks to cause image capture.
    ul = session.update_loop
    ul.draw_new_frame(session)

    if limit_frame_rate:
        dt = time() - t0
        gw = session.ui.main_window.graphics_window
        frame_time= gw.redraw_interval / 1000.0	# seconds
        if dt < frame_time:
            sleep(frame_time - dt)


def register_command(session):
    from . import cli
    desc = cli.CmdDesc(
        optional=[('frames', cli.PositiveIntArg)],
        synopsis='suspend command processing for a specified number of frames'
        ' or until finite motions have stopped ')
    cli.register('wait', desc, wait, logger=session.logger)
