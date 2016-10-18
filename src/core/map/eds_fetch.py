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
#
def fetch_eds_map(session, id, type = '2fofc', ignore_cache=False, **kw):
  '''
  Fetch crystallographic density maps from the Upsalla Electron Density Server.
  
   2fofc:    http://eds.bmc.uu.se/eds/sfd/1cbs/1cbs.omap
    fofc:    http://eds.bmc.uu.se/eds/sfd/1cbs/1cbs_diff.omap
    Info:    http://eds.bmc.uu.se/cgi-bin/eds/uusfs?pdbCode=1cbs
  Holdings:  http://eds.bmc.uu.se/eds/eds_holdings.txt
  '''

  site = 'eds.bmc.uu.se'
  url_pattern = 'http://%s/eds/dfs/%s/%s/%s'

  # Fetch map.
  log = session.logger
  log.status('Fetching %s from web site %s...' % (id,site))

  if type == 'fofc':
    map_name = id + '_diff.omap'
  elif type == '2fofc':
    map_name = id + '.omap'
  map_url = url_pattern % (site, id[1:3], id, map_name)

  from ..fetch import fetch_file
  filename = fetch_file(session, map_url, 'map %s' % id, map_name, 'EDS',
                        ignore_cache=ignore_cache)

  from .. import io
  models, status = io.open_data(session, filename, format = 'dsn6', name = id, **kw)
  return models, status

# -----------------------------------------------------------------------------
# Register to fetch EMDB maps with open command.
#
def register_eds_fetch():
    from .. import fetch
    fetch.register_fetch('eds', fetch_eds_map, 'dsn6', prefixes = ['eds'])
#    reg('EDS', fetch_eds_map, '1A0M', 'eds.bmc.uu.se/eds',
#        'http://eds.bmc.uu.se/cgi-bin/eds/uusfs?pdbCode=%s', session)
