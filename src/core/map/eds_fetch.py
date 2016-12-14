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
  Fetch crystallographic density maps from PDBe (formerly the Upsalla Electron Density Server).

  2fofc: http://www.ebi.ac.uk/pdbe/coordinates/files/1cbs.ccp4
   fofc: http://www.ebi.ac.uk/pdbe/coordinates/files/1cbs_diff.ccp4
  '''

  url_pattern = 'http://www.ebi.ac.uk/pdbe/coordinates/files/%s'
  
  # Fetch map.
  log = session.logger
  log.status('Fetching %s from PDBe...' % (id,))

  if type == 'fofc':
    map_name = id.lower() + '_diff.ccp4'
  elif type == '2fofc':
    map_name = id.lower() + '.ccp4'
  map_url = url_pattern % map_name

  from ..fetch import fetch_file
  filename = fetch_file(session, map_url, 'map %s' % id, map_name, 'EDS',
                        ignore_cache=ignore_cache)

  from .. import io
  models, status = io.open_data(session, filename, format = 'ccp4', name = id,
                                polar_values = (type == 'fofc'), **kw)
  for v in models:
    v.set_representation('mesh')
    v.show()
    
  return models, status

# -----------------------------------------------------------------------------
#
def fetch_edsdiff_map(session, id, ignore_cache=False, **kw):
  return fetch_eds_map(session, id, type = 'fofc', ignore_cache=False, **kw)

# -----------------------------------------------------------------------------
# Register to fetch EMDB maps with open command.
#
def register_eds_fetch():
    from .. import fetch
    fetch.register_fetch('eds', fetch_eds_map, 'ccp4', prefixes = ['eds'])
    fetch.register_fetch('edsdiff', fetch_edsdiff_map, 'ccp4', prefixes = ['edsdiff'])
