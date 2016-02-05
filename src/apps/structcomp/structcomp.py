# vim:set et sw=4:
import os
import sys

session = session  # noqa

PDB_DIR = '/databases/mol/pdb'
MMCIF_DIR = '/databases/mol/mmCIF'


class FlushFile:
    def __init__(self, fd):
        self.fd = fd

    def write(self, x):
        ret = self.fd.write(x)
        self.fd.flush()
        return ret

    def writelines(self, lines):
        ret = self.writelines(lines)
        self.fd.flush()
        return ret

    def flush(self):
        return self.fd.flush()

    def close(self):
        return self.fd.close()

    def fileno(self):
        return self.fd.fileno()

# always flush stdout, even if output is to a file
sys.stdout = FlushFile(sys.stdout)


def bonds(atoms):
    a0s, a1s = atoms
    bonds = set()
    for rstr0, aname0, rstr1, aname1 in zip(
            a0s.residues.strs, a0s.names, a1s.residues.strs,
            a1s.names):
        id0 = (rstr0, aname0)
        id1 = (rstr1, aname1)
        if id0 > id1:
            id0, id1 = id1, id0
        bonds.add("%s@%s/%s@%s" % (id0 + id1))
    return bonds


def compare(session, pdb_id, pdb_path, mmcif_path):
    # return True if they differ
    session.logger.info('Comparing %s' % pdb_id)
    from chimerax.core import io, fetch
    try:
        if os.path.isabs(pdb_path):
            pdb_models = io.open_data(session, pdb_path)[0]
        else:
            pdb_models = fetch.fetch_from_database(
                session, 'pdb', pdb_path, format='pdb')[0]
    except Exception as e:
        session.logger.error("error: %s: unable to open pdb file: %s" % (pdb_id, e))
        return True
    try:
        if os.path.isabs(mmcif_path):
            mmcif_models = io.open_data(session, mmcif_path)[0]
        else:
            mmcif_models = fetch.fetch_from_database(
                session, 'pdb', mmcif_path, format='mmcif')[0]
    except Exception as e:
        session.logger.error("error: %s: unable to open mmcif file: %s" % (pdb_id, e))
        return

    if len(pdb_models) != len(mmcif_models):
        session.logger.error("error: %s: pdb version has %d model(s), mmcif version has %d" % (
            pdb_id, len(pdb_models), len(mmcif_models)))
        all_same = False
    else:
        all_same = True
        save_pdb_id = pdb_id
        i = 0
        for p, m in zip(pdb_models, mmcif_models):
            same = True

            if len(pdb_models) > 1:
                i += 1
                pdb_id = "%s#%d" % (save_pdb_id, i)
            # atoms
            pdb_atoms = ['%s@%s' % (r, a) for r, a in zip(
                p.atoms.residues.strs, p.atoms.names)]
            mmcif_atoms = ['%s@%s' % (r, a) for r, a in zip(
                m.atoms.residues.strs, m.atoms.names)]
            pdb_tmp = set(pdb_atoms)
            mmcif_tmp = set(mmcif_atoms)
            common = pdb_tmp & mmcif_tmp
            extra = pdb_tmp - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra pdb atom(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            extra = mmcif_tmp - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra mmcif atom(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False

            # bonds
            pdb_bonds = bonds(p.bonds.atoms)
            mmcif_bonds = bonds(m.bonds.atoms)
            common = pdb_bonds & mmcif_bonds
            extra = pdb_bonds - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra pdb bond(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            extra = mmcif_bonds - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra mmcif bond(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False

            # residues
            pdb_residues = set(p.residues.strs)
            mmcif_residues = set(m.residues.strs)
            common = pdb_residues & mmcif_residues
            extra = pdb_residues - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra pdb residue(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            extra = mmcif_residues - common
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra mmcif residue(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            pdb_helix = p.residues.is_helix
            mmcif_helix = m.residues.is_helix
            for (pr, ph, mr, mh) in zip(pdb_residues, pdb_helix, mmcif_residues, mmcif_helix):
                if ph == mh:
                    continue
                same = False
                if ph:
                    session.logger.info('pdb %s is a helix, and mmcif %s is not' % (pr, mr))
                if mh:
                    session.logger.info('mmcif %s is a helix, and pdb %s is not' % (mr, pr))

            pdb_sheet = p.residues.is_sheet
            mmcif_sheet = m.residues.is_sheet
            for (pr, ps, mr, ms) in zip(pdb_residues, pdb_sheet, mmcif_residues, mmcif_sheet):
                if ps == ms:
                    continue
                same = False
                if ps:
                    session.logger.info('pdb %s is a sheet, and mmcif %s is not' % (pr, mr))
                if ms:
                    session.logger.info('mmcif %s is a sheet, and pdb %s is not' % (mr, pr))

            # chains
            diff = p.num_chains - m.num_chains
            if diff != 0:
                session.logger.error("error: %s: pdb has %d chain(s) vs. %d mmcif chain(s)" % (
                    pdb_id, p.num_chains, m.num_chains))
                same = False

            # coord_sets
            diff = p.num_coord_sets - m.num_coord_sets
            if diff != 0:
                session.logger.error("error: %s: pdb has %d coord_set(s) %s than mmcif" % (
                    pdb_id, abs(diff), "fewer" if diff < 0 else "more"))
                same = False

            # pseudobonds
            pdb_pbg_map = p.pbg_map
            pdb_pbgs = set(pdb_pbg_map.keys())
            mmcif_pbg_map = m.pbg_map
            mmcif_pbgs = set(mmcif_pbg_map.keys())
            common_pbgs = pdb_pbgs & mmcif_pbgs
            extra = pdb_pbgs - common_pbgs
            if extra:
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra pdb pseudobond group(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            pbg = 'missing structure'
            if pbg in extra:
                pdb_bonds = bonds(pdb_pbg_map[pbg].pseudobonds.atoms)
                pdb_bonds = list(pdb_bonds)
                pdb_bonds.sort()
                session.logger.error('error: %s: %d extra pdb %s bond(s): %s' % (
                    pdb_id, len(extra), pbg, pdb_bonds))
            extra = mmcif_pbgs - common_pbgs
            if extra and not (len(extra) == 1 and 'hydrogen bonds' in extra):
                # mmcif file having hydrogen bonds and pdb not is normal
                extra = list(extra)
                extra.sort()
                session.logger.error('error: %s: %d extra mmcif pseudobond group(s): %s' % (
                    pdb_id, len(extra), extra))
                same = False
            pbg = 'missing structure'
            if pbg in extra:
                mmcif_bonds = bonds(mmcif_pbg_map[pbg].pseudobonds.atoms)
                mmcif_bonds = list(mmcif_bonds)
                mmcif_bonds.sort()
                session.logger.error('error: %s: %d extra mmcif %s bond(s): %s' % (
                    pdb_id, len(extra), pbg, mmcif_bonds))
            for pbg in common_pbgs:
                pdb_bonds = bonds(pdb_pbg_map[pbg].pseudobonds.atoms)
                mmcif_bonds = bonds(mmcif_pbg_map[pbg].pseudobonds.atoms)
                common = pdb_bonds & mmcif_bonds
                extra = pdb_bonds - common
                if extra:
                    extra = list(extra)
                    extra.sort()
                    session.logger.error('error: %s: %d extra pdb %s bond(s): %s' % (
                        pdb_id, len(extra), pbg, pdb_bonds))
                    same = False
                extra = mmcif_bonds - common
                if extra:
                    extra = list(extra)
                    extra.sort()
                    session.logger.error('error: %s: %d extra mmcif %s bond(s): %s' % (
                        pdb_id, len(extra), pbg, extra))
                    same = False

            all_same = all_same and same
        pdb_id = save_pdb_id
    if all_same:
        session.logger.info('same: %s' % pdb_id)
    for m in pdb_models:
        m.delete()
    for m in mmcif_models:
        m.delete()
    return all_same


def file_gen(dir):
    # generate files in 2-character subdirectories of dir
    for root, dirs, files in os.walk(dir):
        if root == dir:
            root = ''
        else:
            root = root[len(dir) + 1:]
        # dirs = [d for d in dirs if len(d) == 2]
        dirs.sort()
        if not root:
            files = []
        else:
            files.sort()
        for f in files:
            yield root, f


def next_info(gen):
    try:
        return next(gen)
    except StopIteration:
        return None


def pdb_id(pdb_file):
    if len(pdb_file) != 11:
        return None
    n, ext = os.path.splitext(pdb_file)
    if ext != '.ent' or not n.startswith('pdb'):
        return None
    return n[3:]


def mmcif_id(mmcif_file):
    if len(mmcif_file) != 8:
        return None
    n, ext = os.path.splitext(mmcif_file)
    if ext != '.cif':
        return None
    return n


def compare_all(session):
    from datetime import datetime, timedelta
    start_time = datetime.now()
    pdb_files = file_gen(PDB_DIR)
    mmcif_files = file_gen(MMCIF_DIR)

    pdb_info = next_info(pdb_files)
    mmcif_info = next_info(mmcif_files)

    all_same = True
    while pdb_info and mmcif_info:
        pdb_dir, pdb_file = pdb_info
        pid = pdb_id(pdb_file)
        mmcif_dir, mmcif_file = mmcif_info
        mid = mmcif_id(mmcif_file)
        if (pdb_dir < mmcif_dir or
                (pdb_dir == mmcif_dir and
                 (pid is not None and mid is not None and pid < mid) or
                 ((pid is None or mid is None) and pdb_file < mmcif_file))):
            session.logger.info('Skipping pdb: %s' % os.path.join(pdb_dir, pdb_file))
            pdb_info = next_info(pdb_files)
            continue
        if (mmcif_dir < pdb_dir or
                (pdb_dir == mmcif_dir and
                 (pid is not None and mid is not None and mid < pid) or
                 ((pid is None or mid is None) and mmcif_file < pdb_file))):
            session.logger.info('Skipping mmcif: %s' % os.path.join(mmcif_dir, mmcif_file))
            mmcif_info = next_info(mmcif_files)
            continue
        assert(pid == mid)
        same = compare(session, pid, os.path.join(PDB_DIR, pdb_dir, pdb_file),
                       os.path.join(MMCIF_DIR, mmcif_dir, mmcif_file))
        all_same = all_same and same
        pdb_info = next_info(pdb_files)
        mmcif_info = next_info(mmcif_files)
    end_time = datetime.now()
    delta = end_time - start_time
    days = delta // timedelta(days=1)
    delta -= days * timedelta(days=1)
    hours = delta // timedelta(hours=1)
    delta -= hours * timedelta(hours=1)
    minutes = delta // timedelta(minutes=1)
    delta -= minutes * timedelta(minutes=1)
    seconds = delta // timedelta(seconds=1)
    delta -= seconds * timedelta(seconds=1)
    microseconds = delta // timedelta(microseconds=1)
    session.logger.info('Total time: %d days, %d hours, %d minutes, %d seconds, %d microseconds' % (days, hours, minutes, seconds, microseconds))
    session.logger.info('Total time: %s' % (end_time - start_time))
    raise SystemExit(os.EX_OK if all_same else os.EX_DATAERR)


def compare_id(session, pdb_id):
    if len(pdb_id) != 4:
        session.logger.error('PDB ids should be 4 characters long')
        raise SystemExit(os.EX_DATAERR)
    pdb_id = pdb_id.lower()
    if os.path.exists(PDB_DIR):
        pdb_path = os.path.join(PDB_DIR, pdb_id[1:3], 'pdb%s.ent' % pdb_id)
    else:
        pdb_path = pdb_id
    if os.path.exists(MMCIF_DIR):
        mmcif_path = os.path.join(MMCIF_DIR, pdb_id[1:3], '%s.cif' % pdb_id)
    else:
        mmcif_path = pdb_id
    return compare(session, pdb_id, pdb_path, mmcif_path)


def usage():
    import sys
    session.logger.warning(
        'usage: %s: [-a] [-h] [--all] [--help] [pdb_id(s)]' % sys.argv[0])
    session.logger.warning(
        '''Compare the structures produced by the PDB and mmCIF readers.
        Give one or more pdb identifiers to compare just those structures.
        Or give --all (-a) option to compare all of the strctures in
        /databases/mol/{pdb,mmCIF}/.''')


def main():
    import getopt
    import sys
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "ahi:", ["all", "help"])
    except getopt.GetoptError as err:
        session.logger.error(err)
        usage()
        raise SystemExit(os.EX_USAGE)
    all = False
    for opt, arg in opts:
        if opt in ('-a', '--all'):
            all = True
        elif opt in ('-h', '--help'):
            usage()
            raise SystemExit(os.EX_OK)
    if not all and not args:
        usage()
        raise SystemExit(os.EX_USAGE)
    if all:
        if not os.path.exists(PDB_DIR) or not os.path.exists(MMCIF_DIR):
            session.logger.error("pdb and/or mmCIF databases missing")
            raise SystemExit(os.EX_DATAERR)
        compare_all(session)
    same = True
    for pdb_id in args:
        is_same = compare_id(session, pdb_id)
        if not is_same:
            same = False
    raise SystemExit(os.EX_OK if same else os.EX_DATAERR)


if __name__ == '__main__':
    main()
