#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys
import os
import os.path
import re
import pyparsing as pg
import numpy as np
from collections import OrderedDict
from pmx import Model, ffparser, forcefield2
from pmx.ndx import IndexGroup, IndexFile


class QMsystem(object):

    GMXLIB = '/usr/share/gromacs/top'
    LA = 'LA'  # Type of linking atoms

    igro = None
    ogro = None

    itop = None
    otop = None

    indx = None
    ondx = None

    system = dict()
    la_ind = list()
    group = IndexGroup()

    def __init__(self, *args, **kwargs):
        """Initialize class instance"""
        self.args = args

        try:
            self.igro = Model(kwargs['igro'])
        except:
            sys.stderr.write('Unable to open input gro %s' % kwargs['igro'])
            sys.exit(1)

        try:
            self.itop = forcefield2.Topology(kwargs['itop'])
        except:
            sys.stderr.write('Unable to open input top %s' % kwargs['itop'])
            sys.exit(1)

        if kwargs['indx']:
            self.indx = IndexFile(kwargs['indx'])

        if 'GMXLIB' in os.environ:
            self.GMXLIB = os.environ['GMXLIB']

    def process(self):
        """Main function"""
        if 'qmsystem' in self.args:
            self.add_residues(self.args['qmsystem'])
        else:
            sys.stderr.write('Empty qmsystem')
            sys.exit(1)

        self.process_topology()

        self.write_index(self.args.output_index)

        self.add_la_coordinates()
        self.gro.write(self.args.output_coordinates)

    def __add_residues(self, residues):
        for i in residues.items():
            resi, atoms = i
            self.__add_residue(resi, **atoms)

    def __missing_atom(at, resi, resn):
        sys.stderr.write(
            'There are no atom %s in residue %d%s') % \
            (at, resi, resn)
        sys.exit(1)

    def __add_residue(self, resi, include=None, exclude=None):
        if include and exclude:
            sys.stderr.write(
                'Can not use include and exclude simultaneously')
            sys.exit(1)

        atoms = list()

        res = self.itop.residues[resi - 1]
        tatoms = np.array(res.atoms)
        tnames = map(lambda x: x.name, tatoms)

        tmask = np.zeros(len(tatoms), dtype=bool)
        tmask.fill(True)

        if include:
            tmask = np.zeros(len(tatoms), dtype=bool)

            for at in include:
                try:
                    aind = tnames.index(at)
                    tmask[aind] = True

                except ValueError:
                    self.__missing_atom(at, res.id, res.resname)

        elif exclude:

            for at in exclude:
                try:
                    aind = tnames.index(at)
                    tmask[aind] = False

                except ValueError:
                    self.__missing_atom(at, res.id, res.resname)

        self.system[resi] = atoms[tmask]

    def add_la_coordinates(self):

        atoms = list()
        num = len(self.bonds2break)

        for i in range(num):

            atom = list()
            top = self.la_atoms[i]
            atom.append(top[2])  # resi
            atom.append(top[3])  # resn
            atom.append(top[4])  # aname
            atom.append(top[0])  # num

            bb = self.bonds2break[i]
            ain, ajn = bb[:2]
            ai = np.array(self.gro.atoms[ain - 1][4:])
            aj = np.array(self.gro.atoms[ajn - 1][4:])

            la = ai + (aj - ai) * self.get_ratio(*bb[:4])

            atom.extend(la)
            atoms.append(atom)
        for i in range(num):
            self.gro.atoms.insert(self.la_ind[i] - 1, atoms[i])
        for i in range(self.la_ind[-1], len(self.gro.atoms)):
            print i
            self.gro.atoms[i][0] += num
#        self.gro.atoms.extend(atoms)

    def process_topology(self, fi=None, fo=None):

        self.check_forcefield()

        self.process_bonds(fi, fo)

        self.add_section_virtual_sites2(fo)

        self.process_angles(fi, fo)

        self.process_dihedrals(fi, fo)

    def check_forcefield(self, f):
        self.ff_name = self.detect_ff(f)

        lpath = self.ff_name
        gpath = os.path.join(self.GMXLIB, self.ff_name)

        if os.path.isdir(lpath):
            self.ff_path = lpath
        elif os.path.isdir(gpath):
            self.ff_path = gpath
        else:
            sys.stderr.write(
                'Unable to locate directory for {} forcefield'.format(
                    self.ff_name))
            sys.exit(1)

        print 'Forcefield {0} in {1} directory'.format(
            self.ff_name, self.ffpath)

        atpfn = os.path.join(self.ff_path, 'atomtypes.atp')
        atypes = ffparser.ATPParser(atpfn)

        if self.LA not in atypes.dic:
            sys.stderr.write(
                'Missing {} atomtype in atomtypes.atp'.format(self.LA))
            sys.stderr.write(
                'Please, add following string to {}'.format(atpfn))
            sys.stderr.write(
                'LA         0.00000  ; QMMM')
            sys.exit(1)

        nbfn = os.path.join(self.ff_path, 'ffnonbonded.itp')
        nonbonded = ffparser.NBParser(nbfn)

        if self.LA not in nonbonded.atomtypes:
            sys.stderr.write(
                'Missing {} atomtype in ffnonbonded.itp'.format(self.LA))
            sys.stderr.write(
                'Please, add following string to {}'.format(nbfn))
            sys.stderr.write(
                'LA           0       0.0000  0.0000  A   \
                0.00000e+00  0.00000e+00')
            sys.exit(1)

        sys.stderr.write(
            'Found {} atomtype in forcefield'.format(self.LA))

        return True

    def add_section_virtual_sites2(self, fo):
        """Add section virtual_sites
        [ virtual_sites2 ]
         LA QMatom MMatom 1 X.XXX
        """

        fmt = ['{:5d}', '{:5d}', '{:5d}', '    1', '{:8.3f}', '\n']

        fo.write('\n')
        fo.write('; Linkage atoms for QM/MM\n')
        fo.write('[ virtual_sites2 ]\n')
        for i in range(len(self.la_ind)):
            site = list()
            site.append(self.la_ind[i])
            bb = self.bonds2break[i]
            site.extend(bb[:2])
            site.append(self.get_ratio(*bb[:4]))
            fo.write(' '.join(fmt).format(*site))
        fo.write('\n')

        return True

    def adjust_index(self):
        """Adjust indexes after adding of linking atoms"""
        def adjust(start, offset, v):
            if v > start:
                return v + offset
            else:
                return v
        for i in self.groups.keys():
            self.groups[i] = map(lambda x: adjust(
                self.la_ind[0], len(self.la_ind), x), self.groups[i])

        self.groups[self.group].extend(self.la_ind)
        self.groups['System'].extend(self.la_ind)

    def process_dihedrals(self, fi, fo):
        """Comment out dihedrals of QM system"""

        group = self.groups[self.group]

        p = fi.tell()

        end = re.compile('\[')
        com = re.compile('(;|#)')

        limit = 3

        it = pg.Word(pg.nums).setParseAction(lambda v: int(v[0]))

        ai = it
        aj = it
        ak = it
        al = it
        func = it
        comment = pg.Optional(pg.Literal(';') + pg.restOfLine)

        parser = ai + aj + ak + al + func + comment

        l = fi.readline()
        while l:
            c = 0
            if end.match(l.lstrip()):
                fi.seek(p)
                break
            elif com.match(l.strip()):
                pass
            elif l.strip():
                atoms = parser.parseString(l).asList()[0:3]
                for i in atoms:
                    if i in group:
                        c += 1
                if c >= limit:
                    l = ';' + l

            fo.write(l)
            p = fi.tell()
            l = fi.readline()

        return True
        pass

    def process_angles(self, fi, fo):
        """Comment out angles of QM system"""
        group = self.groups[self.group]
        p = fi.tell()

        end = re.compile('\[')
        com = re.compile('(;|#)')

        limit = 2

        it = pg.Word(pg.nums).setParseAction(lambda v: int(v[0]))

        ai = it
        aj = it
        ak = it
        func = it
        comment = pg.Optional(pg.Literal(';') + pg.restOfLine)

        parser = ai + aj + ak + func + comment

        l = fi.readline()
        while l:
            c = 0
            if end.match(l.lstrip()):
                fi.seek(p)
                break
            elif com.match(l.strip()):
                pass
            elif l.strip():
                atoms = parser.parseString(l).asList()[0:3]
                for i in atoms:
                    if i in group:
                        c += 1
                if c >= limit:
                    l = ';' + l

            fo.write(l)
            p = fi.tell()
            l = fi.readline()

        return True

    def process_bonds(self, fi, fo):
        """Changes bondtype to 5 for all QMatoms"""
        group = self.groups[self.group]
        p = fi.tell()

        end = re.compile('\[')
        com = re.compile('(;|#)')

        fmt = ['{:5d}', '{:5d}', '{:5d}']
        limit = 2

        it = pg.Word(pg.nums).setParseAction(lambda v: int(v[0]))

        ai = it
        aj = it
        func = it
        comment = pg.Optional(pg.Literal(';') + pg.restOfLine)

        parser = ai + aj + func + comment

        l = fi.readline()
        while l:
            c = 0
            if end.match(l.lstrip()):
                fi.seek(p)
                break
            elif com.match(l.strip()):
                pass
            elif l.strip():
                atoms = parser.parseString(l).asList()[0:2]
                for i in atoms:
                    if i in group:
                        c += 1
                if c >= limit:
                    atoms.append(5)
                    l = ' '.join(fmt).format(*atoms) + ' ; QM system\n'

            fo.write(l)
            p = fi.tell()
            l = fi.readline()

        return True

    def add_la_atoms(self, f, atype=None, astart=None, rstart=None, num=None):
        """Add linking atoms to [ atoms ] section in topology"""
        if astart is None:
            astart = len(self.atoms)

        if num is None:
            num = len(self.bonds2break)

        if atype is None:
            atype = self.LA

        if rstart is None:
            rstart = self.atoms[-1][1]

        # Copied from GromacsWrapper
        fmt = ["{:6d}", "{:>10s}", "{:-6d}", "{:>6s}",
               "{:>6s}", "{:>6d}",
               "{:-10.4f}", "{:-10.3f}",
               ]
        atoms = list()
        la = list()
        f.write(';Linking atoms for QM/MM\n')
        for i in range(num):
            c = i + 1
            la.append(astart + c)
            atom = [astart + c, atype, rstart + c,
                    'XXX', atype, astart + c, 0.0, 0.0,
                    ]
            atoms.append(atom)
            astr = ' '.join(fmt).format(*atom)
            f.write(astr + '\n')

        return(la, atoms)

    def check_section_name(self, fi, fo, name=None):
        l = fi.readline()
        if re.match('\[ {} \]'.format(name), l):
            fo.write(l)
        else:
            raise RuntimeError(
                'Wrong section. {} is not [ {} ]'.format(l, name))

        return True

    def detect_ff(self, f):
        """Reads topology and guesses forcefield path"""
        inc = re.compile('#include')
        for l in self.itop.header:
            if inc.match(l.lstrip()):
                ff = re.search('"(.*\.ff)', l).group(1)
                break
        if not ff:
            raise RuntimeError('Unable to determine forcefield from topology')
        return ff

    def select_group(self, group=None, groups=None):

        if group is None:
            group = self.group

        if groups is None:
            groups = self.groups

        def print_groups():
            print '\n'
            print 'Groups in index file:\n'
            for i, v in enumerate(groups.keys()):
                print "Group %6d (%16s) has %5d elements" % (
                    i, v, len(groups[v]))

        def resolve_group(groups, i):
            return groups[groups.keys()[i]]

        while (
                group is None or
                group < 0 or
                group > len(self.groups) or
                len(resolve_group(groups, group)) == 0):
            print_groups()
            group = int(raw_input("Select QM/MM group: "))

        return group.keys()[group]

    def read_la_list(self, f=None):
        """
        Reads LA list in format
        QM MM ; Optional comment
        """
        if f is None:
            f = self.args.la_list

        atype = pg.Word(pg.alphas)
        it = pg.Word(pg.nums).setParseAction(lambda v: int(v[0]))

        ai = it
        aj = it
        aitype = atype
        htype = atype

        comment = pg.Optional(pg.Literal(';') + pg.restOfLine)

        parser = ai + aj + aitype + htype + comment

        bonds = list()
        for line in f:
            bonds.append(parser.parseString(line).asList())

        if len(bonds) > 0:
            print 'Found following bonds to breake:'
            for i in bonds:
                print '{0} - {1}'.format(*i)
        else:
            raise RuntimeError('Empty linking atom list')

        return(bonds)


if __name__ == '__main__':
    import glob
    import argparse as ag
    import shutil

    def backup_output(fn):
        try:
            if os.stat(fn).st_size:
                path = os.path.dirname(fn)
                name = os.path.basename(fn)
                guess = sorted(map(
                    lambda x: int(re.sub(r'#.*\.(\d+)#', r'\1', x)),
                    glob.glob(path + "#" + name + ".*#")))
                if not guess:
                    i = 1
                else:
                    i = guess[-1]
                    if i < 99:
                        i += 1
                    else:
                        raise Exception("Too many backups")
                backup = "#" + name + '.' + str(i) + "#"
                backup = os.path.join(path, backup)
                shutil.copy(fn, backup)
                print "Back Off! I just backed up {0} to {1}".format(fn, backup)
        except OSError:
            pass

        return open(fn, 'w')

    parser = ag.ArgumentParser(
        description="""Processes GROMACS topology file (.top)
        and adds linking atoms (LA) for QM/MM.""")

    parser.add_argument('-l', '--la-list', type=ag.FileType('r'),
                        help='List of bonds to breake', required=True)
    parser.add_argument('-p', '--input-topology', type=ag.FileType('r'),
                        help='Input topology', required=True)
    parser.add_argument('-n', '--input-index', type=ag.FileType('r'),
                        help='Index file with QM/MM group', required=True)
    parser.add_argument('-c', '--input-coordinates', type=ag.FileType('r'),
                        help='Index file with QM/MM group', required=True)
    parser.add_argument('-o', '--output-topology', type=backup_output,
                        help='Output PDB file', required=True)
    parser.add_argument('-x', '--output-index', type=backup_output,
                        help='Output index file', required=True)
    parser.add_argument('-y', '--output-coordinates', type=backup_output,
                        help='Output coordinates file', required=True)
    get_args = parser.parse_args()

    add = QMsystem(get_args)
    process = add.process()
