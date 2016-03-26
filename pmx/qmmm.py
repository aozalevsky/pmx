#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pmx
import copy
import numpy as np
from pmx import Model, forcefield2
from pmx.ndx import IndexGroup, IndexFile
from pmx.futil import backup_output, Error


class QMsystem(object):

    GMXLIB = '/usr/share/gromacs/top'
    LA = 'LA'  # Type of linking atoms
    breakable_bonds = [
        ('N', 'CA'),  # break backbone at beggining
        ('N', 'CD')  # break backbone at beggining in PRO
        ('CA', 'CB'),  # break sidechain
        ('CA', 'C')  # break backbone near end
    ]

    iqm = None
    iqmfn = None

    itop = None
    itopfn = None
    otop = None
    otopfn = None

    igro = None
    igrofn = None
    ogro = None
    ogrofn = None

    indx = None
    indxfn = None
    ondx = None
    ondxfn = None

    system = list()
    bonds2break = list()
    vsites2 = list()
    group = IndexGroup()

    def __init__(
            self,
            iqmfn=None,
            itopfn=None,
            otopfn=None,
            igrofn=None,
            ogrofn=None,
            indxfn=None,
            ondxfn=None,
            *args, **kwargs):
        """Initialize class instance"""
        self.iqmfn = iqmfn

        self.itopfn = itopfn
        self.otopfn = otopfn

        self.igrofn = igrofn
        self.ogrofn = ogrofn

        self.indxfn = indxfn
        self.ondxfn = ondxfn

    def open_inputs(self):
        try:
            with open(self.iqmfn, 'r') as f:
                iqm = eval(f.read())
            if len(iqm) == 0:
                Error('Empty qmsystem')
            self.iqm = iqm
        except:
            Error('Unable to qm input %s' % self.iqmfn)

        try:
            self.itop = forcefield2.Topology(self.itopfn)
        except:
            Error('Unable to open input top %s' % self.itopfn)

        if not self.otopfn:
            self.otopfn = self.igrofn.replace('.top', '_qm.top')

        try:
            self.igro = Model(self.igrofn)
        except:
            Error('Unable to open input gro %s' % self.igrofn)

        if not self.ogrofn:
            self.ogrofn = self.igrofn.replace('.gro', '_qm.gro')

        if self.indxfn:
            try:
                self.indx = IndexFile(self.indxfn)
            except:
                Error('Unable to open input ndx %s' % self.indxfn)

        else:
            self.indx = IndexFile()
            group = IndexGroup(name='freeze', atoms=self.itop.atoms)
            self.indx.add_group(group)

        if not self.ondxfn:
            if self.indxfn:
                self.ondxfn = self.indxfn.replace('.ndx', '_qm.ndx')
            else:
                self.ondxfn = 'qm.ndx'

    def process(self):
        """Main function"""

        self.read_inputs()

        self.process_all()

        self.write_outputs()

    def read_inputs(self):
        self.open_inputs()
        self.add_residues(self.iqm)

    def process_all(self):

        self.process_topology()
        self.process_coordinates()
        self.process_index()

    def write_outputs(self):

        self.otop.write(backup_output(self.otopfn))
        self.ogro.write(backup_output(self.ogrofn))
        self.ondx.write(backup_output(self.ondxfn))

    def add_residues(self, residues):
        for i in residues.items():
            resi, atoms = i
            self.add_residue(resi, **atoms)

    def __missing_atom(self, at, res):
        Error(
            'There is no atom %s in residue %d%s' %
            (at, res.id, res.resname))

    def add_residue(self, resi, include=None, exclude=None):
        if include and exclude:
            Error(
                'Can not use include and exclude directives '
                'simultaneously in residue %d in %s' %
                (resi, self.iqmfn))

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
                    self.__missing_atom(at, res)

        elif exclude:

            for at in exclude:
                try:
                    aind = tnames.index(at)
                    tmask[aind] = False

                except ValueError:
                    self.__missing_atom(at, res)

        self.system = list(set(self.system.extend(tatoms[tmask])))

    def process_coordinates(self):
        self.ogro = copy.copy(self.igro)

        for vs in self.vsites2:
            aLA, qm, mm, t, ratio = vs
            ai = np.array(self.ogro.atoms[qm.id - 1].x)
            aj = np.array(self.ogro.atoms[mm.id - 1].x)
            aLA.x = ai + (aj - ai) * ratio
            self.ogro.atoms.append(aLA)

    def process_topology(self):
        self.otop = copy.copy(self.itop)

        self.check_forcefield()

        self.process_bonds()

        self.process_angles()

        self.process_dihedrals()

        self.add_virtual_sites2()

    def check_forcefield(self):
        if self.LA not in self.itop.NBParams.atomtypes:
            Error(
                """Missing QMMM parameters in forcefield.
Create file qmmm.itp with following content:
; New atom and atomtype for QMMM
LA         0.00000  ; QMMM

[ atomtypes ]
LA           0       0.0000  0.0000  A                0.00000e+00  0.00000e+00

and add it to topology like:
#include "qmmm.itp" ; QMMM parameters
""")

        return True

    def add_virtual_sites2(self):
        """Add section virtual_sites
        [ virtual_sites2 ]
         LA QMatom MMatom 1 X.XXX
        """
        vsites2 = list()

        for bond in self.bonds2break:
            site = list()
            ai, aj = bond

            if (ai in self.system and aj not in self.system):
                qm, mm = ai, aj
            elif (aj in self.system and ai not in self.system):
                qm, mm = aj, ai
            else:
                Error('Something went wrong with breaking bonds')

            rLA = pmx.Molecule()
            rLA.id = len(self.otop.residues) + 1
            rLA.resname = 'XXX'
            rLA.charge = 0

            aLA = pmx.Atom()
            aLA.id = len(self.otop.atoms) + 1
            aLA.atomtype = self.LA
            aLA.atomtypeB = None
            aLA.resnr = aLA.id
            aLA.resname = rLA.resname
            aLA.name = self.LA
            aLA.cgnr = aLA.resnr
            aLA.q = 0.0
            aLA.m = 0.0

            rLA.atoms.append(aLA)

            self.otop.residues.append(rLA)
            self.otop.atoms.append(aLA)

            ratio = self.__get_ratio(qm, mm)

            site = [aLA, qm, mm, 1, ratio]

            vsites2.append(site)

        if len(vsites2) > 0:
            self.vsites2 = vsites2
            self.otop.virtual_sites2.extend(vsites2)
            self.otop.has_vsites2 = True

        return True

    def process_index(self):
        self.ondx = copy.copy(self.indx)

        atoms = list()
        atoms.extend(self.system)
        latoms = map(lambda x: x[0], self.vsites2)
        atoms.extend(latoms)

        group = IndexGroup(name='QM', atoms=atoms)

        self.ondx.add_group(group)
        self.ondx.dic['freeze'].ids.extend(map(lambda x: x.id, latoms))

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

    def process_angles(self, limit=2):
        """Comment out angles of QM system"""
        angles = list()

        for angle in self.itop.angles:
            c = 0
            atoms = angle[:3]

            for a in atoms:
                if a in self.system:
                    c += 1
            if c < limit:
                angles.append(angle)

        self.otop.angles = angles

    def process_dihedrals(self, limit=3):
        """Comment out dihedrals of QM system"""
        angles = list()

        for angle in self.itop.dihedrals:
            c = 0
            atoms = angle[:4]

            for a in atoms:
                if a in self.system:
                    c += 1
            if c < limit:
                angles.append(angle)

        self.otop.dihedrals = angles

    def process_bonds(self):
        """Changes bondtype to 5 for all QMatoms"""

        bonds = list()
        bonds2break = list()

        for bond in self.itop.bonds:
            ai, aj = bond[:2]
            t = bond[2]

            if self.__is_qm_bond(ai, aj):
                t = 5

            elif self.__is_qmmm_bond(ai, aj):
                if self.__is_breakable_bond(ai.name, aj.name):
                    t = 5
                    bonds2break.append((ai, aj))
                else:
                    self.__unbreakable_bond(ai, aj)

            bond[2] = t
            bonds.append(bond)

        self.otop.bonds = bonds
        self.bonds2break = bonds2break

    def __is_qm_bond(self, ai, aj):
        result = False

        if ai in self.system:
            if aj in self.system:
                result = True

        return result

    def __is_qmmm_bond(self, ai, aj):
        result = False

        if (ai in self.system and aj not in self.system) or \
                (aj in self.system and ai not in self.system):

            result = True

        return result

    def __is_breakable_bond(self, ai, aj):
        result = False

        if (ai, aj) in self.breakable_bonds:
            result = True
        elif (aj, ai) in self.breakable_bonds:
            result = True

        return result

    def __unbreakable_bond(self, ai, aj):
        Error(
            'Bond %d (%s) - %d (%s) not in the list of breakable bonds' %
            (ai.id, ai.name, aj.id, aj.name))

    def __get_ratio(self, ai, aj):

        ratio = None

        abond = (ai.atomtype, aj.atomtype)

        for bond in self.itop.BondedParams.bondtypes:
            bondt1 = (bond[0], bond[1])
            bondt2 = (bond[1], bond[0])

            if abond == bondt1 or abond == bondt2:
                    ratio = 0.1 / bond[3]  # 0.1 because of nm in params

        if ratio is None:
            Error(
                'Did not find parameters for breaking bond %d (%s) - %d (%s)' %
                (ai.id, ai.atomtype, aj.id, aj.atomtype))

        return ratio


def run():
    import argparse as ag

    parser = ag.ArgumentParser(
        description="""Processes GROMACS topology file (.top)
        and adds linking atoms (LA) for QM/MM.""")

    parser.add_argument('-q', '--qm-system',
                        type=str,
                        help='Input file with description of qm system',
                        dest='iqmfn',
                        required=True)
    parser.add_argument('-p', '--input-topology',
                        type=str,
                        help='Input topology',
                        dest='itopfn',
                        required=True)
    parser.add_argument('-c', '--input-coordinates',
                        type=str,
                        dest='igrofn',
                        help='Index file with QM/MM group',
                        required=True)
    parser.add_argument('-n', '--input-index',
                        type=str,
                        dest='indxfn',
                        help='Index file')
    parser.add_argument('-o', '--output-topology',
                        dest='otopfn',
                        type=str,
                        help='Output topology file')
    parser.add_argument('-y', '--output-coordinates',
                        type=str,
                        dest='ogrofn',
                        help='Output coordinates file')
    parser.add_argument('-x', '--output-index',
                        type=str,
                        dest='ondxfn',
                        help='Output index file')
    args = vars(parser.parse_args())

    QMMM = QMsystem(**args)
    QMMM.process()

if __name__ == '__main__':
    run()
