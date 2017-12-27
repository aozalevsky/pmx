#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pmx
import copy
import numpy as np
from pmx import forcefield2
from .model import Model

from pmx.ndx import IndexGroup, IndexFile
from pmx.futil import backup_output, Error
from .atomselection import Atomselection
from collections import OrderedDict as OD


class QMsystem(object):

    GMXLIB = '/usr/share/gromacs/top'
    LA = 'LA'  # Type of linking atoms
    breakable_bonds = [
        ('N', 'CA'),  # break backbone at beggining
        ('N', 'CD'),  # break backbone at beggining in PRO
        ('CA', 'CB'),  # break sidechain
        ('CA', 'C')  # break backbone near end
    ]

    sol_ions = [
        'HOH',
        'SOL',
        'NA',
        'CL',
        'K',
        'MG',
        'CRW',
    ]

    sol_ions_type = {
        'OW': {'type': 'OW', 'mass': 16.00000, 'charge': -0.834},
        'HW1': {'type': 'HW', 'mass': 1.00800, 'charge': 0.417},
        'HW2': {'type': 'HW', 'mass': 1.00800, 'charge': 0.417},
        'NA': {'type': 'Na', 'mass': 22.99, 'charge': 1.0000},
        'K': {'type': 'K', 'mass': 39.10, 'charge': 1.0000},
        'Mg': {'type': 'MG', 'mass': 24.305, 'charge': 2.0000},
        'MG': {'type': 'MG', 'mass': 24.305, 'charge': 2.0000},
        'CL': {'type': 'Cl', 'mass': 35.45, 'charge': -1.0000},
    }

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

    refine = False

    def __init__(
            self,
            iqmfn=None,
            itopfn=None,
            otopfn=None,
            igrofn=None,
            ogrofn=None,
            indxfn=None,
            ondxfn=None,
            refine=False,
            *args, **kwargs):
        """Initialize class instance"""
        self.iqmfn = iqmfn

        self.itopfn = itopfn
        self.otopfn = otopfn

        self.igrofn = igrofn
        self.ogrofn = ogrofn

        self.indxfn = indxfn
        self.ondxfn = ondxfn

        self.refine = refine

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
        except Exception:
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

        if not self.ondxfn:
            if self.indxfn:
                self.ondxfn = self.indxfn.replace('.ndx', '_qm.ndx')
            else:
                self.ondxfn = 'qm.ndx'

    def process(self, refine=False):
        """Main function"""

        self.read_inputs()

        if refine or self.refine:
            self.refine_system()

        self.main_job()

    def main_job(self):
        self.init_outputs()

        self.process_topology()

        self.write_outputs()

    def read_inputs(self):
        self.open_inputs()
        self.add_residues(self.iqm)

    def init_outputs(self):
        self.otop = self.itop
        self.ogro = self.igro
        self.ondx = self.indx

    def write_outputs(self):

        self.otop.write(backup_output(self.otopfn))
        self.ogro.write(backup_output(self.ogrofn))
        self.ondx.write(backup_output(self.ondxfn))

    def add_residues(self, residues):
        for i in residues.items():
            resi, atoms = i
#            resn = self.igro.residues[resi - 1].resname
#            if resn not in self.sol_ions:
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

        res = self.igro.residues[resi - 1]
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

        self.system.extend(tatoms[tmask])
        self.system = list(set(self.system))

    def check_qm_system(self):

        print(len(self.system))
        print(len(set(self.system)))

    def get_system_charge(self, atoms=None):
        if not atoms:
            atoms = self.system

        satoms = pmx.atomselection.Atomselection(atoms=atoms)

        atoms = satoms.fetch_atoms(satoms.backbone, inv=True)

        satoms = pmx.atomselection.Atomselection(atoms=atoms)
        atoms = satoms.expand_byres(self.igro.atoms)

        atoms = self.swap_gro2top_ids(atoms, self.itop.atoms)

        fcharge = sum((map(lambda x: x.q, atoms)))
        icharge = int(round(fcharge))
        if np.isclose(fcharge, icharge, rtol=1e-03):
            return icharge
        else:
            print(fcharge, icharge)
            print('QM charge %.6f is not roundable' % fcharge)
            return icharge

    def process_topology(self):

        self.check_forcefield()
        self.swap_waters()
        self.charge = self.get_system_charge(self.system)
        self.process_bonds()
        self.process_angles()
        self.process_dihedrals()
        self.add_virtual_sites2()
        self.process_index()

    def refine_system(self):

        self.system = self.swap_gro2top_ids(self.system, self.itop.atoms)

        ssystem = pmx.atomselection.Atomselection(atoms=self.system)
        self.system = ssystem.expand_byres_minimal(self.igro.atoms)
        self.system = self.swap_gro2top_ids(self.system, self.itop.atoms)

        self.system_raw = copy.deepcopy(self.system)

        self.system = self.expand_ends(self.system)
        self.system = self.swap_gro2top_ids(self.system, self.itop.atoms)

        self.charge = self.get_system_charge(self.system)

        print('Total size of QM system: %d; Charge: %.1f' % (
            len(self.system), self.charge))

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

    def strip_vsites2(self):
        # if self.otop.has_vsites2:
        #    pass
        # else:
        #    return

        syss = Atomselection(atoms=self.system)
        gros = Atomselection(atoms=self.ogro.atoms)
        tops = Atomselection(atoms=self.otop.atoms)
        tlacount = len(tops.fetch_atoms('LA'))
        glacount = len(gros.fetch_atoms('LA'))

        gro = gros.fetch_atoms('LA', inv=True)
        new = Model(atoms=gro)
        self.ogro.atoms = new.atoms

        top = tops.fetch_atoms('LA', inv=True)
        new = Model(atoms=top)
        self.otop.atoms = new.atoms

        sys = syss.fetch_atoms('LA', inv=True)
        self.system = sys

        self.otop.virtual_sites2 = []
        self.otop.has_vsites2 = False

        print('Stripped previous LA atoms: %d' % max(tlacount, glacount))

    def add_virtual_sites2(self):
        """Add section virtual_sites
        [ virtual_sites2 ]
         LA QMatom MMatom 1 X.XXX
        """
        self.strip_vsites2()

        vsites2 = list()

        count = 0
        for bond in self.bonds2break:
            site = list()
            ai, aj = bond

            if (ai in self.system and aj not in self.system):
                qm, mm = ai, aj
            elif (aj in self.system and ai not in self.system):
                qm, mm = aj, ai
            else:
                Error('Something went wrong with breaking bonds')

            aLA = pmx.Atom()
            aLA.id = len(self.otop.atoms) + count
            aLA.cgnr = aLA.id + 1
            aLA.atomtype = self.LA
            aLA.atomtypeB = None
            aLA.name = self.LA
            aLA.q = 0.0
            aLA.m = 0.0
            aLA.resname = 'XXX'
            aLA.resnr = len(self.otop.residues) + count + 1
            aLA.symbol = 'H'

            ratio = self.__get_ratio(qm, mm)

            ai = np.array(self.ogro.atoms[qm.id - 1].x)
            aj = np.array(self.ogro.atoms[mm.id - 1].x)
            aLA.x = ai + (aj - ai) * ratio

            site = [aLA, qm, mm, 1, ratio, '; qmmm']

            vsites2.append(site)

            count += 1

        print('Total new virtual sites: ', len(vsites2))

        if len(vsites2) > 0:
            self.vsites2 = vsites2

            # strip old vsites2
            # self.otop.virtual_sites2.extend(vsites2)
            self.otop.virtual_sites2 = vsites2
            self.otop.has_vsites2 = True

            self.update_atoms(map(lambda x: x[0], vsites2))

        return True

    def update_atoms(self, atoms):

        start = len(self.otop.atoms)
        l = len(atoms)

        newtop_atoms = self.otop.atoms
        newtop_atoms.extend(atoms)
        newtop = Model(atoms=newtop_atoms)

        newgro_atoms = self.ogro.atoms[:start]
        newgro_atoms.extend(atoms)
        newgro_atoms.extend(self.ogro.atoms[start:])
        newgro = Model(atoms=newgro_atoms)

        self.otop.atoms = newtop.atoms
        self.ogro.atoms = newgro.atoms

        for a in newgro.atoms:
            a.cgnr = a.id

        return range(start + 1, start + 1 + l)

    def process_index(self):
        atoms = list()
        atoms.extend(self.system)
        latoms = map(lambda x: x[0], self.vsites2)
        atoms.extend(latoms)

        group_qm = IndexGroup(name='QM', atoms=atoms)

        if 'QM' in self.ondx.dic:
            self.ondx.delete_group('QM')

        self.ondx.add_group(group_qm)

        group_freeze = IndexGroup(name='freeze', atoms=self.otop.atoms)

        if 'freeze' in self.ondx.dic:
            self.ondx.delete_group('freeze')

        self.ondx.add_group(group_freeze)

        group_sol_and_ions = IndexGroup(
            name='Water_and_ions',
            atoms=self.ogro.atoms[len(self.otop.atoms):])

        if 'Water_and_ions' in self.ondx.dic:
            self.ondx.delete_group('Water_and_ions')

        self.ondx.add_group(group_sol_and_ions)

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

        sys = set(self.system)

        for angle in self.itop.angles:
            atoms = angle[:3]
            c = len(set.intersection(sys, atoms))

            if c < limit:
                angles.append(angle)

        self.otop.angles = angles

    def process_dihedrals(self, limit=3):
        """Comment out dihedrals of QM system"""
        angles = list()

        sys = set(self.system)

        for angle in self.itop.dihedrals:
            atoms = angle[:4]
            c = len(set.intersection(sys, atoms))

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

    def expand_ends(self, atoms):
        result = [i for i in atoms]

        resnrs = list(set(map(lambda x: x.resnr, atoms)))

        for i in resnrs:
            tatoms = Atomselection(atoms=atoms).fetch_atoms(i, how='byresnr')
            tresult = self.expand_terminus(tatoms)
            result.extend(tresult)

        result = list(set(result))
        result = sorted(result, key=lambda x: x.id)

        return result

    def expand_terminus(self, res):
        result = res

        anames = set(map(lambda x: x.name, res))

        if 'N' in anames:
            tresult = self.expand_n_terminus(res)
            if tresult:
                result.extend(tresult)

        if 'C' in anames:
            tresult = self.expand_c_terminus(res)
            if tresult:
                result.extend(tresult)

        result = list(set(result))
        result = sorted(result, key=lambda x: x.id)

        return result

    def expand_n_terminus(self, res):
        result = list()

        resnr = res[0].resnr

        if resnr == 1:
            return result

        mmsys = Atomselection(atoms=self.itop.atoms)
        preres = Atomselection(
            atoms=mmsys.fetch_atoms(resnr - 1, how='byresnr'))

        N = Atomselection(atoms=res).fetch_atoms(['N'])[0]

        try:
            C = preres.fetch_atoms(['C'])[0]
        except:
            Error('Bad residue %d %s' % (res[0].resnr, res[0].resname))

        if self.itop.is_bond(N, C):
            result = preres.fetch_atoms(['C', 'O'])

        result.extend(res)

        return result

    def expand_c_terminus(self, res):
        C, N = None, None
        result = [i for i in res]

        resnr = res[-1].resnr

        C = Atomselection(atoms=res).fetch_atoms(['C'])[-1]

        mmsys = Atomselection(atoms=self.itop.atoms)

        try:
            preres = Atomselection(
                atoms=mmsys.fetch_atoms(resnr + 1, how='byresnr'))
            N = preres.fetch_atoms(['N'])[0]
        except:
            print('Warning: May be terminal residue?')

        if (C and N) and self.itop.is_bond(C, N):
            result.extend(preres.fetch_atoms(['N', 'H']))

        return result

    def strip_sol_ions(self):
        molecules = OD(self.otop.molecules)

        top = self.otop.atoms
        tops = Atomselection(atoms=top)
        start = len(self.otop.atoms)

        # gro = self.ogro.atoms
        # gros = Atomselection(atoms=gro)

        print("Searching sol_ions")

        for r in self.sol_ions:
            rtop = tops.fetch_atoms(r, how='byresname')
            if rtop:
                nr = len(set(map(lambda x: x.resnr, rtop)))
                try:
                    molecules[r] += nr
                except KeyError:
                    molecules[r] = nr

                top = tops.fetch_atoms(r, how='byresname', inv=True)
                tops = Atomselection(atoms=top)

                # gro = gros.fetch_atoms(r, how='byresname', inv=True)
                # gros = Atomselection(atoms=gro)

        new = Model(atoms=top)
        self.otop.atoms = new.atoms
        end = len(self.otop.atoms)

        self.otop.molecules = map(list, molecules.items())

        # new = Model(atoms=gro)
        # self.ogro.atoms = new.atoms

        print('Stripped %d sol and ions atoms from topology' % (start - end))

    def swap_waters(self):
        self.strip_sol_ions()

        res = OD(map(lambda x: (x.resnr, x.resname), self.system))

        if len(set.intersection(set(self.sol_ions), set(res.values()))) > 0:
            pass
        else:
            return

        sogro = Atomselection(atoms=self.ogro.atoms)

        res2top = []

        molecules = OD(self.otop.molecules)

        for k, v in res.items():
            if v in self.sol_ions:
                res2top.append(k)
                try:
                    molecules[v] -= 1
                except KeyError:
                    pass

        self.otop.molecules = map(list, molecules.items())

        a2top = sogro.fetch_atoms(res2top, how='byresnr')

        for i in range(len(a2top)):
            a2top[i] = self.assign_sol_ions_type(a2top[i])

        a2top_inv = sogro.fetch_atoms(res2top, how='byresnr', inv=True)

        self.ogro.atoms = a2top_inv
        self.update_atoms(a2top)

        print('Found water and/or ions in QM system. Total: %d' % len(a2top))

    def assign_sol_ions_type(self, a):

        a.atomtype = self.sol_ions_type[a.name]['type']
        a.cgnr = a.id
        a.atomtypeB = None
        a.q = self.sol_ions_type[a.name]['charge']
        a.m = self.sol_ions_type[a.name]['mass']

        return a

    @staticmethod
    def swap_gro2top_ids(gro, top):
        tops = Atomselection(atoms=top)
        gros = Atomselection(atoms=gro)

        swap = tops.fetch_atoms(
            map(lambda x: x.id, gro), how='byid')

        add_swap = gros.fetch_atoms(
            map(lambda x: x.id, swap), how='byid', inv=True)

        if add_swap:
            swap.extend(add_swap)

        return swap

    def write_system(self, fname):
        res = OD(map(lambda x: (x.resnr, x.resname), self.system))
        syss = Atomselection(atoms=self.system)

        line = "# QM system\n"
        line += "# CHARGE: %d\n" % self.charge
        line += "{\n"
        for i in res.keys():
            atoms = syss.fetch_atoms(i, how='byresnr')
            aname = map(lambda x: x.name, atoms)
            line += "    # %d %s\n" % (i, res[i])
            line += "    %d: {'include': %s},\n\n" % (i, aname.__repr__())
        line += '}'
        with open(fname, 'w') as f:
            f.write(line)


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
    parser.add_argument('-r', '--refine',
                        action='store_true',
                        dest='refine',
                        help='Refine system: expand byres, terminal NH and CO')
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
