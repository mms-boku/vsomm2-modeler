#!/usr/bin/env python3

from vsomm_modeler.system import System
from vsomm_modeler.pdb import PDB

import argparse
from copy import deepcopy
from collections import Counter, OrderedDict

import itertools
from random import randint
import numpy as np
import random


class System2(System):
    def gen_coordinates(self):
        """Generate the PDB coordinates."""

        for m, mol in enumerate(self.molecules):  # no take in consideration the start_end_group
            molpdb = PDB()
            # PDB

            tmpvec = np.array([0., 0., 0.])
            for idx, building_block in enumerate(mol[1:-1]):

                HA = deepcopy(self.HAs[building_block])

                for i in range(len(HA.resseqs)):
                    HA.resseqs[i] = idx + 1

                # rotate
                angle = np.deg2rad(random.randint(0, 360))
                #angle = np.deg2rad(20 * idx)
                rotmat = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])

                HA.coords = HA.coords - HA.coords[0]

                tmp_coords = []
                for coord in HA.coords:
                    tmp_coords.append(np.matmul(rotmat, coord))
                HA.coords = np.array(tmp_coords)

                #HA.coords = HA.coords - HA.coords[0] + tmpvec + np.array([5, 0, 0])
                HA.coords = HA.coords + tmpvec + np.array([5, 0, 0])

                molpdb += HA

                tmpvec = molpdb.coords[-1]

            with open("%s/mol_%s.pdb" % (self.workdir, m), 'w') as f:
                f.writelines(molpdb.output())

            # G96
            input_pdb2g96 = {}
            input_pdb2g96['pdb2g96'] = self.gromospp_bin_dir + "pdb2g96"
            input_pdb2g96['topo'] = self.workdir + "/mol_" + str(m) + ".top"
            input_pdb2g96['pdb'] = "%s/mol_%s.pdb" % (self.workdir, m)
            input_pdb2g96['lib'] = self.vsomm_building_blocks_dir + self.forcefield + "/gromos_lib/pdb2g96.lib"
            input_pdb2g96['output'] = "%s/pdb2g96_mol_%s.cnf" % (self.workdir, m)

            command = "{pdb2g96} @topo {topo} @pdb {pdb} @lib {lib} > {output}".format(**input_pdb2g96)
            self.run_command(command)

            if not self.debug:
                os.remove(input_pdb2g96["pdb"])


    def mol_assembly(self):
        """Molecular assembly in accordance to the bond a angle information of building blocks."""

        for m, mol in enumerate(self.molecules):  # no take in consideration the start_end_group
            # calcular numero de atomos en current_state
            list_gca = []
            atomnum = 0

            for idx, x in enumerate(mol[:-1]):
                if idx == 0:
                    atomnum += self.start_end_groups[x]['atom_count']
                elif idx == 1:
                    atomnum += self.building_blocks[x]['atom_count']
                else:
                    i, j, k = atomnum, atomnum + 1, atomnum + 2
                    bond_angle = self.building_blocks[x]['bond_angle']
                    bond_length = self.building_blocks[x]['bond_length']
                    list_gca.append("d%1:{i},{j}%{length} ".format(i=i, j=j, length=0.5))
                    #list_gca.append("d%1:{i},{j}%{length} ".format(i=i, j=j, length=bond_length * 4 / 10))
                    #list_gca.append("a%1:{i},{j},{k}%{angle} ".format(i=i, j=j, k=k, angle=bond_angle))
                    atomnum += self.building_blocks[x]['atom_count']

            # GCA
            input_gca = {}
            input_gca['gca'] = self.gromospp_bin_dir + "gca"
            input_gca['topo'] = self.workdir + "/mol_" + str(m) + ".top"
            input_gca['prop'] = " ".join(list_gca)
            input_gca['traj'] = "%s/pdb2g96_mol_%s.cnf" % (self.workdir, m)
            input_gca['output'] = "%s/gca_mol_%s.cnf" % (self.workdir, m)

            command = "{gca} @topo {topo} @prop {prop} @traj {traj} @outformat cnf @pbc v > {output}".format(
                **input_gca)
            self.run_command(command)

            if not self.debug:
                os.remove(input_gca["traj"])


def main():

    system = System(multiprocessing_enabled=True, forcefield="54a8", debug=True)


    building_blocks = list(system.building_blocks.keys())
    #building_blocks = [bb for bb in building_blocks if bb.startswith("HS37")]

    building_blocks = np.array(building_blocks)


    molecules = []
    for i in range(10000):

        idxs = np.random.randint(low=0, high=len(building_blocks), size=2)
        l = list(building_blocks[idxs])

        start = system.building_blocks[l[0]]['start_group']
        l.insert(0, start)

        end = system.building_blocks[l[-1]]['end_group']
        l.append(end)
        
        molecules.append(l)


    for idx, i in enumerate(molecules):
        print(idx, " ".join(i))


    system.molecules = molecules

    system.gen_topology()
    system.gen_coordinates()
    #system.mol_assembly()
    system.fix_hydrogens()
    system.minimize()
    system.equilibrate()
    #system.statistics_status()


if __name__ == "__main__":
    main()
