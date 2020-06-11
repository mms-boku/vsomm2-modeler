#!/usr/bin/env python3

from vsomm_modeler.system import System

import argparse
from copy import deepcopy
from collections import Counter, OrderedDict

import itertools
from random import randint

def main():
    system = System(multiprocessing_enabled=True, forcefield="54a8", debug=True)


    building_blocks = list(system.building_blocks.keys())
    #building_blocks = [bb for bb in building_blocks if bb.startswith("HS16") or bb.startswith("HS12")  or bb.startswith("HS17")  ]


    print(building_blocks)
    building_blocks.sort()
    molecules = []
    l = []
    for i in building_blocks:
        l.append(system.building_blocks[i]['start_group'])
        x = [i] * 2
        l.extend(x)
        l.append(system.building_blocks[i]['end_group'])
        molecules.append(l)
        l = []


    for i in itertools.combinations(building_blocks, 2):
        #print(i)
        l.append(system.building_blocks[i[0]]['start_group'])
        x = list(i) * 2
        l.extend(x)
        l.append(system.building_blocks[i[1]]['end_group'])
        #print(l)
        molecules.append(l)
        l = []


    for idx, i in enumerate(molecules):
        print(idx, " ".join(i))


    system.molecules = molecules

    system.gen_topology()
    ##system.gen_coordinates()
    ##system.mol_assembly()
    ##system.fix_hydrogens()
    ##system.minimize()
    ##system.equilibrate()
    ###system.statistics_status()


if __name__ == "__main__":
    main()
