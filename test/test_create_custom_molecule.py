#!/usr/bin/env python3

import sys
sys.path.append('/home/yescalona/Development/VSOMM/vsomm_git/modeler')
from modeler import *

import argparse
from copy import deepcopy
from collections import Counter, OrderedDict

import itertools
from random import randint

def main(args):
    #mol = ['HA19', 'HA28', 'HA23', 'HA14D', 'HA12', 'HA21', 'HA11', 'HA32' ]
    mol = ['HA29', 'HA18', 'HA24D', 'HA31D', 'HA34D', 'HA35D']
    mol = mol*33 + ['HA29', 'HA18']
    args.BBsmol = 5 #len(mol) # number of BBs per mol
    
    system = System(args)
    system.set_molecule(mol)

    system.gen_topology()
    system.gen_coordinates()
    system.mol_assembly()
    system.fix_hydrogens()
    system.minimize()
    system.equilibrate()

    #system.statistics()

    system.combine_topologies()
    system.set_counterions()
    system.solvation()
    system.ionize()

    # md of the system
    system.minimize_system()
    system.equilibrate_system()
    system.production_system()

    #system.statistics()
    system.get_system_pdb()
    system.retrieve_system()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    class Range(object):
        "To require an argument to be a float between 0.0-1.0 using argparse"
        def __init__(self, start, end):
            self.start = start
            self.end = end
        def __eq__(self, other):
            return self.start <= other <= self.end
        def __repr__(self):
            return '{0}-{1}'.format(self.start, self.end)

    # General configuration
    parser.add_argument('--model', default="molecule", help="Model")
    #parser.add_argument('--type', default="advanced",  default=0.0, choices=["basic", "advanced", "random"], help="Type of input")
    parser.add_argument('--Waters', type=int, default=400, help="Number of water molecules", required=False)
    #parser.add_argument('--Density', type=float, default=0, help="Density")
    parser.add_argument('--ion', choices=["Ca2+", "Na+"], default="Na+", help="Ion", required=False)
    parser.add_argument('--pH', type=float, default=7.0, choices=[Range(1,14)], help="pH", required=False)
    parser.add_argument('--NumBBs', type=int, default=0, help="Numer of Building Blocks", required=False)
    parser.add_argument('--BBsmol', type=int, default=0, help="Building Blocks per molecule", required=False)

    # Elemental composition
    parser.add_argument('--CarbonFraction', type=float, default=0.0, choices=[Range(0.3, 1.0)], help="Carbon fraction", required=False)
    parser.add_argument('--HydrogenFraction', type=float, default=0.0, help="Hydrogen fraction", required=False)
    parser.add_argument('--OxygenFraction', type=float, default=0.0, help="Oxygen fraction", required=False)
    parser.add_argument('--NitrogenFraction', type=float, choices=[Range(0.0, 1.0)], help="Nitrogen fraction", required=False)
    parser.add_argument('--SulphurFraction', type=float, default=0.0, help="Sulphur fraction", required=False)
    
    
    subparser = parser.add_subparsers(dest="assignments", help='commands')
    orgcomp_IHSS = subparser.add_parser('IHSS', help='List contents')
    orgcomp_IHSS.add_argument('--CarbonylFraction', type=float, choices=[Range(0.0, 1.0)], help="Ketone/Quinone fraction", required=False)
    orgcomp_IHSS.add_argument('--CarboxylFraction', type=float, choices=[Range(0.0, 1.0)], help="Carboxyl fraction", required=False)
    orgcomp_IHSS.add_argument('--AromaticFraction', type=float, choices=[Range(0.0, 1.0)], help="Aromatic fraction", required=False)
    orgcomp_IHSS.add_argument('--AcetalFraction', type=float, choices=[Range(0.0, 1.0)], help="Acetal/aromatic fraction", required=False)
    orgcomp_IHSS.add_argument('--HeteroaliphaticFraction', type=float, choices=[Range(0.0, 1.0)], help="Heteroaliphatic fraction", required=False)
    orgcomp_IHSS.add_argument('--AliphaticFraction', type=float, choices=[Range(0.0, 1.0)], help="Aliphatic fraction", required=False)

    
    orgcomp_NMR = subparser.add_parser('NMR', help='List contents')
    orgcomp_NMR.add_argument('--CarbonylC_Fraction', type=float, choices=[Range(0.0, 1.0)], help="Carbonyl-C fraction", required=False)
    orgcomp_NMR.add_argument('--CarboxylC_Fraction', type=float, choices=[Range(0.0, 1.0)], help="Carboxyl-C fraction", required=False)
    orgcomp_NMR.add_argument('--O_ArylFraction', type=float, choices=[Range(0.0, 1.0)], help="O-Aryl fraction", required=False)
    orgcomp_NMR.add_argument('--ArylFraction', type=float, choices=[Range(0.0, 1.0)], help="Aryl fraction", required=False)
    orgcomp_NMR.add_argument('--Di_O_AlkylFraction', type=float, choices=[Range(0.0, 1.0)], help="Di-O-Alkyl fraction", required=False)
    orgcomp_NMR.add_argument('--O_AlkylFraction', type=float, choices=[Range(0.0, 1.0)], help="O-Alkyl fraction", required=False)
    orgcomp_NMR.add_argument('--MethoxylFraction', type=float, choices=[Range(0.0, 1.0)], help="Methoxyl fraction", required=False)
    orgcomp_NMR.add_argument('--Alkyl_C_Fraction', type=float, choices=[Range(0.0, 1.0)], help="Alkyl-C fraction", required=False)

    # Path configuration
    parser.add_argument('--seed', type=int, default=randint(0, 999999), help="seed")
    #parser.add_argument('--seed', type=int, default=123456, help="seed")
    parser.add_argument('--WORKDIR_PATH', default=tempfile.mkdtemp(), help="Working directory")
    parser.add_argument('--GROMOS_BIN_PATH', default="/home/yescalona/Development/VSOMM/vsomm_git/bin/", help="GROMOS binaries")
    parser.add_argument('--GROMOS_LIB_PATH', default="/home/yescalona/Development/VSOMM/vsomm_git/gromos_lib/", help="GROMOS libraries")
    parser.add_argument('--FORCEFIELD_PATH', default="/home/yescalona/Development/VSOMM/vsomm_git/forcefield/", help="GROMOS forcefield")

    args = parser.parse_args()

    print(args)
    main(args)
