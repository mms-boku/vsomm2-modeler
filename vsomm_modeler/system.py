import os
import sys
import subprocess
import json
import billiard as multiprocessing
import random
import tempfile
import numpy as np
import collections

from copy import deepcopy
from collections import Counter, OrderedDict

from vsomm_modeler.pdb import PDB, get_mol_coordinates
from vsomm_modeler.state import State_IHSS, State_SSNMR, PKA_COOH
from vsomm_modeler import imds


class System():
    """
    System class

    """

    def __init__(self, model=None, number_of_building_blocks=200,
                 building_blocks_per_molecule=5, water_molecules=None,
                 counterion="Na+", ionrandom=False,
                 pH=7.0, logK1=4.2, initial_density=900, boxsize=None,
                 seed=random.randint(0, 999999), seed_gen_algorithm=None,
                 assignment="IHSS",
                 multiprocessing_enabled=False, num_threads=8,
                 vsomm_building_blocks_dir=os.environ.get("VSOMM_BUILDING_BLOCKS"),
                 forcefield="54a7",
                 gromosxx_bin_dir=os.environ.get("GROMOSXX_BIN"),
                 gromospp_bin_dir=os.environ.get("GROMOSPP_BIN"),
                 gromacs_bin_dir=os.environ.get("GROMACS_BIN"),
                 workdir=tempfile.mkdtemp(), debug=False,
                 run_equilibration=False, togromacs=False,
                 **kwargs):
        """ Initialization of the system."""

        self.model = model
        self.number_of_building_blocks = number_of_building_blocks
        self.building_blocks_per_molecule = building_blocks_per_molecule
        self.water_molecules = water_molecules
        self.kwargs = kwargs
        self.counterion = counterion
        self.counterions = 0
        self.ionrandom = ionrandom

        # include these properties to the input_state
        self.kwargs['pH'] = pH
        self.deprot_prot_fraction = 10**(pH - logK1)
        self.kwargs['deprot_prot_fraction'] = self.deprot_prot_fraction

        self.initial_density = initial_density
        self.boxsize = boxsize
        self.seed = seed
        self.seed_gen_algorithm = seed_gen_algorithm
        self.multiprocessing_enabled = multiprocessing_enabled
        self.num_threads = num_threads
        self.assignment = assignment
        self.vsomm_building_blocks_dir = vsomm_building_blocks_dir
        self.forcefield = forcefield
        self.gromosxx_bin_dir = gromosxx_bin_dir
        self.gromospp_bin_dir = gromospp_bin_dir
        self.gromacs_bin_dir = gromacs_bin_dir
        self.workdir = workdir
        self.debug = debug
        self.run_equilibration = run_equilibration
        self.togromacs = togromacs

        self.states = {
            'IHSS': State_IHSS,
            'SSNMR': State_SSNMR
        }
        # Definition of the input state of the system
        self.input_state = self.states[self.assignment]()
        self.set_input_state()

        # Definition of the current state of the system
        self.current_state = self.states[self.assignment]()
        self.building_blocks_start_end_groups = []

        # Get information from the json databases
        self.building_blocks = json.load(open(self.vsomm_building_blocks_dir + self.forcefield + "/properties/Building_Blocks.json", 'r'), object_pairs_hook=OrderedDict)
        self.start_end_groups = json.load(open(self.vsomm_building_blocks_dir + self.forcefield + "/properties/Start_End_Groups.json", 'r'), object_pairs_hook=OrderedDict)

        # Creating a subset of building blocks without protonated states
        # also a dictionary to see the relation between deprotonated / protonated states
        self.protonated_state = {}
        self.subset_building_blocks = collections.OrderedDict()
        for building_block, properties in self.building_blocks.items():
            if "p" in building_block:
                self.protonated_state[building_block[:-1]] = building_block
            else:
                self.subset_building_blocks[building_block] = properties

        # Get coordinates of the building blocks
        self.HAs = get_mol_coordinates(self.vsomm_building_blocks_dir + self.forcefield + "/coordinates/HSs.pdb")

    def set_input_state(self):
        """Set input state to compare with current state."""

        for frac in self.input_state.fractions:
            try:
                setattr(self.input_state, frac, self.kwargs[frac])
            except:
                print("WARNING setting input fraction " + frac + ": 0.0", file=sys.stderr)
                setattr(self.input_state, frac, 0.0)

        return self.input_state

    def set_molecule(self, molecule):
        """Update the list of properties of the molecule."""

        for building_block in molecule:
            self.current_state.append(self.building_blocks[building_block])

        self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)
        self.gen_molecules()

    def eval_func_ga(self, chromosome):
        """Evaluation of the individual chromosome (list of building blocks)."""

        mol = self.states[self.assignment]()
        for i, c in enumerate(chromosome):
            # TODO: es necesario que esta ordenado
            building_block = list(self.subset_building_blocks)[c]

            if i % self.building_blocks_per_molecule == 0:
                start_group = self.subset_building_blocks[building_block]['start_group']
                mol.append(self.start_end_groups[start_group])

            mol.append(self.subset_building_blocks[building_block])

            if (i + 1) % self.building_blocks_per_molecule == 0:
                end_group = self.subset_building_blocks[building_block]['end_group']
                mol.append(self.start_end_groups[end_group])

        # mol.set_list(chromosome_building_blocks)
        mol.set_status()
        # print(mol.list_building_blocks)
        score = self.diff_status_ga(mol)
        return score

    def diff_status_ga(self, temp_state):
        """Calculate the difference between the current and input state for the GA"""

        diff = []

        for frac in temp_state.fractions:
            # TODO: consider these fractions
            if frac is "hydrogen_fraction" or frac is "oxygen_fraction" or frac is "sulphur_fraction":
                continue
            #elif frac is "carbon_fraction" or frac is "nitrogen_fraction":
            #    continue
            elif frac is "deprot_prot_fraction":
                continue
            #    #print(vars(temp_state))
            #    if self.input_state.deprot_prot_fraction > temp_state.carboxyl_c_count:
            #        if temp_state.protonated >= 1:
            #            diff.append(1)
            #    else:
            #         diff.append(0.2 * ((temp_state.deprot_prot_fraction - self.input_state.deprot_prot_fraction) / self.input_state.deprot_prot_fraction)**2)
            else:
                try:
                    if frac is "carbon_fraction" or frac is "nitrogen_fraction":
                        # permit a 20% percentage error 
                        # (0.05*0.05)/(0.2*0.2) = 0.0625
                        diff.append(0.0625 * ((np.abs(getattr(temp_state, frac) - getattr(self.input_state, frac)) / getattr(self.input_state, frac))**2))
                    else:
                        diff.append((np.abs(getattr(temp_state, frac) - getattr(self.input_state, frac)) / getattr(self.input_state, frac))**2)
                except:
                    diff.append(1)

        diff = np.array(diff)
        diff = diff.max()
        return diff

    def gen_system_ga(self):
        """Selection of the building blocks of the system using a Genetic Algorithm."""

        from pyevolve import G1DList, GSimpleGA, Selectors, Statistics, Scaling
        from pyevolve import Initializators, Mutators, Consts, DBAdapters, Crossovers

        ngen = 1000
        condition = True

        # in order to mantain the seed for the creation of molecules   
        if self.seed_gen_algorithm:
            seed = self.seed_gen_algorithm
        else:
            seed = self.seed

        attempt = 0

        #while condition and attempt < 1:
        genome = G1DList.G1DList(self.number_of_building_blocks)
        genome.evaluator.set(self.eval_func_ga)
        genome.setParams(rangemin=0, rangemax=len(self.subset_building_blocks) - 1,
                         bestrawscore=0.0025)

        ga = GSimpleGA.GSimpleGA(genome, seed=seed)
        ga.setMinimax(Consts.minimaxType["minimize"])
        ga.terminationCriteria.set(GSimpleGA.RawScoreCriteria)
        ga.setGenerations(ngen)
        ga.setPopulationSize(500)

        # TODO: make GA billiard compatible
        #if self.multiprocessing_enabled:
        #    ga.setMultiProcessing(True, max_processes=self.num_threads)

        ga.evolve(freq_stats=50)

        #    if ga.getCurrentGeneration() < ngen:
        #        condition = False
        #    else:
        #        seed += 1
        #        attempt += 1

        #if condition is False:
        #    raise RuntimeError("GA function is not able to satisfy the experimental input")
        #    sys.exit(1)

        print(ga.bestIndividual(), file=sys.stderr)

        lista = ga.bestIndividual()
        lista = lista.getInternalList()

        best = [self.subset_building_blocks[list(self.subset_building_blocks)[i]] for i in lista]
        for i in best:
            self.current_state.append(i)

        self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)

        #self.gen_molecules()

    def eval_func_deap(self, chromosome):
        """Evaluation of the individual chromosome (list of building blocks)."""

        mol = self.states[self.assignment]()
        for i, c in enumerate(chromosome):
            building_block = list(self.building_blocks)[c]

            if i % self.building_blocks_per_molecule == 0:
                start_group = self.building_blocks[building_block]['start_group']
                mol.append(self.start_end_groups[start_group])

            mol.append(self.building_blocks[building_block])

            if (i + 1) % self.building_blocks_per_molecule == 0:
                end_group = self.building_blocks[building_block]['end_group']
                mol.append(self.start_end_groups[end_group])

        # mol.set_list(chromosome_building_blocks)
        mol.set_status()
        # print(mol.list_building_blocks)
        score = self.diff_status_deap(mol)
        return (score,)

    def diff_status_deap(self, temp_state):
        """Calculate the difference between the current and input state for the GA"""

        diff = []

        for frac in temp_state.fractions:
            if frac is "hydrogen_fraction" or frac is "oxygen_fraction" or frac is "sulphur_fraction":
                continue
            if frac in ["carbon_fraction", "nitrogen_fraction"]:
                continue
            elif frac is "pH":
                frac_deprot_input_state = (10**(-PKA_COOH)) / (10**(-1 * getattr(self.input_state, frac)))
                frac_deprot_temp_state = (10**(-PKA_COOH)) / (10**(-1 * getattr(temp_state, frac)))
                diff.append(((frac_deprot_temp_state - frac_deprot_input_state) / frac_deprot_input_state)**2)
            else:
                try:
                    diff.append((np.abs(getattr(temp_state, frac) - getattr(self.input_state, frac)) /
                                 getattr(self.input_state, frac))**2)
                except:
                    diff.append(1)

        diff = np.array(diff)
        diff = diff.max()
        return diff

    def gen_system_deap(self):
        """Selection of the building blocks of the system using a Genetic Algorithm."""

        from deap import algorithms, base, creator, tools

        random.seed(self.seed)

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)

        toolbox = base.Toolbox()

        if self.multiprocessing_enabled:
            pool = multiprocessing.Pool(processes=self.num_threads)
            toolbox.register("map", pool.map)

        toolbox.register("attr_int", random.randint, 0, len(self.building_blocks) - 1)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, n=self.number_of_building_blocks)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", self.eval_func_deap)
        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutUniformInt, low=0, up=len(self.building_blocks) - 1, indpb=(1/len(self.building_blocks)))

        #toolbox.register("select", tools.selBest, k=2)
        #toolbox.register("select", tools.selRandom, k=300)
        toolbox.register("select", tools.selTournament, tournsize=3)
        #toolbox.register("evaluate", benchmarks.rastrigin)

        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        #stats.register("Min", min)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)

        #population = toolbox.population(n=400)
        #checkpoint = tools.Checkpoint(population=population)

        gen, ngen, cxpb, mutpb = 0, 1000, 0.9, 0.02

        # create an initial population of 300 individuals (where
        # each individual is a list of integers)
        population = toolbox.population(n=100)

        print("Start of evolution")

        # Evaluate the entire population
        #fitnesses = list(map(toolbox.evaluate, population))
        #for ind, fit in zip(population, fitnesses):
        #    ind.fitness.values = fit

        #print("  Evaluated %i individuals" % len(population))

        # Extracting all the fitnesses of
        #fits = [ind.fitness.values[0] for ind in population]

        # Variable keeping track of the number of generations
        g = 0

        bestIndScore = 1000
        # Begin the evolution
        while bestIndScore > 0.1 and g < ngen:
            # A new generation
            g += 1
            #print("-- Generation %i --" % g)

            # select the next generation of individuals
            #offspring = toolbox.select(population, k=20)
            #offspring += toolbox.population(n=80)
            offspring = toolbox.select(population, k=len(population))

            # Apply crossover and mutation on the offspring
            offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            #print("  Evaluated %i individuals" % len(invalid_ind))

            # The population is entirely replaced by the offspring
            population[:] = offspring

            #fits = [ind.fitness.values[0] for ind in population]
            #print("fits", fits)
            record = stats.compile(population)
            bestIndScore = record['min']

            if g % 20 == 0:
                #fitnesses = toolbox.map(toolbox.evaluate, population)

                print("Gen. %s: Min/Max/Avg/Std %.2f/%.2f/%.2f/%.2f" %
                      (g, record['min'], record['max'], record['avg'], record['std']))

        print("Gen. %s: Min/Max/Avg/Std %.2f/%.2f/%.2f/%.2f" %
              (g, record['min'], record['max'], record['avg'], record['std']))

        lista = tools.selBest(population, k=1)[0]

        best = [self.building_blocks[list(self.building_blocks)[i]] for i in lista]
        for i in best:
            self.current_state.append(i)

        self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)

        #self.gen_molecules()

    def diff_status(self, temp_state):
        """Calculate the difference between the current and input state for the simple algorithm """

        diff = 0.

        for frac in temp_state.fractions:
            if frac is "hydrogen_fraction" or frac is "oxygen_fraction" or frac is "sulphur_fraction":
                continue
            elif frac is "pH":
                frac_deprot_input_state = (10**(-4.7)) / (10**(-1 * getattr(self.input_state, frac)))
                frac_deprot_temp_state = (10**(-4.7)) / (10**(-1 * getattr(temp_state, frac)))
                diff += ((frac_deprot_temp_state - frac_deprot_input_state) / frac_deprot_input_state)**2 * 0.1
            #elif frac is "carbon_fraction" or frac is "nitrogen_fraction":
            #    continue
            #    diff += (0.01 * ((np.abs(getattr(temp_state, frac) -
            #                               getattr(self.input_state, frac)) /
            #                        getattr(self.input_state, frac)) ** 2))
            else:
                diff += ((np.abs(getattr(temp_state, frac) - getattr(self.input_state, frac))) /
                         getattr(self.input_state, frac))**2

        return diff

    def gen_system(self):
        """"Selection of the building blocks of the system using a simple algorithm."""

        new_building_block = {}
        for n in range(self.number_of_building_blocks):
            min_diff = 1000
            for building_block in self.subset_building_blocks.keys():
                # temporal state = current state with a new building block
                temp_state = deepcopy(self.current_state)
                temp_state.append(self.building_blocks[building_block])
                temp_state.set_status()
                temp_diff = self.diff_status(temp_state)

                # if reduce de difference between the input and current state ...
                if temp_diff < min_diff:
                    new_building_block = building_block
                    min_diff = temp_diff

            # add new building block
            self.current_state.append(self.building_blocks[new_building_block])
            self.current_state.set_status()

            # TODO: consider start and end groups at this level

        self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)
        self.gen_molecules()

    def gen_molecules(self):
        """Generate list of molecules."""

        # create a list of molecules with specific building_blocks
        self.molecules = [self.current_state.list_building_blocks[i:i + self.building_blocks_per_molecule]
                          for i in range(0, len(self.current_state.list_building_blocks), self.building_blocks_per_molecule)]

        # add the start and end group for each molecule and update current state
        for mol in self.molecules:
            start_group = self.building_blocks[mol[0]]['start_group']
            mol.insert(0, start_group)
            self.current_state.append(self.start_end_groups[start_group])

            end_group = self.building_blocks[mol[-1]]['end_group']
            mol.append(end_group)
            self.current_state.append(self.start_end_groups[end_group])

        self.current_state.set_status()

        return self.current_state

    def gen_topology(self):
        """Generate the topology of the building blocks"""

        for m, mol in enumerate(self.molecules):

            input_make_top = {}
            input_make_top['make_top'] = self.gromospp_bin_dir + "make_top"
            input_make_top['mtb'] = self.vsomm_building_blocks_dir + self.forcefield + "/forcefield/" + self.forcefield + "_vsomm.mtb"
            input_make_top['ifp'] = self.vsomm_building_blocks_dir + self.forcefield + "/forcefield/" + self.forcefield + ".ifp"
            input_make_top['seq'] = " ".join(mol)
            input_make_top['output'] = self.workdir + \
                "/mol_" + str(m) + ".top"

            command = "{make_top} @build {mtb} @param {ifp} @seq {seq} @solv H2O > {output}".format(**input_make_top)
            self.run_command(command)

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
                #angle = np.deg2rad(20 * idx)
                angle = np.deg2rad(random.randint(0, 360))
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
                self.remove_files([input_pdb2g96["pdb"]])


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
                    list_gca.append("d%1:{i},{j}%{length} ".format(i=i, j=j, length=bond_length / 10))
                    list_gca.append("a%1:{i},{j},{k}%{angle} ".format(i=i, j=j, k=k, angle=bond_angle))
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
                self.remove_files([input_gca["traj"]])

    def fix_hydrogens(self):
        """Add or fix the position of hydrogens in the molecules."""

        for m, mol in enumerate(self.molecules):  # no take in consideration the start_end_group

            # GCH
            input_gch = {}
            input_gch['gch'] = self.gromospp_bin_dir + "gch"
            input_gch['topo'] = self.workdir + "/mol_" + str(m) + ".top"

            #input_gch['cnf'] = "%s/pdb2g96_mol_%s.cnf" % (self.workdir, m)
            input_gch['cnf'] = "%s/pdb2g96_mol_%s.cnf" % (self.workdir, m)

            input_gch['output'] = "%s/mol_%s.cnf" % (self.workdir, m)

            command = "{gch} @topo {topo} @pos {cnf} @tol 0.1 @pbc v > {output}".format(**input_gch)
            self.run_command(command)

            if not self.debug:
                self.remove_files([input_gch['cnf']])

    def minimize(self):
        """Minimize each molecule."""



        ###############################################
        jobs = [
            ["",  "min_", 50000, 0, 0],
            ["min_", "min2_", 50000, 1, 0],
            ["min2_", "min3_", 50000, 1, 3]
        ]

        to_clean = []
        for prefix_in, prefix_out, steps, charge, shake in jobs:

            if self.multiprocessing_enabled:
                pool = multiprocessing.Pool(processes=self.num_threads)

            for m, mol in enumerate(self.molecules):

                # calcular numero de atomos en current_state
                atomnum = 0
                atomnum += self.start_end_groups[mol[0]]['atom_count']
                atomnum += self.start_end_groups[mol[-1]]['atom_count']
                for building_block in mol[1:-1]:
                    atomnum += self.building_blocks[building_block]['atom_count']

                with open(self.workdir + "/" + prefix_out + "mol_" + str(m) + ".imd", "w") as imdfile:
                    imdfile.write(imds.minimize_molecule.format(seed=self.seed,
                                                                atomnum=atomnum,
                                                                charge=charge,
                                                                shake=shake,
                                                                steps=steps))

                input_min = {}

                input_min['md'] = self.gromosxx_bin_dir + "/md"
                input_min['topo'] = self.workdir + "/mol_" + str(m) + ".top"
                input_min['conf'] = self.workdir + "/" + prefix_in + "mol_" + str(m) + ".cnf"
                input_min['fin'] = self.workdir + "/" + prefix_out + "mol_" + str(m) + ".cnf"
                input_min['input'] = self.workdir + "/" + prefix_out + "mol_" + str(m) + ".imd"
                input_min['output'] = self.workdir + "/" + prefix_out + "mol_" + str(m) + ".omd"
                command = "export OMP_NUM_THREADS=1 && {md} @topo {topo} @conf {conf} @fin {fin} @input {input} > {output}".format(**input_min)

                # saving the trajectory and energy of the system
                # input_min['trc'] = self.workdir + "/" + prefix_out+"mol_" + str(m) + ".trc"
                # input_min['tre'] = self.workdir + "/" + prefix_out+"mol_" + str(m) + ".tre"
                # command = "{md} @topo {topo} @conf {conf} @fin {fin} @input {input} @trc {trc} @tre {tre}  > {output}".format(**input_min)

                if self.multiprocessing_enabled:
                    pool.apply_async(self.run_command, args=(command,))
                else:
                    self.run_command(command)


                if not self.debug:
                        to_clean.extend([input_min['conf'], input_min['input'], input_min['output']])

            if self.multiprocessing_enabled:
                pool.close()
                pool.join()


            if not self.debug:
                self.remove_files([c for c in to_clean])
                to_clean = []



    def equilibrate(self):
        """Equilibration of each molecule."""

        to_clean = []
        if self.multiprocessing_enabled:
            pool = multiprocessing.Pool(processes=self.num_threads)

        for m, mol in enumerate(self.molecules):



            # calcular numero de atomos en current_state
            atomnum = 0
            atomnum += self.start_end_groups[mol[0]]['atom_count']
            atomnum += self.start_end_groups[mol[-1]]['atom_count']
            for building_block in mol[1:-1]:
                atomnum += self.building_blocks[building_block]['atom_count']

            with open(self.workdir + "/eq_mol_" + str(m) + ".imd", "w") as imdfile:
                imdfile.write(imds.equilibrate_molecule.format(seed=self.seed,
                                                               atomnum=atomnum))

            input_eq = {}

            input_eq['md'] = self.gromosxx_bin_dir + "/md"
            input_eq['topo'] = self.workdir + "/mol_" + str(m) + ".top"
            input_eq['conf'] = self.workdir + "/min3_mol_" + str(m) + ".cnf"
            input_eq['fin'] = self.workdir + "/eq_mol_" + str(m) + ".cnf"
            input_eq['input'] = self.workdir + "/eq_mol_" + str(m) + ".imd"
            input_eq['output'] = self.workdir + "/eq_mol_" + str(m) + ".omd"
            command = "export OMP_NUM_THREADS=1 && {md} @topo {topo} @conf {conf} @fin {fin} " \
                      "@input {input} > {output}".format(**input_eq)

            # saving the trajectory and energy of the system
            # input_eq['trc'] = self.workdir + "/eq_mol_" + str(m) + ".trc"
            # input_eq['tre'] = self.workdir + "/eq_mol_" + str(m) + ".tre"
            # command = "{md} @topo {topo} @conf {conf} @fin {fin} @input {input} \
            #        @trc {trc} @tre {tre}  > {output}".format(**input_eq)

            if self.multiprocessing_enabled:
                pool.apply_async(self.run_command, args=(command,))
            else:
                self.run_command(command)

            if not self.debug:
                to_clean.extend([input_eq['conf'], input_eq['input'], input_eq['output']])

        if self.multiprocessing_enabled:
            pool.close()
            pool.join()

        if not self.debug:
            self.remove_files([c for c in to_clean])
            to_clean = []




    def combine_topologies(self):
        """Combine the topologies of the different molecules."""

        input_com_top = {}
        input_com_top['com_top'] = self.gromospp_bin_dir + "/com_top"
        input_com_top['topo'] = " ".join([self.workdir + "/mol_" + str(m) +
                                          ".top" for m in range(len(self.molecules))])
        input_com_top['output'] = self.workdir + "/system.top"

        command = "{com_top} @topo {topo} @param 1 @solv 1 > {output}".format(**input_com_top)
        self.run_command(command)

        if not self.debug:
            self.remove_files([self.workdir + "/mol_" + str(m) + ".top" for m in range(len(self.molecules))])

    def set_counterions(self):
        valence = {'Na+': 1, 'Ca2+': 2}
        topofiles = {'Na+': "forcefield/na.top", 'Ca2+': "forcefield/ca.top"}

        # neutralize with a odd number of charge
        # protonate the first deprotonated building block
        #print(self.current_state.list_building_blocks)
        if self.counterion == "Ca2+" and self.current_state.charge % 2 != 0:
            #print("DENTRO DEL IF")
            for position, building_block in enumerate(self.current_state.list_building_blocks):
                # protonated state, and charge is -1
                #print(position, building_block, self.building_blocks[building_block]['charge'])
                if building_block in self.protonated_state and self.building_blocks[building_block]['charge'] == -1:
                    #print("ENCONTRO UN BUILDING BLOCK")
                    new_building_block = self.protonated_state[building_block]
                    self.current_state.modify(self.building_blocks[building_block], position, self.building_blocks[new_building_block])
                    self.current_state.set_status()
                    self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)
                    #print(self.current_state.list_building_blocks)
                    break

        self.topo_counterion = topofiles[self.counterion]
        self.counterions = int(-1 * self.current_state.charge / valence[self.counterion])

    def set_pH(self):
        for position, building_block in enumerate(self.current_state.list_building_blocks):
            if building_block in self.protonated_state and self.building_blocks[building_block]['charge'] == -1:
                # throw a dice and look if you get 1
                if 1 == random.randint(1, int(self.deprot_prot_fraction)):
                    new_building_block = self.protonated_state[building_block]
                    self.current_state.modify(self.building_blocks[building_block], position, self.building_blocks[new_building_block])
                    self.current_state.set_status()
                    self.current_state.list_building_blocks_to_count = deepcopy(self.current_state.list_building_blocks)


    def solvation_boxsize(self):
        """Solvation of the system."""

        # TODO: use ran_solvation if water_molecules is None

        # self.counterions = int(-1 * self.current_state.charge / 2) ### ions will replace this water molecules

        input_ran_box = {}
        input_ran_box['ran_box'] = self.gromospp_bin_dir + "/ran_box"
        input_ran_box['topo'] = " ".join([self.workdir + "/mol_" + str(m) + ".top" for m in range(len(self.molecules))])
        input_ran_box['pos'] = " ".join([self.workdir + "/eq_mol_" + str(m) + ".cnf" for m in range(len(self.molecules))])
        input_ran_box['nsm'] = " ".join(["1"] * len(self.molecules))
        input_ran_box['boxsize'] = " ".join([str(b) for b in self.boxsize])
        input_ran_box['seed'] = self.seed
        input_ran_box['output'] = self.workdir + "/system_solute.cnf"

        command = "{ran_box} @topo {topo} @pos {pos} @nsm {nsm} @boxsize {boxsize} @seed {seed} @pbc r> {output}".format(**input_ran_box)
        self.run_command(command)

        if not self.debug:
            self.remove_files([self.workdir + "/eq_mol_" + str(m) + ".cnf" for m in range(len(self.molecules))])

        input_sim_box = {}
        input_sim_box['sim_box'] = self.gromospp_bin_dir + "/sim_box"
        input_sim_box['topo'] = self.workdir + "/system.top"
        input_sim_box['pos'] = self.workdir + "/system_solute.cnf"
        input_sim_box['solvent'] = self.vsomm_building_blocks_dir + self.forcefield + "/coordinates/spc.cnf"
        input_sim_box['boxsize'] = " ".join([str(b) for b in self.boxsize])
        input_sim_box['output'] = self.workdir + "/system.cnf"

        command = "{sim_box} @topo {topo} @pbc r @pos {pos} @solvent {solvent} @boxsize {boxsize} > {output}".format(**input_sim_box)
        self.run_command(command)


        with open(self.workdir + "/system.cnf") as file:
            for line in file:
                if line.startswith("Added"):
                    self.water_molecules = int(line.split()[1]) - self.counterions
                    break

        # TODO: for proteins, check system charge using check_top


    def solvation(self):
        """Solvation of the system."""

        # TODO: use ran_solvation if water_molecules is None

        # self.counterions = int(-1 * self.current_state.charge / 2) ### ions will replace this water molecules

        input_ran_box = {}
        input_ran_box['ran_box'] = self.gromospp_bin_dir + "/ran_box"
        input_ran_box['topo'] = " ".join([self.workdir + "/mol_" + str(m) +
                                          ".top" for m in range(len(self.molecules))] + [self.vsomm_building_blocks_dir + self.forcefield + "/forcefield/spc.top"])
        input_ran_box['pos'] = " ".join([self.workdir + "/eq_mol_" + str(m) +
                                         ".cnf" for m in range(len(self.molecules))] + [self.vsomm_building_blocks_dir + self.forcefield + "/forcefield/spc.dat"])

        input_ran_box['nsm'] = " ".join(["1"] * len(self.molecules) + [str(self.water_molecules + self.counterions)])
        ##input_ran_box['nsm'] = " ".join(["1"] * len(self.molecules) + [str(self.counterions)])

        input_ran_box['dens'] = self.initial_density
        input_ran_box['seed'] = self.seed
        input_ran_box['output'] = self.workdir + "/system.cnf"

        command = "{ran_box} @topo {topo} @pos {pos} @nsm {nsm} @dens {dens} @seed {seed} @pbc r> {output}".format(**input_ran_box)
        self.run_command(command)

        if not self.debug:
            self.remove_files([self.workdir + "/eq_mol_" + str(m) + ".cnf" for m in range(len(self.molecules))])


    def ionize(self):
        """Ionization of the system"""

        #self.counterions = int(-1 * self.current_state.charge / 2)

        input_com_top = {}
        input_com_top['com_top'] = self.gromospp_bin_dir + "/com_top"
        input_com_top['topo'] = self.workdir + "/system.top" + " " + \
            str(self.counterions) + ":" + self.vsomm_building_blocks_dir + self.forcefield + "/" + self.topo_counterion
        input_com_top['output'] = self.workdir + "/system_ions.top"

        command = "{com_top} @topo {topo} @param 1 @solv 1 > {output}".format(**input_com_top)
        self.run_command(command)

        input_ion = {}
        input_ion['ion'] = self.gromospp_bin_dir + "/ion"
        input_ion['topo'] = self.workdir + "/system.top"
        input_ion['pos'] = self.workdir + "/system.cnf"
        input_ion['charge'] = self.counterions
        input_ion['chargename'] = self.counterion
        input_ion['seed'] = self.seed
        input_ion['output'] = self.workdir + "/system_ions.cnf"

        if self.ionrandom:
            command = "{ion} @topo {topo} @pos {pos} @pbc r @positive {charge} {chargename} \
                       @random {seed} > {output}".format(**input_ion)
        else:
            command = "{ion} @topo {topo} @pos {pos} @pbc r @positive {charge} {chargename} \
                       @potential 1.4 @mindist 0.25 > {output}".format(**input_ion)

        self.run_command(command)

        if not self.debug:
            self.remove_files([[input_ion['topo'], input_ion['pos']])

    def get_system_pdb(cnf_filename, pdb_filename, self):
        """Obtain PDB of the system."""

        input_frameout = {}
        input_frameout['frameout'] = self.gromospp_bin_dir + "/frameout"
        input_frameout['topo'] = self.workdir + "/system_ions.top"
        input_frameout['traj'] = self.workdir + cnf_filename
        # input_frameout['output'] = self.workdir + "/system_ions.pdb"

        # fixme: FRAMEOUT DOESNT HAVE AN OUTPUT FLAG :@
        command = "cd {workdir} && {frameout} @topo {topo} @traj {traj} @pbc r cog @outformat pdb \
                @include ALL @time 0 0.5".format(**input_frameout, workdir=self.workdir)
        self.run_command(command)
        command = "mv {workdir}/FRAME_00001.pdb {workdir}/{pdb_filename}".format(workdir=self.workdir, pdb_filename=pdb_filename)
        self.run_command(command)

    def statistics(self):
        """Statistics of the system."""

        dict_args = vars(self)
        dict_input_state = vars(self.input_state)
        dict_current_state = vars(self.current_state)

        # Calculate the percentage error
        dict_diff = {}
        for frac in self.input_state.fractions:
            if frac is "hydrogen_fraction" or frac is "oxygen_fraction" or frac is "sulphur_fraction":
                dict_diff["diff_" + frac] = 0
                continue
            #elif frac is "pH":
            #    frac_deprot_input_state = ((10**(-PKA_COOH)) /
            #                               (10**(-1 * dict_input_state[frac])))
            #    dict_diff["input_frac_deprot"] = frac_deprot_input_state
            #    frac_deprot_temp_state = ((10**(-PKA_COOH)) /
            #                              (10**(-1 * dict_current_state[frac])))
            #    dict_diff["frac_deprot"] = frac_deprot_temp_state
            #    dict_diff["diff_frac_deprot"] = (np.abs(frac_deprot_temp_state -
            #                                     frac_deprot_input_state) /
            #                                     frac_deprot_input_state)

            dict_diff["diff_" + frac] = (np.abs(dict_current_state[frac] -
                                         dict_input_state[frac]) /
                                         dict_input_state[frac])

        # Merging all the variables in one dictionary
        dict_input_state = dict(("input_" + key, value) for (key, value) in dict_input_state.items())
        dict_all = {**dict_args, **dict_input_state, **dict_current_state, **dict_diff}

        # Count the number of building blocks in the system
        dict_all['counter'] = "\n".join(["{:>5}: {:5}".format(key, value)
                                         for key, value in OrderedDict(
                Counter(self.current_state.list_building_blocks_to_count).most_common()).items()])

        counter = Counter(self.current_state.list_building_blocks_to_count).most_common()
        #pH              {input_pH:6.3f} {pH:6.3f}

        fmt_SSNMR = """
Model:                     {model}
Random Seed:               {seed}
Number of building blocks: {number_of_building_blocks}
Builing blocks per mol.:   {building_blocks_per_molecule}
Number of Water Molecules: {water_molecules}
Number of ions:            {counterions} {counterion}

Fraction         Input Output  Score
Carbon          {input_carbon_fraction:6.3f} {carbon_fraction:6.3f} {diff_carbon_fraction:6.3f}
Hydrogen           - {hydrogen_fraction:6.3f}
Oxygen             - {oxygen_fraction:6.3f}
Nitrogen        {input_nitrogen_fraction:6.3f} {nitrogen_fraction:6.3f} {diff_nitrogen_fraction:6.3f}
Sulphur            - {sulphur_fraction:6.3f}

CarbonylC       {input_carbonyl_c_fraction:6.3f} {carbonyl_c_fraction:6.3f} {diff_carbonyl_c_fraction:6.3f}
CarboxylC       {input_carboxyl_c_fraction:6.3f} {carboxyl_c_fraction:6.3f} {diff_carboxyl_c_fraction:6.3f}
O_Aryl          {input_o_aryl_fraction:6.3f} {o_aryl_fraction:6.3f} {diff_o_aryl_fraction:6.3f}
Aryl            {input_aryl_fraction:6.3f} {aryl_fraction:6.3f} {diff_aryl_fraction:6.3f}
Di_O_Alkyl      {input_di_o_alkyl_fraction:6.3f} {di_o_alkyl_fraction:6.3f} {diff_di_o_alkyl_fraction:6.3f}
O_Alkyl         {input_o_alkyl_fraction:6.3f} {o_alkyl_fraction:6.3f} {diff_o_alkyl_fraction:6.3f}
Methoxyl        {input_methoxyl_fraction:6.3f} {methoxyl_fraction:6.3f} {diff_methoxyl_fraction:6.3f}
Alkyl_C         {input_alkyl_c_fraction:6.3f} {alkyl_c_fraction:6.3f} {diff_alkyl_c_fraction:6.3f}

                Input Output
[COO-]/[COOH]   {input_deprot_prot_fraction:6.1f} {deprot_prot_fraction:6.1f} {diff_deprot_prot_fraction:6.1f}

Counter
Carbon:         {carbon_count:>10}
Nitrogen:       {nitrogen_count:>10}
CarbonylC:      {carbonyl_c_count:>10}
CarboxylC:      {carboxyl_c_count:>10}
O_Aryl:         {o_aryl_count:>10}
Aryl:           {aryl_count:>10}
Di_O_Alkyl:     {di_o_alkyl_count:>10}
O_Alkyl:        {o_alkyl_count:>10}
Methoxyl:       {methoxyl_count:>10}
Alkyl_C:        {alkyl_c_count:>10}

Building Blocks used
{counter}

charge and masses
Protonated:     {protonated:>6}
Deprotonated:   {deprotonated:>6}
charge:         {charge:>6}
mass:           {mass:.3f}

        """

        fmt_IHSS = """
Model:                     {model}
Random Seed:               {seed}
Number of building blocks: {number_of_building_blocks}
Builing blocks per mol.:   {building_blocks_per_molecule}
Number of Water Molecules: {water_molecules}
Number of ions:            {counterions} {counterion}

Fraction         Input Output  Score
Carbon          {input_carbon_fraction:6.3f} {carbon_fraction:6.3f} {diff_carbon_fraction:6.3f}
Hydrogen           - {hydrogen_fraction:6.3f}
Oxygen             - {oxygen_fraction:6.3f}
Nitrogen        {input_nitrogen_fraction:6.3f} {nitrogen_fraction:6.3f} {diff_nitrogen_fraction:6.3f}
Sulphur            - {sulphur_fraction:6.3f}

CarbonylC       {input_carbonyl_fraction:6.3f} {carbonyl_fraction:6.3f} {diff_carbonyl_fraction:6.3f}
CarboxylC       {input_carboxyl_fraction:6.3f} {carboxyl_fraction:6.3f} {diff_carboxyl_fraction:6.3f}
Aromatic        {input_aromatic_fraction:6.3f} {aromatic_fraction:6.3f} {diff_aromatic_fraction:6.3f}
Acetal          {input_acetal_fraction:6.3f} {acetal_fraction:6.3f} {diff_acetal_fraction:6.3f}
Heteroaliphatic {input_heteroaliphatic_fraction:6.3f} {heteroaliphatic_fraction:6.3f} {diff_heteroaliphatic_fraction:6.3f}
Aliphatic       {input_aliphatic_fraction:6.3f} {aliphatic_fraction:6.3f} {diff_aliphatic_fraction:6.3f}

                 Input Output
[COO-]/[COOH]   {input_deprot_prot_fraction:6.1f} {deprot_prot_fraction:6.1f} {diff_deprot_prot_fraction:6.1f}

Counter
Carbon          {carbon_count:>10}
Nitrogen        {nitrogen_count:>10}
CarbonylC       {carbonyl_count:>10}
CarboxylC       {carboxyl_count:>10}
Aromatic        {aromatic_count:>10}
Acetal          {acetal_count:>10}
Heteroaliphatic {heteroaliphatic_count:>10}
Aliphatic       {aliphatic_count:>10}

Building Blocks used
{counter}

charge and masses
Protonated:     {protonated:>6}
Deprotonated:   {deprotonated:>6}
charge:         {charge:>6}
mass:           {mass:.3f}

        """
        formats = {
            "IHSS": fmt_IHSS,
            "SSNMR": fmt_SSNMR
        }

        print(formats[self.assignment].format(**dict_all))

        with open(self.workdir + "/stats.txt", "w") as outfile:
            outfile.writelines(formats[self.assignment].format(**dict_all))

        return dict_diff, counter

    def minimize_system(self):
        """Steps of minimization of the system."""

        mol = self.number_of_building_blocks + self.counterions
        watermol = self.water_molecules  # - self.counterions
        cationatom = int(self.current_state.atom_count + self.counterions)
        lastatom = int(self.current_state.atom_count + self.counterions +
                       3 * self.water_molecules)  # - self.counterions))

        with open(self.workdir + "/min_system.imd", "w") as imdfile:
            imdfile.write(imds.minimize_system.format(seed=self.seed, mol=mol, watermol=watermol, bbsatom=int(
                self.current_state.atom_count), cationatom=cationatom, lastatom=lastatom))

        input_min = {}
        input_min['md'] = self.gromosxx_bin_dir + "/md"
        input_min['topo'] = self.workdir + "/system_ions.top"
        input_min['conf'] = self.workdir + "/system_ions.cnf"
        input_min['fin'] = self.workdir + "/min_system.cnf"
        input_min['input'] = self.workdir + "/min_system.imd"
        input_min['output'] = self.workdir + "/min_system.omd"
        #input_min['trc'] = self.workdir + "/min_system.trc"
        input_min['tre'] = self.workdir + "/min_system.tre"
        #command = "export OMP_NUM_THREADS=4 && {md} @topo {topo} @conf {conf} @fin {fin} @input {input} @trc {trc} @tre {tre}  > {output}".format(**input_min)
        command = "export OMP_NUM_THREADS=4 && {md} @topo {topo} @conf {conf} @fin {fin} @input {input} @tre {tre}  > {output}".format(**input_min)
        self.run_command(command)

        if not self.debug:
            self.remove_files([input_min['conf'],input_min['input'],input_min['output'],input_min['tre']])

    def equilibrate_system(self):
        """Equilibration of the system."""

        mol = self.number_of_building_blocks + self.counterions
        watermol = self.water_molecules  # - self.counterions
        cationatom = int(self.current_state.atom_count + self.counterions)
        lastatom = int(self.current_state.atom_count + self.counterions +
                       3 * self.water_molecules)  # - self.counterions))

        # with open(self.workdir + "/eq_system.imd", "w") as imdfile:
        #    imdfile.write(imd.format(seed=self.seed, mol=mol, watermol=watermol, bbsatom=int(self.current_state.atom_count), cationatom=cationatom, lastatom=lastatom))

        #input_eq = {}

        #input_eq['md'] = self.gromosxx_bin_dir + "/md"
        #input_eq['topo'] = self.workdir + "/system_ions.top"
        #input_eq['conf'] = self.workdir + "/min_system.cnf"
        #input_eq['fin'] = self.workdir + "/eq_system.cnf"
        #input_eq['input'] = self.workdir + "/eq_system.imd"
        #input_eq['output'] = self.workdir + "/eq_system.omd"
        #input_eq['trc'] = self.workdir + "/eq_system.trc"
        #input_eq['tre'] = self.workdir + "/eq_system.tre"
        #command = "{md} @topo {topo} @conf {conf} @fin {fin} @input {input} @trc {trc} @tre {tre}  > {output}".format(**input_eq)

        # self.run_command(command)

        jobs = [
            #input output  steps  vel  temp
            ["min",  "eq", 2000000, 1, 420],
            ["eq",  "eq2",  500000, 0, 360],
            ["eq2", "eq3",  500000, 0, 300],
        ]

        for prefix_in, prefix_out, steps, ntivel, temp in jobs:
            with open(self.workdir + "/" + prefix_out + "_system.imd", "w") as imdfile:
                imdfile.write(imds.equilibration_system.format(seed=self.seed, mol=mol, watermol=watermol, bbsatom=int(
                    self.current_state.atom_count), cationatom=cationatom, lastatom=lastatom, temp=temp, ntivel=ntivel, steps=steps))

            input_eq = {}

            input_eq['md'] = self.gromosxx_bin_dir + "/md"
            input_eq['topo'] = self.workdir + "/system_ions.top"
            input_eq['conf'] = self.workdir + "/" + prefix_in + "_system.cnf"
            input_eq['fin'] = self.workdir + "/" + prefix_out + "_system.cnf"
            input_eq['input'] = self.workdir + "/" + prefix_out + "_system.imd"
            input_eq['output'] = self.workdir + "/" + prefix_out + "_system.omd"
            #input_eq['trc'] = self.workdir + "/" + prefix_out + "_system.trc"
            input_eq['tre'] = self.workdir + "/" + prefix_out + "_system.tre"
            #command = "export OMP_NUM_THREADS=4 && {md} @topo {topo} @conf {conf} @fin {fin} @input {input} @trc {trc} @tre {tre}  > {output}".format(
            command = "export OMP_NUM_THREADS=4 && {md} @topo {topo} @conf {conf} @fin {fin} @input {input} @tre {tre}  > {output}".format(
                **input_eq)

            if self.run_equilibration:
                self.run_command(command)

            if not self.debug:
                self.remove_files([input_eq['conf'], input_eq['input'], input_eq['output'], input_eq['tre']])

    def production_system(self):
        """Write production file."""

        mol = self.number_of_building_blocks + self.counterions
        watermol = self.water_molecules  # - self.counterions
        cationatom = int(self.current_state.atom_count + self.counterions)
        lastatom = int(self.current_state.atom_count + self.counterions +
                       3 * self.water_molecules)  # - self.counterions))

        with open(self.workdir + "/md_system.imd", "w") as imdfile:
            imdfile.write(imds.production_system.format(seed=self.seed, mol=mol, watermol=watermol, bbsatom=int(
                self.current_state.atom_count), cationatom=cationatom, lastatom=lastatom))

        #input_md = {}

        #input_md['md'] = self.gromosxx_bin_dir + "/md"
        #input_md['topo'] = self.workdir + "/system_ions.top"
        #input_md['conf'] = self.workdir + "/eq3_system.cnf"
        #input_md['fin'] = self.workdir + "/md_system.cnf"
        #input_md['input'] = self.workdir + "/md_system.imd"
        #input_md['output'] = self.workdir + "/md_system.omd"
        #input_md['trc'] = self.workdir + "/md_system.trc"
        #input_md['tre'] = self.workdir + "/md_system.tre"
        #command = "{md} @topo {topo} @conf {conf} @fin {fin} @input {input} @trc {trc} @tre {tre}  > {output}".format(**input_md)

        # self.run_command(command)


    def retrieve_data(self, files=[], output_path=None):
        """Retrieve data of the system."""

        if output_path is None:
            output_path = "."

        if not files:
            files = ["system_ions.top", "stats.txt", "min_system.cnf", "min_system.pdb"]

            if self.run_equilibration:
                files.extend(["eq_system.cnf", "eq2_system.cnf", "eq3_system.cnf", "eq3_system.pdb", "md_system.imd"])
            else:
                files.extend(["eq_system.imd", "eq2_system.imd", "eq3_system.imd", "md_system.imd"])

            if self.togromacs:
                files.extend(["min_system.gro", "eq3_system.gro", "system_ions_gmx.top", "*.itp", "md_system_gmx.mdp"])

            if self.debug:
                files.extend(["min_system.tre",  "eq_system.tre", "eq2_system.tre", "eq3_system.tre"])

        command = "cp "
        for f in files:
            command += self.workdir + "/" + f + " "
        command += output_path

        self.run_command(command)


    def gen_GROMACS_input_files(self):

        from vsomm_modeler import gromos2gromacs
        gromos2gromacs.gen_GROMACS_topology(self.workdir, "system_ions.top", len(self.molecules), self.counterion, self.counterions, self.water_molecules)
        gromos2gromacs.gen_GROMACS_mdp(self.workdir)
        if self.run_equilibration:
            gromos2gromacs.gen_GROMACS_coordinates(self.workdir + "/eq3_system.cnf", self.gromacs_bin_dir)
        else:
            gromos2gromacs.gen_GROMACS_coordinates(self.workdir + "/min_system.cnf", self.gromacs_bin_dir)

    def remove_files(self, files):

        for f in files:
            try:
                os.remove(f)
            except:
                print("Error while deleting file ", f)


    def run_command(self, command):
        """Function to run the GROMOS commands."""

        try:
            with open(self.workdir + "/log.txt", "a") as outfile:
                outfile.write("\n[+] " + command + "\n")

            output = subprocess.check_output(command,
                                             shell=True,
                                             stderr=subprocess.STDOUT,
                                             # timeout=1800,
                                             universal_newlines=True)
        except subprocess.CalledProcessError as err:
            print('[!] ERROR:', command + "\n" + err.output)
            with open(self.workdir + "/log.txt", "a") as outfile:
                # outfile.writelines(command+"\n")
                outfile.writelines(err.output)
            # sys.exit(1) # comment if you want to debug this tool
        else:
            with open(self.workdir + "/log.txt", "a") as outfile:
                # outfile.writelines(command)
                outfile.writelines(output)


