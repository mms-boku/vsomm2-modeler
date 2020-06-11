import os
from vsomm_modeler.cnf import CNF
from vsomm_modeler.top import TOP


def generate_gromacs_topology(workdir, topology, molecules, counterion=None, counterions=0, water_molecules=0):

    topo = TOP()
    topo.processFile(workdir + "/" + topology)

    itp_lines = '''
; Include forcefield parameters
#include "gromos54a7.ff/forcefield.itp"

; Include water topology
#include "gromos54a7.ff/spc.itp"

; Include HS topologies
'''

    subtopos = [t for t in topo.generate_subtopologies()]

    for m in range(molecules):
        mol_name = "HS_" + str(m+1)

        sub_itp_filename = mol_name + ".itp"
        with open(workdir + "/" + sub_itp_filename, "w") as subitp_file:
            subitp_file.write(subtopos[m].generate_itp(mol_name))

        itp_lines += '#include "' + sub_itp_filename + '"\n'

    if counterion:
        itp_lines += "\n; Include ion topology\n"
        itp_lines += "#include \"" + counterion + ".itp\"\n"

        sub_itp_filename = counterion + ".itp"
        with open(workdir + "/" + sub_itp_filename, "w") as subitp_file:
            subitp_file.write(subtopos[m+1].generate_itp(counterion))

    itp_lines += '''
[ system ]
; Name
system_ions

[ molecules ]
; Compound            #mols
'''

    for m in range(molecules):
        mol_name = "HS_" + str(m+1)
        itp_lines += '{:12}'.format(mol_name) + '    1\n'
    if counterion:
        itp_lines += '{:12}{:5}\n'.format(counterion, counterions)

    if water_molecules:
        itp_lines += '{:12}{:5}\n'.format("SOL", water_molecules)

    itp_output = open(workdir + "/" + "system_ions_gmx.top", 'w')
    itp_output.write(itp_lines)
    itp_output.close()


def generate_gromacs_mdp(workdir):

    from . import mdps

    with open(workdir + "/" + "/md_system_gmx.mdp", "w") as mdpfile:
        mdpfile.write(mdps.production_system)


def generate_gromacs_index(workdir, humicatom, cation, cationatom, lastatom):

    with open(workdir + "/" + "/index.ndx", "w") as index_file:
        index_file.write("[ HS ]\n")
        for i in range(1, humicatom+1):
            if i % 15 == 0: index_file.write("\n")
            index_file.write(str(i) + " ")
        index_file.write("\n[ "  + cation  + " ]\n")
        for i in range(humicatom + 1, cationatom + 1):
            if i % 15 == 0: index_file.write("\n")
            index_file.write(str(i) + " ")
        index_file.write("\n[ SOLV ]\n")
        for i in range(cationatom + 1, lastatom + 1):
            if i % 15 == 0: index_file.write("\n")
            index_file.write(str(i) + " ")
        index_file.write("\n[ solvent ]\n")
        for i in range(humicatom + 1, lastatom + 1):
            if i % 15 == 0: index_file.write("\n")
            index_file.write(str(i) + " ")
        index_file.write("\n")


def generate_gro_file(cnf_filename):

    workdir, filename = os.path.split(cnf_filename)
    basename = os.path.splitext(filename)[0]

    if not workdir:
        workdir = "."

    cnf = CNF(workdir + "/" + filename)
    cnf.output_gro(workdir + "/" + basename + ".gro")
