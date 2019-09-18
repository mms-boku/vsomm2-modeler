
def at_sort(ats):
    al = []
    for a in ats:
        al.append(int(a))
    sal = []
    for a in sorted(al):
        sal.append(str(a))
    return sal


def gen_itp(mol_name, top):
    itp_lines = '''[ moleculetype ]
; Name            nrexcl
''' + mol_name + '''     0

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
'''
    cg_num = 1
    for at_id in top.atoms:
        at = top.atoms[at_id]
        itp_lines += '{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:7d}'.format(
            at.id, at.a_type.name, at.res.id, at.res.name, at.name, cg_num)
        cg_num += at.mk_cg
        itp_lines += '{:>11s}{:>11s}\n'.format(str(at.p_ch), str(at.m))

    itp_lines += '''
[ bonds ]
;  ai    aj funct            c0            c1            c2            c3
'''

    for b in top.bonds:
        itp_lines += '{:>5s}{:>6s}{:>6s}    {:}\n'.format(
            b.atoms[0].id, b.atoms[1].id, '2', 'gb_' + b.p.id)

    itp_lines += '''
[ exclusions ]
;  ai    aj
'''
    for at_id in top.atoms:
        at = top.atoms[at_id]
        at1 = int(at_id)
        for at_excl in at.e_l:
            if at1 < int(at_excl.id):
                itp_lines += '{:>5s}{:>6s}\n'.format(at_id, at_excl.id)
    for at_id in top.atoms:
        at = top.atoms[at_id]
        at1 = int(at_id)
        # p_l: pairlist
        for at_excl in at.p_l:
            if at1 < int(at_excl.id):
                itp_lines += '{:>5s}{:>6s}\n'.format(at_id, at_excl.id)

    itp_lines += '''
[ pairs ]
;  ai    aj funct
'''
    for at_id in top.atoms:
        at = top.atoms[at_id]
        at1 = int(at_id)
        for at_excl in at.p_l:
            if at1 < int(at_excl.id):
                itp_lines += '{:>5s}{:>6s}     1\n'.format(at_id, at_excl.id)

    itp_lines += '''
[ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3
'''
    for b in top.angles:
        itp_lines += '{:>5s}{:>6s}{:>6s}{:>6s}    {:}\n'.format(
            b.atoms[0].id, b.atoms[1].id, b.atoms[2].id, '2', 'ga_' + b.p.id)

    itp_lines += '''
[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
'''
    for b in top.dihedrals:
        itp_lines += '{:>5s}{:>6s}{:>6s}{:>6s}{:>6s}    {:}\n'.format(
            b.atoms[0].id, b.atoms[1].id, b.atoms[2].id, b.atoms[3].id, '1', 'gd_' + b.p.id)

    itp_lines += '''
[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
'''
    for b in top.impropers:
        itp_lines += '{:>5s}{:>6s}{:>6s}{:>6s}{:>6s}    {:}\n'.format(
            b.atoms[0].id, b.atoms[1].id, b.atoms[2].id, b.atoms[3].id, '2', 'gi_' + b.p.id)

    itp_lines += '\n'

    return itp_lines


def gen_GROMACS_topology(workdir, topo, molecules, counterion=None, counterions=0, water_molecules=0):

    import SMArt.md.gromos as gr
    #from gromos2gromacs import gen_itp
    #from gromos2gromacs import at_sort

    topo_full = gr.parse_top(workdir + "/" + topo)

    itp_lines = '''
; Include forcefield parameters
#include "gromos54a7.ff/forcefield.itp"

; Include water topology
#include "gromos54a7.ff/spc.itp"

; Include HS topologies
'''

    mols = topo_full.get_molecules()

    for m in range(molecules):
        mol_name = "HS_" + str(m+1)

        sub_topo = topo_full.reduce_top(at_sort(mols[m]))
        sub_topo.renumber()

        sub_itp_filename = mol_name + ".itp"
        sub_itp_output = open(workdir + "/" + sub_itp_filename, "w")
        sub_itp_output.write(gen_itp(mol_name, sub_topo))
        sub_itp_output.close()

        itp_lines += '#include "' + sub_itp_filename + '"\n'

    if counterion:
        itp_lines += "\n; Include ion topology\n"
        itp_lines += "#include \"" + counterion + ".itp\"\n"

        # search for the topology of the next "molecule" that is a cation
        sub_topo = topo_full.reduce_top(at_sort(mols[molecules+1]))
        sub_topo.renumber()

        sub_itp_filename = counterion + ".itp"
        sub_itp_output = open(workdir + "/" + sub_itp_filename, "w")
        sub_itp_output.write(gen_itp(counterion, sub_topo))
        sub_itp_output.close()

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


def gen_GROMACS_mdp(workdir):

    from . import mdps

    with open(workdir + "/" + "/md_system_gmx.mdp", "w") as mdpfile:
        mdpfile.write(mdps.production_system)


def gen_GROMACS_coordinates(cnf):
    import os
    from SMArt.md.gromos import parse_cnf

    workdir, filename = os.path.split(cnf)
    basename = os.path.splitext(filename)[0]

    if not workdir:
        workdir = "."

    os.system("cp %s/%s %s/%s.g96" % (workdir, filename, workdir, basename))

    smartcnf = parse_cnf(cnf)
    boxsize = " ".join(str(i) for i in smartcnf.box.abc)
    print(boxsize)

    # TODO: use self.gromacs_bin_dir variable
    os.system("gmx editconf -f %s/%s.g96 -o %s/%s.gro -box %s" % (workdir, basename, workdir, basename, boxsize))
