

minimize_molecule = """TITLE
steepest descent energy minimization of the complex in vacuum
END
# minimize for 20 cycles
# do 2000 steps
ENERGYMIN
#     NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM
#         1       0      0.1   0.01    0.05   2000    0.0
         1       1      0.1   0.01    0.05   {steps}    0.0    # problem with gromosXX 2018
END
# we have 1 solute and no solvent
SYSTEM
#      NPM      NSM
         1        0
END
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0     {seed}    300.0
END
STEP
#    NSTLIM     T      DT
       {steps}   0.0   0.0005
END
# do it in vacuum
BOUNDCOND
#      NTB      NDFMIN
         0           1
END
# every 10 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR      NTPP
        100         0
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)    NTCS      NTCS0(1)
    {shake}          1    0.00010       1      0.00010
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
#  0  0     1  1      1  1     1  1     1  1  WARNING In_Parameter : Old FORCE block used
 1     1     1      1       {charge}        1
# 1     1     1      1       0        0 # all dummy atoms
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
#
  1  {atomnum}
END
PAIRLIST
#       algorithm: standard (0) (gromos96 like pairlist)
#                      grid (1) (XX grid pairlist)
#       SIZE:      grid cell size (or auto = 0.5 * RCUTP)
#       TYPE:      chargegoup (0) (chargegroup based cutoff)
#                      atomic (1) (atom based cutoff)
#
#   ALGORITHM    NSNB     RCUTP     RCUTL     SIZE   TYPE
            0       5       0.8       1.4      0.4      0
END
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF
       1      0.0       1.4         1 1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
      -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
#WRITETRAJ
# NTWX       controls writing of coordinate trajectory
#       0: no coordinate trajectory is written (default)
#      >0: write solute and solvent coordinates every NTWX steps
#      <0: write solute coordinates every |NTWX| steps
# NTWSE >= 0 selection criteria for coordinate trajectory writing
#       0: write normal coordinate trajectory
#      >0: write minimum-energy coordinate and energy trajectory (based on the
#          energy entry selected by NTWSE and as blocks of length NTWX)
#          (see configuration/energy.cc or ene_ana library for indices)
# NTWV       controls writing of velocity trajectory
#       0: no velocity trajectory is written (default)
#      >0: write solute and solvent velocities every NTWV steps
#      <0: write solute velocities every |NTWV| steps
# NTWF       controls writing of force trajectory
#       0: no force trajectory is written (default)
#      >0: write solute and solvent forces every NTWF steps
#      <0: write solute forces every |NTWF| steps
# NTWE >= 0 controls writing of energy trajectory
#       0: no energy trajectory is written (default)
#      >0: write energy trajectory every NTWE steps
# NTWG >= 0 controls writing of free energy trajectory
#       0: no free energy trajectory is written (default)
#      >0: write free energy trajectory every NTWG steps
# NTWB >= 0 controls writing of block-averaged energy trajectory
#       0: no block averaged energy trajectory is written (default)
#      >0: write block-averaged energy variables every |NTWB| steps
#          (and free energies if NTWG > 0) trajectory
#
#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
#      100          0         0         0       100         0       0
#END
"""

equilibrate_molecule = """TITLE
equilibration of the Building Block in vacuum
END
# we have 1 solute and 910 solvent molecules
SYSTEM
#      NPM      NSM
         1        0
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         1         0        0         0         1       0          0        0  {seed}     300.0
END
# do 10000 steps
STEP
#   NSTLIM         T        DT
    1000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         0         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         1
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
        300     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         1
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       {atomnum}         1         1
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    1000
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 100 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
       100        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     1      {atomnum}
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
# use grid based pairlist generation to speed up
PAIRLIST
#    algorithm: standard(0) (gromos96 like pairlist)
#                     grid(1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup(0) (chargegroup based cutoff)
#                     atomic(1) (atom based cutoff)
#
#    algorithm      NSNB    RCUTP    RCUTL      SIZE    TYPE
            0           5      0.8      1.4       0.4       0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       3       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
# every 100 steps write the energy and coordinates to the
# trajectory
#WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
#       100         0         0         0     100         0         0
#END
"""

minimize_system = """TITLE
equilibration of the system
END
# we have 1 solute and 910 solvent molecules
SYSTEM
#      NPM      NSM
        1   {watermol}
END
ENERGYMIN
#     NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM
#         1       0      0.1   0.01    0.05   2000    0.0
         1       1      0.1   0.01    0.05  10000    0.0    # problem with gromosXX 2018
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         1         0        0         0         1       0          0        0  {seed}     300.0
END
# do 10000 steps
STEP
#   NSTLIM         T        DT
      10000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 100 steps write the energy and coordinates to the
# trajectory
WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
       100         0         0         0     100         0         0
END
# every 100 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
       100        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     3      {bbsatom}   {cationatom}       {lastatom}
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
# use grid based pairlist generation to speed up
PAIRLIST
#    algorithm: standard(0) (gromos96 like pairlist)
#                     grid(1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup(0) (chargegroup based cutoff)
#                     atomic(1) (atom based cutoff)
#
#    algorithm      NSNB    RCUTP    RCUTL      SIZE    TYPE
                0            5      0.8      1.4       0.4       0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       3       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END"""

equilibration_system = """TITLE
equilibration of the system
END
# we have 1 solute and 910 solvent molecules
SYSTEM
#      NPM      NSM
        1   {watermol}
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         {ntivel}         0        0         0         1       0          0        0  {seed}     300.0
END
# do 10000 steps
STEP
#   NSTLIM         T        DT
     {steps}       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
        {temp}     0.1      {temp}     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       {bbsatom}         1         1       {lastatom}        2         2
END
PRESSURESCALE
#       COUPLE: off(0), calc(1), scale(2)
#       SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4)
#       VIRIAL: none(0), atomic(1), group(2)
#
#   COUPLE  SCALE   COMP        TAUP    VIRIAL
    2       2     4.575E-4      2.0     2
#   SEMI (semianisotropic couplings: X, Y, Z)
#       e.g. 1 1 2: x and y jointly coupled and z separately coupled
#       e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath
    1 1 1
#   reference pressure
    0.06102     0.00000     0.00000
    0.00000     0.06102     0.00000
    0.00000     0.00000     0.06102
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    1000
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 100 steps write the energy and coordinates to the
# trajectory
WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
       100         0         0         0     100         0         0
END
# every 100 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
       100        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     3      {bbsatom}    {cationatom}   {lastatom}
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
# use grid based pairlist generation to speed up
PAIRLIST
#    algorithm: standard(0) (gromos96 like pairlist)
#                     grid(1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup(0) (chargegroup based cutoff)
#                     atomic(1) (atom based cutoff)
#
#    algorithm      NSNB    RCUTP    RCUTL      SIZE    TYPE
            0           5      0.8      1.4       0.4       0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       3       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END"""

production_system = """TITLE
equilibration of the system
END
# we have 1 solute and 910 solvent molecules
SYSTEM
#      NPM      NSM
        1   {watermol}
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0  {seed}     300.0
END
# do 10000 steps
STEP
#   NSTLIM         T        DT
    500000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
        300     0.1      300     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       {bbsatom}         1         1       {lastatom}        2         2
END
PRESSURESCALE
#       COUPLE: off(0), calc(1), scale(2)
#       SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4)
#       VIRIAL: none(0), atomic(1), group(2)
#
#   COUPLE  SCALE   COMP        TAUP    VIRIAL
    2       2     4.575E-4      2.0     2          # PREGUNTAR??????s
#   SEMI (semianisotropic couplings: X, Y, Z)
#       e.g. 1 1 2: x and y jointly coupled and z separately coupled
#       e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath
    1 1 1
#   reference pressure
    0.06102     0.00000     0.00000
    0.00000     0.06102     0.00000
    0.00000     0.00000     0.06102
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    1000
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 100 steps write the energy and coordinates to the
# trajectory
WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
       100         0         0         0     100         0         0
END
# every 100 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
       100        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     3      {bbsatom}    {cationatom}   {lastatom}
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
# use grid based pairlist generation to speed up
PAIRLIST
#    algorithm: standard(0) (gromos96 like pairlist)
#                     grid(1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup(0) (chargegroup based cutoff)
#                     atomic(1) (atom based cutoff)
#
#    algorithm      NSNB    RCUTP    RCUTL      SIZE    TYPE
            0           5      0.8      1.4       0.4       0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       3       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END"""