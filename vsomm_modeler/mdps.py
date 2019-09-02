production_system="""
;Run control: A leap-frog algorithm for integrating Newton's equations
integrator               = md
;Total simulation time:
:time step in femtoseconds
dt                       = 0.002
;number of steps
nsteps                   = 50000000
;frequency to write coordinates to output trajectory file
nstxout                  = 500000
;frequency to write velocities to output trajectory file
nstvout                  = 500000
;frequency to write energies to log file
nstlog                   = 5000
;frequency to write energies to energy file
nstenergy                = 5000
;frequency to write coordinates to xtc trajectory
nstxtcout                = 5000
;group(s) to write to xtc trajectory
xtc_grps                 = non-Water
;group(s) to write to energy file
energygrps               = non-Water
;Frequency to update the neighbor list (and the long-range forces,
;when using twin-range cut-off's).
nstlist                  = 5
;Make a grid in the box and only check atoms in neighboring grid cells
;when constructing a new neighbor list every nstlist steps.
ns_type                  = grid
;treatment of electrostatic interactions
coulombtype              = reaction-field
epsilon_rf               = 61.0
rcoulomb                 = 1.4
;treatment of van der waals interactions
rvdw                     = 1.4
; Periodic boudary conditions in all the directions
pbc                      = xyz
;Temperature coupling
tcoupl                   = berendsen
tc-grps                  = non-Water Water
tau_t                    = 0.1 0.1
ref_t                    = 300 300
;Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = isotropic
tau_p                    = 2.0
compressibility          = 4.5e-5 4.5e-5
ref_p                    = 1.0 1.0
;Velocity generation
;gen_vel                 = yes
gen_temp                 = 300
gen_seed                 = -1
;Constrain all bonds
constraints              = all-bonds
constraint-algorithm     = Lincs
lincs-order              = 4
lincs-iter               = 4
"""
