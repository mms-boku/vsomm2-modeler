production_system="""
title                   = HS system
; Run parameters
integrator              = md            ; leap-frog integrator
nsteps                  = 50000000      ; 100 ns
dt                      = 0.002         ; 2 fs
; Output control
nstxout 		= 500000	;frequency to write coordinates to output trajectory file
nstvout 		= 500000	;frequency to write velocities to output trajectory file
nstlog  		= 5000		;frequency to write energies to log file
nstenergy		= 5000		;frequency to write energies to energy file
nstxtcout		= 5000		;frequency to write coordinates to xtc trajectory 
;xtc_grps		= 	        ;group(s) to write to xtc trajectory 
;energygrps		= 	        ;group(s) to write to energy file
;Frequency to update the neighbor list (and the long-range forces, 
;when using twin-range cut-off's).
; Bond parameters
continuation            = yes       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = all-bonds   ; bonds involving H are constrained
lincs-order		= 4
lincs-iter		= 4
; Neighborsearching
cutoff-scheme           = Verlet                   ; Buffered neighbor searching
verlet-buffer-tolerance = 0.0001
ns_type                 = grid                     ; search neighboring grid cells
;nstlist                 = 5                        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.4                      ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.4                      ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype		= reaction-field
epsilon_rf		= 61.0
; Temperature coupling is on
tcoupl                  = berendsen                ; modified Berendsen thermostat
tc-grps                 = HS solvent   	   ; two coupling groups - more accurate
tau_t                   = 0.1    0.1       ; time constant, in ps
ref_t                   = 300    300       ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = berendsen	           ; Pressure coupling on in NPT
pcoupltype              = isotropic                ; uniform scaling of box vectors
tau_p                   = 0.5                      ; time constant, in ps
ref_p                   = 1.0    1.0        ; reference pressure, in bar
compressibility         = 4.5e-5 4.5e-5     ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz                      ; 3-D PBC
; Velocity generation
gen_vel                 = no                       ; Velocity generation is off 
gen_temp		= 300
gen_seed		= -1
"""