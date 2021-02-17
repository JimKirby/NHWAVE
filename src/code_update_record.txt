--------------------------------------------------------------------------------------

                               NHWAVE update records

--------------------------------------------------------------------------------------

UPDATE_00001
by Cheng Zhang
11/23/2016

■Bug Fixes
·Fixed two bugs that leaded to unstability since two variables were not initialized.

■Additions to Functions
·Added a function to record max runup of wave.
·Added a function to record varying bathymetry.
·Added a function to hot start simulation.
·Added a function to limit max maximum Froude number.
·Added a function to generate irregular waves with a internal wave maker.
·Added a function to impose nesting data as boundary conditions.

■Changes
·Changed Sommerfeld open boundary condition to Orlanski open boundary condition.
·Changed subroutine flux_left_right_boundary to impose velocity at inflow boundary.
·Changed subroutine probes to output velocity components at each layer when setting 
 z=-1.0 in stat.txt
 
 
UPDATE_00002
by Cheng Zhang
03/08/2017

■Additions to Functions
·Added a subroutine for fluid slide simulation.
·Added a subroutine for comprehensive landslide bathymetry discription(from Fengyan).
·Added a function to impose wall friction on east, west, south and/or north walls.
·Added a function to identify 2D and 3D solid slide in input file.
·Added a function for 3D solid triangular slide.

■Changes
·Changed long sentences to short ones to avoid compile problems.


UPDATE_00003
by Cheng Zhang
03/23/2017

■Additions to Functions
·Added a function to hot start simulation with moving boundary.

■Bug Fixes
·Fixed a bug that can cause parallelization problem.


UPDATE_00004
by Cheng Zhang
09/05/2017

■Changes
·Changed the old viscous slide model with a new one with TVD scheme, removing
 artificial diffusion from mass conservation and double counted bottom friction 
 from momentum equations.

■Bug Fixes
·Fixed a bug that can cause parallelization problem.
·Fixed two bugs in the granular slide model.


UPDATE_00005
by Cheng Zhang
06/02/2017

■Changes
·Changed influx boundary condition to impose fixed flux.


UPDATE_00006
by Cheng Zhang
04/12/2017

■Bug Fixes
·Fixed a bug in subroutine update_mask to ensure a wet-dry front moves correctly 
 indepedent on the direction.
·Fixed a bug in subroutine eval_duvw that can lead non-zero DU/DV/DW in grid with 
 Mask(i,j)=0.
·Fixed a bug in wall friction boundary condition.

■Additions to Functions
·Added a function to hot start simulation with viscous landslide.

■Changes
·Changed delxH and delyxH at the edge of wet-dry boundary by using Mask9 to avoid 
 extremely height slope at that place.
 
 
UPDATE_00007
by Cheng Zhang
02/04/2018

■Additions to Functions
·Added a function to simulate slump case of rigid landslide.

■Changes
·Changed the IBM code to be applicable to the hydrostatic case. 

UPDATE_00008
J. Kirby
02/08/2019

Collected all known instances (so far of solutions of linear dispersion relation and replaced by a call to new subroutine wave_dispersion.  Created new file waves.F that contains wave_dispersion plus boundary or initial wave generation.  
subroutine source_terms remains to be modified (to extract portions that only need to be executed once.)  Subroutines in waves.F all extracted from nhwave.F except for teh new wave_dispersion.

Update_00009
J. Kirby
02/11/2019

Corrected error made in wave_dispersion subroutine.  Changed some of the copyright statements (gnu -> bsd)

Update_00010
J. Kirby
02/20/2019

Corrected two typos in waves.F
Moved subroutines related to Poisson solver from nhwave.F to pressure.F

UPDATE_00011
by Cheng Zhang
05/28/2019

■Additions to Functions
·Added subroutines to simulate deformable landslide for granular and viscous flow.

■Changes
·Changed the flux boundary condition at wet-dry boundary to avoid numerical oscillation.

UPDATE_00012
by Cheng Zhang
02/16/2021

■Update landslide model
·Change surface gradient method (SGM) in TVD scheme to upwind method for source terms to avoiding negative depth at cell face and ensure well-balance property;
·Changed the flux boundary condition at wet-dry boundary to avoid numerical oscillation;
·Change cell-center based method of non-hydrostatic model to cell-corner based method to avoid numerical osscillation of non-hydrostatic pressure;
·Add a dispersion correction term into governing equations to improve vertical effect of motion.
·Add immersed boundary method (IBM) for slide-structureinteractions.
