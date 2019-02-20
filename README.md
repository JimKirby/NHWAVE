# NHWAVE
NHWAVE Developer Group and Software Distribution

Purpose:

     Manage code modifications and updates for NHWAVE, and provide code access to end users.
     
Development team members:

                           github user           email                 affiliation                 role
     Jim Kirby             JimKirby              kirby@udel.edu        University of Delaware      Group owner
     Gangfeng Ma           gangfengma            gma@odu.edu           Old Dominion University
     Fengyan Shi           fengyanshi            fyshi@udel.edu        University of Delaware    
     Morteza Derakhti      derakhti              derakhti@uw.edu       University of Washington
     Cheng Zhang           chzhangudel           chzhang@udel.edu      University of Delaware
     
Change history:

11/07/16 - Kirby extracts initialize.F and two_layer_slide.F from nhwave.F

11/23/16 - Cheng 
 Bug Fixes
·Fixed two bugs that leaded to unstability since two variables were not initialized.

 Additions to Functions
·Added a function to record max runup of wave.
·Added a function to record varying bathymetry.
·Added a function to hot start simulation.
·Added a function to limit max maximum Froude number.
·Added a function to generate irregular waves with a internal wave maker.
·Added a function to impose nesting data as boundary conditions.

 Changes
·Changed Sommerfeld open boundary condition to Orlanski open boundary condition.
·Changed subroutine flux_left_right_boundary to impose velocity at inflow boundary.
·Changed subroutine probes to output velocity components at each layer when setting 
 z=-1.0 in stat.txt
 
08/03/2017 - Cheng Zhang

 Additions to Functions
·Added a subroutine for fluid slide simulation.
·Added a subroutine for comprehensive landslide bathymetry discription(from Fengyan).
·Added a function to impose wall friction on east, west, south and/or north walls.
·Added a function to identify 2D and 3D solid slide in input file.
·Added a function for 3D solid triangular slide.

 Changes
·Changed long sentences to short ones to avoid compile problems.

23/03/2017 - Cheng Zhang

 Additions to Functions
·Added a function to hot start simulation with moving boundary.

 Bug Fixes
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
02/06/2017

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
04/02/2018

■Additions to Functions
·Added a function to simulate slump case of rigid landslide.

■Changes
·Changed the IBM code to be applicable to the hydrostatic case. 


02/08/2019 - Kirby

  Collected all known instances (so far of solutions of linear dispersion relation and replaced by a call to new subroutine   
  wave_dispersion.  Created new file waves.F that contains wave_dispersion plus boundary or initial wave generation.  
  subroutine source_terms remains to be modified (to extract portions that only need to be executed once.)  Subroutines in waves.F 
  all extracted from nhwave.F except for the new wave_dispersion.

02/11/2019 - Kirby

   Corrected error made in wave_dispersion subroutine.  Changed some of the copyright statements (gnu -> bsd)

02/20/2019 - Kirby

   Corrected two typos in waves.F
   Moved subroutines related to Poisson solver from nhwave.F to pressure.F




