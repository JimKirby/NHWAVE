!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!
     program NHWAVE
!
!-----------------------------------------------------------------------------------------
!
!    NHWAVE: Nonhydrostatic WAVE dynamics
!
!    Code developer: Gangfeng Ma, University of Delaware
!    Last update: 14/04/2011
!    Last update: 14/08/2011, parallel implementation

!    NHWAVE V2.0_kirby
!
!   This file is part of NHWAVE.

!    NHWAVE is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    NHWAVE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with NHWAVE.  If not, see <http://www.gnu.org/licenses/>.!
!
!----------------------------------------------------------------------------------------
!
     use global
     implicit none
     integer :: j,Istage
     real(SP) :: tbegin,tend

# if defined (PARALLEL)
     call MPI_INIT(ier)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ier)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,NumP,ier)
# endif
!
!-----------------------------------------------------------------------------------------
!
!    Model configuration
!
!    These subroutines are in initialize.F   (Kirby, 6/27/16)
!
!-----------------------------------------------------------------------------------------
!
! record wall time
!
     call wall_time_secs(tbegin)
!
!    read input data
!
     call read_input
!
!    work index (what is a work index?  Kirby 11/7/16
!
     call index
!
!    allocate variables
!
     call allocate_variables
!
!    generate grids
!
     call generate_grid
!
!    read bathymetry
!
# if defined(LANDSLIDE_COMPREHENSIVE)
     call read_bathymetry_comprehensive
# else
     call read_bathymetry
# endif
!
!    initialize model run
!
     call initial
!
!    read nesting data (added by Cheng)
!
# if defined (COUPLING)
     call read_nesting_file
     TIME=TIME_COUPLING_1
# endif
!
!-----------------------------------------------------------------------------------------
!
!    Model time stepping
!
!-----------------------------------------------------------------------------------------
!
     do while (TIME<TOTAL_TIME.and.RUN_STEP<SIM_STEPS)
	 
	   ! nesting (added by Cheng)
# if defined (COUPLING)
       call OneWayCoupling
# endif

       ! time step     
       call estimate_dt

# if defined (LANDSLIDE)
       ! run landslide-generated tsunami
       if(SlideType(1:5)=='RIGID') then ! modified by Cheng to identify 2d and 3d landslide
         call update_bathymetry
       endif
# endif

!added by Cheng for viscous slide
# if defined (FLUIDSLIDE)
       call update_bathymetry
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
       call update_bathymetry_comprehensive
# endif
! end landslide comprehensive

# if defined (POROUSMEDIA)
       ! update porosity
       call read_porosity
# endif

       ! update boundary conditions       
       call update_wave_bc

       ! update mask
       call update_mask

# if defined (TWOLAYERSLIDE)
       ! update mask for lower layer
       call update_maska
# endif

       ! update wind
       call update_wind

       ! update vars
       call update_vars

       ! SSP Runge-Kutta time stepping
       do Istage = 1,It_Order

# if defined (OBSTACLE)
         ! obstacle velocity 
         call set_obsvel

         ! set obstacle flag                                                                                                   
         call set_obsflag
# endif

         ! well-balanced source terms
         call source_terms

         ! fluxes at cell faces
         call fluxes

         ! update all variables
         call eval_duvw(Istage)

# if defined (TWOLAYERSLIDE)
         ! adjust gravitational acceleration
         call adjust_grav

         ! fluxes for lower layer
         call fluxes_ll

         ! source terms for lower layer
         call source_terms_ll

         ! update variables for lower layer
         call eval_huv_ll(Istage) 

         ! update the thickness of upper layer
         call update_hc_ul
# endif

         ! sponge layer
         if(SPONGE_ON) then
           call sponge_damping
         endif

         ! turbulence model
         if(VISCOUS_FLOW) call eval_turb(Istage)

# if defined (SALINITY)
         ! update salinity
         call eval_sali(Istage)

         ! update density
         call eval_dens  
# endif

# if defined (BUBBLE)
         if(TIME>=TIM_B) then
           ! bubble rise velocity
           call bslip_velocity

           ! update bubble concentration
           call eval_bub(Istage)
         endif
# endif


# if defined (SEDIMENT)
         if(TIME>=TIM_Sedi) then

           ! sediment settling velocity
           call settling_velocity

           ! update sediment concentration
           call eval_sedi(Istage)
   
           ! update mixture density
           if(COUPLE_FS) call eval_dens

           ! bed-load sediment transport
           if(BED_LOAD) call eval_bedload

           ! update bed elevation
           if(BED_CHANGE) call update_bed(Istage)          

         endif
# endif

         ! nesting (added by Cheng)
# if defined (COUPLING)
         call OneWayCoupling
# endif

       enddo

# if defined (BALANCE2D)
       ! evaluate momentum balance in cross-shore
       call eval_balance
# endif

       ! wave average quantities
       if(WAVE_AVERAGE_ON) then
         call wave_average
       endif

       ! screen output
       Screen_Count = Screen_Count+dt
       if(Screen_Count>=Screen_Intv) then
         Screen_Count = Screen_Count-Screen_Intv
         call statistics
       endif
	   
	   ! added by Cheng for recording Hmax
	   call max_min_property

       ! probe output to files
       if(NSTAT>0) then
         Plot_Count_Stat = Plot_Count_Stat+dt
         if(Plot_Count_Stat>=Plot_Intv_Stat) then
           Plot_Count_Stat=Plot_Count_Stat-Plot_Intv_Stat
           call probes
         endif
       endif

       ! field output to files
       if(TIME>=Plot_Start) then
         Plot_Count = Plot_Count+dt
         if(Plot_Count>=Plot_Intv) then
           Plot_Count=Plot_Count-Plot_Intv
           call preview
         endif
       endif
      
     end do
	 
     ! close nesting file (added by Cheng)
# if defined (COUPLING)
     CLOSE(11)
# endif

     ! write out wave height and setup
     if(WAVE_AVERAGE_ON) then
       call print_wh_setup
     endif

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Normal Termination!'
     if(myid.eq.0) write(3,*) 'Normal Termination!'
# else
     write(*,*) 'Normal Termination!'
     write(3,*) 'Normal Termination!'
# endif

     ! wall time at the end
     call wall_time_secs(tend)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
# else
     write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     write(3,*) 'Simulation takes',tend-tbegin,'seconds'
# endif

# if defined (PARALLEL)
     call MPI_FINALIZE(ier)
# endif

     end

# if defined (BALANCE2D)
     subroutine eval_balance
!-------------------------------------------                                                  
!    Evaluate momentum balance in y-dir
!    Called by                         
!       main
!    Last update: 25/12/2012, Gangfeng Ma
!-------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
     real(SP) :: Dz1,Cdrag,Umag

     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+(P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+(P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))   
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     ! hydrostatic pressure gradient
     DEDX2D = zero
     DEDY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         DEDX2D(i,j) = -grav*D(i,j)*(Eta(i+1,j)-Eta(i-1,j))/(2.*dx)
         DEDY2D(i,j) = -grav*D(i,j)*(Eta(i,j+1)-Eta(i,j-1))/(2.*dy)
       endif
     enddo
     enddo

     ! dynamic pressure contribution
     DPDX2D = zero
     DPDY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend
           DPDX2D(i,j) = DPDX2D(i,j)-D(i,j)/Rho0*  &
                (DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))*dsig(k)
           DPDY2D(i,j) = DPDY2D(i,j)-D(i,j)/Rho0*  &
                (DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))*dsig(k)
         enddo
       endif
     enddo
     enddo

     ! turbulent horizontal diffusion
     DIFFX2D = zero
     DIFFY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend
           DIFFX2D(i,j) = DIFFX2D(i,j)+(Diffxx(i,j,k)+Diffxy(i,j,k))*dsig(k)
           DIFFY2D(i,j) = DIFFY2D(i,j)+(Diffyx(i,j,k)+Diffyy(i,j,k))*dsig(k)
         enddo
       endif
     enddo
     enddo

# if defined (VEGETATION)
     FVEGX2D = zero
     FVEGY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
         if(sigc(k)*D(i,j)<=VegH) then
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           FVEGX2D(i,j) = FVEGX2D(i,j)-0.5*StemD*VegDens*VegDrag*Umag*DU(i,j,k)*dsig(k)
           FVEGY2D(i,j) = FVEGY2D(i,j)-0.5*StemD*VegDens*VegDrag*Umag*DV(i,j,k)*dsig(k)
         endif
       endif
     enddo
     enddo
# endif

     ! bottom friction contribution
     TAUBX2D = zero
     TAUBY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         Dz1 = 0.5*D(i,j)*dsig(Kbeg)
         if(ibot==1) then
           if(hc(i,j)>0.1) then
             Cdrag = Cd0
           else
             Cdrag = 0.2
           endif
         else
           Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
         endif
         TAUBX2D(i,j) = -Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*U(i,j,Kbeg)
         TAUBY2D(i,j) = -Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*V(i,j,Kbeg)
       endif
     enddo
     enddo

     ! acceleration
     DUDT2D = zero
     DVDT2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend 
           DVDT2D(i,j) = DVDT2D(i,j)+(DV(i,j,k)-DV0(i,j,k))/dt*dsig(k)
         enddo
       endif      
     enddo
     enddo

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     end subroutine eval_balance
# endif

     subroutine projection_corrector
!-------------------------------------------
!    Correct the velocity field
!    Called by
!       eval_duvw
!    Last update: 25/03/2011, Gangfeng Ma
!-------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
 
     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+  &
              (P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+  &
              (P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

# if defined (POROUSMEDIA)
       DU(i,j,k) = DU(i,j,k)-D(i,j)*dt/Rho0*(DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))/(1.0+Cp_Por(i,j,k))
       DV(i,j,k) = DV(i,j,k)-D(i,j)*dt/Rho0*(DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))/(1.0+Cp_Por(i,j,k))               
       DW(i,j,k) = DW(i,j,k)-dt/Rho0*DelzP(i,j,k)/(1.0+Cp_Por(i,j,k))
# else
       DU(i,j,k) = DU(i,j,k)-D(i,j)*dt/Rho0*(DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))
       DV(i,j,k) = DV(i,j,k)-D(i,j)*dt/Rho0*(DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))
       DW(i,j,k) = DW(i,j,k)-dt/Rho0*DelzP(i,j,k)
# endif
     enddo
     enddo
     enddo

# if defined (OBSTACLE)
     ! update velocities for calculating IB forces
     call get_UVW     

     ! calculate forcing at obstacle boundary                                                                                    
     call imm_obs

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       if(set_flag(i,j,k)==1) then
         DU(i,j,k) = D(i,j)*obs_u
         DV(i,j,k) = D(i,j)*obs_v
         DW(i,j,k) = D(i,j)*obs_w
       else
         DU(i,j,k) = DU(i,j,k)+dt*ObsForceX(i,j,k)
         DV(i,j,k) = DV(i,j,k)+dt*ObsForceY(i,j,k)
         DW(i,j,k) = DW(i,j,k)+dt*ObsForceZ(i,j,k)
       endif
     enddo
     enddo
     enddo
# endif

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     return
     end subroutine projection_corrector
     

     subroutine poisson_solver
!--------------------------------------------
!    Solve poisson equation for dynamic pressure
!    Called by
!       eval_duvw
!    Last update: 24/03/2011, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer :: i,j,k,imask
# if !defined (PARALLEL)
     ! variables for serial computation
     real(SP), dimension(:), allocatable :: Wksp
     integer,  dimension(:), allocatable :: IWksp
     real(SP), dimension(neqns) :: Phi
     real(SP) :: RPARM(30),Pbar(1)
     integer :: IPARM(30),S(1),IS(1),nwksp,inwksp,Ndim,Mdim,N,Maxnz,ierr,neq
     external :: MIC3,IC3,SOR3,GMRES,CG,BCGS
# endif

     ! generate coefficient matrix and rhs
     call generate_coef_rhs

# if defined (PARALLEL)
     ! use HYPRE package for parallel computation
     call hypre_pres_solver
# else
     ! use NSPCG package for serial computation
     call dfault(IPARM,RPARM)

     ! reset default values
     IPARM(2) = itmax
     IPARM(3) = 3
     IPARM(4) = 33
     RPARM(1) = tol

     Ndim = 5*neqns
     Mdim = 5*15
     N = neqns
     Maxnz = 15
     nwksp = 30*neqns
     inwksp = 10*neqns

     allocate(Wksp(nwksp))
     allocate(Iwksp(inwksp))
   
     ! initial guess
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       Phi(neq) = P(i,j,k)
     enddo
     enddo
     enddo
     
     ! solve poisson equation using preconditioned CG or GMRES
     if(isolver==1) then
       call nspcg(MIC3,CG,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     elseif(isolver==2) then
       call nspcg(IC3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)  
     elseif(isolver==3) then
       call nspcg(SOR3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     endif   

     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       P(i,j,k) = Phi(neq)
     enddo
     enddo
     enddo

     deallocate(Wksp)
     deallocate(Iwksp)
# endif
  
     ! fyshi gave boundary condition for dry cells
     ! set zero for dry set is inaccurate
     ! dry cells
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) then
         P(i,j,k) = Zero
         
         ! south boundary 
         if(Mask(i,j+1)==1)then
           do imask=1,Nghost
             P(i,j-imask+1,k)=P(i,j+imask,k)
           enddo
         ! north boundary
         elseif(Mask(i,j-1)==1)then
           do imask=1,Nghost
             P(i,j+imask-1,k)=P(i,j-imask,k)
           enddo
         ! west boundary
         elseif(Mask(i+1,j)==1)then
           do imask=1,Nghost
             P(i-imask+1,j,k)=P(i+imask,j,k)
           enddo
         ! east boundary
         elseif(Mask(i-1,j)==1)then
           do imask=1,Nghost
             P(i+imask-1,j,k)=P(i-imask,j,k)
           enddo
         endif
       endif 
     enddo
     enddo
     enddo

     ! collect into ghost cells
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       do i = 1,Nghost
         P(Ibeg-i,j,k) = P(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend 
       do i = 1,Nghost     
         P(Iend+i,j,k) = P(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jbeg-j,k) = P(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jend+j,k) = P(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_3D_exch(P)
# endif

     end subroutine poisson_solver


     subroutine hypre_pres_solver
!---------------------------------------------
!    solve pressure using hypre package
!    called by
!       poisson_solver
!    Last update: 22/08/2011, Gangfeng Ma
!---------------------------------------------
# if defined (PARALLEL)
     use global
     implicit none
     integer, parameter :: ndim=3
     integer, parameter :: nentries=15
     integer :: i,j,k,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(Mloc),j_glob(Nloc),k_glob(Kloc)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries), &
                periodic_shift(ndim)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,0,1,0,0,0,1,0,0,-1,1,-1,0,1,  &
             0,0,1,1,0,1,0,1,1,-1,0,0,0,-1,0,  &
             0,1,-1,1,0,-1,0,0,-1,-1,0,-1,0,-1,-1/

     ! set up a three dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
       k_glob(k) = k-Nghost
     enddo
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     ilower(3) = k_glob(Kbeg)
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)
     iupper(3) = k_glob(Kend)

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)

     if(PERIODIC_X.or.PERIODIC_Y) then
       if(PERIODIC_X) then
         periodic_shift(1) = Mglob
       else
         periodic_shift(1) = 0
       endif
       if(PERIODIC_Y) then
         periodic_shift(2) = Nglob
       else
         periodic_shift(2) = 0
       endif
       periodic_shift(3) = 0
       call HYPRE_StructGridSetPeriodic(grid,periodic_shift,ierr)
     endif

     call HYPRE_StructGridAssemble(grid,ierr)

     ! define the discretization stencil
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)

     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

     ! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)

     call HYPRE_StructMatrixInitialize(matrix,ierr)

     ! set the matrix coefficient
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(Mloc,Nloc,Kloc))
 
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)*nentries
     allocate(values(nvalues))

     ivalues = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coef(indx(i,j,k),n)
       enddo
     enddo
     enddo
     enddo

     call HYPRE_StructMatrixSetBoxValues(matrix,ilower,iupper,nentries,  &
                                  stencil_indices,values,ierr) 
     call HYPRE_StructMatrixAssemble(matrix,ierr)
     !call HYPRE_StructMatrixPrint(matrix,zero,ierr)

     ! set up struct vectors for b and x
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_b,ierr)
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_x,ierr)

     call HYPRE_StructVectorInitialize(vec_b,ierr)
     call HYPRE_StructVectorInitialize(vec_x,ierr)

     ! set the vector coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqns))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j,k)) = P(i,j,k)
     enddo
     enddo
     enddo
     
     call HYPRE_StructVectorSetBoxValues(vec_x,ilower,iupper,Phi,ierr)
     call HYPRE_StructVectorAssemble(vec_x,ierr)
     !call HYPRE_StructVectorPrint(vec_x,zero,ierr)

     ! set up and use a solver
     call HYPRE_StructGMRESCreate(MPI_COMM_WORLD,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,itmax,ierr)
     call HYPRE_StructGMRESSetTol(solver,tol,ierr)
     call HYPRE_StructGMRESSetPrintLevel(solver,0,ierr)
     call HYPRE_StructGMRESSetLogging(solver,0,ierr)

     ! use symmetric SMG as preconditioner
     n_pre = 1; n_post = 1
     call HYPRE_StructSMGCreate(MPI_COMM_WORLD,precond,ierr)
     call HYPRE_StructSMGSetMemoryUse(precond,0,ierr)
     call HYPRE_StructSMGSetMaxIter(precond,1,ierr)
     call HYPRE_StructSMGSetTol(precond,0.0,ierr)
     call HYPRE_StructSMGSetNumPreRelax(precond,n_pre,ierr)
     call HYPRE_StructSMGSetNumPostRelax(precond,n_post,ierr)
     call HYPRE_StructSMGSetLogging(precond,0,ierr)

     ! set up preconditioner
     precond_id = 0
     call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     
     ! do the setup
     call HYPRE_StructGMRESSetup(solver,matrix,vec_b,vec_x,ierr)
 
     ! do the solve
     call HYPRE_StructGMRESSolve(solver,matrix,vec_b,vec_x,ierr)

     ! get results
     call HYPRE_StructVectorGetBoxValues(vec_x,ilower,iupper,Phi,ierr)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       P(i,j,k) = Phi(indx(i,j,k))
     enddo
     enddo
     enddo

     ! get some info
     !call HYPRE_StructGMRESGetFinalRelati(solver,final_res_norm,ierr)
     !call HYPRE_StructGMRESGetNumIteratio(solver,num_iterations,ierr);
     !
     !if(myid.eq.0) then
     !  write(*,*)'Iterations = ',num_iterations
     !  write(*,*)'Final Relative Residual Norm = ',final_res_norm
     !endif

     ! free memory
     call HYPRE_StructGridDestroy(grid,ierr)
     call HYPRE_StructStencilDestroy(stencil,ierr)
     call HYPRE_StructMatrixDestroy(matrix,ierr)
     call HYPRE_StructVectorDestroy(vec_b,ierr)
     call HYPRE_StructVectorDestroy(vec_x,ierr)
     call HYPRE_StructGMRESDestroy(solver,ierr)
     call HYPRE_StructSMGDestroy(precond,ierr)

     deallocate(indx)
     deallocate(values)
     deallocate(Phi)

# endif
     return
     end subroutine hypre_pres_solver


     subroutine generate_coef_rhs
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    Called by 
!       poisson_solver
!    Last update: 24/03/2011, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k,neq,n,ic
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
     integer,  dimension(:,:,:), allocatable :: indx

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))

     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero
     do k = Kbeg,Kend1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
     enddo
   
     ! generate coefficient matrix
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo 
     enddo

     ! generate source term 
     Rhs = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
# if defined (POROUSMEDIA)
       Rhs(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt* &
              (1+Cp_Por(i,j,k))
# else
       Rhs(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
# endif
     enddo
     enddo
     enddo

     Coef = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Coef(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))+  &
                A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1)))
       Coef(indx(i,j,k),2) = -1./(dx*dx)
       Coef(indx(i,j,k),3) = -1./(dy*dy)
       Coef(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))   
       Coef(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))
       Coef(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),9) = -1./(dx*dx)
       Coef(indx(i,j,k),10) = -1./(dy*dy)
       Coef(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1))
       Coef(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
     enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks face 02/15/2013
     do i = Ibeg+1,Iend-1
     do j = Jbeg+1,Jend-1
     do k = Kbeg,Kend
       if(mask(i,j)==0) then
         ! left 
         if(mask(i+1,j)==1) then
           ic = indx(I+1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
           Coef(ic,9) = Zero
           Coef(ic,5) = Zero
           Coef(ic,14) = Zero
         ! right 
         elseif(mask(i-1,j)==1) then
           ic = indx(I-1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
           Coef(ic,2) = Zero
           Coef(ic,7) = Zero
           Coef(ic,12) = Zero
         ! south
         elseif(mask(i,j+1)==1) then
           ic = indx(i,J+1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
           Coef(ic,10) = Zero
           Coef(ic,4) = Zero
           Coef(ic,15) = Zero
         ! north
         elseif(mask(i,j-1)==1) then
           ic = indx(i,J-1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
           Coef(ic,3) = Zero
           Coef(ic,8) = Zero
           Coef(ic,11) = Zero
         endif ! end mask+1=1 
       endif ! end mask=0
     enddo
     enddo
     enddo

# if defined (OBSTACLE)
     do i = Ibeg+1,Iend-1
     do j = Jbeg+1,Jend-1
     do k = Kbeg+1,Kend-1
       if(set_flag(i,j,k)==1) then
        ! left 
         if(set_flag(i+1,j,k)==0) then
           ic = indx(I+1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,9) = Zero
         ! right 
         elseif(set_flag(i-1,j,k)==0) then
           ic = indx(I-1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,2) = Zero
         ! south
         elseif(set_flag(i,j+1,k)==0) then
           ic = indx(i,J+1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,10) = Zero
         ! north
         elseif(set_flag(i,j-1,k)==0) then
           ic = indx(i,J-1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,3) = Zero
         ! bottom
         elseif(set_flag(i,j,k+1)==0) then
           ic = indx(i,j,k+1)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,13)
           Coef(ic,13) = Zero
         ! upper
         elseif(set_flag(i,j,k-1)==0) then
           ic = indx(i,j,k-1)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,6)
           Coef(ic,6) = Zero
         endif  
       endif
     enddo
     enddo
     enddo
# endif

     ! boundary conditions
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by cheng for nesting, search (COUPLING) to find rest in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     i = Ibeg
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
       Coef(ic,9) = Zero
       Coef(ic,5) = Zero
       Coef(ic,14) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! right side
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     i = Iend
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
       Coef(ic,2) = Zero
       Coef(ic,7) = Zero
       Coef(ic,12) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! front side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     j = Jbeg
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)         
       Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
       Coef(ic,10) = Zero
       Coef(ic,4) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! back side
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     j = Jend
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
       Coef(ic,3) = Zero
       Coef(ic,8) = Zero
       Coef(ic,11) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! bottom side
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
# if defined (LANDSLIDE)
       if(SlideType(1:5)=='RIGID') then ! modified by Cheng to identify 2d and 3d landslide
         Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
       endif
# endif

!added by Cheng for viscous slide
# if defined (FLUIDSLIDE)
       Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
       Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
# endif
! end landslide comprehense

# if defined (TWOLAYERSLIDE)
       if(D(i,j)>0.05) then
         Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &             
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
       endif
# endif

       Coef(ic,6) = Coef(ic,6)+Coef(ic,13)
       Coef(ic,7) = Coef(ic,7)+Coef(ic,12)
       Coef(ic,8) = Coef(ic,8)+Coef(ic,11)
       Coef(ic,5) = Coef(ic,5)+Coef(ic,14)
       Coef(ic,4) = Coef(ic,4)+Coef(ic,15)
       Coef(ic,13) = Zero
       Coef(ic,12) = Zero
       Coef(ic,11) = Zero
       Coef(ic,14) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo

     ! top side (Dirichlet boundary)
     k = Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,4) = Zero
       Coef(ic,5) = Zero
       Coef(ic,6) = Zero
       Coef(ic,7) = Zero
       Coef(ic,8) = Zero
     enddo
     enddo

     ! take (i=2,j=2,k=2) to obtain the diagonal information
     JCoef(1) = indx(Ibeg+1,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k)
     JCoef(2) = indx(Ibeg+2,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k) 
     JCoef(3) = indx(Ibeg+1,Jbeg+2,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k)
     JCoef(4) = indx(Ibeg+1,Jbeg,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i,j-1,k+1)
     JCoef(5) = indx(Ibeg,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k+1)
     JCoef(6) = indx(Ibeg+1,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k+1)
     JCoef(7) = indx(Ibeg+2,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k+1)
     JCoef(8) = indx(Ibeg+1,Jbeg+2,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k+1)
     JCoef(9) = indx(Ibeg,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k)
     JCoef(10) = indx(Ibeg+1,Jbeg,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j-1,k)
     JCoef(11) = indx(Ibeg+1,Jbeg+2,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j+1,k-1)
     JCoef(12) = indx(Ibeg+2,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i+1,j,k-1)
     JCoef(13) = indx(Ibeg+1,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j,k-1)
     JCoef(14) = indx(Ibeg,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i-1,j,k-1)
     JCoef(15) = indx(Ibeg+1,Jbeg,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i,j-1,k-1)

     deallocate(DelxS)
     deallocate(DelyS)
     deallocate(DelzS)
     deallocate(A1) 
     deallocate(indx)

     return
     end subroutine generate_coef_rhs


     subroutine eval_duvw(ISTEP)
!-----------------------------------------------
!    Update all variables D,U,V,W,Omega
!    Called by
!       main
!    Last update: 25/12/2010, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
     real(SP),dimension(:,:),allocatable :: R1,qz
     real(SP),dimension(:,:,:),allocatable :: R2,R3,R4
     real(SP) :: dedt,Umag,Dz1,Cdrag,Ustar2,Wtop,Wbot
     integer :: i,j,k,n,Ista,Nlen
     ! added by Cheng for limiting the maximum Froude number
# if defined(FROUDE_CAP)
     REAL(SP) :: FroudeU,DUU,Dangle
# endif
	 
     Nlen = Kend-Kbeg+1

     allocate(qz(Mloc,Nloc))
     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc,Kloc))
     allocate(R3(Mloc,Nloc,Kloc))
     allocate(R4(Mloc,Nloc,Kloc))
     allocate(Acoef(Nlen))
     allocate(Bcoef(Nlen))
     allocate(Ccoef(Nlen))
     allocate(Xsol(Nlen))
     allocate(Rhs0(Nlen))

     ! calculate baroclinic pressure gradient
     if(.not.BAROTROPIC) call baropg_z

     ! estimate horizontal diffusion terms
     if(VISCOUS_FLOW) call diffusion

     ! external forcing
     if(EXTERNAL_FORCING) call driver(ExtForceX,ExtForceY)

# if defined (VEGETATION)
     if(trim(Veg_Type)=='FLEXIBLE') then
       ! update flexible vegetation height
       qz = 0.0
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         ! find depth-averaged load                                                            
         do k = Kbeg,Kend
           if(sigc(k)*D(i,j)<=(2./3.*FVegH(i,j)+1./3.*VegH)) then
             Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
             qz(i,j) = qz(i,j)+0.5*Rho0*VegDrag*foliage(i,j)*VegH/FVegH(i,j)*  &
                 Umag**2*StemD*dsig(k)*D(i,j)
           endif
         enddo
         qz(i,j) = qz(i,j)/(2./3.*FVegH(i,j)+1./3.*VegH)

         ! non-dimensionalize load
         qz(i,j) = qz(i,j)*VegH**3/EI
       enddo
       enddo

       ! estimate height of flexible vegetation                                         
       call veg_height(qz)

       ! estimate effects of foliage                                                         
       call foli(qz)
     endif
# endif

     ! solve total water depth D
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R1(i,j) = R1(i,j)-1.0/dx*(Ex(i+1,j,k)-Ex(i,j,k))*dsig(k)  &
                        -1.0/dy*(Ey(i,j+1,k)-Ey(i,j,k))*dsig(k)
       enddo
       ! internal wavemaker
       R1(i,j) = R1(i,j)+D(i,j)*SourceC(i,j)
       D(i,j) = ALPHA(ISTEP)*D0(i,j)+BETA(ISTEP)*(D(i,j)+dt*R1(i,j)) 
     enddo
     enddo

     ! update D and Eta          
     D = max(D,MinDep)
     call wl_bc
     Eta = D-Hc

     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)

     ! sigma transformation coefficient                                                  
     call sigma_transform

     ! prepare right-hand side terms
     R2 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle  
       do k = Kbeg,Kend
         R2(i,j,k) = -1.0/dx*(Fx(i+1,j,k)-Fx(i,j,k))-1.0/dy*(Fy(i,j+1,k)-Fy(i,j,k)) &
                -1.0/dsig(k)*(Fz(i,j,k+1)-Fz(i,j,k))+fcor*DV(i,j,k)+DRhoX(i,j,k)+  &
                SourceX(i,j)+Diffxx(i,j,k)+Diffxy(i,j,k)+ExtForceX(i,j,k)
# if defined (CORALREEF)
         if(Hc(i,j)<0.12) then
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           R2(i,j,k) = R2(i,j,k)-0.5*Creef*Umag*U(i,j,k)
         endif
# endif
       enddo
     enddo
     enddo

     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R3(i,j,k) = -1.0/dx*(Gx(i+1,j,k)-Gx(i,j,k))-1.0/dy*(Gy(i,j+1,k)-Gy(i,j,k)) &   
                        -1.0/dsig(k)*(Gz(i,j,k+1)-Gz(i,j,k))-fcor*DU(i,j,k)+DRhoY(i,j,k)  &
                        +SourceY(i,j)+Diffyx(i,j,k)+Diffyy(i,j,k)+ExtForceY(i,j,k)
# if defined (CORALREEF)
         if(Hc(i,j)<0.12) then
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           R3(i,j,k) = R3(i,j,k)-0.5*Creef*Umag*V(i,j,k)
         endif
# endif
       enddo
     enddo
     enddo

     R4 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R4(i,j,k) = -1.0/dx*(Hx(i+1,j,k)-Hx(i,j,k))-1.0/dy*(Hy(i,j+1,k)-Hy(i,j,k)) &  
                        -1.0/dsig(k)*(Hz(i,j,k+1)-Hz(i,j,k))+Diffzx(i,j,k)+Diffzy(i,j,k)
# if defined (CORALREEF)
         if(Hc(i,j)<0.12) then
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           R4(i,j,k) = R4(i,j,k)-0.5*Creef*Umag*W(i,j,k)
         endif
# endif
       enddo
     enddo
     enddo

     ! solve DU
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+  &
                     Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                     0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                     (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+  &
                     CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                     0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                     (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
           endif

# if defined (POROUSMEDIA)
           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip  
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Acoef(Nlen)-Ccoef(Nlen)
           endif
 
           Umag = dsqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(U(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           Bcoef(Nlen) = Bcoef(Nlen)+dt*(Ap_Por(i,j,k)+Bp_Por(i,j,k)*(1+7.5/Xkc)*Umag)
# else        

           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+  &
                CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif
# endif

# if defined (VEGETATION)
           ! account for drag force
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation 
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*Umag+  &
                                   VegVM*(0.25*pi*StemD**2)*VegDens
               endif
             else  ! flexible vegetation    
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*foliage(i,j)*  &       
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*StemD**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(k==Kbeg.and.Bc_Z0==5) then  ! friction law                                                                        
             if(ibot==1) then
               Cdrag = Cd0
             else
               if(D(i,j)>Dfric_Min) then
                 Dz1 = 0.5*D(i,j)*dsig(Kbeg)
               else
                 Dz1 = 0.5*Dfric_Min*dsig(Kbeg)
               endif
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*U(i,j,k)
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DU(i,j,k)+dt*R2(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DU(i,j,k)+dt*R2(i,j,k)+dt*Wsx(i,j)/dsig(k)
           else
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DU(i,j,k)+dt*R2(i,j,k)
           endif
           Rhs0(Nlen) = Rhs0(Nlen)+dt*Cp_Por(i,j,k)*U(i,j,k)*R1(i,j)
# else
           if(k==Kbeg.and.Bc_Z0==5) then  ! friction law
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
# if defined (SEDIMENT)
               Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2
# else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*U(i,j,k)
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)+dt*Wsx(i,j)/dsig(k) 
           else
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)
           endif
# endif

# if defined (VEGETATION)
           ! account for virtual mass force 
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then       
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                  
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DU(i,j,k)
               endif
             else  ! flexible vegetation                                          
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DU(i,j,k)             
               endif
             endif
           endif
# endif
         enddo
      
         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DU(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
# if defined (OBSTACLE)
           if(set_flag(i,j,k)==1) then
             DU(i,j,k) = D(i,j)*obs_u
           else
             DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
           endif
# else

# if defined (POROUSMEDIA)
           DU(i,j,k) = DU(i,j,k)*(1.0+Cp_Por(i,j,k))+  &
                       dt*R2(i,j,k)+dt*Cp_Por(i,j,k)*U(i,j,k)*R1(i,j)

           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(U(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           DU(i,j,k) = DU(i,j,k)/((1.0+Cp_Por(i,j,k))+dt*(Ap_Por(i,j,k)+  &
                       Bp_Por(i,j,k)*(1.0+7.5/Xkc)*Umag))
# else
           DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
# endif
# endif
         enddo
       endif
     enddo
     enddo

     ! solve DV
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+  &
                  CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k-1)))  
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+  &
                  Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                  0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k+1)))   
           endif

# if defined (POROUSMEDIA)
           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip 
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+  &
                CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Acoef(Nlen)-Ccoef(Nlen)
           endif

           Umag = dsqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(V(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           Bcoef(Nlen) = Bcoef(Nlen)+dt*(Ap_Por(i,j,k)+Bp_Por(i,j,k)*(1.0+7.5/Xkc)*Umag)

# else
           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip                                             
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+  &
                CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif
# endif

# if defined (VEGETATION)
           ! account for drag force       
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then         
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                             
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*Umag+  &
                                 VegVM*(0.25*pi*StemD**2)*VegDens                                   
               endif
             else  ! flexible vegetation                                                           
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*foliage(i,j)*  &            
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*StemD**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(k==Kbeg.and.Bc_Z0==5) then
             if(ibot==1) then
               Cdrag = Cd0
             else
               if(D(i,j)>Dfric_Min) then
                 Dz1 = 0.5*D(i,j)*dsig(Kbeg)
               else
                 Dz1 = 0.5*Dfric_Min*dsig(Kbeg)
               endif
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*V(i,j,k)

             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DV(i,j,k)+dt*R3(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DV(i,j,k)+dt*R3(i,j,k)+dt*Wsy(i,j)/dsig(k)
           else
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DV(i,j,k)+dt*R3(i,j,k)
           endif
           Rhs0(Nlen) = Rhs0(Nlen)+dt*Cp_Por(i,j,k)*V(i,j,k)*R1(i,j)
# else
           if(k==Kbeg.and.Bc_Z0==5) then
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
# if defined (SEDIMENT)
               Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2
# else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*V(i,j,k)

             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)+dt*Wsy(i,j)/dsig(k)
           else
             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)
           endif
# endif

# if defined (VEGETATION)
           ! account for virtual mass force                                                       
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then      
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                          
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DV(i,j,k)           
               endif
             else  ! flexible vegetation                                    
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DV(i,j,k)          
               endif
             endif
           endif
# endif

         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DV(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
# if defined (OBSTACLE)
           if(set_flag(i,j,k)==1) then
             DV(i,j,k) = D(i,j)*obs_v
           else
             DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
           endif
# else

# if defined (POROUSMEDIA)
           DV(i,j,k) = (1.0+Cp_Por(i,j,k))*DV(i,j,k)+  &
                       dt*R3(i,j,k)+dt*Cp_Por(i,j,k)*V(i,j,k)*R1(i,j)
 
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(V(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           DV(i,j,k) = DV(i,j,k)/((1.0+Cp_Por(i,j,k))+dt*(Ap_Por(i,j,k)+  &
                       Bp_Por(i,j,k)*(1.0+7.5/Xkc)*Umag))
# else 
           DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
# endif
# endif
         enddo
       endif
     enddo
     enddo

     ! solve DW
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+  &
                  Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k-1))) 
           endif
 
           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+  &
                 Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k+1)))    
           endif

# if defined (POROUSMEDIA)
           if(k==Kbeg) then
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Ccoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &    
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/(0.25*dsig(k)*(dsig(k)+dsig(k-1)))
           else
             Bcoef(Nlen) = (1.0+Cp_Por(i,j,k))-Acoef(Nlen)-Ccoef(Nlen)
           endif

           Umag = dsqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(W(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           Bcoef(Nlen) = Bcoef(Nlen)+  &
                         dt*(Ap_Por(i,j,k)+Bp_Por(i,j,k)*(1.0+7.5/Xkc)*Umag)
# else
           if(k==Kbeg) then
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/(0.25*dsig(k)*(dsig(k)+dsig(k-1)))
           elseif(k==Kend) then
             Bcoef(Nlen) = 1.0-Acoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/(0.25*dsig(k)*(dsig(k)+dsig(k+1)))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif
# endif

# if defined (VEGETATION)
           ! account for drag force
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then    
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                              
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*Umag+  &
                                  VegVM*(0.25*pi*StemD**2)*VegDens                               
               endif
             else  ! flexible vegetation                                                     
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*StemD*VegDens*VegDrag*foliage(i,j)*  &         
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*StemD**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(k==Kbeg) then
             Wbot = -U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j)                                                       
             Rhs0(Nlen) =(1.0+Cp_Por(i,j,k))*DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+  &                                                  
                CmuR(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &                                               
                (0.25*dsig(k)*(dsig(k)+dsig(k-1)))*Wbot
           elseif(k==Kend) then
             Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j)                         
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+CmuR(i,j,k)+  &
                CmuR(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  & 
                (0.25*dsig(k)*(dsig(k)+dsig(k+1)))*Wtop
           else
             Rhs0(Nlen) = (1.0+Cp_Por(i,j,k))*DW(i,j,k)+dt*R4(i,j,k)
           endif
           Rhs0(Nlen) = Rhs0(Nlen)+dt*Cp_Por(i,j,k)*W(i,j,k)*R1(i,j)
# else
           if(k==Kbeg) then
             Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+  &
                CmuR(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k-1)))*Wbot
           elseif(k==Kend) then
             Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+CmuR(i,j,k)+  &
                CmuR(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k+1)))*Wtop
           else
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)
           endif
# endif

# if defined (VEGETATION)
           ! account for virtual mass force  
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then      
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                               
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DW(i,j,k)             
               endif
             else  ! flexible vegetation                                                       
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*StemD**2)*VegDens*DW(i,j,k)              
               endif
             endif
           endif
# endif

         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DW(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
# if defined (OBSTACLE)
           if(set_flag(i,j,k)==1) then
             DW(i,j,k) = D(i,j)*obs_w
           else
             DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
           endif
# else

# if defined (POROUSMEDIA)
           DW(i,j,k) = (1.0+Cp_Por(i,j,k))*DW(i,j,k)+  &
                       dt*R4(i,j,k)+dt*Cp_Por(i,j,k)*W(i,j,k)*R1(i,j)

           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           Xkc = max(abs(W(i,j,k)),0.01)*max(1.0,Per_Wave)/D50_Por
           DW(i,j,k) = DW(i,j,k)/((1.0+Cp_Por(i,j,k))+dt*(Ap_Por(i,j,k)+  &
                       Bp_Por(i,j,k)*(1.0+7.5/Xkc)*Umag))
# else
           DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
# endif
# endif
         enddo 
       endif
     enddo
     enddo

     ! run non-hydrostatic simulation  
     if(NON_HYDRO) then
       ! obtain hydrostatic velocity
       call get_UVW

       ! interpolate velocity into vertical faces
       call interpolate_velocity_to_faces

       ! solve dynamic pressure 
       call poisson_solver

       ! correct velocity field  
       call projection_corrector
     endif

     ! SSP Runge-Kutta time stepping
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DU(i,j,k) = ALPHA(ISTEP)*DU0(i,j,k)+BETA(ISTEP)*DU(i,j,k)
       DV(i,j,k) = ALPHA(ISTEP)*DV0(i,j,k)+BETA(ISTEP)*DV(i,j,k)
       DW(i,j,k) = ALPHA(ISTEP)*DW0(i,j,k)+BETA(ISTEP)*DW(i,j,k)

       if(Mask(i,j)==0) then
         DU(i,j,k) = Zero
         DV(i,j,k) = Zero
         DW(i,j,k) = Zero
       endif
     enddo
     enddo
     enddo
	 
	 ! added by Cheng for limiting the maximum Froude number
# if defined(FROUDE_CAP)
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(I,J)>0)THEN
         FroudeU=SQRT(grav*D(I,J))*FROUDECAP*D(I,J)
         DUU=SQRT(DU(I,J,K)**2+DV(I,J,K)**2)
         IF(DUU>FroudeU)THEN
          Dangle=atan2(DV(I,J,K),DU(I,J,K))
          DU(I,J,K)=FroudeU*COS(Dangle)
          DV(I,J,K)=FroudeU*SIN(Dangle)
         ENDIF
         IF(ABS(DW(I,J,K))>FroudeU)THEN
          DW(I,J,K)=SIGN(FroudeU,DW(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO
# endif

     ! boundary conditions and final velocity
     call get_UVW  

     ! update Omega
     call get_Omega(R1)

     ! if running hydrostatic mode, replace vertical velocity
     ! in fact, W is useless. Only for output
     if(.not.NON_HYDRO) then
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         W(i,j,k) = 0.5*(Omega(i,j,k)+Omega(i,j,k+1))-  &
              DeltH(i,j)+sigc(k)*R1(i,j)-U(i,j,k)*  &
              ((1.0-sigc(k))*DelxH(i,j)+sigc(k)*DelxEta(i,j))-  &
              V(i,j,k)*((1.0-sigc(k))*DelyH(i,j)+sigc(k)*DelyEta(i,j))
         if(Mask(i,j)==0) W(i,j,k) = Zero
         DW(i,j,k) = D(i,j)*W(i,j,k)         
       enddo
       enddo
       enddo

       ! update velocity field
       call get_UVW
     endif

     deallocate(qz)
     deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(R4)
     deallocate(Acoef)
     deallocate(Bcoef)
     deallocate(Ccoef)
     deallocate(Xsol)
     deallocate(Rhs0)

     end subroutine eval_duvw


     subroutine driver(ExtForceX,ExtForceY)
!--------------------------------------------------------------------------
!    Specify external forcing
!    Called by 
!       eval_duvw
!    Last Update: 18/07/2012, Gangfeng Ma
!--------------------------------------------------------------------------
     use global, only: SP,Zero,Mloc,Nloc,Kloc,Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,D, &
                       Pgrad0
     implicit none
     real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: ExtForceX,ExtForceY
     integer :: i,j,k
    
     ExtForceX = Zero
     ExtForceY = Zero

     ! specify energy slope for open channel flow
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ExtForceX(i,j,k) = Pgrad0*D(i,j)
     enddo
     enddo
     enddo

     return
     end subroutine driver


     subroutine interpolate_velocity_to_faces
!------------------------------------------------                                        
!    Interpolate U,V,W to vertical faces                                                 
!    Called by                                                                           
!       main                                                                             
!    Last Update: 19/03/2011, Gangfeng Ma                                                
!------------------------------------------------                                        
     use global, only: SP,U,V,W,Uf,Vf,Wf,dsig,sig,sigc,  &
                       Mloc,Nloc,Kloc,Kloc1,Nghost
     implicit none
     integer  :: i,j,k
     real(SP) :: Int_factor,Int_factor1,Int_factor2,Int_factor3
     logical  :: Linear_Interp=.True.

     if(Linear_Interp) then
       ! first-order linear interpolation
       do k = 2,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor = dsig(k)/(dsig(k)+dsig(k-1))
         Uf(i,j,k) = (1.0-Int_factor)*U(i,j,k)+Int_factor*U(i,j,k-1)
         Vf(i,j,k) = (1.0-Int_factor)*V(i,j,k)+Int_factor*V(i,j,k-1)
         Wf(i,j,k) = (1.0-Int_factor)*W(i,j,k)+Int_factor*W(i,j,k-1)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     else
       ! second-order lagrange interpolation
       do k = 3,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(k)-sigc(k-1))*(sig(k)-sigc(k))/  &
             ((sigc(k-2)-sigc(k-1))*(sigc(k-2)-sigc(k)))
         Int_factor2 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k))/  &
             ((sigc(k-1)-sigc(k-2))*(sigc(k-1)-sigc(k)))
         Int_factor3 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k-1))/  &
             ((sigc(k)-sigc(k-2))*(sigc(k)-sigc(k-1)))
         Uf(i,j,k) = Int_factor1*U(i,j,k-2)+Int_factor2*U(i,j,k-1)+Int_factor3*U(i,j,k)
         Vf(i,j,k) = Int_factor1*V(i,j,k-2)+Int_factor2*V(i,j,k-1)+Int_factor3*V(i,j,k)
         Wf(i,j,k) = Int_factor1*W(i,j,k-2)+Int_factor2*W(i,j,k-1)+Int_factor3*W(i,j,k)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(2)-sigc(2))*(sig(2)-sigc(3))/  &
             ((sigc(1)-sigc(2))*(sigc(1)-sigc(2)))
         Int_factor2 = (sig(2)-sigc(1))*(sig(2)-sigc(3))/  &
             ((sigc(2)-sigc(1))*(sigc(2)-sigc(3)))
         Int_factor3 = (sig(2)-sigc(1))*(sig(2)-sigc(2))/  &
             ((sigc(3)-sigc(1))*(sigc(3)-sigc(2)))
         Uf(i,j,2) = Int_factor1*U(i,j,1)+Int_factor2*U(i,j,2)+Int_factor3*U(i,j,3)
         Vf(i,j,2) = Int_factor1*V(i,j,1)+Int_factor2*V(i,j,2)+Int_factor3*V(i,j,3)
         Wf(i,j,2) = Int_factor1*W(i,j,1)+Int_factor2*W(i,j,2)+Int_factor3*W(i,j,3)
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     endif

     end subroutine interpolate_velocity_to_faces


     subroutine get_Omega(R1)
!-----------------------------------------------
!    Obtain vertical velocity in sigma-coord.
!    Called by 
!       eval_duvw
!    Last update: 30/08/2013, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(Mloc,Nloc), intent(in) :: R1
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! reconstruct flux using new velocities
     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j) 
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
     enddo
     enddo
     enddo     

     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     call wave_speed

     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)

     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Ex(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
         Ex(Iend1,j,k) = Zero
       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front and back side  
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
         Ey(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
         Ey(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! update Omega
     Omega = zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg+1,Kend1
       Omega(i,j,k) = Omega(i,j,k-1)-dsig(k-1)*  &
          (R1(i,j)+(Ex(i+1,j,k-1)-Ex(i,j,k-1))/dx+(Ey(i,j+1,k-1)-Ey(i,j,k-1))/dy)                      
       if(Mask(i,j)==0) Omega(i,j,k) = zero
     enddo
     enddo
     enddo

!     ! adjust omega to make omega=zero at free surface
!     do i = Ibeg,Iend
!     do j = Jbeg,Jend
!       if(abs(Omega(i,j,Kend1))>1.e-8) then
!         do k = Kbeg+1,Kend1
!           Omega(i,j,k) = Omega(i,j,k)-  &
!              float(k-Kbeg)/float(Kend-Kbeg+1)*Omega(i,j,Kend1)                                        
!         enddo
!       endif
!     enddo
!     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine get_Omega


     subroutine get_UVW
!------------------------------------------------
!    Obtain U,V,W
!    Called by
!       eval_duvw
!    Last update: 25/12/2010, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     do j = Jbeg,Jend
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       U(i,j,k) = DU(i,j,k)/D(i,j)
       V(i,j,k) = DV(i,j,k)/D(i,j)
       W(i,j,k) = DW(i,j,k)/D(i,j)
     enddo
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
     call phi_3D_exch(DU)
     call phi_3D_exch(DV)
     call phi_3D_exch(DW)
# endif

     end subroutine get_UVW


     subroutine fluxes
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none

     ! second order construction
     call delxyzFun
     call construction  

     ! calculate wave speed
     call wave_speed

     ! calculate fluxes at faces
     if(ADV_HLLC) then
       call fluxes_at_faces_HLLC
     else
       call fluxes_at_faces_HLL
     endif

     ! impose boundary conditions
     call flux_bc

     end subroutine fluxes


     subroutine construction
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes
!    Last update: 04/01/2011, Gangfeng Ma
!-----------------------------------------
     use global
     implicit none
     integer :: i,j,k

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_x(W,DelxW,WxL,WxR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_x(DW,DelxDW,DWxL,DWxR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
       FxL(i,j,k) = DUxL(i,j,k)*UxL(i,j,k)+0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))
       FxR(i,j,k) = DUxR(i,j,k)*UxR(i,j,k)+0.5*Grav*(EtaxR(i,j)*EtaxR(i,j)+2.0*EtaxR(i,j)*hfx(i,j))
       GxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*VxL(i,j,k)
       GxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*VxR(i,j,k)
       HxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*WxL(i,j,k)
       HxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*WxR(i,j,k)
     enddo
     enddo
     enddo

     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_y(W,DelyW,WyL,WyR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)
     call construct_3D_y(DW,DelyDW,DWyL,DWyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j)
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
       FyL(i,j,k) = DyL(i,j)*UyL(i,j,k)*VyL(i,j,k)
       FyR(i,j,k) = DyR(i,j)*UyR(i,j,k)*VyR(i,j,k)
       GyL(i,j,k) = DVyL(i,j,k)*VyL(i,j,k)+0.5*Grav*(EtayL(i,j)*EtayL(i,j)+2.0*EtayL(i,j)*hfy(i,j))
       GyR(i,j,k) = DVyR(i,j,k)*VyR(i,j,k)+0.5*Grav*(EtayR(i,j)*EtayR(i,j)+2.0*EtayR(i,j)*hfy(i,j))
       HyL(i,j,k) = DyL(i,j)*VyL(i,j,k)*WyL(i,j,k)
       HyR(i,j,k) = DyR(i,j)*VyR(i,j,k)*WyR(i,j,k) 
     enddo
     enddo
     enddo

     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)

     call Kbc_Surface
     call Kbc_Bottom

     end subroutine construction


     subroutine construct_2D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in x-direction
!    Called by
!       construction
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
       OutL(i,j) = Vin(i-1,j)+0.5*dx*Din(i-1,j)
       OutR(i,j) = Vin(i,j)-0.5*dx*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_x


     subroutine construct_3D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in x-direction
!    Called by 
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Kloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i-1,j,k)+0.5*dx*Din(i-1,j,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dx*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_x


     subroutine construct_2D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in y-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------ 
     use global, only: SP,Zero,dy,Mloc,Nloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
       OutL(i,j) = Vin(i,j-1)+0.5*dy*Din(i,j-1)
       OutR(i,j) = Vin(i,j)-0.5*dy*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_y


     subroutine construct_3D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in y-direction 
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dy,Mloc,Nloc,Kloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i,j-1,k)+0.5*dy*Din(i,j-1,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dy*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_y


     subroutine construct_3D_z(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in z-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dsig,Mloc,Nloc,Kloc,Kloc1, &
                       Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,Kend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend1
       OutL(i,j,k) = Vin(i,j,k-1)+0.5*dsig(k-1)*Din(i,j,k-1)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dsig(k)*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_z


     subroutine delxyzFun 
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global
     implicit none
     integer :: i
     
     call delxFun_2D(Eta,DelxEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delxFun_3D(W,DelxW)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delxFun_3D(DW,DelxDW)

     call delyFun_2D(Eta,DelyEta)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delyFun_3D(W,DelyW)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)
     call delyFun_3D(DW,DelyDW)

     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)

     end subroutine delxyzFun

     
     subroutine delxFun_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Mask,Brks
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Mask(i,j)==0) then
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_2D


     subroutine delxFun_3D(Din,Dout)
!------------------------------------------
!    Second-order derivative in x
!    Called by 
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 2,Mloc-1
     do j = 1,Nloc
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i+1,j,k)-Din(i,j,k))/dx
         TMP2 = (Din(i,j,k)-Din(i-1,j,k))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do j = 1,Nloc
     do k = 1,Kloc
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

     return
     end subroutine delxFun_3D

     subroutine delxFun1_3D(Din,Dout)
!------------------------------------------                                                                                          
!    Second-order derivative in x                                                                                                                
!    Called by                                                                                                                                  
!       delxyzFun                                                                                                                               
!    Last update: 04/01/2011, Gangfeng Ma                                                                                                      
!------------------------------------------                                                                                                     
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Kloc1,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc1)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 2,Mloc-1
     do j = 1,Nloc
     do k = 1,Kloc1
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i+1,j,k)-Din(i,j,k))/dx
         TMP2 = (Din(i,j,k)-Din(i-1,j,k))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do j = 1,Nloc
     do k = 1,Kloc1
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

     return
     end subroutine delxFun1_3D


     subroutine delyFun_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y
!    Called by 
!       delxyzFun  
!    Last update: 04/01/2011, Gangfeng Ma  
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Mask(i,j)==0) then 
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_2D


     subroutine delyFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in y
!    Called by
!       delxyzFun 
!    Last update: 04/01/2011, Gangfeng Ma 
!-------------------------------------------
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 2,Nloc-1
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i,j+1,k)-Din(i,j,k))/dy
         TMP2 = (Din(i,j,k)-Din(i,j-1,k))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do k = 1,Kloc
       Dout(i,1,k) = (Din(i,2,k)-Din(i,1,k))/dy
       Dout(i,Nloc,k) = (Din(i,Nloc,k)-Din(i,Nloc-1,k))/dy
     enddo
     enddo

     return 
     end subroutine delyFun_3D


     subroutine delyFun1_3D(Din,Dout)
!------------------------------------------- 
!    Second-order derivative in y       
!    Called by                   
!       delxyzFun                                                                                                                                  
!    Last update: 04/01/2011, Gangfeng Ma                                                                                                          
!-------------------------------------------                                                                                                       
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Kloc1,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc1)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 2,Nloc-1
     do k = 1,Kloc1
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i,j+1,k)-Din(i,j,k))/dy
         TMP2 = (Din(i,j,k)-Din(i,j-1,k))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do k = 1,Kloc1
       Dout(i,1,k) = (Din(i,2,k)-Din(i,1,k))/dy
       Dout(i,Nloc,k) = (Din(i,Nloc,k)-Din(i,Nloc-1,k))/dy
     enddo
     enddo

     return
     end subroutine delyFun1_3D
     
     subroutine delzFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in z
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!-------------------------------------------
     use global, only: SP,Small,Zero,dsig,sigc,Mloc,Nloc,Kloc
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 2,Kloc-1
       TMP1 = (Din(i,j,k+1)-Din(i,j,k))/(sigc(k+1)-sigc(k))
       TMP2 = (Din(i,j,k)-Din(i,j,k-1))/(sigc(k)-sigc(k-1))

       if((abs(TMP1)+abs(TMP2))<Small) then
         Dout(i,j,k) = Zero
       else
         Dout(i,j,k) = LIMITER(TMP1,TMP2)
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do j = 1,Nloc
       Dout(i,j,1) = (Din(i,j,2)-Din(i,j,1))/(0.5*(dsig(1)+dsig(2)))
       Dout(i,j,Kloc) = (Din(i,j,Kloc)-Din(i,j,Kloc-1))/(0.5*(dsig(Kloc-1)+dsig(Kloc)))
     enddo
     enddo

     return
     end subroutine delzFun_3D


     subroutine flux_bc
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    Called by
!       fluxes
!    Last update: 25/12/2010, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by cheng for nesting. Please search (COUPLING) for others in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
	 do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Ex(Ibeg,j,k) = Zero
         Fx(Ibeg,j,k) = 0.5*Grav*(EtaxR(Ibeg,j)*EtaxR(Ibeg,j)+2.0*EtaxR(Ibeg,j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Zero
         Hx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
         Fx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Uin_X0(j,k)+  &
                  0.5*Grav*(Ein_X0(j)*Ein_X0(j)+2.0*Ein_X0(j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Vin_X0(j,k)
         Hx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Win_X0(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
         Ex(Iend1,j,k) = Zero
         Fx(Iend1,j,k) = 0.5*Grav*(EtaxL(Iend1,j)*EtaxL(Iend1,j)+  &
                  2.0*EtaxL(Iend1,j)*hfx(Iend1,j))
         Gx(Iend1,j,k) = Zero
         Hx(Iend1,j,k) = Zero
       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
         Fx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Uin_Xn(j,k)+  &
                 0.5*Grav*(Ein_Xn(j)*Ein_Xn(j)+2.0*Ein_Xn(j)*hfx(Iend1,j)) 
         Gx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Vin_Xn(j,k)
         Hx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Win_Xn(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
         Ey(i,Jbeg,k) = Zero
         Fy(i,Jbeg,k) = Zero
         Gy(i,Jbeg,k) = 0.5*Grav*(EtayR(i,Jbeg)*EtayR(i,Jbeg)+  &
                2.0*EtayR(i,Jbeg)*hfy(i,Jbeg))
         Hy(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
         Ey(i,Jend1,k) = Zero
         Fy(i,Jend1,k) = Zero
         Gy(i,Jend1,k) = 0.5*Grav*(EtayL(i,Jend1)*EtayL(i,Jend1)+  &
                 2.0*EtayL(i,Jend1)*hfy(i,Jend1))
         Hy(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! upper and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,Kbeg) = Zero
       Gz(i,j,Kbeg) = Zero
       Hz(i,j,Kbeg) = Zero
       Fz(i,j,Kend1) = Zero
       Gz(i,j,Kend1) = Zero
       Hz(i,j,Kend1) = Zero
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg-1,Jend+1
     do i = Ibeg-1,Iend+1
       if(Mask(i,j)==0) then
         Ex(i,j,k) = Zero
# if defined (PARALLEL)
         if(i==Ibeg.and.n_west.eq.MPI_PROC_NULL) then
# else
         if(i==Ibeg) then
# endif
           Fx(i,j,k) = Zero
         else
           Fx(i,j,k) = 0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+  &
                2.0*EtaxL(i,j)*hfx(i,j))*Mask(i-1,j)
         endif
         Gx(i,j,k) = Zero
         Hx(i,j,k) = Zero

         Ex(i+1,j,k) = Zero
# if defined (PARALLEL)
         if(i==Iend.and.n_east.eq.MPI_PROC_NULL) then
# else
         if(i==Iend) then
# endif
           Fx(i+1,j,k) = Zero
         else
           Fx(i+1,j,k) = 0.5*Grav*(EtaxR(i+1,j)*EtaxR(i+1,j)+  &
                  2.0*EtaxR(i+1,j)*hfx(i+1,j))*Mask(i+1,j)
         endif
         Gx(i+1,j,k) = Zero
         Hx(i+1,j,k) = Zero

         Ey(i,j,k) = Zero
         Fy(i,j,k) = Zero
# if defined (PARALLEL)
         if(j==Jbeg.and.n_suth.eq.MPI_PROC_NULL) then
# else
         if(j==Jbeg) then
# endif   
           Gy(i,j,k) = Zero
         else
           Gy(i,j,k) = 0.5*Grav*(EtayL(i,j)*EtayL(i,j)+  &
                2.0*EtayL(i,j)*hfy(i,j))*Mask(i,j-1)
         endif
         Hy(i,j,k) = Zero

         Ey(i,j+1,k) = Zero
         Fy(i,j+1,k) = Zero
# if defined (PARALLEL)
         if(j==Jend.and.n_nrth.eq.MPI_PROC_NULL) then
# else
         if(j==Jend) then
# endif
           Gy(i,j+1,k) = Zero
         else
           Gy(i,j+1,k) = 0.5*Grav*(EtayR(i,j+1)*EtayR(i,j+1)+  &
                2.0*EtayR(i,j+1)*hfy(i,j+1))*Mask(i,j+1)
         endif
         Hy(i,j+1,k) = Zero
       endif
     enddo
     enddo
     enddo

     end subroutine flux_bc


     subroutine fluxes_at_faces_HLLC
!---------------------------------------------
!    Fluxes at cell faces estimated by HLLC approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR,D3xLS,D3xRS,  &
                 D3yLS,D3yRS,DUxLS,DUxRS,DVxLS,DVxRS,DWxLS,DWxRS,DUyLS,DUyRS,DVyLS, &
                 DVyRS,DWyLS,DWyRS

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     allocate(D3xLS(Mloc1,Nloc,Kloc))
     allocate(D3xRS(Mloc1,Nloc,Kloc))
     allocate(DUxLS(Mloc1,Nloc,Kloc))
     allocate(DUxRS(Mloc1,Nloc,Kloc))
     allocate(DVxLS(Mloc1,Nloc,Kloc))
     allocate(DVxRS(Mloc1,Nloc,Kloc))
     allocate(DWxLS(Mloc1,Nloc,Kloc))
     allocate(DWxRS(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = DxL(i,j)
       D3xR(i,j,k) = DxR(i,j)
       D3xLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)
       D3xRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)       
       DUxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DUxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DVxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*VxL(i,j,k)
       DVxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*VxR(i,j,k)
       DWxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*WxL(i,j,k)
       DWxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*WxR(i,j,k)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     allocate(D3yLS(Mloc,Nloc1,Kloc))
     allocate(D3yRS(Mloc,Nloc1,Kloc))
     allocate(DUyLS(Mloc,Nloc1,Kloc))
     allocate(DUyRS(Mloc,Nloc1,Kloc))
     allocate(DVyLS(Mloc,Nloc1,Kloc))
     allocate(DVyRS(Mloc,Nloc1,Kloc))
     allocate(DWyLS(Mloc,Nloc1,Kloc))
     allocate(DWyRS(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = DyL(i,j)
       D3yR(i,j,k) = DyR(i,j)
       D3yLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)
       D3yRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)
       DUyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*UyL(i,j,k)
       DUyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*UyR(i,j,k)
       DVyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DVyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DWyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*WyL(i,j,k)
       DWyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*WyR(i,j,k)
     enddo
     enddo
     enddo

     ! horizontal fluxes
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,ExL,ExR,D3xL,D3xLS,D3xR,D3xRS,Ex)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,EyL,EyR,D3yL,D3yLS,D3yR,D3yRS,Ey)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,FxL,FxR,DUxL,DUxLS,DUxR,DUxRS,Fx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,FyL,FyR,DUyL,DUyLS,DUyR,DUyRS,Fy)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,GxL,GxR,DVxL,DVxLS,DVxR,DVxRS,Gx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,GyL,GyR,DVyL,DVyLS,DVyR,DVyRS,Gy)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,HxL,HxR,DWxL,DWxLS,DWxR,DWxRS,Hx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,HyL,HyR,DWyL,DWyLS,DWyR,DWyRS,Hy)     

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
!       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*U(i,j,k))/(dsig(k)+dsig(k-1))
!       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*V(i,j,k))/(dsig(k)+dsig(k-1))                        
!       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)
     deallocate(D3xLS)
     deallocate(D3xRS)
     deallocate(D3yLS)
     deallocate(D3yRS)
     deallocate(DUxLS)
     deallocate(DUxRS)
     deallocate(DUyLS)
     deallocate(DUyRS)
     deallocate(DVxLS)
     deallocate(DVxRS)
     deallocate(DVyLS)
     deallocate(DVyRS)
     deallocate(DWxLS)
     deallocate(DWxRS)
     deallocate(DWyLS)
     deallocate(DWyRS)

     return
     end subroutine fluxes_at_faces_HLLC


     subroutine fluxes_at_faces_HLL
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     ! horizontal fluxes
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,FxL,FxR,DUxL,DUxR,Fx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,FyL,FyR,DUyL,DUyR,Fy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,GxL,GxR,DVxL,DVxR,Gx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,GyL,GyR,DVyL,DVyR,Gy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,HxL,HxR,DWxL,DWxR,Hx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,HyL,HyR,DWyL,DWyR,Hy)     

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
!       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*U(i,j,k))/(dsig(k)+dsig(k-1))
!       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*V(i,j,k))/(dsig(k)+dsig(k-1))
!       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine fluxes_at_faces_HLL


     subroutine HLL(M,N,L,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLL
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSE
         FOUT(I,J,K) = SR(I,J,K)*FL(I,J,K)-SL(I,J,K)*FR(I,J,K)+  &
               SL(I,J,K)*SR(I,J,K)*(UR(I,J,K)-UL(I,J,K))
         IF((ABS(SR(I,J,K)-SL(I,J,K)))<SMALL)THEN
           FOUT(I,J,K) = FOUT(I,J,K)/SMALL
         ELSE
           FOUT(I,J,K) = FOUT(I,J,K)/(SR(I,J,K)-SL(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLL

   
     subroutine HLLC(M,N,L,SL,SR,SS,FL,FR,UL,ULS,UR,URS,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLLC
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,SS,FL,FR,UL,ULS,UR,URS
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSEIF(SS(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)+SL(I,J,K)*(ULS(I,J,K)-UL(I,J,K))
       ELSE
         FOUT(I,J,K) = FR(I,J,K)+SR(I,J,K)*(URS(I,J,K)-UR(I,J,K))
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLLC


     subroutine wave_speed
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    Called by
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!    Last update: 12/04/2011, Gangfeng Ma, wetting-drying
!-----------------------------------------------
     use global, only: SP,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend, &
                       DxL,DxR,DyL,DyR,UxL,UxR,VyL,VyR, &
                       SxL,SxR,SxS,SyL,SyR,SyS,Grav,Mask
     implicit none
     integer  :: i,j,k
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S
     
     ! x-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Mask(i-1,j)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UxL(i,j,k)-UxR(i,j,k))
         U_S = 0.5*(UxL(i,j,k)+UxR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SxL(i,j,k) = min(UxL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SxR(i,j,k) = max(UxR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SxS(i,j,k) = U_S
       elseif(Mask(i-1,j)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SxL(i,j,k) = UxR(i,j,k)-2.0*SQR_PHI_R
         SxR(i,j,k) = UxR(i,j,k)+SQR_PHI_R
         SxS(i,j,k) = SxL(i,j,k)
       elseif(Mask(i-1,j)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SxL(i,j,k) = UxL(i,j,k)-SQR_PHI_L
         SxR(i,j,k) = UxL(i,j,k)+2.0*SQR_PHI_L
         SxS(i,j,k) = SxR(i,j,k)
       endif
     enddo
     enddo
     enddo

     ! y-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Mask(i,j-1)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VyL(i,j,k)-VyR(i,j,k))
         U_S = 0.5*(VyL(i,j,k)+VyR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SyL(i,j,k) = min(VyL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SyR(i,j,k) = max(VyR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SyS(i,j,k) = U_S
       elseif(Mask(i,j-1)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SyL(i,j,k) = VyR(i,j,k)-2.0*SQR_PHI_R
         SyR(i,j,k) = VyR(i,j,k)+SQR_PHI_R
         SyS(i,j,k) = SyL(i,j,k)
       elseif(Mask(i,j-1)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SyL(i,j,k) = VyL(i,j,k)-SQR_PHI_L
         SyR(i,j,k) = VyL(i,j,k)+2.0*SQR_PHI_L
         SyS(i,j,k) = SyR(i,j,k)
       endif
     enddo
     enddo
     enddo

     end subroutine wave_speed


     FUNCTION LIMITER(A,B)
     use global, only: SP,Zero,One,Small
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: LIMITER

!     ! minmod limiter
!     LIMITER=max(Zero,min(A,B))

!     ! van Leer limiter
     LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))

!     ! superbee limiter
!     LIMITER=SIGN(One,B)*MAX(Zero,MIN(2.0*ABS(B),SIGN(One,B)*A),  &
!          MIN(ABS(B),2.0*SIGN(One,B)*A))

     RETURN
     END FUNCTION LIMITER


     subroutine source_terms
!------------------------------------------------
!    This subroutine is used to evaluate source
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,Iter,nn,ndir,nfreq,nk
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Source_Area,myvar, &
                 WnumX,WnumY,Phs_lag,dfreq,ddir,Angle,tmp1,tmp2,tmp3,tmp4,Umag,Ubar,Vbar
     real(SP) :: Ytrough,Mod1,Zup,Zlow,Zmid,Xstart,Zero1,cnoidal_cn,cnoidal_ck,Atmp

     ! internal wavemaker for linear wave
     if(WaveMaker(1:7)=='INT_LIN') then
       ! Find wave number for linear wave (Newton-Ralphson Method)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length
     
       Iter = 0
 55    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
       Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
           Grav*tanh(Wave_Number*Dep_Wave) 
       Wave_Number = Wave_Number-Fk/Fkdif
       Iter = Iter+1
       goto 55
 65    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave
       WnumX = Wave_Number*cos(Theta_Wave*pi/180.)
       WnumY = Wave_Number*sin(Theta_Wave*pi/180.)       

       Source_Area = 0.0    
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           Phs_lag = (j-Jbeg)*dy*WnumY
           SourceC(i,j) = Celerity*Amp_Wave/Source_Area*cos(pi/2-Segma*time+Phs_lag)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for random waves 
     if(WaveMaker(1:7)=='INT_SPC') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
           do ndir = 1,NumDir
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
   75        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 75
   85        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! adjust wave direction for periodic bc
             Angle = Dire(ndir)*pi/180.
             if(Angle>zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3<Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = 0.5*pi-small
                 else
                   tmp3 = asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             else
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3>Angle)
                 nk = nk+1
                tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = -0.5*pi+small
                 else
                   tmp3 = -asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

             ! calculate root-mean-squre wave height for each component
             if(nfreq==1) then
               dfreq = Freq(2)-Freq(1)
             elseif(nfreq==NumFreq) then
               dfreq = Freq(NumFreq)-Freq(NumFreq-1)
             else
               dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
             endif
             dfreq = abs(dfreq)

             if(ndir==1) then
               ddir = Dire(2)-Dire(1)
             elseif(ndir==NumDir) then
               ddir = Dire(NumDir)-Dire(NumDir-1)
             else
               ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
             endif
             ddir = abs(ddir)
         
             Amp_Wave = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*ddir*dfreq)

             Phs_lag = (dy/2.0+(j-Jbeg)*dy)*WnumY
             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
           enddo
           enddo
         endif
       enddo
       enddo
     endif
	 
	 ! Added by Cheng Zhang, internal wavemaker for irregular waves 
     if(WaveMaker(1:7)=='INT_IRR') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
           do ndir = 1,NumDir
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
   95        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 105
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 95
   105        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! adjust wave direction for periodic bc
             Angle = Dire(ndir)*pi/180.
             if(Angle>zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3<Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = 0.5*pi-small
                 else
                   tmp3 = asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             elseif (Angle<zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3>Angle)
                 nk = nk+1
                tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = -0.5*pi+small
                 else
                   tmp3 = -asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
			 else
			   Angle = 0.0
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

             ! wave height for each component
         
             Amp_Wave = 2.0*Wave_Spc2d(ndir,nfreq)
             Phs_lag = (dy/2.0+(j-Jbeg)*dy)*WnumY
			 
			 if(nfreq==1) then
			   SourceC(i,j) = SourceC(i,j)+2*Celerity/Source_Area*  &
                  (Amp_Wave/2*cos(pi/2-Segma*time+Phs_lag+  &
				  Random_Phs(ndir,nfreq))+eta_mean)
			 else
               SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
			 endif
           enddo
           enddo
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for jonswap spectrum
     if((WaveMaker(1:7)=='INT_JON').or.(WaveMaker(1:7)=='INT_TMA')) then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)
# endif

       dfreq = (Freq_Max-Freq_Min)/float(NumFreq)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
 76          Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 76
 86          continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! root-mean-square wave height
             Amp_Wave = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)

             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                cos(pi/2-Segma*time+RanPhs(nfreq))
           enddo
         endif
       enddo
       enddo
     endif


     ! internal wavemaker for cnoidal wave
     if(WaveMaker(1:7)=='INT_CON') then
       call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)

       ! wave number
       Wave_Number = 2.0*pi/Wave_Length

!# if defined(PARALLEL)
!       if(myid.eq.0) write(*,*) 'Mod=',Mod1,'Ytrough=',Ytrough, &
!            'Wave_Number=',wave_number
!# endif      

       ! find zero start
       Zup = 1.0
       Zlow = 0.0
       Zmid= (Zup+Zlow)/2.0
       nn = 0
 200   nn = nn+1
       Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

       if(abs(Zero1)<=1.0e-6) goto 210
       if(nn>1000) then
         write(*,*)'too many iterations; stop'
         stop
       endif
       if(Zero1<0.0) then
         Zup = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       else
         Zlow = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       endif
 210   continue
       Xstart = Zmid

       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*(Ytrough+Amp_Wave*cnoidal_cn(  &
               Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for solitary wave
     if(WaveMaker(1:7)=='INT_SOL') then
       Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*  &
               Amp_Wave/cosh(Atmp*(Xstart-Celerity*TIME))**2
         endif
       enddo
       enddo
     endif

     ! source terms for momentum eqs.
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SourceX(i,j) = Grav*Eta(i,j)*DelxH(i,j)*Mask(i,j)
       SourceY(i,j) = Grav*Eta(i,j)*DelyH(i,j)*Mask(i,j)
     enddo
     enddo

     end subroutine source_terms
 

     subroutine wave_average
!---------------------------------------------------
!    Estimate wave averaged quantities
!    Called by                        
!       main 
!    Last update: 13/01/2012, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP), dimension(:,:), allocatable :: U_Dep_Ave,V_Dep_Ave,S_Dep_Ave
     integer :: i,j,k,n
     real(SP) :: Tmp,Tmp_0,Dz,Zk,U_zk,V_zk,W_zk,S_zk,Zbeg,Zend,Zn,Zn1,Sinterp

     allocate(U_Dep_Ave(Mloc,Nloc))
     allocate(V_Dep_Ave(Mloc,Nloc))
# if defined (SEDIMENT)
     allocate(S_Dep_Ave(Mloc,Nloc))
# endif

     if(TIME>Wave_Ave_Start.and.TIME<=Wave_Ave_End) then
       ! Lagrangian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Lag_Umean(i,j,k) = Lag_Umean(i,j,k)+U(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Vmean(i,j,k) = Lag_Vmean(i,j,k)+V(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Wmean(i,j,k) = Lag_Wmean(i,j,k)+W(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
       enddo
       enddo
       enddo

       ! Eulerian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Dz = Hc(i,j)/float(Kglob)
         Zk = (k-Kbeg)*Dz+Dz/2.0

         U_zk = Zero; V_zk = Zero; W_zk = Zero
         Zbeg = sigc(Kbeg)*(Eta(i,j)+Hc(i,j))
         Zend = sigc(Kend)*(Eta(i,j)+Hc(i,j))
         if(Zk<=Zbeg) then
           U_zk = U(i,j,Kbeg)
           V_zk = V(i,j,Kbeg)
           W_zk = W(i,j,Kbeg)
# if defined (SEDIMENT)
           S_zk = Conc(i,j,Kbeg)
# endif
         elseif(Zk>=Zend) then
           U_zk = U(i,j,Kend)
           V_zk = V(i,j,Kend)
           W_zk = W(i,j,Kend)
# if defined (SEDIMENT)
           S_zk = Conc(i,j,Kend)
# endif
         else
           do n = Kbeg,Kend-1
             Zn = sigc(n)*(Eta(i,j)+Hc(i,j)) 
             Zn1 = sigc(n+1)*(Eta(i,j)+Hc(i,j))
             if(Zk>=Zn.and.Zk<Zn1) then
               Sinterp = (Zk-Zn)/(Zn1-Zn)
               U_zk = U(i,j,n)*(1.0-Sinterp)+U(i,j,n+1)*Sinterp
               V_zk = V(i,j,n)*(1.0-Sinterp)+V(i,j,n+1)*Sinterp
               W_zk = W(i,j,n)*(1.0-Sinterp)+W(i,j,n+1)*Sinterp
# if defined (SEDIMENT)
               S_zk = Conc(i,j,n)*(1.0-Sinterp)+Conc(i,j,n+1)*Sinterp
# endif
             endif 
           enddo
         endif
         Euler_Umean(i,j,k) = Euler_Umean(i,j,k)+U_zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Vmean(i,j,k) = Euler_Vmean(i,j,k)+V_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Wmean(i,j,k) = Euler_Wmean(i,j,k)+W_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
# if defined (SEDIMENT)
         Euler_Smean(i,j,k) = Euler_Smean(i,j,k)+S_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
# endif        
       enddo
       enddo
       enddo         

       ! depth-averaged velocity
       U_Dep_Ave = Zero
       V_Dep_Ave = Zero
# if defined (SEDIMENT)
       S_Dep_Ave = Zero
# endif
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         do k = Kbeg,Kend
           U_Dep_Ave(i,j) = U_Dep_Ave(i,j)+U(i,j,k)/float(Kend-Kbeg+1)
           V_Dep_Ave(i,j) = V_Dep_Ave(i,j)+V(i,j,k)/float(Kend-Kbeg+1)
# if defined (SEDIMENT)
           S_Dep_Ave(i,j) = S_Dep_Ave(i,j)+Conc(i,j,k)/float(Kend-Kbeg+1)
# endif
         enddo
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Setup(i,j) = Setup(i,j)+Eta(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Umean(i,j) = Umean(i,j)+U_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Vmean(i,j) = Vmean(i,j)+V_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
# if defined (SEDIMENT)
         Smean(i,j) = Smean(i,j)+S_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
# endif

         if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
         if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
         
         Tmp = Eta(i,j)
         Tmp_0 = Eta0(i,j)
         if(Tmp>Tmp_0.and.Tmp*Tmp_0<=Zero) then
           Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
           if(Num_Zero_Up(i,j)>=2) then
             if(WaveheightID==1) then  ! Average wave height
               WaveHeight(i,j) = WaveHeight(i,j)+Emax(i,j)-Emin(i,j)
             elseif(WaveheightID==2) then  ! RMS wave height
               WaveHeight(i,j) = WaveHeight(i,j)+(Emax(i,j)-Emin(i,j))**2
             endif
           endif

           ! reset Emax and Emin to find next wave
           Emax(i,j) = -1000.
           Emin(i,j) = 1000.
         endif  
       enddo
       enddo
     endif

     deallocate(U_Dep_Ave)
     deallocate(V_Dep_Ave)
# if defined (SEDIMENT)
     deallocate(S_Dep_Ave)
# endif

     end subroutine wave_average

     
     subroutine print_wh_setup
!---------------------------------------------------
!    Estimate wave averaged quantities
!    Called by
!       main 
!    Last update: 13/01/2012, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: i,j
     character(len=80) :: FDIR,file

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Num_Zero_Up(i,j)>=2) then
         if(WaveheightID==1) then
           WaveHeight(i,j) = WaveHeight(i,j)/float(Num_Zero_Up(i,j)-1)
         elseif(WaveheightID==2) then
           WaveHeight(i,j) = sqrt(WaveHeight(i,j))/float(Num_Zero_Up(i,j)-1)
         endif
       else
         WaveHeight(i,j) = Zero
       endif
     enddo
     enddo

     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

     file = TRIM(FDIR)//'waveheight'
     call putfile2D(file,WaveHeight)

     file = TRIM(FDIR)//'setup'
     call putfile2D(file,Setup)
     
     file = TRIM(FDIR)//'umean'
     call putfile2D(file,Umean)

     file = TRIM(FDIR)//'vmean'
     call putfile2D(file,Vmean)

# if defined (SEDIMENT)
     file = TRIM(FDIR)//'smean'
     call putfile2D(file,Smean)
# endif

     file = TRIM(FDIR)//'lag_umean'
     call putfile3D(file,Lag_Umean)

     file = TRIM(FDIR)//'lag_vmean'
     call putfile3D(file,Lag_Vmean)

     file = TRIM(FDIR)//'lag_wmean'
     call putfile3D(file,Lag_Wmean)

     file = TRIM(FDIR)//'euler_umean'
     call putfile3D(file,Euler_Umean)

     file = TRIM(FDIR)//'euler_vmean'
     call putfile3D(file,Euler_Vmean)

     file = TRIM(FDIR)//'euler_wmean'
     call putfile3D(file,Euler_Wmean)

# if defined (SEDIMENT)
     file = TRIM(FDIR)//'euler_smean'
     call putfile3D(file,Euler_Smean)
# endif

     end subroutine print_wh_setup


     subroutine statistics
!---------------------------------------------------
!    This subroutine is used to show statistics
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP) :: MassVolume,CellMass,Energy,MaxEta,MinEta,MaxU, &
                 MaxV,MaxW,MaxS,MinS
     integer :: i,j,k
# if defined (PARALLEL)
     real(SP) :: myvar
# endif

     ! Vol = sum(D*dx*dy)
     ! Energy = sum(m*g*h+0.5*m*u^2), reference is at z = 0
     MassVolume = Zero
     Energy = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       MassVolume = MassVolume+D(i,j)*dx*dy
       do k = Kbeg,Kend
         CellMass = Rho0*dsig(k)*D(i,j)*dx*dy
         Energy = Energy+CellMass*Grav*(D(i,j)*sigc(k)-Hc(i,j))+  &
                    0.5*CellMass*(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
       enddo
     enddo
     enddo

     MaxEta = MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta = MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU = MAXVAL(abs(U(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxV = MAXVAL(abs(V(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxW = MAXVAL(abs(W(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend))) 
# if defined (SALINITY)
     MaxS = MAXVAL(abs(Sali(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MinS = MINVAL(abs(Sali(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
# endif


# if defined (PARALLEL)
     call MPI_ALLREDUCE(MassVolume,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)        
     MassVolume = myvar
     call MPI_ALLREDUCE(Energy,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)            
     Energy = myvar
     call MPI_ALLREDUCE(MaxEta,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)            
     MaxEta = myvar
     call MPI_ALLREDUCE(MinEta,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)            
     MinEta = myvar
     call MPI_ALLREDUCE(MaxU,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxU = myvar
     call MPI_ALLREDUCE(MaxV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxV = myvar
     call MPI_ALLREDUCE(MaxW,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxW = myvar
# if defined (SALINITY)
     call MPI_ALLREDUCE(MaxS,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxS = myvar
     call MPI_ALLREDUCE(MinS,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)                                                 
     MinS = myvar
# endif
# endif

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     ! print screen
     WRITE(*,*),'----------------- STATISTICS ----------------'
     WRITE(*,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(*,102) TIME,dt,TRIM(dt_constraint)
     WRITE(*,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(*,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS

     ! print log file 
     WRITE(3,*),'----------------- STATISTICS ----------------'
     WRITE(3,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(3,102) TIME,dt,TRIM(dt_constraint)
     WRITE(3,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(3,101), MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS
# if defined (PARALLEL)
     endif
# endif

101  FORMAT(10E12.4)
102  FORMAT(2E12.4,A8)
103  FORMAT(A97)
 
     end subroutine statistics

 
     subroutine probes
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 16/11/2011, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: n,iu,i,j,k
     character(len=80) :: STAT_FILE,FDIR,FILE_NUM
     real(SP) :: zlev1(Kloc),zlev2(Kloc),ulev,vlev,wlev,plev,  &
                 klev,sinterp

     FDIR = TRIM(RESULT_FOLDER)
     
     do n = 1,NSTAT
       iu = 100+n
       write(FILE_NUM(1:4),'(I4.4)') n
       STAT_FILE = TRIM(FDIR)//'probe_'//TRIM(FILE_NUM)
       open(iu,file=TRIM(STAT_FILE),access='APPEND')
       
       if(zstat(n).eq.0.0) then
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xstat(n)>=x(i).and.xstat(n)<x(i+1).and.  &
              ystat(n)>=y(j).and.ystat(n)<y(j+1)) then
             write(iu,'(100E12.4)') time,eta(i,j)
           endif
         enddo
         enddo
         close(iu)
	   ! added by Cheng to get the velocity at each layer
	   elseif(zstat(n).eq.-1.0) then
	     do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xstat(n)>=x(i).and.xstat(n)<=x(i+1).and.  &
              ystat(n)>=y(j).and.ystat(n)<=y(j+1)) then
             write(iu,'(100E14.6)') time,eta(i,j),(u(i,j,k),v(i,j,k),w(i,j,k),k=Kbeg,Kend)
           endif
         enddo
         enddo
		 close(iu)
       else
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xstat(n)>=x(i).and.xstat(n)<x(i+1).and.  &
              ystat(n)>=y(j).and.ystat(n)<y(j+1)) then

             do k = 1,Kloc
               zlev1(k) = sig(k)*D(i,j)
               zlev2(k) = sigc(k)*D(i,j)
             enddo

             do k = Kbeg-1,Kend
               if(zstat(n)>=zlev1(k).and.zstat(n)<zlev1(k+1)) then
                 sinterp = (zstat(n)-zlev1(k))/(zlev1(k+1)-zlev1(k))
                 plev = (1.0-sinterp)*P(i,j,k)+sinterp*P(i,j,k+1)
               endif

               if(zstat(n)>=zlev2(k).and.zstat(n)<zlev2(k+1)) then
                 sinterp = (zstat(n)-zlev2(k))/(zlev2(k+1)-zlev2(k))
                 ulev = (1.0-sinterp)*U(i,j,k)+sinterp*U(i,j,k+1)
                 vlev = (1.0-sinterp)*V(i,j,k)+sinterp*V(i,j,k+1)
                 wlev = (1.0-sinterp)*W(i,j,k)+sinterp*W(i,j,k+1)
                 klev = (1.0-sinterp)*Tke(i,j,k)+sinterp*Tke(i,j,k+1)
               endif
             enddo

             write(iu,'(100E12.4)') time,eta(i,j),ulev,vlev,wlev,plev,klev
           endif
         enddo
         enddo
         close(iu)
       endif
     enddo

     end subroutine probes


     subroutine preview
!--------------------------------------------------- 
!    This subroutine is used to preview
!    Called by                         
!       main 
!    Last update: 23/12/2010, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,I1,I2,I3,I4,I5
     character(len=80) :: FDIR=''
     character(len=80) :: FILE_NAME=''
     character(len=80) :: file=''

     ! file number
     Icount = Icount+1
   
     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     if(myid.eq.0) write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME     
# else
     write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
# endif
102  FORMAT(A20,I5,A14,F8.3,A2,F8.3)
100  FORMAT(5000E16.6)

     I1 = mod(Icount/10000,10)
     I2 = mod(Icount/1000,10)
     I3 = mod(Icount/100,10)
     I4 = mod(Icount/10,10)
     I5 = mod(Icount,10)

     write(FILE_NAME(1:1),'(I1)') I1
     write(FILE_NAME(2:2),'(I1)') I2
     write(FILE_NAME(3:3),'(I1)') I3
     write(FILE_NAME(4:4),'(I1)') I4
     write(FILE_NAME(5:5),'(I1)') I5

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     open(5,file=TRIM(FDIR)//'time',position="append")
     write(5,*) TIME
     close(5)
# if defined (PARALLEL)
     endif
# endif

     if(Icount==1) then
       if(OUT_H) then
         file=TRIM(FDIR)//'depth'
         call putfile2D(file,Hc0)
       endif
     endif

     if(OUT_E) then
       file=TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
       call putfile2D(file,Eta)
     endif

     if(OUT_U) then
       file=TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
       call putfile3D(file,U)
     endif

     if(OUT_V) then
       file=TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
       call putfile3D(file,V)
     endif

     if(OUT_W) then
       file=TRIM(FDIR)//'w_'//TRIM(FILE_NAME)
       call putfile3D(file,W)
     endif

     if(OUT_P) then
       file=TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
       call putfile3D(file,P)
     endif

     if(OUT_K) then
       file=TRIM(FDIR)//'k_'//TRIM(FILE_NAME)
       call putfile3D(file,Tke)
# if defined (VEGETATION)
       file=TRIM(FDIR)//'vk_'//TRIM(FILE_NAME)
       call putfile3D(file,Tke_w)
# endif
     endif

     if(OUT_D) then
       file=TRIM(FDIR)//'d_'//TRIM(FILE_NAME)
       call putfile3D(file,Eps)
     endif

     if(OUT_S) then
       file=TRIM(FDIR)//'s_'//TRIM(FILE_NAME)
       call putfile3D(file,Prod_s)
     endif

     if(OUT_C) then
       file=TRIM(FDIR)//'c_'//TRIM(FILE_NAME)
       call putfile3D(file,CmuVt)
     endif

     if(OUT_A) then
       file=TRIM(FDIR)//'upwp_'//TRIM(FILE_NAME)
       call putfile3D(file,UpWp)
     endif
	 
	 !added by Cheng for varying depth
     if(OUT_Z) then
       file=TRIM(FDIR)//'depth_'//TRIM(FILE_NAME)
       call putfile2D(file,Hc)
     endif
	 
	 !added by Cheng for recording Hmax
     if(OUT_M) then
       file=TRIM(FDIR)//'hmax_'//TRIM(FILE_NAME)
       call putfile2D(file,HeightMax)
     endif

# if defined (VEGETATION)
     if(TRIM(Veg_Type)=='FLEXIBLE') then
       file=TRIM(FDIR)//'vegh_'//TRIM(FILE_NAME)
       call putfile2D(file,FVegH)
     endif
# endif

# if defined (BUBBLE)
     if(OUT_B) then
       file=TRIM(FDIR)//'b_'//TRIM(FILE_NAME)
       call putfile3D(file,Vbg)
     endif
# endif

# if defined (SEDIMENT)
     if(OUT_F) then
       file=TRIM(FDIR)//'f_'//TRIM(FILE_NAME)
       call putfile3D(file,Conc)
     endif
     if(OUT_T) then
       file=TRIM(FDIR)//'t_'//TRIM(FILE_NAME)
       call putfile2D(file,Taub)
     endif
!     if(OUT_G) then
!       file=TRIM(FDIR)//'g_'//TRIM(FILE_NAME)
!       call putfile2D(file,Bed)
!       file=TRIM(FDIR)//'qbed_'//TRIM(FILE_NAME)
!       call putfile2D(file,Qbedx)
!     endif
# endif

# if defined (SALINITY)
     if(OUT_I) then
       file=TRIM(FDIR)//'i_'//TRIM(FILE_NAME)
       call putfile3D(file,Sali)
     endif
# endif

# if defined (POROUSMEDIA)
     file=TRIM(FDIR)//'por_'//TRIM(FILE_NAME)
     call putfile3D(file,Porosity)
# endif

# if defined (BALANCE2D)
     file=TRIM(FDIR)//'dudt2d_'//TRIM(FILE_NAME)
     call putfile2D(file,DUDT2D)
     file=TRIM(FDIR)//'dedx2d_'//TRIM(FILE_NAME)
     call putfile2D(file,DEDX2D)
     file=TRIM(FDIR)//'dpdx2d_'//TRIM(FILE_NAME)
     call putfile2D(file,DPDX2D)
     file=TRIM(FDIR)//'diffx2d_'//TRIM(FILE_NAME)
     call putfile2D(file,DIFFX2D)
     file=TRIM(FDIR)//'taubx2d_'//TRIM(FILE_NAME)
     call putfile2D(file,TAUBX2D)
# if defined (VEGETATION)
     file=TRIM(FDIR)//'fvegx2d_'//TRIM(FILE_NAME)
     call putfile2D(file,FVEGX2D)
# endif
# endif

# if defined (TWOLAYERSLIDE)
     file=TRIM(FDIR)//'Hc_'//TRIM(FILE_NAME)
     call putfile2D(file,Hc)
     file=TRIM(FDIR)//'Ha_'//TRIM(FILE_NAME)
     call putfile2D(file,Ha)
     file=TRIM(FDIR)//'Ua_'//TRIM(FILE_NAME)
     call putfile2D(file,Ua)
     file=TRIM(FDIR)//'Va_'//TRIM(FILE_NAME)
     call putfile2D(file,Va)
# endif

     end subroutine preview


# if defined (PARALLEL)
    subroutine putfile2D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(Mglob,Nglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,iglob,jglob,len,n

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    ! put the data in master processor into the global var
    if(myid==0) then
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        phiglob(iglob,jglob) = Phi(i,j)
      enddo
      enddo
    endif

    ! collect data from other processors into the master processor
    len = Mloc*Nloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          phiglob(iglob,jglob) = xx(i,j)
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo       

    if(myid==0) then
      open(5,file=TRIM(file))
      do j = 1,Nglob
        write(5,100) (phiglob(i,j),i=1,Mglob)
      enddo
      close(5)
    endif
100 FORMAT(5000f15.6)

    end subroutine putfile2D


    subroutine putfile3D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,philoc
    real(SP),dimension(Mglob,Nglob,Kglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,k,jk,iglob,jglob,kk,n,len,nreq,NKloc

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)

    NKloc = Nloc*Kloc

    ! put the data in master processor into the global var
    if(myid==0) then
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        kk = k-Nghost
        phiglob(iglob,jglob,kk) = Phi(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(philoc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      philoc(i,jk) = phi(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = Kbeg,Kend
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          kk = k-Nghost
          jk = (k-1)*Nloc+j
          phiglob(iglob,jglob,kk) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(philoc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do k = 1,Kglob
      do j = 1,Nglob
        write(5,100) (phiglob(i,j,k),i=1,Mglob)
      enddo
      enddo
      close(5)
    endif
100 FORMAT(5000f15.6)

    deallocate(philoc)
    deallocate(xx)

    end subroutine putfile3D
# else
     subroutine putfile2D(file,phi)
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j

     open(5,file=trim(file))
     do j = Jbeg,Jend
       write(5,100) (phi(i,j),i=Ibeg,Iend)
     enddo
     close(5)
 100 FORMAT(5000f15.6)
     end subroutine putfile2D

     subroutine putfile3D(file,phi)
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j,k

     open(5,file=trim(file))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       write(5,100) (phi(i,j,k),i=Ibeg,Iend)
     enddo
     enddo
     close(5)
 100 FORMAT(5000f15.6)
     end subroutine putfile3D  
# endif

 
     subroutine estimate_dt
!----------------------------------------------------
!    This subroutine is used to estimate dt
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: tmp1,tmp2,dxonu,dyonv,dzonw,dt_growth,dt_courant,dt_viscous
# if defined (PARALLEL)
     real(SP) :: myvar
# endif
     ! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     real(SP) :: dt_landslide
# endif

     ! save previous time step
     dt_old = dt
     dt_growth = 1.05*dt_old     

     tmp2 = Large
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       tmp1 = abs(U(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(V(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv

       tmp1 = max(abs(W(i,j,k)),Small)
       dzonw = dsig(k)*D(i,j)/tmp1
       if(dzonw<tmp2) tmp2=dzonw
     enddo
     enddo
     enddo

# if defined (TWOLAYERSLIDE)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       tmp1 = abs(Ua(i,j))+sqrt(Grav_Lz(i,j)*Ha(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(Va(i,j))+sqrt(Grav_Lz(i,j)*Ha(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv
     enddo
     enddo
# endif

# if defined (PARALLEL)
     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
# endif
     dt_courant = CFL*tmp2
	 
	 ! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     ! time step limit due to the landslide 
     tmp2 = Large
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       tmp1 = abs(SldUs(i,j))+sqrt(Grav*SldDZs(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(SldVs(i,j))+sqrt(Grav*SldDZs(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv
     enddo
     enddo
# if defined (PARALLEL)
     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
# endif
     dt_landslide = 0.1d0*CFL*tmp2
# endif

     ! time step limit due to explicit viscous stress terms
     dt_viscous = Large
     if(VISCOUS_FLOW) then
       tmp2 = Large
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         tmp1 = dx**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1

         tmp1 = dy**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1
       enddo
       enddo
       enddo
# if defined (PARALLEL)
       call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
       tmp2 = myvar
# endif
       dt_viscous = VISCOUS_NUMBER*tmp2
     endif 

     ! get dt    
     dt = min(dt_growth,dt_courant,dt_viscous,dt_max)
     if(dt<dt_min) then
# if defined (PARALLEL)
       if(myid.eq.0) then
         write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
         stop
       endif
# else
       write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
       stop
# endif
     endif
     TIME = TIME+dt
     RUN_STEP = RUN_STEP+1 
# if defined (PARALLEL)
     if(myid.eq.0) write(3,*) RUN_STEP,dt,TIME
# else 
     write(3,*) RUN_STEP,dt,TIME
# endif

     if(dt==dt_growth) then
       dt_constraint = 'GROWTH'
     elseif(dt==dt_courant) then
       dt_constraint = 'COURANT'
     elseif(dt==dt_viscous) then
       dt_constraint = 'VISCOUS'
     elseif(dt==dt_max) then
       dt_constraint = 'MAXIMUM'
	 ! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     elseif(dt==dt_landslide) then
       dt_constraint = 'LANDSLD'	
# endif
     endif  

     end subroutine estimate_dt


     subroutine vel_bc_old
!----------------------------------------------------
!    Boundary conditions for velocity
!    Called by 
!       main and get_UVW
!    Last update: 01/02/2011, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,imask
     real(SP) :: Wtop,Wbot,Cdrag,Phi,Dz1,Cg
     real(SP), dimension(:,:,:), allocatable :: DelxUzL,DelxVzL,DelxWzL,&                               
          DelxUzR,DelxVzR,DelxWzR,DelyUzL,DelyVzL,DelyWzL,DelyUzR,DelyVzR,DelyWzR


     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by Cheng for nesting. Please search for others with (COUPLING) in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_X0==1) then  ! free-slip wall
         do i = 1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==2) then ! no-slip wall
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = -V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = -W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = -DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = -DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==3) then ! inflow and outflow
         if(WaveMaker(1:7)=='LEF_TID') then ! for long-wave
           do i =1,Nghost
             U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
             V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
             W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
             DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
             DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
             DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
           enddo
         else
           do i = 1,Nghost
             U(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)-U(Ibeg+i-1,j,k)
             V(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)-V(Ibeg+i-1,j,k)
             W(Ibeg-i,j,k) = 2.0*Win_X0(j,k)-W(Ibeg+i-1,j,k)
             DU(Ibeg-i,j,k) = 2.0*Din_X0(j)*Uin_X0(j,k)-DU(Ibeg+i-1,j,k)
             DV(Ibeg-i,j,k) = 2.0*Din_X0(j)*Vin_X0(j,k)-DV(Ibeg+i-1,j,k)
             DW(Ibeg-i,j,k) = 2.0*Din_X0(j)*Win_X0(j,k)-DW(Ibeg+i-1,j,k)
           enddo
         endif
       elseif(Bc_X0==4) then
         do i =1,Nghost
           U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dx/2.0/Zob))**2
         endif
         Phi = dx*Cdrag*sqrt(V(Ibeg,j,k)**2+W(Ibeg,j,k)**2)/(Cmu(Ibeg,j,k)+CmuHt(Ibeg,j,k))
         Phi = dmin1(Phi,2.0)
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = (1.0-Phi)*V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = (1.0-Phi)*W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = D(i,j)*U(Ibeg-i,j,k)
           DV(Ibeg-i,j,k) = D(i,j)*V(Ibeg-i,j,k)
           DW(Ibeg-i,j,k) = D(i,j)*W(Ibeg-i,j,k)
         enddo
       elseif(Bc_X0==8) then  ! specify u,v,w at ghost cells
         do i = 1,Nghost
           U(Ibeg-i,j,k) = Uin_X0(j,k)
           V(Ibeg-i,j,k) = Vin_X0(j,k)
           W(Ibeg-i,j,k) = Win_X0(j,k)
           DU(Ibeg-i,j,k) = Din_X0(j)*Uin_X0(j,k)
           DV(Ibeg-i,j,k) = Din_X0(j)*Vin_X0(j,k)
           DW(Ibeg-i,j,k) = Din_X0(j)*Win_X0(j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_Xn==1) then  ! free-slip wall 
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==2) then ! no-slip wall
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = -V(Iend-i+1,j,k)
           W(Iend+i,j,k) = -W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = -DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = -DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==3) then
         do i = 1,Nghost
           U(Iend+i,j,k) = 2.0*Uin_Xn(j,k)-U(Iend-i+1,j,k)
           V(Iend+i,j,k) = 2.0*Vin_Xn(j,k)-V(Iend-i+1,j,k)
           W(Iend+i,j,k) = 2.0*Win_Xn(j,k)-W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = 2.0*Din_Xn(j)*Uin_Xn(j,k)-DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = 2.0*Din_Xn(j)*Vin_Xn(j,k)-DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = 2.0*Din_Xn(j)*Win_Xn(j,k)-DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==4) then 
         do i = 1,Nghost
           U(Iend+i,j,k) = U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dx/2.0/Zob))**2
         endif
         Phi = dx*Cdrag*sqrt(V(Iend,j,k)**2+W(Iend,j,k)**2)/(Cmu(Iend,j,k)+CmuHt(Iend,j,k))
         Phi = dmin1(Phi,2.0)
         do i =1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = (1.0-Phi)*V(Iend-i+1,j,k)
           W(Iend+i,j,k) = (1.0-Phi)*W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = D(i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(i,j)*W(Iend+i,j,k)
         enddo
       elseif(Bc_Xn==6) then
         do i = 1,Nghost
!           Cg = -(U0(Iend+i-1,j,k)-U00(Iend+i-1,j,k)+1.e-16)/  &
!                 (U00(Iend+i-1,j,k)-U00(Iend+i-2,j,k)+1.e-16)
           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Cg = max(min(Cg,1.0),0.0)
           U(Iend+i,j,k) = Cg*U0(Iend+i-1,j,k)+(1.0-Cg)*U0(Iend+i,j,k)

!           Cg =-(V0(Iend+i-1,j,k)-V00(Iend+i-1,j,k)+1.e-16)/  &
!                 (V00(Iend+i-1,j,k)-V00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           V(Iend+i,j,k) = Cg*V0(Iend+i-1,j,k)+(1.0-Cg)*V0(Iend+i,j,k)
           
!           Cg =-(W0(Iend+i-1,j,k)-W00(Iend+i-1,j,k)+1.e-16)/  &
!                 (W00(Iend+i-1,j,k)-W00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           W(Iend+i,j,k) = Cg*W0(Iend+i-1,j,k)+(1.0-Cg)*W0(Iend+i,j,k)

           DU(Iend+i,j,k) = D(Iend+i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(Iend+i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(Iend+i,j)*W(Iend+i,j,k)
         enddo 
       elseif(Bc_Xn==8) then
         do i = 1,Nghost
           U(Iend+i,j,k) = Uin_Xn(j,k)
           V(Iend+i,j,k) = Vin_Xn(j,k)
           W(Iend+i,j,k) = Win_Xn(j,k)
           DU(Iend+i,j,k) = Din_Xn(j)*Uin_Xn(j,k)
           DV(Iend+i,j,k) = Din_Xn(j)*Vin_Xn(j,k)
           DW(Iend+i,j,k) = Din_Xn(j)*Win_Xn(j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Y0==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = -U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = -W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = -DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = -DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==4) then
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dy/2.0/Zob))**2
         endif
         Phi = dy*Cdrag*sqrt(U(i,Jbeg,k)**2+W(i,Jbeg,k)**2)/(Cmu(i,Jbeg,k)+CmuHt(i,Jbeg,k))
         Phi = dmin1(Phi,2.0)
         do j =1,Nghost
           U(i,Jbeg-j,k) = (1.0-Phi)*U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = (1.0-Phi)*W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = D(i,j)*U(i,Jbeg-j,k)
           DV(i,Jbeg-j,k) = D(i,j)*V(i,Jbeg-j,k)
           DW(i,Jbeg-j,k) = D(i,j)*W(i,Jbeg-j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Yn==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = -U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = -W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = -DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = -DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==4) then
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dy/2.0/Zob))**2
         endif
         Phi = dy*Cdrag*sqrt(U(i,Jend,k)**2+W(i,Jend,k)**2)/(Cmu(i,Jend,k)+CmuHt(i,Jend,k))
         Phi = dmin1(Phi,2.0)
         do j =1,Nghost
           U(i,Jend+j,k) = (1.0-Phi)*U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jend+j,k) = (1.0-Phi)*W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = D(i,j)*U(i,Jend+j,k)
           DV(i,Jend+j,k) = D(i,j)*V(i,Jend+j,k)
           DW(i,Jend+j,k) = D(i,j)*W(i,Jend+j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     allocate(DelxUzL(Mloc,Nloc,Kloc1))
     allocate(DelxVzL(Mloc,Nloc,Kloc1))
     allocate(DelxWzL(Mloc,Nloc,Kloc1))
     allocate(DelxUzR(Mloc,Nloc,Kloc1))
     allocate(DelxVzR(Mloc,Nloc,Kloc1))
     allocate(DelxWzR(Mloc,Nloc,Kloc1))
     allocate(DelyUzL(Mloc,Nloc,Kloc1))
     allocate(DelyVzL(Mloc,Nloc,Kloc1))
     allocate(DelyWzL(Mloc,Nloc,Kloc1))
     allocate(DelyUzR(Mloc,Nloc,Kloc1))
     allocate(DelyVzR(Mloc,Nloc,Kloc1))
     allocate(DelyWzR(Mloc,Nloc,Kloc1))

     ! top and bottom
     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)
     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)
     
     call Kbc_Surface
     call Kbc_Bottom

     call delxFun1_3d(UzL,DelxUzL)
     call delxFun1_3d(VzL,DelxVzL)
     call delxFun1_3d(WzL,DelxWzL)
     call delxFun1_3d(UzR,DelxUzR)
     call delxFun1_3d(VzR,DelxVzR)
     call delxFun1_3d(WzR,DelxWzR)
     call delyFun1_3d(UzL,DelyUzL)
     call delyFun1_3d(VzL,DelyVzL)
     call delyFun1_3d(WzL,DelyWzL)
     call delyFun1_3d(UzR,DelyUzR)
     call delyFun1_3d(VzR,DelyVzR)
     call delyFun1_3d(WzR,DelyWzR)

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Dz1 = 0.5*D(i,j)*dsig(Kbeg)
       if(ibot==1) then
         Cdrag = Cd0
       else
# if defined (SEDIMENT)
         Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                          
# else
         Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
       endif
       Phi = 2.0*Dz1*Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)/(Cmu(i,j,Kbeg)+CmuVt(i,j,Kbeg))
       Phi = dmin1(Phi,2.0)

       if(Bc_Z0==1) then  ! free-slip
         !Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*WzR(i,j,Kbeg)-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==2) then  ! no-slip
         !Wbot = -DeltH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = -U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = -V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*WzR(i,j,Kbeg)-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==5) then
         do k = 1,Nghost
           U(i,j,Kbeg-k) = (1.0-Phi)*U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = (1.0-Phi)*V(i,j,Kbeg+k-1)
           !Wbot = -DeltH(i,j)-0.5*(U(i,j,Kbeg)+U(i,j,Kbeg-1))*DelxH(i,j)-  &
           !         0.5*(V(i,j,Kbeg)+V(i,j,Kbeg-1))*DelyH(i,j)
           W(i,j,Kbeg-k) = 2.0*WzR(i,j,Kbeg)-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       endif

       ! at the surface (no stress)
       ! Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j)
       do k = 1,Nghost
         W(i,j,Kend+k) = 2.0*WzL(i,j,Kend1)-W(i,j,Kend-k+1)
         U(i,j,Kend+k) = U(i,j,Kend-k+1)
         V(i,j,Kend+k) = V(i,j,Kend-k+1)
         DU(i,j,Kend+k) = D(i,j)*U(i,j,Kend+k)
         DV(i,j,Kend+k) = D(i,j)*V(i,j,Kend+k)
         DW(i,j,Kend+k) = D(i,j)*W(i,j,Kend+k)
       enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks 02/15/2013
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(i,j)==0) THEN
         ! south boundary 
         IF(Mask(i,j+1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j-imask+1,k) = U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = W(i,j+imask,k)
               DU(i,j-imask+1,k) = DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = DW(i,j+imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j-imask+1,k) = -U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = -W(i,j+imask,k)
               DU(i,j-imask+1,k) = -DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = -DW(i,j+imask,k)
             enddo
           endif
         ! north  
         ELSEIF(Mask(i,j-1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j+imask-1,k) = U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = W(i,j-imask,k)
               DU(i,j+imask-1,k) = DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = DW(i,j-imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j+imask-1,k) = -U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = -W(i,j-imask,k)
               DU(i,j+imask-1,k) = -DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = -DW(i,j-imask,k)
             enddo
           endif
         ! west
         ELSEIF(Mask(i+1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = V(I+imask,j,k)
               W(I-imask+1,j,k) = W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = DV(I+imask,j,k)
               DW(I-imask+1,j,k) = DW(I+imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall
             do imask =1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = -V(I+imask,j,k)
               W(I-imask+1,j,k) = -W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = -DV(I+imask,j,k)
               DW(I-imask+1,j,k) = -DW(I+imask,j,k)
             enddo
           endif
         ! east 
         ELSEIF(Mask(i-1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall  
             do imask = 1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = V(i-imask,j,k)
               W(i+imask-1,j,k) = W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = DV(i-imask,j,k)
               DW(i+imask-1,j,k) = DW(i-imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = -V(i-imask,j,k)
               W(i+imask-1,j,k) = -W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = -DV(i-imask,j,k)
               DW(i+imask-1,j,k) = -DW(i-imask,j,k)
             enddo
           endif
         ENDIF ! end mask+1=1 
       ENDIF ! end mask=0 
     ENDDO
     ENDDO
     ENDDO

     Deallocate(DelxUzL)
     Deallocate(DelxVzL)
     Deallocate(DelxWzL)
     Deallocate(DelxUzR)
     Deallocate(DelxVzR)
     Deallocate(DelxWzR)
     Deallocate(DelyUzL)
     Deallocate(DelyVzL)
     Deallocate(DelyWzL)
     Deallocate(DelyUzR)
     Deallocate(DelyVzR)
     Deallocate(DelyWzR)

     end subroutine vel_bc_old


     subroutine vel_bc
!----------------------------------------------------
!    Boundary conditions for velocity
!    Called by 
!       main and get_UVW
!    Last update: 01/02/2011, Gangfeng Ma
!    
!    The top and bottom boundary conditions are 
!    replaced by M.Derakhti as in Derakhti etal 2015a
!    Last update: 24/09/2015
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,imask,ii
     real(SP) :: Wtop,Wbot,Cdrag,Phi,Dz1,Cg
     !added by M.Derakhti
     real(SP) :: DUmag,Ustar,Z0,TauBot1,TauBot2,Fext_Bot1,Fext_Bot2,Fext_turb1,Fext_turb2, &
                 Fext_Wind1,Fext_Wind2,Atop,Btop,Abot,Bbot,DUprime,DVprime,dUdX,dUdy,dVdx,dVdy,dWdx,dWdy, &
                 DUtop_old,DUbot_old,DVtop_old,DVbot_old,MinDepBC
     real(SP), dimension(:,:,:), allocatable :: DelxUzL,DelxVzL,DelxWzL,&
          DelxUzR,DelxVzR,DelxWzR,DelyUzL,DelyVzL,DelyWzL,DelyUzR,DelyVzR,DelyWzR
                            
     MinDepBC = 0.01

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by Cheng for nesting. Please search for others with (COUPLING) in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_X0==1) then  ! free-slip wall
         do i = 1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==2) then ! no-slip wall
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = -V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = -W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = -DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = -DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==3) then ! inflow and outflow
         do i = 1,Nghost
           U(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)-U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)-V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = 2.0*Win_X0(j,k)-W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = 2.0*Din_X0(j)*Uin_X0(j,k)-DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = 2.0*Din_X0(j)*Vin_X0(j,k)-DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = 2.0*Din_X0(j)*Win_X0(j,k)-DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==4) then
         do i =1,Nghost
           U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dx/2.0/Zob))**2
         endif
         Phi = dx*Cdrag*sqrt(V(Ibeg,j,k)**2+W(Ibeg,j,k)**2)/(Cmu(Ibeg,j,k)+CmuHt(Ibeg,j,k))
         Phi = dmin1(Phi,2.0)
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = (1.0-Phi)*V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = (1.0-Phi)*W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = D(i,j)*U(Ibeg-i,j,k)
           DV(Ibeg-i,j,k) = D(i,j)*V(Ibeg-i,j,k)
           DW(Ibeg-i,j,k) = D(i,j)*W(Ibeg-i,j,k)
         enddo
       elseif(Bc_X0==8) then  ! specify u,v,w at ghost cells
         do i = 1,Nghost
           U(Ibeg-i,j,k) = Uin_X0(j,k)
           V(Ibeg-i,j,k) = Vin_X0(j,k)
           W(Ibeg-i,j,k) = Win_X0(j,k)
           DU(Ibeg-i,j,k) = Din_X0(j)*Uin_X0(j,k)
           DV(Ibeg-i,j,k) = Din_X0(j)*Vin_X0(j,k)
           DW(Ibeg-i,j,k) = Din_X0(j)*Win_X0(j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_Xn==1) then  ! free-slip wall 
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==2) then ! no-slip wall
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = -V(Iend-i+1,j,k)
           W(Iend+i,j,k) = -W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = -DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = -DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==3) then
         do i = 1,Nghost
           U(Iend+i,j,k) = 2.0*Uin_Xn(j,k)-U(Iend-i+1,j,k)
           V(Iend+i,j,k) = 2.0*Vin_Xn(j,k)-V(Iend-i+1,j,k)
           W(Iend+i,j,k) = 2.0*Win_Xn(j,k)-W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = 2.0*Din_Xn(j)*Uin_Xn(j,k)-DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = 2.0*Din_Xn(j)*Vin_Xn(j,k)-DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = 2.0*Din_Xn(j)*Win_Xn(j,k)-DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==4) then ! inflow and outflow
         do i = 1,Nghost
           U(Iend+i,j,k) = U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dx/2.0/Zob))**2
         endif
         Phi = dx*Cdrag*sqrt(V(Iend,j,k)**2+W(Iend,j,k)**2)/(Cmu(Iend,j,k)+CmuHt(Iend,j,k))
         Phi = dmin1(Phi,2.0)
         do i =1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = (1.0-Phi)*V(Iend-i+1,j,k)
           W(Iend+i,j,k) = (1.0-Phi)*W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = D(i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(i,j)*W(Iend+i,j,k)
         enddo
       elseif(Bc_Xn==6) then
         do i = 1,Nghost
		 ! Changed by Cheng from Sommerfeld to Orlanski open boundary condition
           Cg = -(U0(Iend+i-1,j,k)-U00(Iend+i-1,j,k)+1.e-16)/  &
                 (U00(Iend+i-1,j,k)-U00(Iend+i-2,j,k)+1.e-16)
!           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Cg = max(min(Cg,1.0),0.0)
           U(Iend+i,j,k) = Cg*U0(Iend+i-1,j,k)+(1.0-Cg)*U0(Iend+i,j,k)

           Cg =-(V0(Iend+i-1,j,k)-V00(Iend+i-1,j,k)+1.e-16)/  &
                 (V00(Iend+i-1,j,k)-V00(Iend+i-2,j,k)+1.e-16)
           Cg= max(min(Cg,1.0),0.0)
           V(Iend+i,j,k) = Cg*V0(Iend+i-1,j,k)+(1.0-Cg)*V0(Iend+i,j,k)
           
           Cg =-(W0(Iend+i-1,j,k)-W00(Iend+i-1,j,k)+1.e-16)/  &
                 (W00(Iend+i-1,j,k)-W00(Iend+i-2,j,k)+1.e-16)
           Cg= max(min(Cg,1.0),0.0)
           W(Iend+i,j,k) = Cg*W0(Iend+i-1,j,k)+(1.0-Cg)*W0(Iend+i,j,k)

           DU(Iend+i,j,k) = D(Iend+i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(Iend+i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(Iend+i,j)*W(Iend+i,j,k)
         enddo 
       elseif(Bc_Xn==8) then
         do i = 1,Nghost
           U(Iend+i,j,k) = Uin_Xn(j,k)
           V(Iend+i,j,k) = Vin_Xn(j,k)
           W(Iend+i,j,k) = Win_Xn(j,k)
           DU(Iend+i,j,k) = Din_Xn(j)*Uin_Xn(j,k)
           DV(Iend+i,j,k) = Din_Xn(j)*Vin_Xn(j,k)
           DW(Iend+i,j,k) = Din_Xn(j)*Win_Xn(j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Y0==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = -U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = -W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = -DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = -DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==4) then
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dy/2.0/Zob))**2
         endif
         Phi = dy*Cdrag*sqrt(U(i,Jbeg,k)**2+W(i,Jbeg,k)**2)/(Cmu(i,Jbeg,k)+CmuHt(i,Jbeg,k))
         Phi = dmin1(Phi,2.0)
         do j =1,Nghost
           U(i,Jbeg-j,k) = (1.0-Phi)*U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = (1.0-Phi)*W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = D(i,j)*U(i,Jbeg-j,k)
           DV(i,Jbeg-j,k) = D(i,j)*V(i,Jbeg-j,k)
           DW(i,Jbeg-j,k) = D(i,j)*W(i,Jbeg-j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Yn==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = -U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = -W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = -DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = -DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==4) then
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==5) then ! added by Cheng for wall friction
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*dy/2.0/Zob))**2
         endif
         Phi = dy*Cdrag*sqrt(U(i,Jend,k)**2+W(i,Jend,k)**2)/(Cmu(i,Jend,k)+CmuHt(i,Jend,k))
         Phi = dmin1(Phi,2.0)
         do j =1,Nghost
           U(i,Jend+j,k) = (1.0-Phi)*U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jend+j,k) = (1.0-Phi)*W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = D(i,j)*U(i,Jend+j,k)
           DV(i,Jend+j,k) = D(i,j)*V(i,Jend+j,k)
           DW(i,Jend+j,k) = D(i,j)*W(i,Jend+j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif
     
     allocate(DelxUzL(Mloc,Nloc,Kloc1))
     allocate(DelxVzL(Mloc,Nloc,Kloc1))
     allocate(DelxWzL(Mloc,Nloc,Kloc1))
     allocate(DelxUzR(Mloc,Nloc,Kloc1))
     allocate(DelxVzR(Mloc,Nloc,Kloc1))
     allocate(DelxWzR(Mloc,Nloc,Kloc1))
     allocate(DelyUzL(Mloc,Nloc,Kloc1))
     allocate(DelyVzL(Mloc,Nloc,Kloc1))
     allocate(DelyWzL(Mloc,Nloc,Kloc1))
     allocate(DelyUzR(Mloc,Nloc,Kloc1))
     allocate(DelyVzR(Mloc,Nloc,Kloc1))
     allocate(DelyWzR(Mloc,Nloc,Kloc1))

     !added by M.Derakhti
     !here we only need zL(kend+1) and zR(Kbeg)
     !we can later define new subroutines to remove unnecessary calculations
     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)
     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)
     !kinematic boundary conditions
     call KBC_surface
     call KBC_bottom
     call delxFun1_3d(UzL,DelxUzL)
     call delxFun1_3d(VzL,DelxVzL)
     call delxFun1_3d(WzL,DelxWzL)
     call delxFun1_3d(UzR,DelxUzR)
     call delxFun1_3d(VzR,DelxVzR)
     call delxFun1_3d(WzR,DelxWzR)
     call delyFun1_3d(UzL,DelyUzL)
     call delyFun1_3d(VzL,DelyVzL)
     call delyFun1_3d(WzL,DelyWzL)
     call delyFun1_3d(UzR,DelyUzR)
     call delyFun1_3d(VzR,DelyVzR)
     call delyFun1_3d(WzR,DelyWzR)

     ! top and bottom
     ! rewritten by M.Derakhti as in derakhti etal 2015a
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       !these are added to calculate the tangential external stress parallel to the surface and bed
       Atop = dsqrt(1.0+DelxEta(i,j)**2+DelyEta(i,j)**2)
       Btop = dsqrt(1.0+DelxEta(i,j)**2)
       Abot = dsqrt(1.0+DelxH(i,j)**2+DelyH(i,j)**2)
       Bbot = dsqrt(1.0+DelxH(i,j)**2)
       !velocities parallel to the bed
       DUprime = (DU(i,j,Kbeg)-DelxH(i,j)*DW(i,j,Kbeg))/Bbot 
       DVprime = (-DelxH(i,j)*DelyH(i,j)*DU(i,j,Kbeg)+(1.0+DelxH(i,j)**2)*DV(i,j,Kbeg)-DelyH(i,j)*DW(i,j,Kbeg))/Abot/Bbot
       DUmag = dsqrt(DUprime**2+DVprime**2)
       !the distance of the first grid point from the bed
       Dz1 = 0.5*D(i,j)*dsig(Kbeg)/Abot
       !estimate friction velocity
       Z0 = Zob / 30.0
       if(ibot==1) then
         Cdrag = Cd0
       else
# if defined (SEDIMENT)
         Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(Dz1/Z0))**2                                                          
# else
         Cdrag = (Kappa/log(Dz1/Z0))**2
# endif
       endif
       Cdrag = dmin1(Cdrag,(Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))/Dz1/(DUmag+Small)*D(i,j))

       !BSS/rho parallel to the bottom
       TauBot1 = zero
       TauBot2 = zero
       if(Bc_Z0==1) then  ! free-slip
          TauBot1 = zero   
          TauBot2 = zero
       elseif(Bc_Z0==2) then !no slip
          TauBot1 = DUprime / D(i,j) / Dz1 * (Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))   
          TauBot2 = DVprime / D(i,j) / Dz1 * (Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))
       elseif(Bc_Z0==5) then !friction law
          TauBot1 = Cdrag * DUprime * DUmag / D(i,j)**2   
          TauBot2 = Cdrag * DVprime * DUmag / D(i,j)**2
       endif

       !equation (46)
       Fext_Bot1 = Bbot * TauBot1
       Fext_Bot2 = (DelxH(i,j)*DelyH(i,j)*TauBot1+Abot*TauBot2)/Bbot
       !Wind stress paralle to the free surface WS/rho (equation 49)
       Fext_Wind1 = 1.0/Atop*((1-DelxEta(i,j)**2)*Wsx(i,j)-DelxEta(i,j)*DelyEta(i,j)*Wsy(i,j))
       Fext_Wind2 = 1.0/Atop*((1-DelyEta(i,j)**2)*Wsy(i,j)-DelxEta(i,j)*DelyEta(i,j)*Wsx(i,j))
       !Using the exact stress BC as given in Derakhti etal 2015a to get DU and DV at the ghost cells
       !No need to calculate Wbot here, because Kinematic B.C.s are imposed during the construction 
       !step in the z direction, such that Wtop = WzL(kend+1) and Wbot = WzR(Kbeg).
       !DW at ghost cells are obtained using K.B.C such that 
       !dW/DSigma = 2*(Wtop-W(kend-k+1))/(sigc(Kend+k)-sigc(Kbeg-k+1)) or = 2*(W(kbeg+k-1)-Wbot)/(sigc(Kbeg+k-1)-sigc(Kbeg-k))
       !during a strong opossing flow near shorelines the below formulation needs some
       !numerical caps. I have tried several numerical caps to recognize and prevent this problem.
       !The following two criteria gave the best results as in Derakhti etal 2015a and NHWAVE_benchmark cases
       !However, we may find a better solution for it.
       !In all of my cases MinDepBC = MinDep, however if very fine grid used we may need to use a larger value for MinDepBC. 
       do k = 1,Nghost
          !bottom BCs
          DW(i,j,Kbeg-k) = 2.0*D(i,j)*WzR(i,j,Kbeg)-DW(i,j,Kbeg+k-1)
          !equation (47)
          if(Bc_Z0==1) then
            DU(i,j,Kbeg-k)=DU(i,j,Kbeg+k-1)
            DV(i,j,Kbeg-k)=DV(i,j,Kbeg+k-1)
          else
          DU(i,j,Kbeg-k) = DU(i,j,Kbeg+k-1) - (sigc(Kbeg+k-1)-sigc(Kbeg-k)) &
              *(Fext_Bot1*(D(i,j)**2/Abot/(Small+Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))) &
              + DelxH(i,j)*(DW(i,j,Kbeg+k-1)-D(i,j)*WzR(i,j,Kbeg))/((sigc(Kbeg+k-1)-sigc(Kbeg-k))/2.0) &
              + D(i,j)**2/Abot**2 &
              *(-2.0*DelxH(i,j)*DelxUzR(i,j,Kbeg) &
              -(1.0-DelxH(i,j)**2)* DelxWzR(i,j,Kbeg) &
              -DelyH(i,j)*( DelyUzR(i,j,Kbeg)+DelxVzR(i,j,Kbeg)-DelxH(i,j)* DelyWzR(i,j,Kbeg))))
          DUbot_old = 2.0*D(i,j)*UzR(i,j,Kbeg)-DU(i,j,Kbeg+k-1)
          if (DUbot_old>0.)then
             DU(i,j,Kbeg-k) = dmin1(3.0*DUbot_old,dmax1(-1.0*DUbot_old,DU(i,j,Kbeg-k)))
          elseif (DUbot_old<0.)then
             DU(i,j,Kbeg-k) = dmax1(3.0*DUbot_old,dmin1(-1.0*DUbot_old,DU(i,j,Kbeg-k))) 
          endif
          if (D(i,j)<MinDepBC) then
             DU(i,j,Kbeg-k) = DU(i,j,Kbeg+k-1)
          endif 

          ! y-dir
          DV(i,j,Kbeg-k) = DV(i,j,Kbeg+k-1) - (sigc(Kbeg+k-1)-sigc(Kbeg-k)) &
              *(Fext_Bot2*(D(i,j)**2/Abot/(Small+Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))) & 
              +DelyH(i,j)*(DW(i,j,Kbeg+k-1)-D(i,j)*WzR(i,j,Kbeg))/((sigc(Kbeg+k-1)-sigc(Kbeg-k))/2.0) &
              +D(i,j)**2/Abot**2 &
              *(-2.0*DelyH(i,j)*DelyVzR(i,j,Kbeg) &
              -(1.0-DelyH(i,j)**2)*DelyWzR(i,j,Kbeg) &
              -DelxH(i,j)*( DelxVzR(i,j,Kbeg)+DelyUzR(i,j,Kbeg)-DelyH(i,j)*DelxWzR(i,j,Kbeg))))
          DVbot_old = 2.0*D(i,j)*VzR(i,j,Kbeg)-DV(i,j,Kbeg+k-1)
          if (DVbot_old>0.)then
             DV(i,j,Kbeg-k) = dmin1(3.0*DVbot_old,dmax1(-1.0*DVbot_old,DV(i,j,Kbeg-k)))
          elseif (DVbot_old<0.)then
             DV(i,j,Kbeg-k) = dmax1(3.0*DVbot_old,dmin1(-1.0*DVbot_old,DV(i,j,Kbeg-k))) 
          endif
          if (D(i,j)<MinDepBC) then
             DV(i,j,Kbeg-k) = DV(i,j,Kbeg+k-1) 
          endif
          endif 

          ! surface BCs
          DW(i,j,Kend+k) = 2.0*D(i,j)*WzL(i,j,Kend+1)-DW(i,j,Kend-k+1)
          DU(i,j,Kend+k) = DU(i,j,Kend-k+1) + (sigc(Kend+k)-sigc(Kend-k+1)) &
            *(Fext_Wind1*(D(i,j)**2/Atop/(Small+0.5*(Cmu(i,j,Kend)+CmuR(i,j,Kend)+CmuVt(i,j,Kend)+ & 
            Cmu(i,j,Kend+1)+CmuR(i,j,Kend+1)+CmuVt(i,j,Kend+1)))) & 
            - DelxEta(i,j)*(D(i,j)*WzL(i,j,Kend+1)-DW(i,j,Kend-k+1))/((sigc(Kend+k)-sigc(Kend-k+1))/2.0) &
            + D(i,j)**2/Atop**2 &
            *( + 2.0*DelxEta(i,j)*DelxUzL(i,j,Kend+1) &
            -(1.0-DelxEta(i,j)**2)* DelxWzL(i,j,Kend+1) &
            +DelyEta(i,j)*(DelyUzL(i,j,Kend+1)+DelxVzL(i,j,Kend+1)+DelxEta(i,j)*DelyWzL(i,j,Kend+1))))

          DUtop_old = 2.0*D(i,j)*UzL(i,j,Kend+1)-DU(i,j,Kend-k+1)
          if(DUtop_old>0.)then
             DU(i,j,Kend+k) = dmin1(3.0*DUtop_old,dmax1(-1.0*DUtop_old,DU(i,j,Kend+k)))
          elseif(DUtop_old<0.)then
             DU(i,j,Kend+k) = dmax1(3.0*DUtop_old,dmin1(-1.0*DUtop_old,DU(i,j,Kend+k))) 
          endif

          if(D(i,j)<MinDepBC) then
             DU(i,j,Kend+k) = DU(i,j,Kend-k+1)
          endif

          ! y-dir
          DV(i,j,Kend+k) = DV(i,j,Kend-k+1) + (sigc(Kend+k)-sigc(Kend-k+1)) &
              *(Fext_Wind2*(D(i,j)**2/Atop/(Small+0.5*(Cmu(i,j,Kend  )+CmuR(i,j,Kend  )+CmuVt(i,j,Kend  )+ & 
                                                       Cmu(i,j,Kend+1)+CmuR(i,j,Kend+1)+CmuVt(i,j,Kend+1)))) & 
              - DelyEta(i,j)*(D(i,j)*WzL(i,j,Kend+1)-DW(i,j,Kend-k+1))/((sigc(Kend+k)-sigc(Kend-k+1))/2.0) &
              + D(i,j)**2/Atop**2 &
              *( + 2.0*DelyEta(i,j)*DelyVzL(i,j,Kend+1) &
              -(1.0-DelyEta(i,j)**2)*DelyWzL(i,j,Kend+1) &
              +DelxEta(i,j)*(DelxVzL(i,j,Kend+1)+DelyUzL(i,j,Kend+1)+DelyEta(i,j)*DelxWzL(i,j,Kend+1))))

          DVtop_old = 2.0*D(i,j)*VzL(i,j,Kend+1)-DV(i,j,Kend-k+1)
          if (DVtop_old>0.)then
             DV(i,j,Kend+k) = dmin1(3.0*DVtop_old,dmax1(-1.0*DVtop_old,DV(i,j,Kend+k)))
          elseif (DVtop_old<0.)then
             DV(i,j,Kend+k) = dmax1(3.0*DVtop_old,dmin1(-1.0*DVtop_old,DV(i,j,Kend+k))) 
          endif

          if (D(i,j)<MinDepBC) then
             DV(i,j,Kend+k) = DV(i,j,Kend-k+1)
          endif

         U(i,j,Kbeg-k) = DU(i,j,Kbeg-k)/D(i,j)
         V(i,j,Kbeg-k) = DV(i,j,Kbeg-k)/D(i,j)
         W(i,j,Kbeg-k) = DW(i,j,Kbeg-k)/D(i,j)      
         U(i,j,Kend+k) = DU(i,j,Kend+k)/D(i,j)
         V(i,j,Kend+k) = DV(i,j,Kend+k)/D(i,j)
         W(i,j,Kend+k) = DW(i,j,Kend+k)/D(i,j)
       enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks 02/15/2013
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(i,j)==0) THEN
         ! south boundary 
         IF(Mask(i,j+1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j-imask+1,k) = U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = W(i,j+imask,k)
               DU(i,j-imask+1,k) = DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = DW(i,j+imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j-imask+1,k) = -U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = -W(i,j+imask,k)
               DU(i,j-imask+1,k) = -DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = -DW(i,j+imask,k)
             enddo
           endif
         ! north  
         ELSEIF(Mask(i,j-1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j+imask-1,k) = U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = W(i,j-imask,k)
               DU(i,j+imask-1,k) = DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = DW(i,j-imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j+imask-1,k) = -U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = -W(i,j-imask,k)
               DU(i,j+imask-1,k) = -DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = -DW(i,j-imask,k)
             enddo
           endif
         ! west
         ELSEIF(Mask(i+1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = V(I+imask,j,k)
               W(I-imask+1,j,k) = W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = DV(I+imask,j,k)
               DW(I-imask+1,j,k) = DW(I+imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall
             do imask =1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = -V(I+imask,j,k)
               W(I-imask+1,j,k) = -W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = -DV(I+imask,j,k)
               DW(I-imask+1,j,k) = -DW(I+imask,j,k)
             enddo
           endif
         ! east 
         ELSEIF(Mask(i-1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall  
             do imask = 1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = V(i-imask,j,k)
               W(i+imask-1,j,k) = W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = DV(i-imask,j,k)
               DW(i+imask-1,j,k) = DW(i-imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = -V(i-imask,j,k)
               W(i+imask-1,j,k) = -W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = -DV(i-imask,j,k)
               DW(i+imask-1,j,k) = -DW(i-imask,j,k)
             enddo
           endif
         ENDIF ! end mask+1=1 
       ENDIF ! end mask=0 
     ENDDO
     ENDDO
     ENDDO

     Deallocate(DelxUzL)
     Deallocate(DelxVzL)
     Deallocate(DelxWzL)
     Deallocate(DelxUzR)
     Deallocate(DelxVzR)
     Deallocate(DelxWzR)
     Deallocate(DelyUzL)
     Deallocate(DelyVzL)
     Deallocate(DelyWzL)
     Deallocate(DelyUzR)
     Deallocate(DelyVzR)
     Deallocate(DelyWzR)

     end subroutine vel_bc


     subroutine KBC_surface
!-------------------------------------------------
!    Applied kinematic boundary conditions at free surface                                                                        
!    should be Called after any z-construction                                                                                   
!                                                                                                                               
!    Last update: 01/20/2015, Morteza Derakhti                                                                                    
!------------------------------------------------                                                                                
     use global, only: WzL,UzL,VzL, &
                       Ibeg,Iend,Jbeg,Jend,Kend1,&
                       DelxEta,DelyEta,Eta,Eta0,dt
     implicit none
     integer :: i,j,k

     do i = Ibeg,Iend
     do j = Jbeg,Jend
       WzL(i,j,Kend1) = (Eta(i,j)-Eta0(i,j))/dt+DelxEta(i,j)*UzL(i,j,Kend1)  & 
                        + DelyEta(i,j)*VzL(i,j,Kend1)                                                     
     enddo
     enddo

     end subroutine KBC_surface

     subroutine KBC_bottom
!-------------------------------------------------  
!    Applied kinematic boundary conditions at bottom                                                                              
!    should be Called after any z-construction                                                                                   
!       construction                                                                                                              
!    Last update: 01/20/2015, Morteza Derakhti                                                                                   
!------------------------------------------------                                                                                 
     use global, only: WzR,UzR,VzR,UzL,VzL, &
                       Ibeg,Iend,Jbeg,Jend,Kbeg,&
                       DelxH,DelyH,DeltH,Zero,Bc_Z0
     implicit none
     integer :: i,j,k

     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Bc_Z0==2) then !no-slip                                                                                              
         UzR(i,j,Kbeg) = Zero
         VzR(i,j,Kbeg) = Zero
         UzL(i,j,Kbeg) = Zero
         VzL(i,j,Kbeg) = Zero
       endif
       WzR(i,j,Kbeg) = -DeltH(i,j)-DelxH(i,j)*UzR(i,j,Kbeg)  &
                       -DelyH(i,j)*VzR(i,j,Kbeg)
     enddo
     enddo

     end subroutine KBC_bottom

     subroutine wl_bc
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duvw
!    Last update: 14/06/2012, Gangfeng Ma
!-----------------------------------------------------------
     use global
     implicit none
     real(SP) :: Cg
     integer :: i,j

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by Cheng for nesting. Please search for others with (COUPLING) in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
       if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! free/no-slip wall ! added by Cheng for wall friction
         do i = 1,Nghost
           D(Ibeg-i,j) = D(Ibeg+i-1,j)
         enddo
       elseif(Bc_X0==3) then ! inflow
         if((WaveMaker(1:7)=='LEF_TID').or.(WaveMaker(1:7)=='FLUX_LR')) then ! FLUX_LR is added by Cheng
           do i = 0,Nghost
             D(Ibeg-i,j) = Din_X0(j)
           enddo
         else
           do i = 1,Nghost
             D(Ibeg-i,j) = 2.0*Din_X0(j)-D(Ibeg+i-1,j)
           enddo
         endif
       elseif(Bc_X0==4) then ! outflow
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       elseif(Bc_X0==8) then
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       endif
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
       if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
         do i = 1,Nghost
           D(Iend+i,j) = D(Iend-i+1,j)
         enddo
       elseif(Bc_Xn==3) then 
         do i = 1,Nghost
           D(Iend+i,j) = 2.0*Din_Xn(j)-D(Iend-i+1,j)
         enddo
       elseif(Bc_Xn==4) then
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       elseif(Bc_Xn==6) then 
         do i = 1,Nghost
		 ! Changed by Cheng from Sommerfeld to Orlanski open boundary condition 
           Cg =-(Eta0(Iend+i-1,j)-Eta00(Iend+i-1,j)+1.e-16)/  &
                 (Eta00(Iend+i-1,j)-Eta00(Iend+i-2,j)+1.e-16)
!           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Cg= max(min(Cg,1.0),0.0)  
           Eta(Iend+i,j) = Cg*Eta0(Iend+i-1,j)+(1.0-Cg)*Eta0(Iend+i,j)
           D(Iend+i,j) = Hc(Iend+i,j)+Eta(Iend+i,j)
         enddo
       elseif(Bc_Xn==8) then
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       endif
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

! y-direction and corners                                                                                                     
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif 
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif     
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,j) = D(i,Jbeg+Nghost-j)
       enddo
       enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,Jend+j) = D(i,Jend-j+1)
       enddo
       enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(D)
# endif
     
     return
     end subroutine wl_bc

     subroutine phi_2D_coll(phi)
!-----------------------------------------------------
!    This subroutine is used to collect data into ghost cells
!    Called by
!       eval_duvw
!    Last update: 22/12/2010, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     real(SP), intent(inout) :: phi(Mloc,Nloc)
     integer :: i,j

     ! x-direction
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(i,j) = phi(Ibeg+Nghost-i,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(Iend+i,j) = phi(Iend-i+1,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif
 
     ! y-direction and corners
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif     
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,j) = phi(i,Jbeg+Nghost-j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,Jend+j) = phi(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(phi)
# endif    

     end subroutine phi_2D_coll


     subroutine update_wind
!--------------------------------------------------------
!    Update wind speed at current time step
!    Called by
!       main
!    Last update: 09/07/2013, Gangfeng Ma
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,iglob,jglob
     real(SP),dimension(Mglob,Nglob) :: WUG,WVG
     real(SP) :: Wds,Cds,Cdsmin,Cdsmax

     ! update wind speed
     if(Iws==1) then
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         WdU(i,j) = WindU
         WdV(i,j) = WindV
       enddo
       enddo
     elseif(Iws==2) then
       open(7,file='wind.txt')

       do j = 1,Nglob
       do i = 1,Mglob
         read(7,*) WUG(i,j),WVG(i,j)
       enddo
       enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         WdU(i,j) = WUG(iglob,jglob)
         WdV(i,j) = WVG(iglob,jglob)
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         WdU(i,j) = WUG(iglob,jglob)
         WdV(i,j) = WVG(iglob,jglob)
       enddo
       enddo
# endif

     endif

     ! wind drag coefficient
     Cdsmin = 1.e-3*(0.61+0.063*6.0)
     Cdsmax = 1.e-3*(0.61+0.063*50.0)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Wds = sqrt(WdU(i,j)**2+WdV(i,j)**2)
       Cds = 1.e-3*(0.61+0.063*Wds)
       Cds = dmin1(dmax1(Cds,Cdsmin),Cdsmax)
       
       Wsx(i,j) = 0.001293*Cds*WdU(i,j)*Wds
       Wsy(i,j) = 0.001293*Cds*WdV(i,j)*Wds
     enddo
     enddo

     return
     end subroutine


     subroutine update_mask
!------------------------------------------------------  
!    This subroutine is used to update mask for wetting-drying
!    Called by                                                
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: Ibeg,Iend,Jbeg,Jend,Eta,Hc,D,MinDep,  &
                       Mask,Mask_Struct,Mask9
     implicit none
     integer :: i,j

     ! Mask at ghost cells keeps no change
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask_Struct(i,j)==0) cycle
       
       ! flooding (dry->wet)
       if(Mask(i,j)==0) then
         if(Mask(i-1,j)==1.and.Eta(i-1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i+1,j)==1.and.Eta(i+1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j-1)==1.and.Eta(i,j-1)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j+1)==1.and.Eta(i,j+1)>Eta(i,j)) Mask(i,j)=1
       else
         ! drying (wet->dry)
         if(abs(D(i,j)-MinDep)<=1.e-6) then
           Mask(i,j) = 0
           Eta(i,j) = MinDep-Hc(i,j)
           D(i,j) = Eta(i,j)+Hc(i,j)           
         endif
       endif
     enddo
     enddo

     Mask = Mask*Mask_Struct

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)    
# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     end subroutine update_mask


     subroutine update_vars
!------------------------------------------------------ 
!    This subroutine is used to save variables at 
!    last time step
!    Called by   
!       main 
!    Last update: 22/12/2010, Gangfeng Ma 
!----------------------------------------------------- 
     use global        
     implicit none

     Eta00 = Eta0
     U00 = U0
     V00 = V0
     W00 = W0

     D0 = D
     Eta0 = Eta
     U0 = U
     V0 = V
     W0 = W
     DU0 = DU
     DV0 = DV
     DW0 = DW
     DTke0 = DTke
     DEps0 = DEps

# if defined (BUBBLE)
     DNbg0 = DNbg
# endif

# if defined (SEDIMENT)
     DConc0 = DConc
     Bed0 = Bed
# endif

# if defined (SALINITY)
     DSali0 = Dsali
# endif

# if defined (TWOLAYERSLIDE)
     Ha0 = Ha
     HUa0 = HUa
     HVa0 = HVa
     Ua0 = Ua
     Va0 = Va
# endif

# if defined (VEGETATION)
     DWke0 = DWke
# endif

     end subroutine update_vars

  
     subroutine update_wave_bc
!------------------------------------------------------
!    This subroutine is used to update boundary conditions
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='LEF_SOL') then
       call solitary_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_LIN') then
       call linear_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_CON') then
       call cnoidal_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_STK') then
       call stokes_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_SPC') then
       call random_wave_left_boundary
     elseif((WaveMaker(1:7)=='LEF_JON').or.(WaveMaker(1:7)=='LEF_TMA')) then
       call jonswap_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_TID') then
       call tidal_wave_left_boundary
     elseif(WaveMaker(1:7)=='FOCUSED') then
       call focused_wave_left_boundary
     elseif(WaveMaker(1:7)=='WAV_CUR') then
       call wave_current_left_boundary
     endif
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='RIG_LIN') then
       call linear_wave_right_boundary
     endif
# if defined (PARALLEL)
     endif
# endif

     if(WaveMaker(1:7)=='FLUX_LR') then
       call flux_left_right_boundary
     endif

     end subroutine update_wave_bc 


     subroutine flux_left_right_boundary
!-----------------------------------------------------------   
!    This subroutine is used to specify left/right boundary                                         
!    Called by 
!       update_wave_bc 
!    Last update: 14/06/2012, Gangfeng Ma
!-----------------------------------------------------------
     use global
     implicit none
     integer :: j,k,n
     real(SP) :: Zlev1,Zlev2,Uavg_Left,Uavg_Right,Ufric,sintep,UU,  &
                 FluxL,FluxR,myvar,Ramp
# if defined (SEDIMENT)
     real(SP) :: Sinterp
     real(SP) :: Tim0(8),Cin0(8)
# endif

     Uavg_Left = 0.39
     Uavg_Right = 0.39

     FluxL = 0.0
     FluxR = 0.0

     if(TRamp>0.0) then
       Ramp = tanh(TIME/TRamp)
     else
       Ramp = 1.0
     endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       ! Ein_X0(j) = 1.5*Eta(Ibeg,j)-0.5*Eta(Ibeg+1,j)
	   Ein_X0(j) = 0.0 ! modified by Cheng to lock the influx boundary condition
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! log-profile of velocity
     Ufric = Uavg_Left*Kappa/(log(30.*Din_X0(Jbeg)/Zob)-1) 

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev1 = sigc(k)*Din_X0(j)
	   ! changed by Cheng to impose velocity rather than flux
       Uin_X0(j,k) = Ufric/Kappa*log(30.*Zlev1/Zob)*Ramp
       FluxL = FluxL+dsig(k)*Din_X0(j)*dy*Uin_X0(j,k)
!       Uin_X0(j,k) = FluxL/Din_X0(j)
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0

# if defined (SALINITY)
       Sin_X0(j,k) = 0.0
# endif

# if defined (SEDIMENT)
       open(41,file='cin0.txt')
       do n = 1,8
         read(41,*) Tim0(n),Cin0(n)
         Tim0(n) = Tim0(n)*3600.
         Cin0(n) = Cin0(n)/Srho
       enddo
       close(41)       

       if(TIME<=Tim0(1)) then
         Sed_X0(j,k) = Cin0(1)
       elseif(TIME>=Tim0(8)) then
         Sed_X0(j,k) = Cin0(8)
       else
         do n = 2,8
           if(TIME>Tim0(n-1).and.TIME<=Tim0(n)) then
             Sinterp = (TIME-Tim0(n-1))/(Tim0(n)-Tim0(n-1))
             Sed_X0(j,k) = Cin0(n-1)*(1.0-Sinterp)+Cin0(n)*Sinterp
           endif
         enddo
       endif
# endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)     
     call MPI_ALLREDUCE(FluxL,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     FluxL = myvar
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif

!     if(TIME>TRamp) then
!       call linear_wave_right_boundary
!     endif

     do j = 1,Nloc
!       Ein_Xn(j) = 1.5*Eta(Iend,j)-0.5*Eta(Iend-1,j)
       Ein_Xn(j) = 0.0
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend+1,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev2 = sigc(k)*Din_Xn(j)
	   ! changed by Cheng to impose velocity rather than flux
       Uin_Xn(j,k) = FluxL/(Din_Xn(Jbeg)*dy*float(Nglob))
!       Uin_Xn(j,k) = FluxR/Din_Xn(j)
       Win_Xn(j,k) = 0.0
       Vin_Xn(j,k) = 0.0

# if defined (SALINITY)
       Sin_Xn(j,k) = 30.0
# endif

# if defined (SEDIMENT)
       Sed_Xn(j,k) = Conc(Iend,j,k)
# endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
     
     return
     end subroutine flux_left_right_boundary


     subroutine focused_wave_left_boundary
     !-----------------------------------------------------
     !    This subroutine is used to specify left boundary 
     !    Called by
     !       update_wave_bc 
     !    Last update: 03/15/2013, Derakhti Morteza (derakhti@udel.edu)
     !-----------------------------------------------------
     use global, only: Amp_i,k_i,f_i,nwave,x_breaking,t_breaking,depth_comp,sd_return, &
      SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
      Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,myid
# if defined (INTEL)
     use IFPORT, only: rand
# endif
     implicit none
     integer  :: j,k,Iter,nw
     real(SP),dimension(12000) :: Segma_i,c_i,eta_i,Stokes_Drift_i
     real(SP) :: Zlev,tc,RandU,RandW

     tc = TIME-int(TIME/20.0)*20.0

     Segma_i = Zero
     c_i = Zero
     eta_i = Zero
     Stokes_Drift_i = Zero
     do nw = 1 , nwave
       Segma_i(nw) = 2. * pi * f_i(nw)
       c_i(nw) = segma_i(nw) / k_i(nw) 
       eta_i(nw) = Amp_i(nw) * &
         dcos( -k_i(nw) * x_breaking - segma_i(nw) * (tc - t_breaking) )
       Stokes_Drift_i(nw) = Grav * Amp_i(nw) ** 2 / c_i(nw) / depth_comp / 2.
     end do

     do j = 1,Nloc
       Ein_X0(j) = Zero
       do nw = 1 , nwave
         Ein_X0(j) = Ein_X0(j) + eta_i(nw)
       enddo
       Din_X0(j) = Ein_X0(j) + Hfx(Ibeg,j) 
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = Zero
       Win_X0(j,k) = Zero
       do nw = 1 , nwave
         Uin_X0(j,k) = Uin_X0(j,k) + Grav * Amp_i(nw) * k_i(nw) / Segma_i(nw) * cosh( k_i(nw)*Zlev ) /  &
           cosh( k_i(nw)*(depth_comp+Ein_X0(j)) ) * dcos(-k_i(nw)*x_breaking - Segma_i(nw)*(tc-t_breaking)) &
           - Stokes_Drift_i(nw) * (2*k_i(nw)*depth_comp) / sinh(2*k_i(nw)*(depth_comp+Ein_X0(j))) * &
           cosh(2*k_i(nw)*Zlev) * sd_return
         Win_X0(j,k) = Win_X0(j,k) + Grav * Amp_i(nw) * k_i(nw) / Segma_i(nw) * sinh( k_i(nw)*Zlev ) / &
           cosh( k_i(nw)*(depth_comp+Ein_X0(j)) ) * dsin(-k_i(nw)*x_breaking - Segma_i(nw)*(tc-t_breaking))
       enddo
       Vin_X0(j,k) = Zero
     enddo         
     enddo

     end subroutine focused_wave_left_boundary


     subroutine stokes_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 26/04/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     
     Wave_Length = 2.0*pi/Wave_Number
     Celerity = Wave_Length/Per_Wave
     
     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME)+  &
                   Amp_Wave**2*Wave_Number/16.*cosh(Wave_Number*Dep_Wave)/  &
                   sinh(Wave_Number*Dep_Wave)**3*(2.0+cosh(2.*Wave_Number*Dep_Wave))* &
                   cos(2.*(pi/2-Segma*TIME))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*cosh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*cos(2.*(pi/2-Segma*TIME))
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*sinh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*sin(2.*(pi/2-Segma*TIME))
       Vin_X0(j,k) = 0.0
     enddo
     enddo

     end subroutine stokes_wave_left_boundary


     subroutine tidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 06/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), parameter, dimension(8) :: &
                 !  s2      m2    n2     k2     k1     p1     o1     q1
        period = (/43200.,44712.,45570.,43082.,86164.,86637.,92950.,96726./)

     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(2.0*pi/period(1)*TIME-pi/2.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Uin_X0(j,k) = sqrt(Grav/Hfx(Ibeg,j))*Ein_X0(j)
       Vin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
     enddo
     enddo

# if defined (SEDIMENT)
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Sed_X0(j,k) = 0.5/2650.
     enddo
     enddo
# endif

     end subroutine tidal_wave_left_boundary


     subroutine cnoidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer :: j,k,nn
     real(SP) :: Wave_Length,Celerity,Ytrough,Mod1,Wave_Number,Zup,Zlow,Zmid,Xstart, &
                 Zero1,Zlev,Xcn,Xsn,Xdn,cnoidal_cn,cnoidal_ck,Stokes_Drift,Fact

     call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)
     
     ! wave number
     Wave_Number = 2.0*pi/Wave_Length

     ! find zero start
     Zup = 1.0
     Zlow = 0.0
     Zmid= (Zup+Zlow)/2.0
     nn = 0
200  nn = nn+1
     Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

     if(abs(Zero1)<=1.0e-6) goto 210
     if(nn>1000) then
       write(*,*)'too many iterations; stop'
       stop
     endif
     if(Zero1<0.0) then
       Zup = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     else
       Zlow = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     endif
 210 continue
     Xstart = Zmid
!     write(*,*) "Mod=",Mod1,"Ytrough=",Ytrough,"Wave_Number=",wave_number,  &    
!                   "Xstart=",Xstart

     do j = 1,Nloc
       Ein_X0(j) = Ytrough+Amp_Wave*cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+  &
             2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! mean mass transport
     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)

       Xcn = cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)
       Xsn = sqrt(1.0-Xcn**2)
       Xdn = sqrt(1.0-Mod1*(1.0-Xcn**2)) 

       Uin_X0(j,k) = sqrt(Grav*Dep_Wave)*(-5./4.+3.*(Dep_Wave+Ytrough)/2.0/Dep_Wave-  &
            (Dep_Wave+Ytrough)**2/4./Dep_Wave**2+(3.*Amp_Wave/2./Dep_Wave-  &                       
            (Dep_Wave+Ytrough)*Amp_Wave/2./Dep_Wave**2)*Xcn**2-Amp_Wave**2/4./  &                      
            Dep_Wave**2*Xcn**4-8.0*Amp_Wave*cnoidal_ck(Mod1)**2/Wave_Length**2*(Dep_Wave/3-  &            
            Zlev**2/2/Dep_Wave)*(-Mod1**2*Xsn**2*Xcn**2+Xcn**2*Xdn**2-Xsn**2*Xdn**2))-  &
            ! substract mean mass transport
            Stokes_Drift*Fact
       Win_X0(j,k) = sqrt(Grav*Dep_Wave)*Zlev*2*Amp_Wave*cnoidal_ck(Mod1)/  &
            Wave_Length/Dep_Wave*(1+(Dep_Wave+Ytrough)/Dep_Wave+Amp_Wave/  &                        
            Dep_Wave*Xcn**2+32*cnoidal_ck(Mod1)**2/3/Wave_Length**2*(Dep_Wave**2-  &
            Zlev**2/2)*(Mod1**2*Xsn**2-Mod1**2*Xcn**2-Xdn**2))*Xsn*Xcn*Xdn
       Vin_X0(j,k) =0.0
     enddo
     enddo

     end subroutine cnoidal_wave_left_boundary


     subroutine cnoidal(height,depth,period,l,c,x2,mod)
!----------------------------------------------------------------
!    Cnoidal Function 
!    Called by
!       cnoidal_wave_left_boundary 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------------------
     use global, only: SP,Grav
     implicit none

     real(SP), intent(in) :: height, depth, period
     real(SP), intent(out) :: l, c, x2, mod
     real(SP) :: xa,xb,xtemp,cnoidal_cn,cnoidal_ck,cnoidal_ce
     integer :: n

     mod=0.99999999d0
     n=0
 40  n=n+1

     xa=mod*depth+2.0d0*height-mod*height-3.0d0*height*  &
        cnoidal_ce(mod)/cnoidal_ck(mod)-16.0d0*depth**3*mod**2*cnoidal_ck(mod)**2/  &
        3.0d0/Grav/height/period**2
     if(abs(xa).le.1.0e-8.or.n.gt.1000) goto 50
  
     xb=depth-height+3.0d0*height/2.0d0/mod/(1.0d0-mod)/cnoidal_ck(mod)**2*  &
        ((1.0d0-mod)*cnoidal_ck(mod)**2+cnoidal_ce(mod)**2-2.0d0*(1.0d0-mod)* &
        cnoidal_ck(mod)*cnoidal_ce(mod))-16.0d0*depth**3*mod*cnoidal_ck(mod)/3.0d0/Grav/ &
        (1.0d0-mod)/height/period**2*((1.0d0-mod)*cnoidal_ck(mod)+cnoidal_ce(mod))                          

     mod=mod-xa/xb
     goto 40
 50  continue

     ! sobet el at (1987, J. Waterway)                        
     l=4.0*cnoidal_ck(mod)*depth*sqrt(mod*depth/height/3.0)
     ! mei (1983) or simply c=L/T 
     xtemp=-mod+2.0-3.0*cnoidal_ce(mod)/cnoidal_ck(mod)
     c=sqrt(Grav*depth*(1.0+height/depth/mod*xtemp))
     x2=height/mod*(1.0-mod-cnoidal_ce(mod)/cnoidal_ck(mod))
 
     return
     end subroutine cnoidal


     function cnoidal_cn(u,mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod,u
     real(SP) ::  mod1,a0,a1,b0,b1,c0,c1,c(1000),a(1000),y(1000),cnoidal_cn
     integer :: n,i

     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     a(n) = a0
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       a(n) = a0
       c(n) = c0

       goto 15
     endif

 30  y(n) = 2.0**(n-1)*a(n)*u

     do i = n-1,1,-1
       y(i) = (y(i+1)+asin(c(i+1)/a(i+1)*sin(y(i+1))))/2.0
     enddo

     cnoidal_cn = cos(y(1))

     return 
     end function cnoidal_cn

     function cnoidal_ck(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: mod1,a0,a1,b0,b1,c0,c1,cnoidal_ck
     integer :: n

     mod1 = 1.0-mod
     cnoidal_ck = 0.0
     a0 = 1.0
     b0 = dsqrt(mod1)
     c0 = dsqrt(mod)
     n = 1

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       goto 15
     endif

 30  cnoidal_ck = 3.1415926535897932384626/2.0/a0

     return
     end function cnoidal_ck

     function cnoidal_ce(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: sum,mod1,a0,a1,b0,b1,c0,c1,ck,c(1000),cnoidal_ce
     integer :: n,k1

     cnoidal_ce = 0.0
     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1
       c(n) = c0

       goto 15
     endif

 30  ck = 3.1415926535897932384626/2.0/a0

     sum = 0.0
     do k1 = 1,n
       sum = sum+2.0**(k1-2)*c(k1)**2
     enddo

     cnoidal_ce = ck*(1.0-sum)

     return
     end function cnoidal_ce


     subroutine random_wave_left_boundary
!-----------------------------------------------------    
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma  
!-----------------------------------------------------  
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
                       MaxNumFreq,MaxNumDir,Nglob,Dep_Wave,Random_Phs
# if defined (PARALLEL)
     use global, only: myid
# endif
     implicit none
     integer :: nfreq,ndir,Iter,nk,i,j,k
     real(SP) :: Per_Wave,Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Angle,tmp1,tmp2,tmp3,tmp4,Fact, &
                 Wtheta(MaxNumDir,MaxNumFreq),Wh(MaxNumDir,MaxNumFreq),Stokes_Drift(MaxNumDir,MaxNumFreq), &
                 Wnum(MaxNumDir,MaxNumFreq),Phs_Lag,dfreq,ddir,Zlev

     ! find the right angle for periodic bc
     do nfreq = 1,NumFreq
     do ndir = 1,NumDir
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 75    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85    continue
       Wave_Length = 2.0*pi/Wave_Number

       Angle = Dire(ndir)*pi/180.
       goto 100
       if(Angle>zero) then
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3<Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = 0.5*pi-small
           else
             tmp3 = asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer
         tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       else
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3>Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = -0.5*pi+small
           else
             tmp3 = -asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer 
         tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       endif
 100   continue
       Wtheta(ndir,nfreq) = Angle
 
       if(nfreq==1) then
         dfreq = Freq(2)-Freq(1)
       elseif(nfreq==NumFreq) then
         dfreq = Freq(NumFreq)-Freq(NumFreq-1)
       else
         dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
       endif
       dfreq = abs(dfreq)

       if(ndir==1) then
         ddir = Dire(2)-Dire(1)
       elseif(ndir==NumDir) then
         ddir = Dire(NumDir)-Dire(NumDir-1)
       else
         ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
       endif
       ddir = abs(ddir)
       
       ! save wave number and wave height
       Wnum(ndir,nfreq) = Wave_Number
       Wh(ndir,nfreq) = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*dfreq*ddir)

       ! Stokes Drift
       Stokes_Drift(ndir,nfreq) = Grav*Wh(ndir,nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
!       if(myid.eq.0)write(12,'(8E20.10)') Freq(nfreq),Dire(ndir),dfreq,ddir,Wtheta(ndir,nfreq)*180./pi,  &
!               Wnum(ndir,nfreq),Wh(ndir,nfreq),Stokes_Drift(ndir,nfreq)
     enddo
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(ndir,nfreq)*cos(pi/2-Segma*TIME+Phs_Lag)
       enddo
       enddo
       Ein_X0(j)=Ein_X0(j)*tanh(TIME/20.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
       !if(myid.eq.0) write(12,*) TIME,Ein_X0(Jbeg+Nglob/2)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  &
             cos(Wtheta(ndir,nfreq)) 
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(ndir,nfreq)*Segma*sinh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag)
         Vin_X0(j,k) = Vin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  & 
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  & 
             sin(Wtheta(ndir,nfreq))
       enddo
       enddo
     enddo
     enddo

     end subroutine random_wave_left_boundary


     subroutine jonswap_wave_left_boundary
!-----------------------------------------------------                  
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,Jon_Spc, &
                       Freq,Freq_Max,Freq_Min,Dep_Wave,Per_Wave,Hfx,RanPhs, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     real(SP), dimension(NumFreq) :: Wh,Wnum,Stokes_Drift
     real(SP) :: dfreq,Segma,Celerity,Wave_Length,Wave_NUmber,Fk,Fkdif,Fact,Zlev
     integer :: nfreq,Iter,i,j,k

     dfreq = (Freq_Max-Freq_Min)/float(NumFreq)
     do nfreq = 1,NumFreq
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 77    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 87
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                 
              Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 77
 87    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave

       ! root-mean-square wave height                                                                   
       Wh(nfreq) = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)
       Wnum(nfreq) = Wave_Number
       Stokes_Drift(nfreq) = Grav*Wh(nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(nfreq)*cos(pi/2-Segma*TIME+RanPhs(nfreq))                           
       enddo
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+0.5*Wh(nfreq)*Segma*cosh(Wnum(nfreq)*Zlev)/  &                         
             sinh(Wnum(nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+RanPhs(nfreq))-Stokes_Drift(nfreq)*Fact          
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(nfreq)*Segma*sinh(Wnum(nfreq)*Zlev)/  &                            
             sinh(Wnum(nfreq)*Dep_Wave)*sin(pi/2-Segma*TIME+RanPhs(nfreq))
         Vin_X0(j,k) = 0.0
       enddo
     enddo
     enddo

     end subroutine jonswap_wave_left_boundary


     subroutine wave_current_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,  &
                 Stokes_Drift,Fact,Phs_Lag(Nloc),Ucur,Ramp

     ! Current speed (ramp up flow speed)
     if(TRamp>0.0) then
       Ramp = tanh(TIME/TRamp)
     else
       Ramp = 1.0
     endif
     Ucur = Ramp*Cur_Wave

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
 55  Fk = (Segma-Wave_Number*Ucur)**2-Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = 2*(Segma-Wave_Number*Ucur)*(-Ucur)-Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)-  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
 65  continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = (Ucur+0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_X0(j,k) = (Ucur+0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

# if defined (SALINITY)
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Sin_X0(j,k) = Sali(Ibeg,j,k)
     enddo
     enddo
# endif

     end subroutine wave_current_left_boundary


     subroutine linear_wave_right_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Iend1,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
 55  Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
 65  continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_Xn(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend1,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_Xn(j)
       Uin_Xn(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_Xn(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_Xn(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

     end subroutine linear_wave_right_boundary


     subroutine linear_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

# if defined (SALINITY)
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Sin_X0(j,k) = Sali(Ibeg,j,k)
     enddo
     enddo
# endif

     end subroutine linear_wave_left_boundary


     subroutine solitary_wave_left_boundary
!------------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_bc
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: SP,Zero,pi,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Kbeg,Kend,sigc,BC_X0
     implicit none
     integer  :: j,k
     real(SP) :: Celerity,Atmp,Xstart,C2,D1,D2,D3,Zlev,Xc

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     
     Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
     do j = 1,Nloc
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       Xc = Xstart-Celerity*TIME
       Ein_X0(j) = Amp_Wave/cosh(Atmp*Xc)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
 
       C2 = sqrt(Grav*Dep_Wave)
       D2 = 2.0*Amp_Wave*Atmp**2*(2.0*cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**4
       D1 = -2.0*Amp_Wave*sinh(Atmp*Xc)*Atmp/cosh(Atmp*Xc)**3
       D3 = -8.0*Amp_Wave*sinh(Atmp*Xc)*Atmp**3*(cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**5
       do k = Kbeg,Kend
         Zlev = sigc(k)*Din_X0(j)
         Uin_X0(j,k) = C2*Ein_X0(j)/Dep_Wave*(1.0-1.0/4.0*Ein_X0(j)/Dep_Wave+Dep_Wave/3.0*(Dep_Wave/Ein_X0(j))*  &
                           (1.0-3.0/2.0*Zlev**2/Dep_Wave**2)*D2)
         Win_X0(j,k) = -C2*Zlev/Dep_Wave*((1.0-1.0/2.0*Ein_X0(j)/Dep_Wave)*D1+1.0/3.0*Dep_Wave**2*  &
                           (1.0-1.0/2.0*Zlev**2/Dep_Wave**2)*D3)
         Vin_X0(j,k) = 0.0
       enddo
     enddo

     end subroutine solitary_wave_left_boundary

	 ! changed by Cheng for fluid slide
# if defined (LANDSLIDE) || defined (FLUIDSLIDE)
     subroutine update_bathymetry
!------------------------------------------------------
!    This subroutine is used to update time-varying bathymetry 
!    Called by 
!       main 
!    Last update: 12/05/2011, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter
# if defined (LANDSLIDE)
     real(SP) :: alpha0,L0,T,bl,wl,e,kb,kw,x0,y0,xt,yt,zt,ut,t0,a0, &
           s0,st,slope0
!     ! added by Cheng for analytical dh/dt (need to comment code for numerical dh/dt)
!     real(SP) :: DeltS,Delt2S
     ! added by Cheng for 3D triangle landslide
!     real(SP) :: coef_a,coef_b,coef_c,xgb,xgn,xg,yg,ztm,alpha1
# endif
! added by Cheng for visvous slide
# if defined (FLUIDSLIDE)
     integer :: RKIT
# endif
     
     ! save old bathymetry
     Ho = Hc

# if defined (LANDSLIDE)
     if(SlideType(1:8)=='RIGID_2D') then
	 ! 2D landslide (added by Cheng for 2d landslide)
       e = 0.75
       slope0 = SlopeAngle*pi/180. 
       kb = 2.0*acosh(1.0/sqrt(e))/SlideL
       ut = SlideUt
       a0 = SlideA0
       t0 = ut/a0
       s0 = ut**2/a0
       st = s0*log(cosh(time/t0))*cos(slope0)
       x0 = SlideX0+st
       y0 = SlideY0
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         xt = xc(i)-x0
         zt = SlideT/(1-e)*((1.0/cosh(kb*xt))**2-e)
         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
       enddo
       enddo
       ! cosine shape (added by Cheng to smooth the angle out)
!       slope0 = SlopeAngle*pi/180. 
!       ut = SlideUt
!       a0 = SlideA0
!       t0 = ut/a0
!       s0 = ut**2/a0
!       st = s0*log(cosh(time/t0))*cos(slope0)
!       x0 = SlideX0+st
!       y0 = SlideY0
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         xt = xc(i)-x0
!		 zt = zero
!		 if (xt<=SlideL/2.0 .and. xt>=-SlideL/2.0) then
!           zt = SlideT/2.0*(1.0-cos(2.0*pi*(xt-SlideL/2.0)/SlideL))
!		 endif
!         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
!       enddo
!       enddo
     elseif(SlideType(1:8)=='RIGID_3D') then
	 ! 3D landslide
!	   ! modified by Cheng for analytical dh/dt
!	   DeltH = zero
!	   Delt2H = zero
       e = 0.717
       alpha0 = SlideAngle*pi/180.
       slope0 = SlopeAngle*pi/180. 
       kb = 2.0*acosh(1.0/e)/SlideL
       kw = 2.0*acosh(1.0/e)/SlideW
       ut = SlideUt
       a0 = SlideA0
       t0 = ut/a0
       s0 = ut**2/a0
       st = s0*log(cosh(time/t0))*cos(slope0)
!	   ! modified by Cheng for analytical dh/dt
!	   DeltS=s0/t0*tanh(time/t0)*cos(slope0)
!	   Delt2s=s0/(t0**2)/cosh(time/t0)**2*cos(slope0)
       x0 = SlideX0+st*cos(alpha0)
       y0 = SlideY0+st*sin(alpha0)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         xt = (xc(i)-x0)*cos(alpha0)+(yc(j)-y0)*sin(alpha0)
         yt = -(xc(i)-x0)*sin(alpha0)+(yc(j)-y0)*cos(alpha0)
         zt = SlideT/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yt)-e)
         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
!		 ! modified by Cheng for analytical dh/dt
!		 if (zt>0.0) then
!		   DeltH(i,j)=-(kb*(zt+e*SlideT/(1-e))*DeltS*tanh(kb*xt))
!		   Delt2H(i,j)=-(kb*(zt+e*SlideT/(1-e))*(Delt2S*tanh(kb*xt) &
!		               -kb*DeltS**2*(1-2.0*tanh(kb*xt)**2)))
!		 endif 
       enddo
       enddo
	 ! 3D triangle landslide (added by Cheng)
!       slope0 = SlopeAngle*pi/180. 
!	   ! subaerial case
!	   coef_a = -0.097588
!	   coef_b = 0.759361
!	   coef_c = 0.078776
!	   ! submerged case
!!	   coef_a = -0.085808
!!	   coef_b = 0.734798
!!	   coef_c = -0.034346
!	   st = (coef_a*time**3+coef_b*time**2+coef_c*time)*cos(slope0)
!       x0 = SlideX0+st
!       y0 = SlideY0
!!      xgb = SlideL*cos(slope0)
!!	   xgn = SlideL/cos(slope0)
!       xgb = SlideL-0.2
!	   xgn = SlideL+0.2
!	   xg = xgn-xgb
!	   ztm = xgb*tan(slope0)
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         if(xc(i)>=x0.and.xc(i)<=x0+xgb) then
!           zt = (xc(i)-x0)*tan(slope0)
!!		   alpha1 = slope0/2.0*(1.0+cos(pi*(xc(i)-x0-xgb)/xgb)) !cosine curve
!!		   yg = zt/2.0*tan(alpha1)
!           yg = 0.1
!		   if (yc(j)<=SlideW/2.0-yg) then
!		     zt = zt
!		   elseif (yc(j)>SlideW/2.0-yg.and.yc(j)<=SlideW/2.0+yg) then
!!		     zt= zt-(yc(j)-(SlideW/2.0-yg))/tan(alpha1)  ! line
!			 zt= zt/2.0*(1.0+cos(pi*(yc(j)-(SlideW/2.0-yg))/(2.0*yg)))  !cosine curve
!		   else
!             zt = 0.0
!		   endif
!         elseif(xc(i)>x0+xgb .and.xc(i)<=x0+xgn) then
!!		   zt = (xg-(xc(i)-x0-xgb))/tan(slope0)
!		   zt = ztm/2.0*(1.0+cos(pi*(xc(i)-x0-xgb)/xg))
!!		   alpha1 = slope0/2.0*(1.0+cos(pi*(xc(i)-x0-xgb)/xg))
!!		   yg = zt/2.0*tan(alpha1)
!           yg = 0.1
!		   if (yc(j)<=SlideW/2.0-yg) then
!		     zt = zt
!		   elseif (yc(j)>SlideW/2.0-yg.and.yc(j)<=SlideW/2.0+yg) then
!!		     zt= zt-(yc(j)-(SlideW/2.0-yg))/tan(alpha1)
!			 zt= zt/2.0*(1.0+cos(pi*(yc(j)-(SlideW/2.0-yg))/(2.0*yg)))
!		   else
!             zt = 0.0
!		   endif
!         else
!           zt = 0.0
!         endif
!         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
!       enddo
!       enddo
     endif
# endif

! added by Cheng for visvous slide
# if defined (FLUIDSLIDE)
     SldWZ = Eta !D-Hc
	 
     SldZs=SldZ0;   SldDZs=SldDZ0
     SldNs=SldN0;   SldMs=SldM0
     DO RKIT=1,3
       IF (RKIT.EQ.1) THEN
	 SLIDE_RKC1=1.0d0; SLIDE_RKC2=0.0d0; SLIDE_RKDT=dt
       ELSEIF (RKIT.EQ.2) THEN
	 SLIDE_RKC1=3.0d0/4.0d0; SLIDE_RKC2=1.0d0/4.0d0; SLIDE_RKDT=1.0d0/4.0d0*dt
       ELSEIF (RKIT.EQ.3) THEN
	 SLIDE_RKC1=1.0d0/3.0d0; SLIDE_RKC2=2.0d0/3.0d0; SLIDE_RKDT=2.0d0/3.0d0*dt
       ENDIF

       call nlmassF90();      call nlmmntF90()

       SldZs=SldZ1;   SldDZs=SldDZ1
       SldNs=SldN1;   SldMs=SldM1
     ENDDO
     
     SldZ0=SldZ1;   SldDZ0=SldDZ1
     SldN0=SldN1;   SldM0=SldM1

     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       write(*,*) 'DEPTH_TYPE has to be cell_center'
     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then

       ! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         Hc(i,j) = Hc0(i,j)-SldDZ1(i,j)
       enddo
       enddo
     endif ! end cell center
# endif

     ! ghost cells
     call phi_2D_coll(Hc)
	 ! modified by Cheng for analytical dh/dt
!	 call phi_2D_coll(DeltH)
!	 call phi_2D_coll(Delt2H)

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i+1,j))
     enddo
     Hfx(1,j) = Hfx(2,j)
     Hfx(Mloc1,j) = Hfx(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j+1))
     enddo
     Hfy(i,1) = Hfy(i,2)
     Hfy(i,Nloc1) = Hfy(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     ! time derivative of water depth
     DeltHo = DeltH

     DeltH = zero
     do j = 1,Nloc
     do i = 1,Mloc
       DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
     enddo
     enddo

     ! second-order time derivative
     if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

     end subroutine update_bathymetry
# endif
 
# if defined (LANDSLIDE_COMPREHENSIVE)

     subroutine update_bathymetry_comprehensive
!------------------------------------------------------
!    This subroutine is used to update time-varying bathymetry 
!    Called by 
!       main 
!    update: 12/05/2011, Gangfeng Ma 
!    update: 04/13/2012, Fengyan Shi
!    fyshi make a standard slide application    
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter,kslide
     real(SP) :: Hd,alpha0,L0,T,bl,wl,e,kb,kw,x0,x1,x2,xr,Hb(Mloc,Nloc), &
          xt,xt1,zt,zt1,ut,t0,a0,ht,s0,st,yt,Hslide(Mloc1,Nloc1)
     REAL(SP) :: edge_width
     real(SP) :: S_Slide, tmax_slide, del_slide, Sdel_slide, ks_slide !CO
     real(SP) :: CD_form,Cn_Coulomb,gamma_rho,Cf_skin,Cm,FF,F1,F2,F3
              
     ! save old bathymetry
     Ho = Hc

     ut = term_v(1)
     ut_slide = term_v(1)

   IF(ANA_BATHY)THEN
!      if use analytical slide, 1)need iteration 2) ac is time-dependent
     Hd=1.5  
     T = T_slide(1)
     bl = b_slide(1)
     wl = w_slide(1)
     e = e_slide(1)
     alpha0 = slope_slide(1)*3.1415926/180.0
     kb = 2.0*acosh(1.0/e)/bl
     kw = 2.0*acosh(1.0/e)/wl
     x0=x0_slide(1)
     xr = 0.0

     a0 = acceleration_lab
     t0 = ut/a0
     s0 = ut**2/a0
     st = s0*log(cosh(time/t0))
     L0 = x0/cos(alpha0)-T*tan(alpha0)+st
     x1 = (L0-bl/2.)*cos(alpha0)+xr
     x2 = (L0+bl/2.)*cos(alpha0)+xr

   ELSE

     IF(SLIDE_SHAPE_BOX) THEN  ! more conservative, usually for box-shape application

! ---

      DO kslide=1,NumSlides

       IF(TIME > SlideStartTime(Kslide))THEN

       CD_form=0.3   ! form drag
       Cn_Coulomb=tan(Coulomb_phi(kslide)*3.1415926/180.0)
       gamma_rho=2650.0/1000.0
       Cf_skin=0.005
       Cm=3.1415926*T_slide(kslide)/2.0/b_slide(kslide)

        FF=b_slide(kslide)*0.5-subaer_center2water(kslide)
        if (FF<0.0)FF=0.0
        if (FF>b_slide(kslide))FF=b_slide(kslide)
        F1=b_slide(kslide)*gamma_rho+Cm*FF
        F2=-0.5*(Cf_skin*FF/T_slide(kslide)+Cd_form)
        F3=(b_slide(kslide)*gamma_rho-FF)*grav*  &
           (sin(slope_slide(kslide)*pi/180.0)-  &
             Cn_Coulomb*cos(slope_slide(kslide)*pi/180.0))

        uslide(kslide)=uslide(kslide)+dt*1.0/F1*(F2*uslide(kslide)*uslide(kslide)+F3)
        subaer_center2water(kslide)=subaer_center2water(kslide)-uslide(kslide)*dt

        if(abs(subaer_center2water(kslide)- &
           subaer_center2water0(kslide))>Smax_slide(kslide)) uslide(kslide)=0.0
       
! ---

       if (uslide(kslide) > term_v(kslide)) uslide(kslide) = term_v(kslide)
  
       x0_slide(kslide) = x0_slide(kslide)+dt*uslide(kslide)*cos(alpha_slide(kslide))
       y0_slide(kslide) = y0_slide(kslide)+dt*uslide(kslide)*sin(alpha_slide(kslide))

       ENDIF ! end slide start time

      ENDDO ! end kslide

     ELSE ! Grilli slide

if (Slide) then !Slide CO
   if (SlideStop) then !CO Add Slide Stopping from Grilli and Watts (2005) appendix
      tmax_slide = ut_slide/ac_slide(1)*acosh(exp(Smax_slide(1)/(ut_slide**2/ac_slide(1))))
      del_slide = tmax_slide/(ut_slide/ac_slide(1))-tmax_slide/(ut_slide/ac_slide(1))*0.1
      Sdel_slide = ut_slide**2/ac_slide(1)*log(cosh(del_slide))
      ks_slide = (ut_slide*(tanh(del_slide)))/(Smax_slide(1)-Sdel_slide)

      S_Slide = ut_slide**2/ac_slide(1)*log(cosh(time*ac_slide(1)/ut_slide))*COS(slope_slide(1)*pi/180)
      if (time > del_slide*ut_slide/ac_slide(1)) then 
      S_Slide = Sdel_slide+(Smax_slide(1)-Sdel_slide)*  &
          (1-exp(-ks_slide*(time-del_slide*ut_slide/ac_slide(1))))  &
          *COS(slope_slide(1)*pi/180)
      endif

   else ! Continuous slide movement
      S_Slide = ut_slide**2/ac_slide(1)*log(cosh(time*ac_slide(1)/ut_slide))*COS(slope_slide(1)*pi/180)
   endif


else !Slump CO
!      uslide = 3.1415926/2.0*sf_slide/tf_slide*sin(time*3.1415926/tf_slide)
      S_Slide = sf_slide/2*(1-cos(time*3.1415926/tf_slide))*COS(slope_slide(1)*pi/180)
     
      if (time > tf_slide)then
      S_Slide = sf_slide*COS(slope_slide(1)*pi/180)
      else
      endif
endif
     x0_slide(1) = XX0_Slide+S_Slide*cos(alpha_slide(1)) !CO
     y0_slide(1) = YY0_Slide+S_Slide*sin(alpha_slide(1)) !CO


     ENDIF ! conservative vs. grilli

   ENDIF

     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then

# if defined (PARALLEL)
       if(myid.eq.0) write(*,*) 'DEPTH_TYPE has to be cell_center'
# else
       write(*,*) 'DEPTH_TYPE has to be cell_center'
# endif       

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! base bathymetry 
   IF(ANA_BATHY)THEN   ! lab analytical slide

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)<=xr) then
           Hb(i,j) = -(xr-xc(i))*tan(alpha0)
         elseif(xc(i)<=Hd/tan(alpha0)+xr) then
           Hb(i,j) = (xc(i)-xr)*tan(alpha0)
         else
           Hb(i,j) = Hd
         endif
         ! temporarily no runup
       !  Hb(i,j) = max(0.01,Hb(i,j))
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)<=x1.or.xc(i)>=x2.or.yc(j)>=wl/2.) then
           Hc(i,j) = Hb(i,j)
         else
           iter = 1
           xt = (xc(i)-xr)/cos(alpha0)-L0
           zt = T/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yc(j))-e)
      60   xt1 = ((xc(i)-xr)/cos(alpha0)-zt*tan(alpha0))-L0
           zt1 = T/(1-e)*(1.0/cosh(kb*xt1)/cosh(kw*yc(j))-e)
           if(abs(zt1-zt)/abs(zt)>1.e-8) then
             iter = iter+1
             if(iter>20) write(*,*) 'too many iterations!'
             zt = zt1
             goto 60
           endif
           Hc(i,j) = Hb(i,j)-max(0.0,zt1)/cos(alpha0)  
         endif
       enddo
       enddo




   ELSEIF(SLIDE_SHAPE_BOX) THEN  ! box shape
       Hb=DepC0
       Hc=DepC0

       ! add landslide
! to keep consistent with slide bathy at grid, assuming slide bathy at grid point

       DO kslide = 1, NumSlides

       edge_width=T_slide(KSLIDE)
       Hslide=0.0_SP

       do j = Jbeg,Jend1
       do i = Ibeg,Iend1
         xt = (x(i)-x0_slide(KSLIDE))*cos(alpha_slide(KSLIDE))+(y(j)-y0_slide(KSLIDE))*sin(alpha_slide(KSLIDE))
         yt = -(x(i)-x0_slide(KSLIDE))*sin(alpha_slide(KSLIDE))+(y(j)-y0_slide(KSLIDE))*cos(alpha_slide(KSLIDE))

         IF(abs(xt)<=0.5*b_slide(KSLIDE)-0.5*edge_width.AND.abs(yt)<=0.5*w_slide(KSLIDE)  &
                 +SlideErrorRange)THEN
            ht=T_slide(KSLIDE)
         ELSEIF(abs(xt)>0.5*b_slide(KSLIDE)-0.5*edge_width  &
            .AND.abs(xt)<0.5*b_slide(KSLIDE)+0.5*edge_width.AND.abs(yt)<=0.5*w_slide(KSLIDE)  &
                 +SlideErrorRange)THEN
            ht=T_slide(KSLIDE)-T_slide(KSLIDE)*(abs(xt)-0.5*b_slide(KSLIDE)+0.5*edge_width)/edge_width
         ELSE
            ht=0.0_SP
         ENDIF
         Hslide(I,J) = max(0.0,ht)
       enddo
       enddo

! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         ht=0.25*(Hslide(I,J)+Hslide(I+1,J)+Hslide(I,J+1) &
                            +Hslide(I+1,J+1))
!         Hc(i,j) = Hb(i,j)-ht
         Hc(i,j) = Hc(i,j)-ht
       enddo
       enddo

      ENDDO ! end Kslide



   ELSE  ! ellipse shape (original in version 1.1)
       Hb=DepC0

       ! add landslide
! to keep consistent with slide bathy at grid, assuming slide bathy at grid point

       do j = Jbeg,Jend1
       do i = Ibeg,Iend1
         xt = (x(i)-x0_slide(1))*cos(alpha_slide(1))+(y(j)-y0_slide(1))*sin(alpha_slide(1))
         yt = -(x(i)-x0_slide(1))*sin(alpha_slide(1))+(y(j)-y0_slide(1))*cos(alpha_slide(1))
         ht = T_slide(1)/(1-e_slide(1))*(1./cosh(kb_slide*xt)/cosh(kw_slide*yt)-e_slide(1))
         Hslide(I,J) = max(0.0,ht)
       enddo
       enddo

! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         ht=0.25*(Hslide(I,J)+Hslide(I+1,J)+Hslide(I,J+1) &
                            +Hslide(I+1,J+1))
         Hc(i,j) = Hb(i,j)-ht
       enddo
       enddo

    ENDIF  ! end analytical slide

     endif ! end cell center

     ! ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hfx(2,j)
     Hfx(Mloc1,j) = Hfx(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hfy(i,2)
     Hfy(i,Nloc1) = Hfy(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     ! time derivative of water depth
     DeltHo = DeltH

     DeltH = zero
     do j = 1,Nloc
     do i = 1,Mloc
       DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
     enddo
     enddo

     ! second-order time derivative
     if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

     end subroutine update_bathymetry_comprehensive

# endif
! end landslide comprehensive
 
     subroutine RandomU
!-------------------------------------------------
!    This subroutine is used to generate random
!    perturbation to the initial velocity field.
!    Called by
!       initial
!    Last update: 16/10/2015, Gangfeng Ma
!-------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: fac,randx,randy,old,seed

     seed = 11.0
     old = seed

     fac = 0.1
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       randx = mod((57*old+1), 256.0)
       randy = mod((57*randx+1), 256.0)
       randy = randy/256.0
       old = randy

       U(i,j,k) = (1.0+0.1*(randy-0.5+0.1*myid))*U(i,j,k)
       V(i,j,k) = (randy-0.5+0.01*myid)*0.01
       W(i,j,k) = (randy-0.5+0.01*myid)*0.01
     enddo
     enddo
     enddo

     return
     end subroutine RandomU

     subroutine calculate_sponge
!-------------------------------------------------
!    Calculate sponge function
!    Called by
!      initial
!    Last update: 12/02/2011, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j

     if(Sponge_West_Width>Zero) then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)<=Sponge_West_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((xc(i)-Sponge_West_Width)/Sponge_West_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_East_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)>=Mglob*dx-Sponge_East_Width) then
           Sponge(i,j) = sqrt(1.0-  &
             min(((xc(i)-(Mglob*dx-Sponge_East_Width))/Sponge_East_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_South_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)<=Sponge_South_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-Sponge_South_Width)/Sponge_South_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_North_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)>=Nglob*dy-Sponge_North_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-(Nglob*dy-Sponge_North_Width))/Sponge_North_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     end subroutine calculate_sponge


     subroutine sponge_damping
!---------------------------------------------------
!    This subroutine is used to damp waves using DHI type
!    sponge layer variables
!    Called by 
!      main
!    Last update: 12/02/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: SP,TIME,Eta,Hc,D,U,V,W,Omega,Sponge,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW,TRamp,Cur_Wave
     implicit none
     integer :: i,j,k
     real(SP) :: Ramp,Ucur

     ! Current speed (ramp up flow speed)                                                                                                                   
     if(TRamp>0.0) then
       Ramp = tanh(TIME/TRamp)
     else
       Ramp = 1.0
     endif
     Ucur = Ramp*Cur_Wave

     do j = 1,Nloc
     do i = 1,Mloc
       if(Mask(i,j)==1) then
         Eta(i,j) = Eta(i,j)*Sponge(i,j)
         D(i,j) = Eta(i,j)+Hc(i,j)

         ! W is estimated from continuity equation
         do k = 1,Kloc
           U(i,j,k) = (U(i,j,k)-Ucur)*Sponge(i,j)+Ucur
           V(i,j,k) = V(i,j,k)*Sponge(i,j)
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
         enddo
       endif
     enddo
     enddo
     
     end subroutine sponge_damping

 
     subroutine sigma_transform
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!      eval_duvw
!    Last update: 29/03/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc, &
                       DelxSl,DelySl,Kloc1,sig
     implicit none
     integer :: i,j,k

     DelxSc = Zero
     DelySc = Zero
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j)
       DelySc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo

     ! added by Morteza
     DelxSl = Zero
     DelySl = Zero
     do k = 1,Kloc1
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSl(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)                                      
       DelySl(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)                                      
     enddo
     enddo
     enddo

     end subroutine sigma_transform


     subroutine eval_turb(ISTEP)
!---------------------------------------------------
!    This subroutine is used to calculate viscosity
!    Called by                                                                                                             
!      main
!    Last update: 21/06/2011, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer, intent(in) :: ISTEP
     integer :: i,j,k,n,m
     real(SP) :: DelsU,DelsV,Strxx,Stryy,Strxy,StrainMag,Smax
     real(SP), dimension(3,3) :: VelGrad,Stress
     real(SP), dimension(Mloc,Nloc,Kloc) :: StressMag

     ! laminar viscosity
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
# if defined (SEDIMENT)
       if(trim(Sed_Type)=='COHESIVE') then
         if(Conc(i,j,k)<0.1) then
           Cmu(i,j,k) = Visc
         else
           Cmu(i,j,k) = Mud_Visc
         endif
       else
         Cmu(i,j,k) = Visc
       endif
# else
       Cmu(i,j,k) = Visc
# endif
     enddo
     enddo
     enddo

# if defined (SEDIMENT)
     ! include rheology effects
     if(RHEOLOGY_ON) then
       StressMag = Zero
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         ! estimate gradient first 
         VelGrad = Zero; Stress = Zero
         VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)            
         VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)          
         VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)              
         VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)               
         VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                  
         VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                  
         VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

         ! strain magnitude
         StrainMag = Zero
         do n = 1,3
         do m = 1,3
           StrainMag = StrainMag+0.5*(0.5*(VelGrad(n,m)+VelGrad(m,n)))**2                                                         
         enddo
         enddo

         ! shear stress
         if(IsMove(i,j,k)==1) then
           do n = 1,3
           do m = 1,3
             Stress(n,m) = 2.0*Rho(i,j,k)*(Cmu(i,j,k)+CmuR(i,j,k)+CmuVt(i,j,k))*(0.5*(VelGrad(n,m)+velGrad(m,n)))
             StressMag(i,j,k) = StressMag(i,j,k)+0.5*Stress(n,m)*Stress(n,m)
           enddo
           enddo
         endif

         if(Conc(i,j,k)>0.1) then
           if(Nglob==1) then  ! 2D simulation
             Smax = dmax1(StressMag(i-1,j,k),StressMag(i+1,j,k),StressMag(i,j,k-1),StressMag(i,j,k+1))
           else
             Smax = dmax1(StressMag(i-1,j,k),StressMag(i+1,j,k),StressMag(i,j,k-1),StressMag(i,j,k+1),  &
                       StressMag(i,j-1,k),StressMag(i,j+1,k))
           endif

           if(IsMove(i,j,k)==0) then  ! originally not moving
             if(sqrt(Smax)<Yield_Stress) then ! still not moving 
               CmuR(i,j,k) = 1.0e+10  ! a large number                                                                              
             else ! start moving
               IsMove(i,j,k) = 1
               CmuR(i,j,k) = Plastic_Visc+Yield_Stress/StrainMag/Rho(i,j,k)                                                          
             endif
           else  ! originally moving
             if(sqrt(StressMag(i,j,k))<Yield_Stress) then  ! stop moving
               IsMove(i,j,k) = 0
               CmuR(i,j,k) = 1.e+10
             else
               CmuR(i,j,k) = Plastic_Visc+Yield_Stress/StrainMag/Rho(i,j,k)
             endif
           endif             
         else
           CmuR(i,j,k) = 0.0
         endif
       enddo
       enddo
       enddo
     endif
# endif

     ! vertical turbulent viscosity
     if(IVturb==1) then
       ! constant vertical viscosity
       CmuVt = Cvs
     elseif(IVturb==2) then
       ! subgrid model
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 2,Kloc-1
         DelsU = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         DelsV = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         CmuVt(i,j,k) = Cvs*D(i,j)*dsig(k)**2*sqrt(DelsU**2+DelsV**2)
       enddo
       CmuVt(i,j,1) = CmuVt(i,j,2)
       CmuVt(i,j,Kloc) = CmuVt(i,j,Kloc-1)
       enddo
       enddo
     elseif(IVturb==3) then
       ! k-epsilon turbulence model
       call kepsilon(ISTEP)
     elseif(IVturb==10) then
       ! 3D turbulence model
# if defined (VEGETATION)
       call kepsilon_veg(ISTEP)
# else
       call kepsilon_3D(ISTEP)
# endif
     elseif(IVturb==20) then
       call les_3D(ISTEP)
     elseif(IVturb==30) then
       call les_dyn(ISTEP)
     endif
     
     ! horizontal turbulent viscosity
     if(IHturb==1) then
       ! constant viscosity
       CmuHt = Chs
     elseif(IHturb==2) then
       ! subgrid model
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         Strxx = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
         Stryy = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
         Strxy = 0.5*((U(i,j+1,j)-U(i,j-1,k))/(2.0*dy)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
               (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
         CmuHt(i,j,k) = Chs*dx*dy*sqrt(Strxx**2+2.0*Strxy**2+Stryy**2)
       enddo
       enddo
       enddo

       ! ghost cell
# if defined (PARALLEL)
       call phi_3D_exch(CmuHt)
# endif

# if defined (PARALLEL)
       if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Ibeg-i,j,k) = CmuHt(Ibeg+i-1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Iend+i,j,k) = CmuHt(Iend-i+1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_suth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jbeg-j,k) = CmuHt(i,Jbeg+j-1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jend+j,k) = CmuHt(i,Jend-j+1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         do k = 1,Nghost
           CmuHt(i,j,Kbeg-k) = CmuHt(i,j,Kbeg+k-1)
         enddo
         do k = 1,Nghost
           CmuHt(i,j,Kend+k) = CmuHt(i,j,Kend-k+1)
         enddo
       enddo
       enddo  
     elseif(IHturb>=10) then
       ! use 3D turbulence model
       ! in this case, the length scales in all directions are
       ! in the same order
       CmuHt = CmuVt
     endif

     end subroutine eval_turb


     subroutine diffusion
!---------------------------------------------------------  
!    This subroutine is used to evaluate diffusion terms
!    in simplified form following FVCOM
!    Called by 
!      eval_duvw 
!    Last update: 20/12/2011, Gangfeng Ma
!
!    M.Derakhti rewrittened this subroutine, including the exact 
!    diffusion terms based on UD,UV,UW as in Derakhti etal 2015a,b.
!    Note that the part of each terms that has only vertical gradient
!    moved to LHS and treated implicitly. Terms with mixed or 
!    purely horizontal gradient are calculated here.
!    The current formulation is based on a non-uniform vertical grid
!
!    Last update: 24/05/2015, Morteza Derakhti
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     !added by M.Derakhti
     real(SP) ::L1top,L2top,L1bot,L2bot,alpha_c,beta_c,gamma_c, &
                dsigck,dsigck1,nuH_top,nuH_bot,nuV_top,nuV_bot
       
     Diffxx = zero; Diffxy = zero; Diffxz = zero
     Diffyx = zero; Diffyy = zero; Diffyz = zero
     Diffzx = zero; Diffzy = zero
     do k = Kbeg,Kend
       !linear interpolation using lagrange polynominal to get
       !values at k+1/2 and k-1/2
       !We have a uniform grid in x and y direction 
       !so the values at the i+1/2,i-1/2
       !and j+1/2, j-1/2 can be obtained using a simple averaging
       L1top = dsig(k+1)/(dsig(k)+dsig(k+1))
       L2top = dsig(k  )/(dsig(k)+dsig(k+1))
       L1bot = dsig(k  )/(dsig(k)+dsig(k-1))
       L2bot = dsig(k-1)/(dsig(k)+dsig(k-1))
       !these are used for vertical gradient in non-uniform grid
       !see Derakhti etal 2015a Appendix B
       dsigck  = (dsig(k  )+dsig(k+1))/2.0 
       dsigck1 = (dsig(k-1)+dsig(k  ))/2.0 
       alpha_c = -dsigck/(dsigck+dsigck1)/dsigck1
       beta_c = (dsigck-dsigck1)/(dsigck*dsigck1)
       gamma_c = dsigck1/(dsigck+dsigck1)/dsigck

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle
       !viscosity at k+1/2 and k-1/2 using linear interpolation method
       nuH_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+  &
                 (L1top*CmuHt(i,j,k)+L2top*CmuHt(i,j,k+1))/Schmidt
       nuV_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+  &
                 (L1top*CmuVt(i,j,k)+L2top*CmuVt(i,j,k+1))/Schmidt
       nuH_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+  &
                 (L1bot*CmuHt(i,j,k-1)+L2bot*CmuHt(i,j,k))/Schmidt
       nuV_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+  &
                 (L1bot*CmuVt(i,j,k-1)+L2bot*CmuVt(i,j,k))/Schmidt

       !Diffxx(i,j,k) = ((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)* &
       !        (D(i+1,j)+D(i,j))*(U(i+1,j,k)-U(i,j,k))-  &
       !        (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))* &
       !        (U(i,j,k)-U(i-1,j,k)))/dx**2
       !Reology effects? CmuR shold be added, I guess.

       Diffxx(i,j,k) = ( 2.0*(0.5*(Cmu(i+1,j,k)+Cmu(i,j,k  ))+  &
                         0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k  ))/Schmidt) &
                         *(DU(i+1,j,k)-DU(i,j,k))/dx &
                        -2.0*(0.5*(Cmu(i,j,k  )+Cmu(i-1,j,k))+  &
                         0.5*(CmuHt(i,j,k  )+CmuHt(i-1,j,k))/Schmidt) &
                         *(DU(i,j,k)-DU(i-1,j,k))/dx &
                       ) /dx & !d/dx[2*nu*d/dx[DU]]
                     + ( 2.0*(Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt) &
                         *(alpha_c*DU(i+1,j,k-1)*DelxSc(i+1,j,k-1)&
                          +beta_c *DU(i+1,j,k  )*DelxSc(i+1,j,k  )&
                          +gamma_c*DU(i+1,j,k+1)*DelxSc(i+1,j,k+1)) & 
                        -2.0*(Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt) &
                         *(alpha_c*DU(i-1,j,k-1)*DelxSc(i-1,j,k-1)&
                          +beta_c *DU(i-1,j,k  )*DelxSc(i-1,j,k  )&
                          +gamma_c*DU(i-1,j,k+1)*DelxSc(i-1,j,k+1)) & 
                       ) /2.0/dx & !d/dx[2*nu*d/dSc[DU*DelxSc]]
                     + ( alpha_c*2.0*(Cmu(i,j,k-1)+  &
                          CmuHt(i,j,k-1)/Schmidt)*DelxSc(i,j,k-1)&
                         *(DU(i+1,j,k-1)-DU(i-1,j,k-1))/2.0/dx &
                        +beta_c *2.0*(Cmu(i,j,k  )+  &
                         CmuHt(i,j,k  )/Schmidt)*DelxSc(i,j,k  )&
                         *(DU(i+1,j,k  )-DU(i-1,j,k  ))/2.0/dx &
                        +gamma_c*2.0*(Cmu(i,j,k+1)+  &
                         CmuHt(i,j,k+1)/Schmidt)*DelxSc(i,j,k+1)&
                         *(DU(i+1,j,k+1)-DU(i-1,j,k+1))/2.0/dx &
                        )  !d/dSc[2*nu*DelxSc*d/dx[DU]]

       !Diffxy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
       !        (D(i,j+1)+D(i,j))*(U(i,j+1,k)-U(i,j,k))- &
       !        (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
       !        (U(i,j,k)-U(i,j-1,k)))/dy**2+  &
       !        ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*D(i,j+1)*  &
       !        (V(i+1,j+1,k)-V(i-1,j+1,k))-(Cmu(i,j-1,k)+  &
       !        CmuHt(i,j-1,k)/Schmidt)*D(i,j-1)*  &
       !        (V(i+1,j-1,k)-V(i-1,j-1,k)))/(4.0*dx*dy)
       Diffxy(i,j,k) = ( (0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
                          0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt) &
                         *(DU(i,j+1,k)-DU(i,j,k))/dy &
                        -(0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
                          0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/Schmidt) &
                         *(DU(i,j,k)-DU(i,j-1,k))/dy &
                        ) /dy & !d/dy[2*nu*d/dy[DU]]
                     + ( (Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt) &
                         *( (DV(i+1,j+1,k)-DV(i-1,j+1,k))/2.0/dx &
                             +alpha_c*(DelxSc(i,j+1,k-1)*DV(i,j+1,k-1)+  &
                              DelySc(i,j+1,k-1)*DU(i,j+1,k-1))&
                             +beta_c *(DelxSc(i,j+1,k  )*DV(i,j+1,k  )+  &
                              DelySc(i,j+1,k  )*DU(i,j+1,k  ))&
                             +gamma_c*(DelxSc(i,j+1,k+1)*DV(i,j+1,k+1)+  &
                              DelySc(i,j+1,k+1)*DU(i,j+1,k+1))&
                          )&
                        -(Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt) &
                         *( (DV(i+1,j-1,k)-DV(i-1,j-1,k))/2.0/dx &
                             +alpha_c*(DelxSc(i,j-1,k-1)*DV(i,j-1,k-1)+  &
                              DelySc(i,j-1,k-1)*DU(i,j-1,k-1))&
                             +beta_c *(DelxSc(i,j-1,k  )*DV(i,j-1,k  )+  &
                              DelySc(i,j-1,k  )*DU(i,j-1,k  ))&
                             +gamma_c*(DelxSc(i,j-1,k+1)*DV(i,j-1,k+1)+  &
                              DelySc(i,j-1,k+1)*DU(i,j-1,k+1))&
                          )&
                        ) /2.0/dy & !d/dy[nu*(d/dx[DV]+d/dSc[DelxSc*DV+DelySc*DU])]
                      + ( alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/Schmidt)*DelySc(i,j,k-1)&
                                 *((DU(i,j+1,k-1)-DU(i,j-1,k-1))/2.0/dy+  &
                                  (DV(i+1,j,k-1)-DV(i-1,j,k-1))/2.0/dx) &
                         +beta_c *(Cmu(i,j,k  )+CmuHt(i,j,k  )/Schmidt)*DelySc(i,j,k  )&
                                 *((DU(i,j+1,k  )-DU(i,j-1,k  ))/2.0/dy+  &
                                  (DV(i+1,j,k  )-DV(i-1,j,k  ))/2.0/dx) &
                         +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/Schmidt)*DelySc(i,j,k+1)&
                                 *((DU(i,j+1,k+1)-DU(i,j-1,k+1))/2.0/dy+  &
                                  (DV(i+1,j,k+1)-DV(i-1,j,k+1))/2.0/dx) &
                         ) & !d/dSc[nu*DelySc*(d\dy[DU]+d\dx[DV])]
                      + ( ( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*alpha_c/dsig(k)&
                           +(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck1/(dsigck+dsigck1)&
                          )*DelxSc(i,j,k-1)*DV(i,j,k-1) &
                         +( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*beta_c /dsig(k)&
                           -(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck1/dsigck&
                          )*DelxSc(i,j,k  )*DV(i,j,k  ) &
                         +( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*gamma_c/dsig(k)&
                           +(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck /(dsigck+dsigck1)&
                          )*DelxSc(i,j,k+1)*DV(i,j,k+1) &
                        ) !d\dSc[nu*DelySc*d\dSc[\DelxSc*DV]]
       Diffxz(i,j,k) =  ( alpha_c*(Cmu(i,j,k-1)+CmuVt(i,j,k-1)/Schmidt)/D(i,j)&
                                 *(DW(i+1,j,k-1)-DW(i-1,j,k-1))/2.0/dx &
                         +beta_c*(Cmu(i,j,k)+CmuVt(i,j,k)/Schmidt)/D(i,j)&
                                 *(DW(i+1,j,k)-DW(i-1,j,k))/2.0/dx &
                         +gamma_c*(Cmu(i,j,k+1)+CmuVt(i,j,k+1)/Schmidt)/D(i,j)&
                                 *(DW(i+1,j,k+1)-DW(i-1,j,k+1))/2.0/dx &
                        ) & !d/dSc[nu/D*d/dx[WD]]
                      + ( ( (nuV_top-nuV_bot)/D(i,j)*alpha_c/dsig(k)&
                           +(nuV_top+nuV_bot)/D(i,j)/dsigck1/(dsigck+dsigck1)&
                          )*DelxSc(i,j,k-1)*DW(i,j,k-1) &
                         +( (nuV_top-nuV_bot)/D(i,j)*beta_c /dsig(k)&
                           -(nuV_top+nuV_bot)/D(i,j)/dsigck1/dsigck&
                          )*DelxSc(i,j,k  )*DW(i,j,k  ) &
                         +( (nuV_top-nuV_bot)/D(i,j)*gamma_c/dsig(k)&
                           +(nuV_top+nuV_bot)/D(i,j)/dsigck /(dsigck+dsigck1)&
                          )*DelxSc(i,j,k+1)*DW(i,j,k+1) &
                        ) !d\dSc[nu/D*d\dSc[\DelxSc*DW]]

       !Diffyx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
       !        (D(i+1,j)+D(i,j))*(V(i+1,j,k)-V(i,j,k))-  & 
       !        (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
       !        (V(i,j,k)-V(i-1,j,k)))/dx**2+  &
       !        ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*D(i+1,j)*  &
       !        (U(i+1,j+1,k)-U(i+1,j-1,k))-(Cmu(i-1,j,k)+  &
       !        CmuHt(i-1,j,k)/Schmidt)*D(i-1,j)*  &
       !        (U(i-1,j+1,k)-U(i-1,j-1,k)))/(4.0*dx*dy)
       Diffyx(i,j,k) = ( (0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt) &
                         *(DV(i+1,j,k)-DV(i,j,k))/dx &
                        -(0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/Schmidt) &
                         *(DV(i,j,k)-DV(i-1,j,k))/dx &
                        ) /dx & !d/dx[2*nu*d/dx[DV]]
                     + ( (Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt) &
                         *( (DU(i+1,j+1,k)-DU(i+1,j-1,k))/2.0/dy &
                            +alpha_c*(DelxSc(i+1,j,k-1)*DV(i+1,j,k-1)+DelySc(i+1,j,k-1)*DU(i+1,j,k-1))&
                            +beta_c *(DelxSc(i+1,j,k  )*DV(i+1,j,k  )+DelySc(i+1,j,k  )*DU(i+1,j,k  ))&
                            +gamma_c*(DelxSc(i+1,j,k+1)*DV(i+1,j,k+1)+DelySc(i+1,j,k+1)*DU(i+1,j,k+1))&
                          )&
                        -(Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt) &
                         *( (DU(i-1,j+1,k)-DU(i-1,j-1,k))/2.0/dy &
                            +alpha_c*(DelxSc(i-1,j,k-1)*DV(i-1,j,k-1)+DelySc(i-1,j,k-1)*DU(i-1,j,k-1))&
                            +beta_c *(DelxSc(i-1,j,k  )*DV(i-1,j,k  )+DelySc(i-1,j,k  )*DU(i-1,j,k  ))&
                            +gamma_c*(DelxSc(i-1,j,k+1)*DV(i-1,j,k+1)+DelySc(i-1,j,k+1)*DU(i-1,j,k+1))&
                          )&
                        ) /2.0/dx & !d/dx[nu*(d/dy[DU]+d/dSc[DelxSc*DV+DelySc*DU])]
                      + ( alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/Schmidt)*DelxSc(i,j,k-1)&
                                 *((DU(i,j+1,k-1)-DU(i,j-1,k-1))/2.0/dy+  &
                                   (DV(i+1,j,k-1)-DV(i-1,j,k-1))/2.0/dx) &
                         +beta_c *(Cmu(i,j,k  )+CmuHt(i,j,k  )/Schmidt)*DelxSc(i,j,k  )&
                                 *((DU(i,j+1,k  )-DU(i,j-1,k  ))/2.0/dy+  &
                                  (DV(i+1,j,k  )-DV(i-1,j,k  ))/2.0/dx) &
                         +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/Schmidt)*DelxSc(i,j,k+1)&
                                 *((DU(i,j+1,k+1)-DU(i,j-1,k+1))/2.0/dy+  &
                                    (DV(i+1,j,k+1)-DV(i-1,j,k+1))/2.0/dx) &
                         ) & !d/dSc[nu*DelxSc*(d\dy[DU]+d\dx[DV])]
                      + ( ( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*alpha_c/dsig(k)&
                           +(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck1/(dsigck+dsigck1)&
                          )*DelySc(i,j,k-1)*DU(i,j,k-1) &
                         +( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*beta_c /dsig(k)&
                           -(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck1/dsigck&
                          )*DelySc(i,j,k  )*DU(i,j,k  ) &
                         +( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*gamma_c/dsig(k)&
                           +(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck /(dsigck+dsigck1)&
                          )*DelySc(i,j,k+1)*DU(i,j,k+1) &
                        )!d\dSc[nu*DelxSc*d\dSc[\DelySc*DU]]

       !Diffyy(i,j,k) = ((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
       !        (D(i,j+1)+D(i,j))*(V(i,j+1,k)-V(i,j,k))-  & 
       !        (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
       !        (V(i,j,k)-V(i,j-1,k)))/dy**2
       Diffyy(i,j,k) = ( 2.0*(0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt) &
                         *(DV(i,j+1,k)-DV(i,j,k))/dy &
                        -2.0*(0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/Schmidt) &
                         *(DV(i,j,k)-DV(i,j-1,k))/dy &
                       ) /dy & !d/dy[2*nu*d/dy[DV]]
                     + ( 2.0*(Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt) &
                         *(alpha_c*DV(i,j+1,k-1)*DelySc(i,j+1,k-1)&
                          +beta_c *DV(i,j+1,k)*DelySc(i,j+1,k)&
                          +gamma_c*DV(i,j+1,k+1)*DelySc(i,j+1,k+1)) & 
                        -2.0*(Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt) &
                         *(alpha_c*DV(i,j-1,k-1)*DelySc(i,j-1,k-1)&
                          +beta_c *DV(i,j-1,k)*DelySc(i,j-1,k)&
                          +gamma_c*DV(i,j-1,k+1)*DelySc(i,j-1,k+1)) & 
                       ) /2.0/dy & !d/dy[2*nu*d/dSc[DV*DelySc]]
                     + (alpha_c*2.0*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/Schmidt)*DelySc(i,j,k-1)&
                         *(DV(i,j+1,k-1)-DV(i,j-1,k-1))/2.0/dy &
                        +beta_c*2.0*(Cmu(i,j,k)+CmuHt(i,j,k)/Schmidt)*DelySc(i,j,k)&
                         *(DV(i,j+1,k)-DV(i,j-1,k))/2.0/dy &
                        +gamma_c*2.0*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/Schmidt)*DelySc(i,j,k+1)&
                         *(DV(i,j+1,k+1)-DV(i,j-1,k+1))/2.0/dy &
                        )  !d/dSc[2*nu*DelySc*d/dy[DV]]

       Diffyz(i,j,k) =  ( alpha_c*(Cmu(i,j,k-1)+CmuVt(i,j,k-1)/Schmidt)/D(i,j)&
                                 *(DW(i,j+1,k-1)-DW(i,j-1,k-1))/2.0/dy &
                         +beta_c*(Cmu(i,j,k)+CmuVt(i,j,k)/Schmidt)/D(i,j)&
                                 *(DW(i,j+1,k)-DW(i,j-1,k))/2.0/dy &
                         +gamma_c*(Cmu(i,j,k+1)+CmuVt(i,j,k+1)/Schmidt)/D(i,j)&
                                 *(DW(i,j+1,k+1)-DW(i,j-1,k+1))/2.0/dy &
                        ) & !d/dSc[nu/D*d/dx[WD]]
                      + ( ( (nuV_top-nuV_bot)/D(i,j)*alpha_c/dsig(k)&
                           +(nuV_top+nuV_bot)/D(i,j)/dsigck1/(dsigck+dsigck1)&
                          )*DelySc(i,j,k-1)*DW(i,j,k-1) &
                         +( (nuV_top-nuV_bot)/D(i,j)*beta_c /dsig(k)&
                           -(nuV_top+nuV_bot)/D(i,j)/dsigck1/dsigck&
                          )*DelySc(i,j,k  )*DW(i,j,k  ) &
                         +( (nuV_top-nuV_bot)/D(i,j)*gamma_c/dsig(k)&
                           +(nuV_top+nuV_bot)/D(i,j)/dsigck /(dsigck+dsigck1)&
                          )*DelySc(i,j,k+1)*DW(i,j,k+1) &
                        ) !d\dSc[nu/D*d\dSc[\DelxSc*DW]]

       !Diffzx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
       !        (D(i+1,j)+D(i,j))*(W(i+1,j,k)-W(i,j,k))-  &
       !        (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
       !        (W(i,j,k)-W(i-1,j,k)))/dx**2+  &
       !        ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*  &
       !        (U(i+1,j,k+1)-U(i+1,j,k-1))/(sigc(k+1)-sigc(k-1))-  &
       !        (Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt)*  &
       !        (U(i-1,j,k+1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dx)
       Diffzx(i,j,k) = ( (0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt) &
                         *(DW(i+1,j,k)-DW(i,j,k))/dx &
                        -(0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/Schmidt) &
                         *(DW(i,j,k)-DW(i-1,j,k))/dx &
                       ) /dx & !d/dx[nu*d/dx[DW]]
                     + (  (Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt) &
                         *(alpha_c*(DW(i+1,j,k-1)*DelxSc(i+1,j,k-1)+DU(i+1,j,k-1)/D(i,j))&
                          +beta_c *(DW(i+1,j,k  )*DelxSc(i+1,j,k  )+DU(i+1,j,k  )/D(i,j))&
                          +gamma_c*(DW(i+1,j,k+1)*DelxSc(i+1,j,k+1)+DU(i+1,j,k+1)/D(i,j))&
                          )&
                        - (Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt) &
                         *(alpha_c*(DW(i-1,j,k-1)*DelxSc(i-1,j,k-1)+DU(i-1,j,k-1)/D(i,j))&
                          +beta_c *(DW(i-1,j,k  )*DelxSc(i-1,j,k  )+DU(i-1,j,k  )/D(i,j))&
                          +gamma_c*(DW(i-1,j,k+1)*DelxSc(i-1,j,k+1)+DU(i-1,j,k+1)/D(i,j))&
                          )&
                       ) /2.0/dx & !d/dx[nu*d/dSc[DW*DelxSc+DU/D]]
                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/Schmidt)*DelxSc(i,j,k-1)&
                         *(DW(i+1,j,k-1)-DW(i-1,j,k-1))/2.0/dx &
                       +beta_c *(Cmu(i,j,k  )+CmuHt(i,j,k  )/Schmidt)*DelxSc(i,j,k  )&
                         *(DW(i+1,j,k  )-DW(i-1,j,k  ))/2.0/dx &
                       +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/Schmidt)*DelxSc(i,j,k+1)&
                         *(DW(i+1,j,k+1)-DW(i-1,j,k+1))/2.0/dx &
                        )&  !d/dSc[nu*DelxSc*d/dx[DW]]
                     + ( ( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*alpha_c/dsig(k)&
                           +(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck1/(dsigck+dsigck1)&
                          )*DU(i,j,k-1)/D(i,j) &
                         +( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*beta_c /dsig(k)&
                           -(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck1/dsigck&
                          )*DU(i,j,k  )/D(i,j) &
                         +( (nuH_top*DelxSl(i,j,k+1)-nuH_bot*DelxSl(i,j,k))*gamma_c/dsig(k)&
                           +(nuH_top*DelxSl(i,j,k+1)+nuH_bot*DelxSl(i,j,k))/dsigck /(dsigck+dsigck1)&
                          )*DU(i,j,k+1)/D(i,j) &
                        ) !d\dSc[nu*DelxSc*d\dSc[DU/D]]   

       !Diffzy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
       !        0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
       !        (D(i,j+1)+D(i,j))*(W(i,j+1,k)-W(i,j,k))-  & 
       !        (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
       !        CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
       !        (W(i,j,k)-W(i,j-1,k)))/dy**2+  &
       !        ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*  &
       !        (V(i,j+1,k+1)-V(i,j+1,k-1))/(sigc(k+1)-sigc(k-1))-  &
       !        (Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt)*  & 
       !        (V(i,j-1,k+1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dy)
       Diffzy(i,j,k) = ( (0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt) &
                         *(DW(i,j+1,k)-DW(i,j,k))/dy &
                        -(0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/Schmidt) &
                         *(DW(i,j,k)-DW(i,j-1,k))/dy &
                       ) /dy & !d/dy[nu*d/dy[DW]]
                     + (  (Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt) &
                         *(alpha_c*(DW(i,j+1,k-1)*DelySc(i,j+1,k-1)+DV(i,j+1,k-1)/D(i,j))&
                          +beta_c *(DW(i,j+1,k  )*DelySc(i,j+1,k  )+DV(i,j+1,k  )/D(i,j))&
                          +gamma_c*(DW(i,j+1,k+1)*DelySc(i,j+1,k+1)+DV(i,j+1,k+1)/D(i,j))&
                          )&
                        - (Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt) &
                         *(alpha_c*(DW(i,j-1,k-1)*DelySc(i,j-1,k-1)+DV(i,j-1,k-1)/D(i,j))&
                          +beta_c *(DW(i,j-1,k  )*DelySc(i,j-1,k  )+DV(i,j-1,k  )/D(i,j))&
                          +gamma_c*(DW(i,j-1,k+1)*DelySc(i,j-1,k+1)+DV(i,j-1,k+1)/D(i,j))&
                          )&
                       ) /2.0/dy & !d/dy[nu*d/dSc[DW*DelySc+DV/D]]
                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/Schmidt)*DelySc(i,j,k-1)&
                         *(DW(i,j+1,k-1)-DW(i,j-1,k-1))/2.0/dy &
                       +beta_c *(Cmu(i,j,k  )+CmuHt(i,j,k  )/Schmidt)*DelySc(i,j,k  )&
                         *(DW(i,j+1,k  )-DW(i,j-1,k  ))/2.0/dy &
                       +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/Schmidt)*DelySc(i,j,k+1)&
                         *(DW(i,j+1,k+1)-DW(i,j-1,k+1))/2.0/dy &
                        )&  !d/dSc[nu*DelySc*d/dy[DW]]
                     + ( ( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*alpha_c/dsig(k)&
                           +(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck1/(dsigck+dsigck1)&
                          )*DV(i,j,k-1)/D(i,j) &
                         +( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*beta_c /dsig(k)&
                           -(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck1/dsigck&
                          )*DV(i,j,k  )/D(i,j) &
                         +( (nuH_top*DelySl(i,j,k+1)-nuH_bot*DelySl(i,j,k))*gamma_c/dsig(k)&
                           +(nuH_top*DelySl(i,j,k+1)+nuH_bot*DelySl(i,j,k))/dsigck /(dsigck+dsigck1)&
                          )*DV(i,j,k+1)/D(i,j) &
                        ) !d/dSc[nu*DelySc*d/dSc[DV/D]]
     enddo
     enddo
     enddo
     !!just 2d later should be deleted
     !Diffxy = zero
     !Diffyx = zero; Diffyy = zero; Diffyz = zero
     !Diffzy = zero 
    end subroutine diffusion

     subroutine diffusion_old
!---------------------------------------------------------  
!    This subroutine is used to evaluate diffusion terms
!    in simplified form following FVCOM
!    Called by 
!      eval_duvw 
!    Last update: 20/12/2011, Gangfeng Ma 
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
   
     Diffxx = zero; Diffxy = zero
     Diffyx = zero; Diffyy = zero
     Diffzx = zero; Diffzy = zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       Diffxx(i,j,k) = ((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)* &
               (D(i+1,j)+D(i,j))*(U(i+1,j,k)-U(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))* &
               (U(i,j,k)-U(i-1,j,k)))/dx**2
       Diffxy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(U(i,j+1,k)-U(i,j,k))- &
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (U(i,j,k)-U(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*D(i,j+1)*  &
               (V(i+1,j+1,k)-V(i-1,j+1,k))-(Cmu(i,j-1,k)+  &
               CmuHt(i,j-1,k)/Schmidt)*D(i,j-1)*  &
               (V(i+1,j-1,k)-V(i-1,j-1,k)))/(4.0*dx*dy)

       Diffyx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i+1,j)+D(i,j))*(V(i+1,j,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
               (V(i,j,k)-V(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*D(i+1,j)*  &
               (U(i+1,j+1,k)-U(i+1,j-1,k))-(Cmu(i-1,j,k)+  &
               CmuHt(i-1,j,k)/Schmidt)*D(i-1,j)*  &
               (U(i-1,j+1,k)-U(i-1,j-1,k)))/(4.0*dx*dy)
       Diffyy(i,j,k) = ((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(V(i,j+1,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (V(i,j,k)-V(i,j-1,k)))/dy**2

       Diffzx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i+1,j)+D(i,j))*(W(i+1,j,k)-W(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
               (W(i,j,k)-W(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*  &
               (U(i+1,j,k+1)-U(i+1,j,k-1))/(sigc(k+1)-sigc(k-1))-  &
               (Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt)*  &
               (U(i-1,j,k+1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dx)
       Diffzy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(W(i,j+1,k)-W(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (W(i,j,k)-W(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*  &
               (V(i,j+1,k+1)-V(i,j+1,k-1))/(sigc(k+1)-sigc(k-1))-  &
               (Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt)*  & 
               (V(i,j-1,k+1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dy)
     enddo
     enddo
     enddo

     end subroutine diffusion_old


# if defined (PARALLEL)
    subroutine phi_2D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    REAL(SP),DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_2D_exch
# endif


# if defined (PARALLEL)
    SUBROUTINE phi_3D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc,Kloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,k,ik,jk,nreq,len
    REAL(SP),DIMENSION(Mloc*Kloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc*Kloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Kloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sWmsg(jk,i) = PHI(Ibeg+i-1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sEmsg(jk,i) = PHI(Iend-i+1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Ibeg-i,j,k) = rWmsg(jk,i)
       enddo
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Iend+i,j,k) = rEmsg(jk,i)
       enddo
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Kloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sSmsg(ik,j) = PHI(i,Jbeg+j-1,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sNmsg(ik,j) = PHI(i,Jend-j+1,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jbeg-j,k) = rSmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jend+j,k) = rNmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_3D_exch
# endif


# if defined(PARALLEL)
    ! Jeff added this subroutine to pass mask 02/14/2011
    SUBROUTINE phi_int_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    INTEGER,DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    INTEGER,DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_INTEGER, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_INTEGER, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_INTEGER, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_INTEGER, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_INTEGER, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_INTEGER, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_INTEGER, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_INTEGER, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif
    END SUBROUTINE phi_int_exch
# endif


    subroutine adv_scalar_hlpa(Flx,Fly,Flz,Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   IVAR: indication of different scalars 
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity 
!    = 4: temperature
!    = 5: bubble number density 
!    = 6: sediment concentration 
!   Last update: Gangfeng Ma, 04/04/2012
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc,Nloc,Kloc),  intent(in) :: Phi
    real(SP), dimension(Mloc1,Nloc,Kloc),  intent(in) :: Flx
    real(SP), dimension(Mloc,Nloc1,Kloc),  intent(in) :: Fly
    real(SP), dimension(Mloc,Nloc,Kloc1),  intent(in) :: Flz
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy
    real(SP) :: DUfs,DVfs,Wfs,Fww,Fw,Fp,Fe,hlpa,SchtH
    real(SP) :: L1top,L2top,L1bot,L2bot,alpha_c,beta_c,gamma_c, &
                dsigck,dsigck1,nuH_top,nuH_bot
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUfs = Flx(i,j,k)
      Fww = Phi(i-2,j,k)
      Fw  = Phi(i-1,j,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i+1,j,k)
      Scalx(i,j,k) = DUfs*hlpa(DUfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend      
      DVfs = Fly(i,j,k)
      Fww = Phi(i,j-2,k)
      Fw  = Phi(i,j-1,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j+1,k)
      Scaly(i,j,k) = DVfs*hlpa(DVfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Wfs = Flz(i,j,k)
      Fww = Phi(i,j,k-2)
      Fw  = Phi(i,j,k-1)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j,k+1)
      Scalz(i,j,k) = Wfs*hlpa(Wfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(IVAR,Scalx,Scaly,Scalz)

    ! Schmidt number (fixed by Cheng for RNG)
    if(IVAR==1) then  ! tke eq.
      if (RNG) then
        SchtH = 0.72
      else
        SchtH = 1.0
      endif
    elseif(IVAR==2) then  ! epsilon eq.
      if (RNG) then
         SchtH = 0.72
      else
         SchtH = 1.3
      endif
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
    elseif(IVAR==6) then ! sediment
      SchtH = 1.0
    else
      SchtH = 1.0
    endif

    ! estimate horizontal diffusion
    Sdiffx = Zero; Sdiffy = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      Sdiffx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH)*  &
           (D(i,j)+D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/dx**2
      Sdiffy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH)*  &
           (D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/dy**2
    enddo
    enddo
    enddo

!    Sdiffx = Zero
!    Sdiffy = Zero
!    do k = Kbeg,Kend
!       !these are used for vertical gradient in non-uniform grid
!       !see Derakhti etal 2015a Appendix B
!       dsigck  = (dsig(k)+dsig(k+1))/2.0 
!       dsigck1 = (dsig(k-1)+dsig(k))/2.0 
!       alpha_c = -dsigck/(dsigck+dsigck1)/dsigck1
!       beta_c = (dsigck-dsigck1)/(dsigck*dsigck1)
!       gamma_c = dsigck1/(dsigck+dsigck1)/dsigck
!    do j = Jbeg,Jend
!    do i = Ibeg,Iend
!      if(Mask(i,j)==0) cycle
!       Sdiffx(i,j,k) = (  (0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)! &
!                         *(Phi(i+1,j,k)*D(i+1,j)-Phi(i,j,k)*D(i,j))/dx &
!                        -(0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH) !&
!                         *(Phi(i,j,k)*D(i,j)-Phi(i-1,j,k)*D(i-1,j))/dx &
!                       ) /dx & 
!                     + ( (Cmu(i+1,j,k)+CmuHt(i+1,j,k)/SchtH) &
!                         *(alpha_c*D(i+1,j)*Phi(i+1,j,k-1)*DelxSc(i+1,j,k-1)&
!                          +beta_c *D(i+1,j)*Phi(i+1,j,k  )*DelxSc(i+1,j,k  )&
!                          +gamma_c*D(i+1,j)*Phi(i+1,j,k+1)*DelxSc(i+1,j,k+1)) & 
!                        -(Cmu(i-1,j,k)+CmuHt(i-1,j,k)/SchtH) &
!                         *(alpha_c*D(i-1,j)*Phi(i-1,j,k-1)*DelxSc(i-1,j,k-1)&
!                          +beta_c *D(i-1,j)*Phi(i-1,j,k  )*DelxSc(i-1,j,k  )&
!                          +gamma_c*D(i-1,j)*Phi(i-1,j,k+1)*DelxSc(i-1,j,k+1)) & 
!                       ) /2.0/dx & 
!                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/SchtH)*DelxSc(i,j,k-1)&
!                         *(D(i+1,j)*Phi(i+1,j,k-1)-D(i-1,j)*Phi(i-1,j,k-1))/2.0/dx &
!                        +beta_c*(Cmu(i,j,k)+CmuHt(i,j,k)/SchtH)*DelxSc(i,j,k)&
!                         *(D(i+1,j)*Phi(i+1,j,k  )-D(i-1,j)*Phi(i-1,j,k  ))/2.0/dx &
!                        +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/SchtH)*DelxSc(i,j,k+1)&
!                         *(D(i+1,j)*Phi(i+1,j,k+1)-D(i-1,j)*Phi(i-1,j,k+1))/2.0/dx &
!                        ) 
!
!      Sdiffy(i,j,k) = (  (0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH) !&
!                         *(Phi(i,j+1,k)*D(i,j+1)-Phi(i,j,k)*D(i,j))/dy &
!                        -(0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH) !&
!                         *(Phi(i,j,k)*D(i,j)-Phi(i,j-1,k)*D(i,j-1))/dy &
!                       ) /dy & 
!                     + ( (Cmu(i,j+1,k)+CmuHt(i,j+1,k)/SchtH) &
!                         *(alpha_c*D(i,j+1)*Phi(i,j+1,k-1)*DelySc(i,j+1,k-1)&
!                          +beta_c *D(i,j+1)*Phi(i,j+1,k  )*DelySc(i,j+1,k  )&
!                          +gamma_c*D(i,j+1)*Phi(i,j+1,k+1)*DelySc(i,j+1,k+1)) & 
!                        -(Cmu(i,j-1,k)+CmuHt(i,j-1,k)/SchtH) &
!                         *(alpha_c*D(i,j-1)*Phi(i,j-1,k-1)*DelySc(i,j-1,k-1)&
!                          +beta_c *D(i,j-1)*Phi(i,j-1,k  )*DelySc(i,j-1,k  )&
!                          +gamma_c*D(i,j-1)*Phi(i,j-1,k+1)*DelySc(i,j-1,k+1)) & 
!                       ) /2.0/dy & 
!                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/SchtH)*DelySc(i,j,k-1)&
!                         *(D(i,j+1)*Phi(i,j+1,k-1)-D(i,j-1)*Phi(i,j-1,k-1))/2.0/dy &
!                        +beta_c*(Cmu(i,j,k)+CmuHt(i,j,k)/SchtH)*DelySc(i,j,k)&
!                         *(D(i,j+1)*Phi(i,j+1,k  )-D(i,j-1)*Phi(i,j-1,k  ))/2.0/dy &
!                        +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/SchtH)*DelySc(i,j,k+1)&
!                         *(D(i,j+1)*Phi(i,j+1,k+1)-D(i,j-1)*Phi(i,j-1,k+1))/2.0/dy &
!                        )  
!    enddo
!    enddo
!    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))  &
                  -1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k)) &
                  +Sdiffx(i,j,k)+Sdiffy(i,j,k)
    enddo
    enddo
    enddo

    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)

    end subroutine adv_scalar_hlpa


    function hlpa(Uw,Fww,Fw,Fp,Fe)
!-------------------------------------------------------
!   HLPA scheme
!-------------------------------------------------------
    use global, only: SP,Zero
    implicit none
    real(SP), intent(in)  :: Uw,Fww,Fw,Fp,Fe
    real(SP) :: hlpa,Alpha_pl,Alpha_mn

    if(Uw>=Zero) then
      if(abs(Fp-2.*Fw+Fww)<abs(Fp-Fww)) then
        Alpha_pl = 1.0
      else
        Alpha_pl = 0.0
      endif

      if(abs(Fp-Fww)<=1.e-16) then
        hlpa = Fw
      else
        hlpa = Fw+Alpha_pl*(Fp-Fw)*(Fw-Fww)/(Fp-Fww)
      endif
    endif

    if(Uw<Zero) then
      if(abs(Fw-2.*Fp+Fe)<abs(Fw-Fe)) then
        Alpha_mn = 1.0
      else
        Alpha_mn = 0.0
      endif

      if(abs(Fw-Fe)<=1.e-16) then
        hlpa = Fp
      else
        hlpa = Fp+Alpha_mn*(Fw-Fp)*(Fp-Fe)/(Fw-Fe)
      endif
    endif

    return
    end function hlpa


    subroutine flux_scalar_bc(IVAR,Scalx,Scaly,Scalz)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
!   Last update: Gangfeng Ma, 09/02/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc1,Nloc,Kloc), intent(inout) :: Scalx
    real(SP), dimension(Mloc,Nloc1,Kloc), intent(inout) :: Scaly
    real(SP), dimension(Mloc,Nloc,Kloc1), intent(inout) :: Scalz
    real(SP), dimension(Nloc,Kloc) :: Scal_X0,Scal_Xn
    integer :: i,j,k

    ! temporarily set it here
	! added by Cheng for initialization
	Scal_X0 = Zero
	Scal_Xn = Zero
# if defined (SEDIMENT)
    Scal_X0 = Sed_X0
    Scal_Xn = Sed_Xn
# endif
# if defined (SALINITY)
    Scal_X0 = Sin_X0
    Scal_Xn = Sin_Xn
# endif

    ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by Cheng for nesting. Please search for others with (COUPLING) in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Scalx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Scalx(Ibeg,j,k) = Ex(Ibeg,j,k)*Scal_X0(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
!       if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Scalx(Iend1,j,k) = Zero
!       elseif(Bc_Xn==3) then
!         Scalx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Scal_Xn(j,k)
!       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif      

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
         Scaly(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
         Scaly(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) then
        Scalx(i,j,k) = Zero
        Scalx(i+1,j,k) = Zero
        Scaly(i,j,k) = Zero
        Scaly(i,j+1,k) = Zero
      endif
    enddo
    enddo
    enddo

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Scalz(i,j,Kbeg) = Zero
      Scalz(i,j,Kend1) = Zero
    enddo
    enddo

    return
    end subroutine flux_scalar_bc

   
    subroutine les_dyn(ISTEP)
!----------------------------------------------------------
!   large eddy simulation (LES) with dynamic subgrid model
!   Last update: Gangfeng Ma, 10/29/2015
!----------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), parameter :: Dmin = 0.02
    integer :: i,j,k,n
    real(SP) :: weit,fact,alpha2
    real(SP), dimension(:,:,:,:), allocatable :: Sij,SabSij,Filtered_Sij,Filtered_SabSij, &
                                                 Mij,Uij,Filtered_Uij,Ui,Filtered_Ui,Lij
    real(SP), dimension(:,:,:), allocatable :: Sab,Filtered_Sab,Smag_Const,Filter_Width,Filtered_Const, &
                                               MijMij,LijMij

    allocate(Sij(Mloc,Nloc,Kloc,6))
    allocate(SabSij(Mloc,Nloc,Kloc,6))
    allocate(Filtered_Sij(Mloc,Nloc,Kloc,6))
    allocate(Filtered_SabSij(Mloc,Nloc,Kloc,6))
    allocate(Mij(Mloc,Nloc,Kloc,6))
    allocate(Uij(Mloc,Nloc,Kloc,6))
    allocate(Filtered_Uij(Mloc,Nloc,Kloc,6))
    allocate(Ui(Mloc,Nloc,Kloc,3))
    allocate(Filtered_Ui(Mloc,Nloc,Kloc,3))
    allocate(Lij(Mloc,Nloc,Kloc,6))
    allocate(MijMij(Mloc,Nloc,Kloc))
    allocate(LijMij(Mloc,Nloc,Kloc))
    allocate(Sab(Mloc,Nloc,Kloc))
    allocate(Filtered_Sab(Mloc,Nloc,Kloc))
    allocate(Smag_Const(Mloc,Nloc,Kloc))
    allocate(Filter_Width(Mloc,Nloc,Kloc))
    allocate(Filtered_Const(Mloc,Nloc,Kloc))
    

!   S:
!   S11 (1) S22 (2) S33 (3)
!   S12 (4) S23 (5) S31 (6)
 
    Sij = 0.0
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask(i,j)==1.and.D(i,j)>Dmin) then
        Sij(i,j,k,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                       (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)       
        Sij(i,j,k,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                       (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)       
        Sij(i,j,k,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
        Sij(i,j,k,4) = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                       (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+ &                             
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-  &
                       V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)) 
        Sij(i,j,k,5) = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                      
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
        Sij(i,j,k,6) = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                      
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-  &
               W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom 
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = 1,Nghost
      do n = 1,6
        Sij(i,j,Kbeg-k,n) = Sij(i,j,Kbeg+k-1,n)
      enddo 
    enddo
    enddo
    enddo

    ! at the free surface                     
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = 1,Nghost
      do n = 1,6
        Sij(i,j,Kend+k,n) = Sij(i,j,Kend-k+1,n)
      enddo
    enddo
    enddo
    enddo

    do n = 1,6
      call phi_3D_exch(Sij(:,:,:,n))
    enddo

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
    do i = 1,Nghost
      do n = 1,6
        Sij(Ibeg-i,j,k,n) = Sij(Ibeg+i-1,j,k,n)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
    do i = 1,Nghost
      do n = 1,6
        Sij(Iend+i,j,k,n) = Sij(Iend-i+1,j,k,n)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
     do j = 1,Nghost
       do n = 1,6
         Sij(i,Jbeg-j,k,n) = Sij(i,Jbeg+j-1,k,n)
       enddo
     enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
     do j = 1,Nghost
       do n = 1,6
         Sij(i,Jend+j,k,6) = Sij(i,Jend-j+1,k,6)
       enddo
     enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif


!   |S|
    Sab = 0.0
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(D(i,j)>Dmin) then
        Sab(i,j,k) = Sij(i,j,k,1)**2+Sij(i,j,k,2)**2+Sij(i,j,k,3)**2  &
               +2.0*(Sij(i,j,k,4)**2+Sij(i,j,k,5)**2+Sij(i,j,k,6)**2)
        Sab(i,j,k) = dsqrt(2.0*Sab(i,j,k))
      endif
    enddo
    enddo
    enddo

!   |S|Sij
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(D(i,j)>Dmin) then
        do n = 1,6
          SabSij(i,j,k,n) = Sab(i,j,k)*Sij(i,j,k,n)
        enddo
      endif
    enddo
    enddo
    enddo

!   <|S|Sij> and <Sij>
    weit = 2.0; fact = 1.0/64.0
    do n = 1,6
      call filter(SabSij(:,:,:,n),Filtered_SabSij(:,:,:,n),weit,fact)
      call filter(Sij(:,:,:,n),Filtered_Sij(:,:,:,n),weit,fact)
    enddo

!   <|S|>
    call filter(Sab,Filtered_Sab,weit,fact)

!   Mij = alpha2*<|S|>*<Sij>-<|S|Sij>
    alpha2 = 4.0
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
    do n = 1,6
      Mij(i,j,k,n) = alpha2*Filtered_Sab(i,j,k)*Filtered_Sij(i,j,k,n)-  &
                Filtered_SabSij(i,j,k,n)
    enddo
    enddo
    enddo
    enddo

!   Uij and <Uij>
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Uij(i,j,k,1) = U(i,j,k)*U(i,j,k)
      Uij(i,j,k,2) = V(i,j,k)*V(i,j,k)
      Uij(i,j,k,3) = W(i,j,k)*W(i,j,k)
      Uij(i,j,k,4) = U(i,j,k)*V(i,j,k)
      Uij(i,j,k,5) = V(i,j,k)*W(i,j,k)
      Uij(i,j,k,6) = U(i,j,k)*W(i,j,k)
    enddo
    enddo
    enddo

    do n = 1,6
      call filter(Uij(:,:,:,n),Filtered_Uij(:,:,:,n),weit,fact)
    enddo

!   Ui and <Ui>
    Ui = 0.0
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Ui(i,j,k,1) = U(i,j,k)
      Ui(i,j,k,2) = V(i,j,k)
      Ui(i,j,k,3) = W(i,j,k)
    enddo
    enddo
    enddo

    do n = 1,3
      call filter(Ui(:,:,:,n),Filtered_Ui(:,:,:,n),weit,fact)
    enddo

!   Lij = <Uij>-<Ui><Uj>
    Lij = 0.0
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Lij(i,j,k,1) = Filtered_Uij(i,j,k,1)-Filtered_Ui(i,j,k,1)**2
      Lij(i,j,k,2) = Filtered_Uij(i,j,k,2)-Filtered_Ui(i,j,k,2)**2
      Lij(i,j,k,3) = Filtered_Uij(i,j,k,3)-Filtered_Ui(i,j,k,3)**2
      Lij(i,j,k,4) = Filtered_Uij(i,j,k,4)-Filtered_Ui(i,j,k,1)*Filtered_Ui(i,j,k,2)
      Lij(i,j,k,5) = Filtered_Uij(i,j,k,5)-Filtered_Ui(i,j,k,2)*Filtered_Ui(i,j,k,3)
      Lij(i,j,k,6) = Filtered_Uij(i,j,k,6)-Filtered_Ui(i,j,k,3)*Filtered_Ui(i,j,k,1)
    enddo
    enddo
    enddo

!   MijMij and LijMij                                                                                                     
    MijMij = 0.0; LijMij = 0.0
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(D(i,j)>Dmin) then
        MijMij(i,j,k) = Mij(i,j,k,1)**2+Mij(i,j,k,2)**2+Mij(i,j,k,3)**2   &                                               
              +2.0*(Mij(i,j,k,4)**2+Mij(i,j,k,5)**2+Mij(i,j,k,6)**2)
        LijMij(i,j,k) = Lij(i,j,k,1)*Mij(i,j,k,1)+Lij(i,j,k,2)*Mij(i,j,k,2)  &                                            
              +Lij(i,j,k,3)*Mij(i,j,k,3)+2.0*(Lij(i,j,k,4)*Mij(i,j,k,4)  &                                                
              +Lij(i,j,k,5)*Mij(i,j,k,5)+Lij(i,j,k,6)*Mij(i,j,k,6))
      endif
    enddo
    enddo
    enddo

!   Smagorinsky constant
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Smag_Const(i,j,k) = 0.5*(-LijMij(i,j,k))/(MijMij(i,j,k)+1.e-16)
      Filter_Width(i,j,k) = (dx*dy*D(i,j)*dsig(k))**(1./3.)
      Smag_Const(i,j,k) = Smag_Const(i,j,k)/Filter_Width(i,j,k)**2
      if(Smag_Const(i,j,k)<0.0) Smag_Const(i,j,k) = 0.0
      if(Smag_Const(i,j,k)>0.04) Smag_Const(i,j,k) = 0.04
    enddo
    enddo
    enddo

    weit = 1.0; fact = 1.0/27.0
    call filter(Smag_Const,Filtered_Const,weit,fact)

    ! eddy viscosity
    CmuVt = 0.0
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(D(i,j)>Dmin) then
        CmuVt(i,j,k) = Filtered_Const(i,j,k)*Filter_Width(i,j,k)**2*Sab(i,j,k)
        CmuVt(i,j,k) = max(CmuVt(i,j,k), 0.0)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(CmuVt)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         ! no turbulence at wave generation region
         CmuVt(Ibeg,j,k) = Zero
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       else       
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         CmuVt(Iend+i,j,k) = CmuVt(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jbeg-j,k) = CmuVt(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jend+j,k) = CmuVt(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif       

     ! no turbulence in the internal wavemaker region   
    if(WaveMaker(1:3)=='INT'.or.WaveMaker(1:3)=='LEF' &
        .or.WaveMaker(1:3)=='FOC') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          CmuVt(i,j,k) = Zero
        endif
      enddo
      enddo
      enddo
    endif

    deallocate(Sij)
    deallocate(SabSij)
    deallocate(Filtered_Sij)
    deallocate(Filtered_SabSij)
    deallocate(Mij)
    deallocate(Uij)
    deallocate(Filtered_Uij)
    deallocate(Ui)
    deallocate(Filtered_Ui)
    deallocate(Lij)
    deallocate(MijMij)
    deallocate(LijMij)
    deallocate(Sab)
    deallocate(Filtered_Sab)
    deallocate(Smag_Const)
    deallocate(Filter_Width)
    deallocate(Filtered_Const)

    return
    end subroutine les_dyn

    subroutine filter(PhiIn,PhiOut,weit,fact)
!-------------------------------------------------------
!   double filtering for dynamic subgrid model
!   by Gangfeng Ma, 30/10/2015
!------------------------------------------------------
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc,Kloc), intent(in)  :: PhiIn
    real(SP), dimension(Mloc,Nloc,Kloc), intent(out) :: PhiOut
    real(SP), dimension(Mloc,Nloc,Kloc) :: PhiS1,PhiS2
    real(SP), intent(in) :: weit,fact
    integer :: i,j,k

    PhiOut = 0.0
   
    ! filtering in x direction
    do k = 1,Kloc
    do j = 1,Nloc
    do i = Ibeg,Iend
      PhiS1(i,j,k) = PhiIn(i-1,j,k)+weit*PhiIn(i,j,k)+PhiIn(i+1,j,k)
    enddo
    enddo
    enddo

    ! filtering in y direction
    do k = 1,Kloc
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      PhiS2(i,j,k) = PhiS1(i,j-1,k)+weit*PhiS1(i,j,k)+PhiS1(i,j+1,k)
    enddo
    enddo
    enddo

    ! filtering in z direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      PhiOut(i,j,k) = PhiS2(i,j,k-1)+weit*PhiS2(i,j,k)+PhiS2(i,j,k+1)
    enddo
    enddo
    enddo

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      PhiOut(i,j,k) = fact*PhiOut(i,j,k)
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        PhiOut(i,j,Kbeg-k) = PhiOut(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        PhiOut(i,j,Kend+k) = PhiOut(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(PhiOut)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         PhiOut(Ibeg-i,j,k) = PhiOut(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         PhiOut(Iend+i,j,k) = PhiOut(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         PhiOut(i,Jbeg-j,k) = PhiOut(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         PhiOut(i,Jend+j,k) = PhiOut(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

    return
    end subroutine filter


    subroutine les_3D(ISTEP)
!------------------------------------------------------
!   large eddy simulation (LES)
!   Last update: Gangfeng Ma, 09/22/2011
!------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), parameter :: Dmin = 0.02
    integer :: i,j,k,Iter
    real(SP) :: S11,S22,S33,S12,S13,S23,SijSij,Filter
    real(SP) :: Umag,Zdis,X0,Xa,Xn,FricU
 
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask9(i,j)==1.and.D(i,j)>Dmin) then
        S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                                      
        S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
              (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                                      
        S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
        S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
              (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &                             
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
              (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)) 
        S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
        S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
        SijSij = S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2)
        Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
        CmuVt(i,j,k) = (Cvs*Filter)**2*sqrt(2.0*SijSij)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        CmuVt(i,j,Kbeg) = Cmu(i,j,Kbeg)
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)
        if(Zob>=0.1) then
          ! rough wall
          FricU = Umag/(1.0/0.41*log(30.0*Zdis/Zob))
        else
          ! smooth wall
          X0 = 0.05
          Iter = 0
       
          Xa = dlog(9.0*Umag*Zdis/Cmu(i,j,Kbeg))
 10       Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
          if(Iter>=20) then
            write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
          endif
          if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
            X0 = Xn
            Iter = Iter+1
            goto 10
          else
            FricU = Xn*Umag
          endif
        endif

        CmuVt(i,j,Kbeg) = 0.41*Zdis*FricU*  &
           (1.0-exp(-Zdis*FricU/Cmu(i,j,Kbeg)/19.0))**2
      endif
 100  continue

      do k = 1,Nghost
        CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(CmuVt)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         ! no turbulence at wave generation region
         CmuVt(Ibeg,j,k) = Zero
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       else       
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         CmuVt(Iend+i,j,k) = CmuVt(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jbeg-j,k) = CmuVt(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jend+j,k) = CmuVt(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif       

     ! no turbulence in the internal wavemaker region   
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          CmuVt(i,j,k) = Zero
        endif
      enddo
      enddo
      enddo
    endif

    ! estimate turbulent dissipation rate (Van den Hengel et al., 2005)
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
      Eps(i,j,k) = 2.0*CmuVt(i,j,k)**3/(Cvs*Filter)**4
    enddo
    enddo
    enddo

    end subroutine les_3D

# if defined (VEGETATION)
    subroutine kepsilon_veg(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.01
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress,Vorticity
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk,c5e,clambda,beta_p,beta_d,ced
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb
    real(SP) :: S11,S22,S33,S12,S13,S23,xlfs,Tkes,Epss
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Vorticity(3,3))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    Sche = 1.3
    Schk = 1.0
    beta_p = 0.2
    beta_d = 1.0
    c5e = 0.0
    clambda = 0.01
    ced = beta_p*1.21**(-1.5)

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                 (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
               (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        ! Notice: if using nonlinear model, the initial seeding of tke and epsilon
        !         cannot be zero in order to generate turbulence production. Check 
        !         tke_min and eps_min in subroutine initial. 
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k) 
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (1./3.)*(1./(3.7+smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 2./3./(123.5+2.0*dmax**2)
          c_2 = -2./3./(39.2+2.0*dmax**2)
          c_3 = 2./3./(246.9+2.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/Eps_Old(i,j,k)*  &
                    (VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                 c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   (VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                     VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm 
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          do n = 1,3
          do m = 1,3
            Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
          enddo
          enddo

          !! no negative production at the surface
          !if(k==Kend.and.Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

          ! Do not allow negative production
          if(Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

        endif
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
      endif
    enddo
    enddo
    enddo

    ! vegetation-induced (wake) turbulence
    Wke_p = 0.0
    Prod_w = 0.0
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                              
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                        
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Wke_p(i,j,k) = beta_d*0.5*StemD*VegDens*VegDrag*Umag*Tke_Old(i,j,k)
            Prod_w(i,j,k) = beta_p*0.5*StemD*VegDens*VegDrag*Umag**3
          else
            Wke_p(i,j,k) = 0.0
            Prod_w(i,j,k) = 0.0
          endif
        endif
      endif
    enddo
    enddo
    enddo

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k)-  &
                c5e*Wke_p(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            

        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = dmax1(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k)-Wke_p(i,j,k))-DEps(i,j,k)      

        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = dmax1(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! solve WKE
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_w,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_w(i,j,k)+Wke_p(i,j,k)-Eps_w(i,j,k))

        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DWke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DWke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1.and.sigc(k)*D(i,j)<=VegH.and. &
         xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        DWke(i,j,k) = ALPHA(ISTEP)*DWke0(i,j,k)+BETA(ISTEP)*DWke(i,j,k)
        DWke(i,j,k) = dmax1(DWke(i,j,k),0.0)
      else
        DWke(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        Tkeb = Tke_min
        Epsb = Eps_min
  
        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)

        X0 = 0.05
        Iter = 0

        Xa = dlog(9.0*Umag*Zdis/Visc)
 10     Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
        if(Iter>=20) then
          write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
        endif
        if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
          X0 = Xn
          Iter = Iter+1
          goto 10
        else
          FricU = Xn*Umag
        endif

!        if(Ibot==1) then
!          FricU = sqrt(Cd0)*Umag
!        else
!          FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
!        endif

        Tkeb = FricU**2/sqrt(cmiu)
        Epsb = FricU**3/(Kappa*Zdis)        

        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      endif

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DWke(i,j,Kbeg-k) = DWke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DWke(i,j,Kend+k) = DWke(i,j,Kend-k+1)
      enddo
   
!      xlfs = 0.5*D(i,j)*dsig(Kend)
!      Tkes = DTke(i,j,Kend)/D(i,j)
!      Epss = cmiu**0.75*Tkes**1.5/xlfs
!      DEps(i,j,Kend) = D(i,j)*Epss
      do k = 1,Nghost
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)                                                                  
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DWke)
    call phi_3D_exch(DEps)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         do i = Ibeg,Ibeg+5
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif
     enddo
     enddo

     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DWke(Ibeg-i,j,k) = DWke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DWke(Iend+i,j,k) = DWke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DWke(i,Jbeg-j,k) = DWke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DWke(i,Jend+j,k) = DWke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

   
    ! update eddy viscosity
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        Tke_w(i,j,k) = DWke(i,j,k)/D(i,j)
        Eps_w(i,j,k) = ced*Tke_w(i,j,k)**(3./2.)/StemD

        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)+  &
             Clambda/ced*Tke_w(i,j,k)**0.5*StemD
 
!        ! If turbulence is too small, assume small eddy viscosity
!        if(Tke(i,j,k)<1.e-6.or.Eps(i,j,k)<1.e-6) then
!          CmuVt(i,j,k) = Cmut_min
!        else
!          CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
!        endif
!
!        if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.  &
!          (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn).and.sigc(k)*D(i,j)<=VegH) then                                         !           
!            CmuVt(i,j,k) = CmuVt(i,j,k)+Clambda/ced*Tke_w(i,j,k)**0.5*StemD
!        endif
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        DWke(i,j,k) = 0.0
        Tke_w(i,j,k) = 0.0
        Eps_w(i,j,k) = 0.0
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
!        VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
!                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
!        VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
!                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
!        UpWp(i,j,k) = CmuVt(i,j,k)*(VelGrad(1,2)+VelGrad(2,1))
         UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j) 
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    end subroutine kepsilon_veg
# endif


    subroutine kepsilon_3D_old(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.01
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress,Vorticity
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb
    real(SP) :: S11,S22,S33,S12,S13,S23,xlfs,Tkes,Epss
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Vorticity(3,3))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    Sche = 1.3
    Schk = 1.0

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        ! Notice: if using nonlinear model, the initial seeding of tke and epsilon
        !         cannot be zero in order to generate turbulence production. Check 
        !         tke_min and eps_min in subroutine initial. 
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first
          VelGrad = Zero

!          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
!                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
!          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
!                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
!          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
!          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
!                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
!          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
!                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k) 
!          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
!          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
!                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
!          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
!                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
!          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          
          VelGrad(1,1) = DelxU(i,j,k) +   &
                      (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = DelyU(i,j,k) +   &
                      (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = DelxV(i,j,k) +   &
                      (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = DelyV(i,j,k) +   &
                      (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = DelxW(i,j,k) +   &
                      (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = DelyW(i,j,k) +   &
                      (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (1./3.)*(1./(3.7+smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 2./3./(123.5+2.0*dmax**2)
          c_2 = -2./3./(39.2+2.0*dmax**2)
          c_3 = 2./3./(246.9+2.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/Eps_Old(i,j,k)*  &
                   (VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   (VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          do n = 1,3
          do m = 1,3
            Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
          enddo
          enddo

          !! no negative production at the surface
          !if(k==Kend.and.Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

          ! Do not allow negative production
          if(Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

        endif
      elseif(ke_model==3) then
        ! Following Mayer and Madsen (2000), instead of determining the production on the
        ! basis of the strain rate, the production is based on the rotation of the velocity field.
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first                              
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! vorticity field and production
          do n = 1,3
          do m = 1,3
            Vorticity(n,m) = VelGrad(n,m)-VelGrad(m,n)
            Prod_s(i,j,k) = Prod_s(i,j,k)+CmuVt(i,j,k)*Vorticity(n,m)*Vorticity(n,m)
          enddo
          enddo
        endif
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
      endif
    enddo
    enddo
    enddo
 
    ! flux Richardson number
    Richf = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Prod_s(i,j,k)>Zero) then
        Richf(i,j,k) = -Prod_b(i,j,k)/(Prod_s(i,j,k)+1.0e-16)
        Richf(i,j,k) = dmax1(0.0,dmin1(0.21,Richf(i,j,k)))
      endif
    enddo
    enddo
    enddo

# if defined (VEGETATION)
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                 
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*StemD*VegDens*VegDrag*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        else  ! flexible vegetation                                                                          
          if(sigc(k)*D(i,j)<=FVegH(i,j)) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*StemD*VegDens*VegDrag*foliage(i,j)*  &                                       
                    VegH/FVegH(i,j)*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        endif
      else
        Prod_v(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k)+  &
                cfe*Prod_v(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            
# else
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+  &
                c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
# endif

# if defined (POROUSMEDIA)
        if(Porosity(i,j,k)<1.0) then
          Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
          Tke_Inf = 3.7*(1.0-Porosity(i,j,k))*Porosity(i,j,k)**1.5*Umag**2
          Eps_Inf = 39.0*(1.0-Porosity(i,j,k))**2.5*Porosity(i,j,k)**2*Umag**3/D50_por
          R5(i,j,k) = R5(i,j,k)+D(i,j)*c2e*Eps_Inf**2/(Tke_Inf+1.e-16)
       endif
# endif

        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = dmax1(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+  &
                    Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)  
# else
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
# endif

# if defined (POROUSMEDIA)
        if(Porosity(i,j,k)<1.0) then
          Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
          Tke_Inf = 3.7*(1.0-Porosity(i,j,k))*Porosity(i,j,k)**1.5*Umag**2                                                    
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+  &
                    Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)  
          R5(i,j,k) = R5(i,j,k)+D(i,j)*Eps_Inf
        endif
# endif

        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = dmax1(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        Tkeb = Tke_min
        Epsb = Eps_min
  
        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)

        X0 = 0.05
        Iter = 0

        Xa = dlog(9.0*Umag*Zdis/Visc)
 10     Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
        if(Iter>=20) then
          write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
        endif
        if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
          X0 = Xn
          Iter = Iter+1
          goto 10
        else
          FricU = Xn*Umag
        endif

!        if(Ibot==1) then
!          FricU = sqrt(Cd0)*Umag
!        else
!          FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
!        endif

        Tkeb = FricU**2/sqrt(cmiu)
        Epsb = FricU**3/(Kappa*Zdis)        

        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      endif

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
   
!      xlfs = 0.5*D(i,j)*dsig(Kend)
!      Tkes = DTke(i,j,Kend)/D(i,j)
!      Epss = cmiu**0.75*Tkes**1.5/xlfs
!      DEps(i,j,Kend) = D(i,j)*Epss
!      do k = 1,Nghost
!        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)                                                                  
!      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='FOC') then
         do i = Ibeg,Ibeg+10
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif
     enddo
     enddo

     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

   
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)

        ! If turbulence is too small, assume small eddy viscosity
        if(Tke(i,j,k)<1.e-6.or.Eps(i,j,k)<1.e-6) then
          CmuVt(i,j,k) = Cmut_min
        else
          CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
        endif
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
!        VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
!                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
!        VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
!                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
!        UpWp(i,j,k) = CmuVt(i,j,k)*(VelGrad(1,2)+VelGrad(2,1))
         UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j) 
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    end subroutine kepsilon_3D_old


    subroutine kepsilon_3D(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!  
!   Modified by M.Derakhti as in derakhti etal 2015a
!   RNG-based k-e model of Yakhot etal 1992 is also added
!   last update: 02/15/2015
!
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    !real(SP), parameter :: Dmin = 0.02
    real(SP) :: Dmin
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress,Vorticity
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,cfk,cfe,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen,ii
    real(SP):: Fk,Fkdif
    !Added by M.Derakhti for gradient BC for Tke and Epsilon at the free surface
    real(SP), dimension(:,:,:), allocatable :: DelzTke,DelzEps,  &
              TkezL,TkezR,EpszL,EpszR,DelxTkezL,DelxEpszL,DelyEpszL,DelyTkezL
    real(SP) :: TkeFlux,EpsFlux
    real(SP) :: L1top,L2top,L1bot,L2bot,alpha_c,beta_c,gamma_c, &
                dsigck,dsigck1,nuH_top,nuH_bot,nuV_top,nuV_bot,DxD,DyD
    real(SP) :: DUmag,Ustar,Z0,Dz1,Cdrag,&
                Atop,Btop,Abot,Bbot,DUprime,DVprime,lambda
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Vorticity(3,3))
    !added b M.Derakhti
    allocate(TkezL(Mloc,Nloc,Kloc1))
    allocate(EpszL(Mloc,Nloc,Kloc1))
    allocate(TkezR(Mloc,Nloc,Kloc1))
    allocate(EpszR(Mloc,Nloc,Kloc1))
    allocate(DelzTke(Mloc,Nloc,Kloc))
    allocate(DelzEps(Mloc,Nloc,Kloc))
    allocate(DelxTkezL(Mloc,Nloc,Kloc1))
    allocate(DelxEpszL(Mloc,Nloc,Kloc1))
    allocate(DelyTkezL(Mloc,Nloc,Kloc1))
    allocate(DelyEpszL(Mloc,Nloc,Kloc1))
    !M.derakhti changed Dmin=0.02
    Dmin = 2.0*MinDep

    ! some parameters
    ! RNG-based coefficients are added by M.derakhti
    if (RNG) then
       c1e = 1.42
       !c2e is later calculated dynamically
!       c2e = 1.92
       c3e = 0.0
       cmiu = 0.085
       cfk = 1.0
       cfe = 1.33
       Sche = 0.72
       Schk = 0.72
    else
       c1e = 1.44
       c2e = 1.92
!     c3e = -1.4
       c3e = 0.0
       cmiu = 0.09
       cfk = 1.0
       cfe = 1.33
       Sche = 1.3
       Schk = 1.0
    endif
    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
              (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
              (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*  &
               (S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        ! Notice: if using nonlinear model, the initial seeding of tke and epsilon
        !         cannot be zero in order to generate turbulence production. Check 
        !         tke_min and eps_min in subroutine initial. 
        ! added by M.Derakhti
        ! replace central differencing by limiter based calculation
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first
          VelGrad = Zero

          !VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
          !        (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,1) = DelxU(i,j,k) + (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)

          !VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
          !        (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,2) = DelyU(i,j,k) + (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)

          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          !VelGrad(1,3) =              + DelzU(i,j,k)/D(i,j)
          
          !VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
          !        (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,1) = DelxV(i,j,k) + (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          
          !VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
          !        (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(2,2) = DelyV(i,j,k) + (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          !VelGrad(2,3) =              + DelzV(i,j,k)/D(i,j)
          
          !VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
          !        (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,1) = DelxW(i,j,k) + (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
         
          !VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
          !        (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,2) = DelyW(i,j,k) + (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          !VelGrad(3,3) =              + DelzW(i,j,k)/D(i,j)
          

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (2./3.)*(1./(7.4+2.*smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 1./(185.2+3.0*dmax**2)
          c_2 = -1./(58.5+2.0*dmax**2)
          c_3 = 1./(370.4+3.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/Eps_Old(i,j,k)*  &
                   (VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   (VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*  &
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm 
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          ! M.Derakhti !user choose the appropriate production estimation
          ! all cases in Derakhti etal (2015a,b) the full 3d version has been used
          if (ProdType==1) then
             !only accounts for the vertical shear by the horizontal velocities
             do n = 1,2
             do m = 3,3
                Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
             enddo
             enddo
          elseif(ProdType==2) then
             !only accounts for the vertical shear
             do n = 1,3
             do m = 3,3
                Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
             enddo
             enddo
          else
             !full 3d
             do n = 1,3
             do m = 1,3
                Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
             enddo
             enddo
          endif
          
          !no negative production at the surface
          if(k==Kend.and.Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

        endif
      elseif(ke_model==3) then
        ! Following Mayer and Madsen (2000), instead of determining the production on the
        ! basis of the strain rate, the production is based on the rotation of the velocity field.
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first                              
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! vorticity field and production
          do n = 1,3
          do m = 1,3
            Vorticity(n,m) = VelGrad(n,m)-VelGrad(m,n)
            Prod_s(i,j,k) = Prod_s(i,j,k)+CmuVt(i,j,k)*Vorticity(n,m)*Vorticity(n,m)
          enddo
          enddo
        endif
      endif
    enddo
    enddo
    enddo
    
    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
      endif
    enddo
    enddo
    enddo
 
    ! flux Richardson number
    Richf = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Prod_s(i,j,k)>Zero) then
        Richf(i,j,k) = -Prod_b(i,j,k)/(Prod_s(i,j,k)+1.0e-16)
        Richf(i,j,k) = dmax1(0.0,dmin1(0.21,Richf(i,j,k)))
      endif
    enddo
    enddo
    enddo

# if defined (POROUSMEDIA)
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Porosity(i,j,k)<0.99) then
        Prod_s(i,j,k) = 0.0
        Prod_b(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

# if defined (VEGETATION)
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                 
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        else  ! flexible vegetation                                                                          
          if(sigc(k)*D(i,j)<=FVegH(i,j)) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                       
                    VegH/FVegH(i,j)*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        endif
      else
        Prod_v(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! added by M.Derakhti
    ! need constructed values of TKE and Epsilon at the free surface
    ! at the bottom the values of Kbeg are used
    DelzTke = Zero
    DelzEps = Zero
    call delzFun_3D(Tke_old,DelzTke)
    call delzFun_3D(Eps_old,DelzEps)
    call construct_3D_z(Tke_old,DelzTke,TkezL,TkezR)
    call construct_3D_z(Eps_old,DelzEps,EpszL,EpszR) 
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc1
      TKezL(i,j,k) = dmax1(TKezL(i,j,k),TKe_min)
      TKezR(i,j,k) = dmax1(TKezR(i,j,k),TKe_min)
      EpszL(i,j,k) = dmax1(EpszL(i,j,k),Eps_min)
      EpszR(i,j,k) = dmax1(EpszR(i,j,k),Eps_min)
    enddo
    enddo
    enddo
    call delxFun1_3D(TkezL,DelxTkezL)
    call delxFun1_3D(EpszL,DelxEpszL)
    call delyFun1_3D(TkezL,DelyTkezL)
    call delyFun1_3D(EpszL,DelyEpszL)

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)
    
    !Modified by M.Derakhti in the same way for DU/DV/DW calculations
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle
      DxD = (DelxH(i,j)+DelxEta(i,j))/D(i,j)
      DyD = (DelyH(i,j)+DelyEta(i,j))/D(i,j)
      Nlen = 0
      do k = Kbeg,Kend
           !linear interpolation using lagrange polynominal to get values at k+1/2 and k-1/2
            L1top = dsig(k+1)  /(dsig(k)+dsig(k+1))
            L2top = dsig(k)    /(dsig(k)+dsig(k+1))
            L1bot = dsig(k)    /(dsig(k)+dsig(k-1))
            L2bot = dsig(k-1)  /(dsig(k)+dsig(k-1))
           !these are used for vertical gradient in non-uniform grid, see Derakhti etal 2015a Appendix B
            dsigck  = (dsig(k)+dsig(k+1))/2.0 
            dsigck1 = (dsig(k-1)+dsig(k))/2.0 
           !viscosity at k+1/2 and k-1/2 using linear interpolation method
           !CmuR is needed here? 
            nuV_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+  &
                (L1top*CmuR(i,j,k)+L2top*CmuR(i,j,k+1))+(L1top*CmuVt(i,j,k)+  &
                L2top*CmuVt(i,j,k+1))/Sche
            nuH_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+(L1top*CmuR(i,j,k)+  &
                L2top*CmuR(i,j,k+1))+(L1top*CmuHt(i,j,k)+L2top*CmuHt(i,j,k+1))/Sche
            nuV_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+(L1bot*CmuR(i,j,k-1)+  &
                L2bot*CmuR(i,j,k))+(L1bot*CmuVt(i,j,k-1)+L2bot*CmuVt(i,j,k))/Sche
            nuH_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+(L1bot*CmuR(i,j,k-1)+  &
                L2bot*CmuR(i,j,k))+(L1bot*CmuHt(i,j,k-1)+L2bot*CmuHt(i,j,k))/Sche
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k)+  &
                cfe*Prod_v(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            
# else
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+  &
                c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
# endif
          if(k==Kbeg) then
           R5(i,j,k) = R5(i,j,k)& ! a plus deleted by cheng
             -( & !(nuH_bot*(DelxH(i,j)**2+DelyH(i,j)**2)+nuV_bot)/D(i,j)**2*  &
                  !(Eps_old(i,j,k)-Eps_old(i,j,k-1))*D(i,j)/dsigck1 &
               -(nuH_bot*(DelxH(i,j)*DxD+DelyH(i,j)*DyD)*Eps_old(i,j,k))&
              )/dsig(k) &
             + (dsig(k)-dsigck)/(dsig(k)**2/2.0*(dsigck+dsig(k)/2.0))&
              *(nuH_top*(DelxSl(i,j,k+1)*DelxH(i,j)+DelySl(i,j,k+1)*DelyH(i,j))+nuV_top/D(i,j))&
              *Eps_old(i,j,k)
          elseif (k==Kend) then
           R5(i,j,k) = R5(i,j,k)& ! a plus deleted by cheng
             +( (nuH_top*(DelxEta(i,j)**2+DelyEta(i,j)**2)+nuV_top)/D(i,j)**2*  &
              (Eps_old(i,j,k+1)-Eps_old(i,j,k))*D(i,j)/dsigck &
               -(nuH_top*(-DelxEta(i,j)*DxD-DelyEta(i,j)*DyD)*EpszL(i,j,k+1))&
              )/dsig(k) &
             + (dsig(k)-dsigck1)/(dsig(k)**2/2.0*(dsigck1+dsig(k)/2.0))&
              *(nuH_bot*(-DelxSl(i,j,k)*DelxEta(i,j)-DelySl(i,j,k)*DelyEta(i,j))+nuV_bot/D(i,j))&
              *EpszL(i,j,k+1)
           endif

        Nlen = Nlen+1
        !Modified by M.derakhti similar to UD/DV/DW calculations
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           elseif(k==Kend) then
             Acoef(Nlen) = -dt*(3.0/2.0*(nuH_bot*(DelxSl(i,j,k)*DelxSc(i,j,k-1)+  &
                                 DelySl(i,j,k)*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                 /(dsigck1*(dsigck1+dsig(k)/2.0)) )
           else
             Acoef(Nlen) = -dt*( &
                                ( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k-1)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k-1))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k-1)+  &
                                  DelySl(i,j,k  )*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                )*alpha_c/dsig(k) &
                               +( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k-1)+  &
                                 DelySl(i,j,k+1)*DelySc(i,j,k-1))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k-1)+  &
                                 DelySl(i,j,k  )*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                )/(dsigck1*(dsigck1+dsigck)) &
                               )
           endif

           if(k==Kbeg) then
            Ccoef(Nlen) = -dt*(3.0/2.0*(nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 /(dsigck*(dsigck+dsig(k)/2.0)) )
           elseif(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt*( &
                               ( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k+1)+  &
                                   DelySl(i,j,k  )*DelySc(i,j,k+1))+nuV_bot/D(i,j)**2) &
                                )*gamma_c/dsig(k) &
                               +( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k+1)+  &
                                   DelySl(i,j,k  )*DelySc(i,j,k+1))+nuV_bot/D(i,j)**2) &
                                )/(dsigck*(dsigck1+dsigck)) &
                               )
           endif
           if (RNG) then
           !calculate c2e dynamically based on Yaghot and Orszag (1986), Yakhot et al 1992
           S11 = DelxU(i,j,k)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
           S22 = DelyV(i,j,k)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
           S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
           S12 = 0.5*(DelyU(i,j,k)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
                      DelxV(i,j,k)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
           S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                   
                      DelxW(i,j,k)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
           S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  & 
                      DelyW(i,j,k)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
           lambda = Tke_old(i,j,k)/Eps_old(i,j,k)*dsqrt(2.0*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2)))
           c2e = 1.68 + cmiu*lambda**3*(1.0-lambda/4.38)/(1.0+0.012*lambda**3)
           c2e = dmax1(0.0,dmin1(c2e,4.0))
           endif

           if(k==Kbeg) then
             Bcoef(Nlen) = 1.0 +dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)&
                          -dt*(-(3.0/2.0-dsigck/dsig(k))*(nuH_top*  &
                      (DelxSl(i,j,k+1)*DelxSc(i,j,k)+DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                               /(dsig(k)/2.0*dsigck) )
           elseif(k==Kend) then
             Bcoef(Nlen) = 1.0 +dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)&
                           -dt*(-(3.0/2.0-dsigck1/dsig(k))*  &
                           (nuH_bot*(DelxSl(i,j,k)*DelxSc(i,j,k)+DelySl(i,j,k)*  &
                           DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                               /(dsig(k)/2.0*dsigck1) )
           else
             Bcoef(Nlen) = 1.0 +dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)&
                           -dt*( &
                                ( (nuH_top*(DelxSl(i,j,k+1)*  &
                    DelxSc(i,j,k)+DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*  &
                    DelxSc(i,j,k)+DelySl(i,j,k  )*DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                                )*beta_c/dsig(k) &
                  -( (nuH_top*(DelxSl(i,j,k+1)*  &
                    DelxSc(i,j,k)+DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*  &
                    DelxSc(i,j,k)+DelySl(i,j,k  )*DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                                )/(dsigck*dsigck1) &
                               )
          endif

        !if(k==Kbeg) then
        !  Acoef(Nlen) = 0.0
        !else
        !  Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
        !      0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        !endif

        !if(k==Kend) then
        !  Ccoef(Nlen) = 0.0
        !else
        !  Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
        !      0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        !endif
        
        !Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = dmax1(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    !Modified by M.Derakhti in the same way for DU/DV/DW calculations
    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
     if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle
      DxD = (DelxH(i,j)+DelxEta(i,j))/D(i,j)
      DyD = (DelyH(i,j)+DelyEta(i,j))/D(i,j)
      Nlen = 0
      do k = Kbeg,Kend
            !linear interpolation using lagrange polynominal to get values at k+1/2 and k-1/2
            L1top = dsig(k+1)/(dsig(k)+dsig(k+1))
            L2top = dsig(k)  /(dsig(k)+dsig(k+1))
            L1bot = dsig(k)/(dsig(k)+dsig(k-1))
            L2bot = dsig(k-1)  /(dsig(k)+dsig(k-1))
           !these are used for vertical gradient in non-uniform grid, see Derakhti etal 2015a Appendix B
            dsigck  = (dsig(k)+dsig(k+1))/2.0 
            dsigck1 = (dsig(k-1)+dsig(k))/2.0 
           !viscosity at k+1/2 and k-1/2 using linear interpolation method
           !CmuR is needed here? 
           !In this routine, kepsilon_3d,we have CmuHt = CmuVt
            nuV_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+  &
       (L1top*CmuR(i,j,k)+L2top*CmuR(i,j,k+1))+  &
       (L1top*CmuVt(i,j,k)+L2top*CmuVt(i,j,k+1))/Schk
            nuH_top = (L1top*Cmu(i,j,k)+L2top*Cmu(i,j,k+1))+  &
       (L1top*CmuR(i,j,k)+L2top*CmuR(i,j,k+1))+  &
       (L1top*CmuHt(i,j,k)+L2top*CmuHt(i,j,k+1))/Schk
            nuV_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+  &
       (L1bot*CmuR(i,j,k-1)+L2bot*CmuR(i,j,k))+  &
       (L1bot*CmuVt(i,j,k-1)+L2bot*CmuVt(i,j,k))/Schk
            nuH_bot = (L1bot*Cmu(i,j,k-1)+L2bot*Cmu(i,j,k))+  &
       (L1bot*CmuR(i,j,k-1)+L2bot*CmuR(i,j,k))+  &
       (L1bot*CmuHt(i,j,k-1)+L2bot*CmuHt(i,j,k))/Schk
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+  &
              Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)      
# else
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                   
                
# endif
        if(k==Kbeg) then
           R5(i,j,k) = R5(i,j,k)& ! a plus deleted by cheng
             -( &!(nuH_bot*(DelxH(i,j)**2+DelyH(i,j)**2)+nuV_bot)/D(i,j)**2*  &
                 !(Tke_old(i,j,k)-Tke_old(i,j,k-1))*D(i,j)/dsigck1 &
               -(nuH_bot*(DelxH(i,j)*DxD+DelyH(i,j)*DyD)*Tke_old(i,j,k))&
              )/dsig(k) &
             + (dsig(k)-dsigck)/(dsig(k)**2/2.0*(dsigck+dsig(k)/2.0))&
              *(nuH_top*(DelxSl(i,j,k+1)*DelxH(i,j)+DelySl(i,j,k+1)*DelyH(i,j))+nuV_top/D(i,j))&
              *Tke_old(i,j,k)
          elseif (k==Kend) then
           R5(i,j,k) = R5(i,j,k)& ! a plus deleted by cheng
             +( (nuH_top*(DelxEta(i,j)**2+DelyEta(i,j)**2)+nuV_top)/D(i,j)**2*  &
                (Tke_old(i,j,k+1)-Tke_old(i,j,k))*D(i,j)/dsigck &
               -(nuH_top*(-DelxEta(i,j)*DxD-DelyEta(i,j)*DyD)*TkezL(i,j,k+1))&
              )/dsig(k) &
             + (dsig(k)-dsigck1)/(dsig(k)**2/2.0*(dsigck1+dsig(k)/2.0))&
              *(nuH_bot*(-DelxSl(i,j,k)*DelxEta(i,j)-DelySl(i,j,k)*DelyEta(i,j))+nuV_bot/D(i,j))&
              *TkezL(i,j,k+1)
        endif

        Nlen = Nlen+1
        !Modified by M.derakhti similar to UD/DV/DW calculations
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           elseif(k==Kend) then
             Acoef(Nlen) = -dt*(3.0/2.0*(nuH_bot*(DelxSl(i,j,k)*DelxSc(i,j,k-1)+  &
                           DelySl(i,j,k)*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                 /(dsigck1*(dsigck1+dsig(k)/2.0)) )
           else
             Acoef(Nlen) = -dt*( &
                                ( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k-1)+  &
                                 DelySl(i,j,k+1)*DelySc(i,j,k-1))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k-1)+  &
                                  DelySl(i,j,k  )*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                )*alpha_c/dsig(k) &
                               +( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k-1)+  &
                                 DelySl(i,j,k+1)*DelySc(i,j,k-1))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k-1)+  &
                                  DelySl(i,j,k  )*DelySc(i,j,k-1))+nuV_bot/D(i,j)**2) &
                                )/(dsigck1*(dsigck1+dsigck)) &
                               )
           endif

           if(k==Kbeg) then
             Ccoef(Nlen) = -dt*(3.0/2.0*(nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 /(dsigck*(dsigck+dsig(k)/2.0)) )
           elseif(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt*( &
                                ( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                 DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k+1)+  &
                                 DelySl(i,j,k  )*DelySc(i,j,k+1))+nuV_bot/D(i,j)**2) &
                                )*gamma_c/dsig(k) &
                               +( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k+1)+  &
                                 DelySl(i,j,k+1)*DelySc(i,j,k+1))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k+1)+  &
                                 DelySl(i,j,k  )*DelySc(i,j,k+1))+nuV_bot/D(i,j)**2) &
                                )/(dsigck*(dsigck1+dsigck)) &
                               )
           endif
           
           if(k==Kbeg) then
             Bcoef(Nlen) = 1.0 &
                          -dt*(-(3.0/2.0-dsigck/dsig(k))*(nuH_top*(DelxSl(i,j,k+1)*  &
                            DelxSc(i,j,k)+DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                               /(dsig(k)/2.0*dsigck) )
           elseif(k==Kend) then
             Bcoef(Nlen) = 1.0 &
                           -dt*(-(3.0/2.0-dsigck1/dsig(k))*(nuH_bot*(DelxSl(i,j,k)*DelxSc(i,j,k)+  &
                              DelySl(i,j,k)*DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                               /(dsig(k)/2.0*dsigck1) )
           else
             Bcoef(Nlen) = 1.0 &
                           -dt*( &
                                ( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                                 -(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k)+  &
                                   DelySl(i,j,k  )*DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                                )*beta_c/dsig(k) &
                               -( (nuH_top*(DelxSl(i,j,k+1)*DelxSc(i,j,k)+  &
                                  DelySl(i,j,k+1)*DelySc(i,j,k))+nuV_top/D(i,j)**2) &
                                 +(nuH_bot*(DelxSl(i,j,k  )*DelxSc(i,j,k)+  &
                                  DelySl(i,j,k  )*DelySc(i,j,k))+nuV_bot/D(i,j)**2) &
                                )/(dsigck*dsigck1) &
                               )
          endif

        !if(k==Kbeg) then
        !  Acoef(Nlen) = 0.0
        !else
        !  Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
        !      0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        !endif
        !
        !if(k==Kend) then
        !  Ccoef(Nlen) = 0.0
        !else
        !  Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
        !      0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        !endif
        !
        !Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = dmax1(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

   ! Modified by M.Derakhti as done in "vel_bc" subroutine
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
       !these are added to calculate the tangential external stress parallel to the surface and bed
       Abot = dsqrt(1.0+DelxH(i,j)**2+DelyH(i,j)**2)
       Bbot = dsqrt(1.0+DelxH(i,j)**2)
       !velocities parallel to the bed
       DUprime = (DU(i,j,Kbeg)-DelxH(i,j)*DW(i,j,Kbeg))/Bbot 
       DVprime = (-DelxH(i,j)*DelyH(i,j)*DU(i,j,Kbeg)+(1.0+DelxH(i,j)**2)*  &
                 DV(i,j,Kbeg)-DelyH(i,j)*DW(i,j,Kbeg))/Abot/Bbot
       DUmag = dsqrt(DUprime**2+DVprime**2)

      if(D(i,j)<Dmin.or.Mask(i,j)==0) cycle
      
       if(DUmag/D(i,j)<1.e-6) then
        Tkeb = Tke_min
        Epsb = Eps_min
       
        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
       else
      ! impose wall function 
      !Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      !Zdis = 0.5*dsig(Kbeg)*D(i,j)
      !
      !if(Ibot==1) then
      !   FricU = sqrt(Cd0)*Umag
      !else
      !!rough wall 
      !   FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
      !endif
      !Tkeb = FricU**2/sqrt(cmiu)
      !Epsb = FricU**3/(Kappa*Zdis)

       !first grid point distance from the bed
       Dz1 = 0.5*D(i,j)*dsig(Kbeg)/Abot
       !estimate friction velocity
       Z0 = Zob / 30.0       
       if(ibot==1) then
         Cdrag = Cd0
       else
# if defined (SEDIMENT)
         Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(Dz1/Z0))**2                                                          
# else
         Cdrag = (Kappa/log(Dz1/Z0))**2
# endif
       endif
       Cdrag = dmin1(Cdrag,(Cmu(i,j,Kbeg)+CmuR(i,j,Kbeg)+CmuVt(i,j,Kbeg))/Dz1/(DUmag+Small)*D(i,j))
       Ustar = dsqrt(Cdrag)*DUmag/D(i,j)

       Tkeb = dmax1(Ustar**2/sqrt(cmiu),Tke_min)
       Epsb = dmax1(Ustar**3/(Kappa*Dz1),Eps_min)
       DTke(i,j,Kbeg) = D(i,j)*Tkeb
       DEps(i,j,Kbeg) = D(i,j)*Epsb
      endif

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg)
      enddo
    enddo
    enddo

    !rewritten by M.Derakhti as in derakhti etal 2015a, equation (60)
    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    if (Mask(i,j)==0.or.D(i,j)<Dmin) cycle
      do k = 1,Nghost
           Atop = dsqrt(1.0+DelxEta(i,j)**2+DelyEta(i,j)**2)
           Btop = dsqrt(1.0+DelxEta(i,j)**2)
            !DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
            DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1) &
           +(sigc(Kend+k)-sigc(Kend-k+1))*D(i,j)**2/Atop**2 * &
            ( DelxEta(i,j)*DelxTkezL(i,j,Kend+1) & 
             +DelyEta(i,j)*DelyTkezL(i,j,Kend+1) &
            )
            !DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
            DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1) &
           +(sigc(Kend+k)-sigc(Kend-k+1))*D(i,j)**2/Atop**2 * &
            ( DelxEta(i,j)*DelxEpszL(i,j,Kend+1) &
             +DelyEta(i,j)*DelyEpszL(i,j,Kend+1) &
            )
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
	   if(Bc_X0==5) then !added by Cheng for wall friction
         Umag = sqrt(V(Ibeg,j,k)**2+W(Ibeg,j,k)**2)
         if(Umag<1.e-6) then
           Tkeb = Tke_min
           Epsb = Eps_min
           DTke(Ibeg,j,k) = D(Ibeg,j)*Tkeb
           DEps(Ibeg,j,k) = D(Ibeg,j)*Epsb
         else
           if(Ibot==1) then
             FricU = sqrt(Cd0)*Umag
           else
             FricU = Umag/(1./Kappa*log(30.*dx/2.0/Zob))
           endif
           Tkeb = FricU**2/sqrt(cmiu)
           Epsb = FricU**3/(Kappa*dx/2.0)
           DTke(Ibeg,j,k) = D(Ibeg,j)*Tkeb
           DEps(Ibeg,j,k) = D(Ibeg,j)*Epsb
		 endif
	   endif
       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
	   if(Bc_Xn==5) then !added by Cheng for wall friction
         Umag = sqrt(V(Iend,j,k)**2+W(Iend,j,k)**2)
         if(Umag<1.e-6) then
           Tkeb = Tke_min
           Epsb = Eps_min
           DTke(Iend,j,k) = D(Iend,j)*Tkeb
           DEps(Iend,j,k) = D(Iend,j)*Epsb
         else
           if(Ibot==1) then
             FricU = sqrt(Cd0)*Umag
           else
             FricU = Umag/(1./Kappa*log(30.*dx/2.0/Zob))
           endif
           Tkeb = FricU**2/sqrt(cmiu)
           Epsb = FricU**3/(Kappa*dx/2.0)
           DTke(Iend,j,k) = D(Iend,j)*Tkeb
           DEps(Iend,j,k) = D(Iend,j)*Epsb
		 endif
	   endif
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
	   if(Bc_Y0==5) then !added by Cheng for wall friction
         Umag = sqrt(U(i,Jbeg,k)**2+W(i,Jbeg,k)**2)
         if(Umag<1.e-6) then
           Tkeb = Tke_min
           Epsb = Eps_min
           DTke(i,Jbeg,k) = D(i,Jbeg)*Tkeb
           DEps(i,Jbeg,k) = D(i,Jbeg)*Epsb
         else
           if(Ibot==1) then
             FricU = sqrt(Cd0)*Umag
           else
             FricU = Umag/(1./Kappa*log(30.*dy/2.0/Zob))
           endif
           Tkeb = FricU**2/sqrt(cmiu)
           Epsb = FricU**3/(Kappa*dy/2.0)
           DTke(i,Jbeg,k) = D(i,Jbeg)*Tkeb
           DEps(i,Jbeg,k) = D(i,Jbeg)*Epsb
		 endif
	   endif
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
	   if(Bc_Yn==5) then !added by Cheng for wall friction
         Umag = sqrt(U(i,Jend,k)**2+W(i,Jend,k)**2)
         if(Umag<1.e-6) then
           Tkeb = Tke_min
           Epsb = Eps_min
           DTke(i,Jend,k) = D(i,Jend)*Tkeb
           DEps(i,Jend,k) = D(i,Jend)*Epsb
         else
           if(Ibot==1) then
             FricU = sqrt(Cd0)*Umag
           else
             FricU = Umag/(1./Kappa*log(30.*dy/2.0/Zob))
           endif
           Tkeb = FricU**2/sqrt(cmiu)
           Epsb = FricU**3/(Kappa*dy/2.0)
           DTke(i,Jend,k) = D(i,Jend)*Tkeb
           DEps(i,Jend,k) = D(i,Jend)*Epsb
		 endif
	   endif
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

   
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)

        ! If turbulence is too small, assume small eddy viscosity (added by Cheng)
        if(Tke(i,j,k)<1.e-6.or.Eps(i,j,k)<1.e-6) then
          CmuVt(i,j,k) = Cmut_min
        else
          CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
        endif
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Nghost
       CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg)
       CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend)
    enddo
    enddo
    enddo

# if defined (POROUSMEDIA)
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Porosity(i,j,k)<0.99) then
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo
# endif

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT'.or.WaveMaker(1:3)=='LEF' &
       .or. WaveMaker(1:3)=='FOC') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
        VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
        UpWp(i,j,k) = CmuVt(i,j,k)*(VelGrad(1,2)+VelGrad(2,1))
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)
    deallocate(TkezL)
    deallocate(EpszL)
    deallocate(TkezR)
    deallocate(EpszR)
    deallocate(DelzTke)
    deallocate(DelzEps)
    deallocate(DelxTkezL)
    deallocate(DelxEpszL)
    deallocate(DelyTkezL)
    deallocate(DelyEpszL)

    end subroutine kepsilon_3D


    subroutine kepsilon(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.01
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb,Xlfs,Epsfs
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    Sche = 1.3
    Schk = 1.0

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzU(i,j,k) = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        DelzV(i,j,k) = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_s(i,j,k) = CmuVt(i,j,k)*(DelzU(i,j,k)**2+DelzV(i,j,k)**2)   
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
     endif
    enddo
    enddo
    enddo

# if defined (VEGETATION)
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                 
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*StemD*VegDens*VegDrag*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        else  ! flexible vegetation                                                                          
          if(sigc(k)*D(i,j)<=FVegH(i,j)) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*StemD*VegDens*VegDrag*foliage(i,j)*  &                                       
                    VegH/FVegH(i,j)*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        endif
      else
        Prod_v(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities                                                                                               
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
       R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k)+  &
                    cfe*Prod_v(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            
# else
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k))*  &
                    Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
# endif
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
               0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo

      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = max(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)                                              
# else
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
# endif
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
               0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = max(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6.or.D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      Zdis = 0.5*dsig(Kbeg)*D(i,j)
      X0 = 0.05
      Iter = 0

      Xa = dlog(9.0*Umag*Zdis/Visc)
 10      Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
      if(Iter>=20) then
        write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
      endif
      if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
        X0 = Xn
        Iter = Iter+1
        goto 10
      else
        FricU = Xn*Umag
      endif

!      if(Ibot==1) then
!        FricU = sqrt(Cd0)*Umag
!      else
!        FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
!      endif

      Tkeb = FricU**2/sqrt(cmiu)
      Epsb = FricU**3/(Kappa*Zdis)

      DTke(i,j,Kbeg) = D(i,j)*Tkeb
      DEps(i,j,Kbeg) = D(i,j)*Epsb

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         do i = Ibeg,Ibeg+5
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif

       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j) 
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    end subroutine kepsilon



# if defined (BUBBLE)
    subroutine bslip_velocity
!----------------------------------------------------------
!   Specify bubble radius and calculate rise velocity
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------
    use global, only: SP,pi,Zero,Rho0,Mg,Rbg,DRbg,Wbg,Entrain,Con_b,Surface_Tension
    implicit none
    integer :: g
    real(SP) :: rlogR0,rlogRn,rlogR1,rlogR2,rlogN0,alpha_b,beta_b,sum_e
    real(SP), dimension(Mg) :: rlogN,specN,binN             

    ! specify bubble radius
    rlogR0 = -1.0
    rlogRn = 1.0
    do g = 1,Mg
      rlogR1 = rlogR0+(rlogRn-rlogR0)*(g-1)/float(Mg)
      rlogR2 = rlogR0+(rlogRn-rlogR0)*g/float(Mg)
      Rbg(g) = (0.5*(10**rlogR1+10**rlogR2))*0.001
      DRbg(g) = (10**rlogR2-10**rlogR1)*0.001
    enddo

    ! slip velocity
    do g = 1,Mg
      if(Rbg(g)<=7.0e-4) then
        Wbg(g) = 4474.*Rbg(g)**1.357
      elseif(Rbg(g)>5.1e-3) then
        Wbg(g) = 4.202*Rbg(g)**0.547
      else
        Wbg(g) = 0.23
      endif
    enddo

    ! read bubble size distribution in Deane and Stokes' paper
    rlogN0 = 4.3
    alpha_b = -3.0/2.0
    beta_b = -10.0/3.0

    rlogN(1) = rlogN0
    do g = 2,Mg/2
      rlogN(g) = rlogN(g-1)+alpha_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo
    do g = Mg/2+1,Mg
      rlogN(g) = rlogN(g-1)+beta_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo

    do g = 1,Mg
      specN(g) = 10**rlogN(g)*1.0e+6  ! convert to 1/(m^4)
      binN(g)  = specN(g)*DRbg(g)
    enddo

    ! entrainment coefficient
    sum_e = zero
    do g = 1,Mg
      sum_e = sum_e+Rbg(g)**2*binN(g)
    enddo

    do g = 1,Mg
      Entrain(g) = Con_b/(4.0*pi)*(Surface_Tension/Rho0)**(-1)*binN(g)/sum_e
    enddo

    end subroutine bslip_velocity

    
    subroutine eval_bub(ISTEP)
!---------------------------------------------------------------------
!   Update bubble concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,Phi,DPhi,bkup_pdf,DUfs,DVfs,Wfs
    real(SP), dimension(:,:), allocatable :: bkup_freq
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: bkup_c,bkup_beta,bkup_f,eps_crit,db_star,db_min,  &
                db_max,db_crit,delta_db,db_l,db_r,db_m,pdf_integral,Schb, &
                Ws
    integer :: i,j,k,g,l,m,n,IVAR,NGRD,IGRD,Nlen

    NGRD = 0
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      NGRD = NGRD+1
    enddo
    enddo
    enddo

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))
    allocate(Phi(Mloc,Nloc,Kloc))
    allocate(DPhi(Mloc,Nloc,Kloc))
    allocate(bkup_freq(NGRD,Mg))
    allocate(bkup_pdf(NGRD,Mg,Mg))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! bubble breakup (Martinez-Bazan et al., 1999, 2010)
    bkup_c = 0.25
    bkup_beta = 8.2

    bkup_freq = 0.0
    bkup_pdf = 0.0
    IGRD = 0
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      IGRD = IGRD+1
      do g = 2,Mg
        bkup_f = bkup_beta*(Eps(i,j,k)*2.0*Rbg(g))**(2.0/3.0)-12.0*Surface_Tension/(Rho0*2.0*Rbg(g))
        if(bkup_f>zero) then
          bkup_freq(IGRD,g) = bkup_c*sqrt(bkup_f)/(2.0*Rbg(g))
        else
          bkup_freq(IGRD,g) = zero
        endif
      enddo

      do n = 2,Mg
      do m = 1,Mg-1
        eps_crit = (12.0*Surface_Tension/(bkup_beta*Rho0))**1.5*  &
             (2.0*Rbg(n))**(-2.5)
        if(Eps(i,j,k)<=eps_crit) then
          bkup_pdf(IGRD,m,n) = zero
        elseif(Rbg(n)<=Rbg(m)) then
          bkup_pdf(IGRD,m,n) = zero
        else
          db_star = Rbg(g)/Rbg(n)
          db_crit = (12.0*Surface_Tension/(bkup_beta*Rho0))**0.6*  &
              Eps(i,j,k)**(-0.4)/(2.0*Rbg(n))
          db_min = (12.0*Surface_Tension/(bkup_beta*Rho0*2.0*Rbg(n)))**1.5*  &
              Eps(i,j,k)**(-1.0)/(2.0*Rbg(n))
          db_max = (1.0-db_min**3)**(1.0/3.0)

          if(db_max<=db_min.or.db_star<db_min.or.db_star>db_max) then
            bkup_pdf(IGRD,m,n) = zero
          else
            pdf_integral = zero
            delta_db = (db_max-db_min)/100.0
            do l = 1,100
              db_l = db_min+(l-1)*delta_db
              db_r = db_min+l*delta_db
              db_m = 0.5*(db_min+db_max)        
              
              pdf_integral = pdf_integral+  &
                 db_m**2.0*(db_m**(2./3.))-db_crit**(5./3.)*  &
                 ((1.0-db_m**3.)**(2./9.)-db_crit**(5./3.))*delta_db
            enddo

            bkup_pdf(IGRD,m,n) = db_star**2.0*(db_star**(2.0/3.0)-db_crit**(5.0/3.0))*  &
                 ((1.0-db_star**3.0)**(2.0/9.0)-db_crit**(5.0/3.0))/pdf_integral
          endif
        endif
      enddo
      enddo
    enddo
    enddo
    enddo

    do g = 1,Mg
      ! temporary arrays
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        Phi(i,j,k) = Nbg(i,j,k,g)
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      DUfs = Ex
      DVfs = Ey
      Wfs = Omega+Wbg(g)

      ! advection and diffusion 
      IVAR = 5     
      call adv_scalar_hlpa(DUfs,DVfs,Wfs,Phi,R5,IVAR)

      ! bubble entrainment at the surface                                                                
      k = Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        if(Eps(i,j,k)>Eps_Cric.and.Vbg(i,j,k)<0.4) then
          R5(i,j,k) = R5(i,j,k)+(1.0-Vbg(i,j,k))*Entrain(g)*D(i,j)*Eps(i,j,k)                                     
        endif
      enddo
      enddo

!      ! account for bubble breakup           
!      IGRD = 0  
!      do k = Kbeg,Kend  
!      do j = Jbeg,Jend  
!      do i = Ibeg,Iend  
!        IGRD = IGRD+1 
!        if(g>=2) then  ! sink   
!          R5(i,j,k) = R5(i,j,k)-bkup_freq(IGRD,g)*DNbg(i,j,k,g)  
!        endif 
!    
!        if(g<=Mg-1) then ! source  
!          do n = 2,Mg 
!            R5(i,j,k) = R5(i,j,k)+  &  
!                2.0*bkup_pdf(IGRD,g,n)/(2.0*Rbg(n))*bkup_freq(IGRD,n)*DNbg(i,j,k,n)*(2.0*DRbg(n)) 
!          enddo 
!        endif 
!      enddo 
!      enddo 
!      enddo

      Schb = 0.7
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
           Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
               0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schb)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schb)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          Rhs0(Nlen) = DPhi(i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DNbg(i,j,k,g) = Xsol(Nlen)
        enddo
      enddo
      enddo

      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DNbg(i,j,k,g) = ALPHA(ISTEP)*DNbg0(i,j,k,g)+BETA(ISTEP)*DNbg(i,j,k,g)
        if(Mask(i,j)==0) DNbg(i,j,k,g) = Zero
      enddo
      enddo
      enddo
    enddo

    ! collect data into ghost cells
# if defined (PARALLEL)
    do g = 1,Mg
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      call phi_3D_exch(DPhi)

      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DNbg(i,j,k,g) = DPhi(i,j,k)
      enddo
      enddo
      enddo
    enddo
# endif          

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Ibeg-i,j,k,g) = DNbg(Ibeg+i-1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Iend+i,j,k,g) = DNbg(Iend-i+1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jbeg-j,k,g) = DNbg(i,Jbeg+j-1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jend+j,k,g) = DNbg(i,Jend-j+1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do g = 1,Mg
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DNbg(i,j,Kbeg-k,g) = DNbg(i,j,Kbeg+k-1,g)
      enddo
      do k = 1,Nghost
        DNbg(i,j,Kend+k,g) = DNbg(i,j,Kend-k+1,g)
      enddo
    enddo
    enddo
    enddo

    ! update bubble number density
    Nbg = Zero
    do g = 1,Mg
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==1) then
        Nbg(i,j,k,g) = DNbg(i,j,k,g)/D(i,j)
      endif
    enddo
    enddo
    enddo
    enddo

    ! calculate void fraction
    Vbg = Zero
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      do g = 1,Mg
        Vbg(i,j,k) = Vbg(i,j,k)+4./3.*pi*Rbg(g)**3*Nbg(i,j,k,g)
        if(abs(Vbg(i,j,k))<1.e-16) Vbg(i,j,k) = 0.0
        Vbg(i,j,k) = min(0.4,Vbg(i,j,k))
      enddo
    enddo
    enddo
    enddo

    ! update cell density
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Rho(i,j,k) = (1.0-Vbg(i,j,k))*Rho0+Vbg(i,j,k)*RhoA
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)
    deallocate(Phi)
    deallocate(DPhi)
    deallocate(bkup_freq)
    deallocate(bkup_pdf)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    end subroutine eval_bub
# endif


# if defined (SEDIMENT)
    subroutine settling_velocity
!----------------------------------------------------------
!   Calculate settling velocity for sediment
!   Last update: Gangfeng Ma, 14/06/2012
!--------------------------------------------------------- 
    use global, only: SP,pi,Zero,Grav,Wset,Sd50,Visc,Rho0,Srho, &
                      Mloc,Nloc,Kloc,Conc,ntyws,Sedi_Ws
    implicit none
    real(SP) :: sr,c1,c2
    integer :: i,j,k

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(ntyws==1) then
        ! pre-defined
        Wset = Sedi_Ws
      elseif(ntyws==2) then
        ! related to d50
        sr = Srho/Rho0
        Wset = sqrt((sr-1.0)*Grav*Sd50)*  &
            (sqrt(2.0/3.0+36.0*Visc**2/((sr-1.0)*Grav*Sd50**3))-  &
             sqrt(36.0*Visc**2/((sr-1.0)*Grav*Sd50**3)))      
      elseif(ntyws==3) then
        sr = Srho/Rho0
        Wset = (sr-1.0)*Grav*Sd50**2/18.0/Visc
      endif
    enddo
    enddo
    enddo

    return
    end subroutine settling_velocity


    subroutine eval_sedi(ISTEP)
!---------------------------------------------------------------------
!   Update sediment concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,IVAR,Nlen

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))   

# if !defined (LANDSLIDE)
    ! sediment entrainment at the bottom
    if(trim(Sed_Type)=='COHESIVE') then
      call SSource
    else
      call SedPickup
    endif
# endif

    DUfs = Ex
    DVfs = Ey
    Wfs = Omega-Wset

    ! advection and horizontal diffusion 
    IVAR = 6     
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Conc,R5,IVAR)

    ! vertical diffusion
    if(VISCOUS_FLOW) then
      SchC = 0.7
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          if(k==Kbeg) then
            Rhs0(Nlen) = DConc(i,j,k)+dt*R5(i,j,k)-dt/D(i,j)*SSour(i,j)/dsig(k)
          else
            Rhs0(Nlen) = DConc(i,j,k)+dt*R5(i,j,k)
          endif
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DConc(i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update sediment concentration
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DConc(i,j,k) = ALPHA(ISTEP)*DConc0(i,j,k)+BETA(ISTEP)*DConc(i,j,k)
        if(Mask(i,j)==0) DConc(i,j,k) = Zero
        DConc(i,j,k) = dmax1(DConc(i,j,k),Zero)
      enddo
      enddo
      enddo
    else
      ! update sediment concentration 
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DConc(i,j,k) = ALPHA(ISTEP)*DConc0(i,j,k)+  &
                  BETA(ISTEP)*(DConc(i,j,k)+dt*R5(i,j,k)) 
        if(Mask(i,j)==0) DConc(i,j,k) = Zero
        DConc(i,j,k) = dmax1(DConc(i,j,k),Zero)
      enddo
      enddo
      enddo
    endif

    ! update sediment concentration
    Conc = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        Conc(i,j,k) = DConc(i,j,k)/D(i,j)
      endif
    enddo
    enddo
    enddo

    ! collect data into ghost cells
    call sedi_bc
# if defined (PARALLEL)
    call phi_3D_exch(Conc)
# endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Conc(i,j,k)>0.1) then
        write(12,*) DConc(i,j,k),D(i,j),Conc(i,j,k)
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    return
    end subroutine eval_sedi

    
    subroutine sedi_bc
! ------------------------------------------------------
!
!   boundary condition for sediment concentration
!
!------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_X0==3) then
          Conc(Ibeg-i,j,k) = Sed_X0(j,k)
        else
          Conc(Ibeg-i,j,k) = Conc(Ibeg+i-1,j,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
!        if(Bc_Xn==3) then
!          Conc(Iend+i,j,k) = Sed_Xn(j,k)
!        else
          Conc(Iend+i,j,k) = Conc(Iend-i+1,j,k)
!        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        Conc(i,Jbeg-j,k) = Conc(i,Jbeg+j-1,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        Conc(i,Jend+j,k) = Conc(i,Jend-j+1,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        Conc(i,j,Kbeg-k) = Conc(i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        Conc(i,j,Kend+k) = Conc(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    end subroutine sedi_bc

    subroutine SSource
!-----------------------------------------------------------
!   Sediment suspension and deposition at the bottom
!   Called by 
!      eval_sedi
!   Last update: 15/06/2012, Gangfeng Ma 
!-----------------------------------------------------------
    use global
    implicit none
    integer, parameter :: ntyss=1
    integer :: i,j,k,k0
    real(SP) :: Dz1,Cdrag,Qero,Qdep,c1,c2,c3,SStar,Umag,Um,Vm,Cm

    SSour = Zero
    Taub = Zero

    k = Kbeg
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

      if(ntyss==1) then
        ! sediment erosion/deposition depends on shear stress
        ! bottom shear stress             
        Dz1 = 0.5*D(i,j)*dsig(k)
        if(ibot==1) then
          Cdrag = Cd0
        else
# if defined (SEDIMENT)
          Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                          
# else
          Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
        endif
        Taub(i,j) = Rho0*Cdrag*(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)

        if(Taub(i,j)>=Tau_ce) then
          Qero = Erate*(Taub(i,j)/Tau_ce-1.0)
        else
          Qero = 0.0
        endif

        if(Taub(i,j)<=Tau_cd) then
          Qdep = Wset*Conc(i,j,k)*(1.0-Taub(i,j)/Tau_cd)
        else
          Qdep = 0.0
        endif

        SSour(i,j) = Qdep-Qero
      elseif(ntyss==2) then
        Um = Zero
        Vm = Zero
        Cm = Zero
        do k0 = Kbeg,Kend
          Um = Um+U(i,j,k0)/float(Kend-Kbeg+1)
          Vm = Vm+V(i,j,k0)/float(Kend-Kbeg+1)
          Cm = Cm+Conc(i,j,k0)/float(Kend-Kbeg+1)
        enddo

        ! use sediment carrying capacity concenpt
        c1 = 0.007
        c2 = 0.92
        c3 = 50.0
        Umag = sqrt(Um**2+Vm**2)
        SStar = c1*(Umag**3/D(i,j)/Wset)**c2
        SSour(i,j) = c3*Wset*(Cm-SStar)
      endif
    enddo
    enddo

    return
    end subroutine SSource

   
    subroutine SedPickup
!-----------------------------------------------------------
!   Calculate sand pickup 
!   Called by
!      eval_sedi
!   Last update: 02/04/2013, Gangfeng Ma
!-----------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: Dz1,Cdrag,sr,Qero,Qdep,Shields,taub1 
 
    SSour = Zero
    Taub = Zero

    k = Kbeg
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

      Dz1 = 0.5*D(i,j)*dsig(Kbeg)
      if(ibot==1) then
        Cdrag = Cd0
      else
        Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2 
      endif
      Taub1 = Cdrag*(U(i,j,k)**2+V(i,j,k)**2)

      ! Shields number
      sr = Srho/Rho0
      Shields = Taub1/(sr-1.0)/Grav/SD50

      ! upper-bound shields number for the start of sheet flow
      Shields = dmin1(Shields,1.0) 

      ! Pickup rate
      if(Shields>Shields_c) then
        Qero = 0.00033*((Shields-Shields_c)/Shields_c)**1.5*  &
             (sr-1.0)**0.6*Grav**0.6*SD50**0.8/Visc**0.2
      else
        Qero = Zero
      endif

      ! Deposition rate
      Qdep = Wset*Conc(i,j,k)
    
      ! Bottom exchange rate
      SSour(i,j) = Qdep-Qero
   
    enddo
    enddo

    end subroutine SedPickup


    subroutine eval_bedload
!---------------------------------------------------------
!   Calculate bed load transport rate
!   Called by 
!      main
!   Last update: 12/07/2012, Gangfeng Ma
!--------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: Taubx,Tauby,Cdrag,Dz1,sr,Shields,  &
                Shields_x,Shields_y,fw,Umag

    k = Kbeg
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

      ! bottom shear stress     
!      Dz1 = 0.5*D(i,j)*dsig(Kbeg)
!      if(ibot==1) then
!        Cdrag = Cd0
!      else
!        Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2     
!      endif
!      Taubx = Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*U(i,j,Kbeg)
!      Tauby = Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*V(i,j,Kbeg)
 
      fw = 0.01
      Taubx = 0.5*fw*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*U(i,j,Kbeg)
      Tauby = 0.5*fw*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*V(i,j,Kbeg)

      ! Shields number
      sr = Srho/Rho0
      Shields_x = Taubx/(sr-1.0)/Grav/SD50
      Shields_y = Tauby/(sr-1.0)/Grav/SD50
      Shields = sqrt(Shields_x**2+Shields_y**2)

      ! Bed-load transport rate                     
      if(Shields>Shields_c) then
        Qbedx(i,j) = 11.0*(Shields-Shields_c)**1.65*Shields_x/Shields*  &
                sqrt((sr-1)*Grav*SD50**3)
        Qbedy(i,j) = 11.0*(Shields-Shields_c)**1.65*Shields_y/Shields*  &
                sqrt((sr-1)*Grav*SD50**3)

        if(x(i)<=1.0) then
          Qbedx(i,j) = 0.0
          Qbedy(i,j) = 0.0
        endif
      else
        Qbedx(i,j) = Zero
        Qbedy(i,j) = Zero
      endif
    enddo
    enddo

    call phi_2D_coll(Qbedx)
    call phi_2D_coll(Qbedy)

    end subroutine eval_bedload


    subroutine update_bed(ISTEP)
!-----------------------------------------------------------
!   Update bed elevation 
!   Called by 
!      main
!   Last update: 15/06/2012, Gangfeng Ma
!-----------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    integer :: i,j,I1,I2
    real(SP) :: rsedi
    character(len=2)  :: FILE_NAME
    character(len=80) :: FDIR,file

    if(ISTEP==2) then

      Update_Bed_T = Update_Bed_T+dt

      do j = Jbeg,Jend
      do i = Ibeg,Iend
        TotDep(i,j) = TotDep(i,j)+SSour(i,j)*Srho*dt
      enddo
      enddo

!      ! assuming equilibrium beach profile
!      ! sediment suspension is not considered
!      do j = Jbeg,Jend
!      do i = Ibeg,Iend
!        Bed(i,j) = Bed(i,j)+dt*(-(Qbedx(i+1,j)-Qbedx(i-1,j))/(2.0*dx)-  &
!            (Qbedy(i,j+1)-Qbedy(i,j-1))/(2.0*dy))/(1.0-Spor)
!         Bed(i,j) = Bed(i,j)+dt*SSour(i,j)/(1.0-Spor)
!      enddo
!      enddo

      ! output bed elevation
      if(OUT_G.and.Update_Bed_T>=300.) then
        nbed = nbed+1

        FDIR = TRIM(RESULT_FOLDER)

        I1 = mod(nbed/10,10)
        I2 = mod(nbed,10)

        write(FILE_NAME(1:1),'(I1)') I1
        write(FILE_NAME(2:2),'(I1)') I2

        file=TRIM(FDIR)//'depos_'//TRIM(FILE_NAME)
        call putfile2D(file,TotDep)

!        file=TRIM(FDIR)//'bed_'//TRIM(FILE_NAME)
!        call putfile2D(file,Bed)

        Update_Bed_T = 0.0
      endif


      ! update water depth
!      do j = Jbeg,Jend
!      do i = Ibeg,Iend 
!        Hc(i,j) = Hc0(i,j)-Bed(i,j)
!      enddo
!      enddo             
!            
!      call phi_2D_coll(Hc)  
!
!      ! reconstruct depth at x-y faces 
!      do j = 1,Nloc 
!      do i = 2,Mloc
!        Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i+1,j)) 
!      enddo 
!      Hfx(1,j) = Hfx(2,j)
!      Hfx(Mloc1,j) = Hfx(Mloc,j) 
!      enddo 
!
!      do i = 1,Mloc 
!      do j = 2,Nloc 
!        Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j+1))  
!      enddo
!      Hfy(i,1) = Hfy(i,2) 
!      Hfy(i,Nloc1) = Hfy(i,Nloc) 
!      enddo
!
!      ! derivatives of water depth at cell center 
!      do j = 1,Nloc 
!      do i = 1,Mloc
!        DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx 
!        DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy 
!      enddo 
!      enddo
    endif

    return
    end subroutine update_bed
# endif


# if defined (SALINITY)
    subroutine eval_sali(ISTEP)
!---------------------------------------------------------------------
!   Update sediment concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,IVAR,Nlen

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! advection and diffusion 
    IVAR = 3 
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Sali,R5,IVAR)

    if(VISCOUS_FLOW) then
      SchC = 1.0
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          Rhs0(Nlen) = DSali(i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DSali(i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update salinity
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DSali(i,j,k) = ALPHA(ISTEP)*DSali0(i,j,k)+BETA(ISTEP)*DSali(i,j,k)
        if(Mask(i,j)==0) DSali(i,j,k) = Zero
      enddo
      enddo
      enddo
    else
      ! update salinity
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DSali(i,j,k) = ALPHA(ISTEP)*DSali0(i,j,k)+BETA(ISTEP)*(DSali(i,j,k)+dt*R5(i,j,k))
        if(Mask(i,j)==0) DSali(i,j,k) = Zero
      enddo
      enddo
      enddo
    endif

    ! update Salinity 
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Sali(i,j,k) = DSali(i,j,k)/D(i,j)
    enddo
    enddo
    enddo

    ! boundary condition and ghost cells
    call sali_bc
# if defined (PARALLEL)
    call phi_3D_exch(Sali)
# endif          

    deallocate(R5)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    return
    end subroutine eval_sali

    subroutine sali_bc
!-------------------------------------------------------------------
!
!   Boundary conditions for salinity
!
!------------------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
          Sali(Ibeg-i,j,k) = Sali(Ibeg+i-1,j,k)
        elseif(Bc_X0==3) then
          Sali(Ibeg-i,j,k) = 2.0*Sin_X0(j,k)-Sali(Ibeg+i-1,j,k)
        elseif(Bc_X0==4) then
          Sali(Ibeg-i,j,k) = Sali(Ibeg+i-1,j,k)
        elseif(Bc_X0==8) then
          Sali(Ibeg-i,j,k) = Sin_X0(j,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
          Sali(Iend+i,j,k) = Sali(Iend-i+1,j,k)
        elseif(Bc_Xn==3) then
          Sali(Iend+i,j,k) = 2.0*Sin_Xn(j,k)-Sali(Iend-i+1,j,k)
        elseif(Bc_Xn==4) then
          Sali(Iend+i,j,k) = Sali(Iend-i+1,j,k)
        elseif(Bc_Xn==8) then
          Sali(Iend+i,j,k) = Sin_Xn(j,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
          Sali(i,Jbeg-j,k) = Sali(i,Jbeg+j-1,k)
        elseif(Bc_Y0==4) then
          Sali(i,Jbeg-j,k) = Sali(i,Jbeg+j-1,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
          Sali(i,Jend+j,k) = Sali(i,Jend-j+1,k)
        elseif(Bc_Yn==4) then
          Sali(i,Jend+j,k) = Sali(i,Jend-j+1,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        Sali(i,j,Kbeg-k) = Sali(i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        Sali(i,j,Kend+k) = Sali(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    end subroutine sali_bc
# endif

    subroutine eval_dens
!---------------------------------------------------------------------
!
!   equation of state
!
!---------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 101
    integer, parameter :: kbbm1 = kbb-1
    real(SP), dimension(kbbm1) :: phy_z,rhoztmp,rhomean
    real(SP), dimension(Mloc,Nloc,kbbm1) :: rhoz
    real(SP), dimension(Kloc) :: zm,rhos
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,rhozloc
    real(SP),dimension(Mglob,Nglob,kbbm1) :: rhozglob
# if defined (PARALLEL)
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer,dimension(NumP) :: npxs,npys
# endif
    real(SP) :: DELTZ,TF,SF,RHOF,HMAX,ETAMAX
    real(SP) :: tmp1,tmp2,myvar
    integer :: I,J,K,isum,jk,iglob,jglob,kk,n,len,nreq,NKloc

    ! calculate density from equation of state
    DO K = 1,KLOC
    DO J = 1,NLOC
    DO I = 1,MLOC
      IF(Mask(I,J)==0) cycle
# if defined (SALINITY)
      TF = Temp(I,J,K)
      SF = Sali(I,J,K)
      RHOF = SF*SF*SF*6.76786136E-6_SP-SF*SF*4.8249614E-4_SP+ &
           SF*8.14876577E-1_SP-0.22584586E0_SP
      RHOF = RHOF*(TF*TF*TF*1.667E-8_SP-TF*TF*8.164E-7_SP+ &
           TF*1.803E-5_SP)
      RHOF = RHOF+1.-TF*TF*TF*1.0843E-6_SP+TF*TF*9.8185E-5_SP-TF*4.786E-3_SP
      RHOF = RHOF*(SF*SF*SF*6.76786136E-6_SP-SF*SF*4.8249614E-4_SP+  &
           SF*8.14876577E-1_SP+3.895414E-2_SP)
      RHOF = RHOF-(TF-3.98_SP)**2*(TF+283.0_SP)/(503.57_SP*(TF+67.26_SP))
      Rho(I,J,K) = Rho0+RHOF
# endif

# if defined (SEDIMENT)
      Rho(I,J,K) = (1.0-Conc(i,j,k))*Rho0+Conc(i,j,k)*SRho
# endif
    ENDDO
    ENDDO
    ENDDO

    ! find maximum water depth
    tmp1 = -large
    tmp2 = -large
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==0) cycle 
      if(hc(i,j)>tmp1) tmp1 = Hc(i,j)  
      if(eta(i,j)>tmp2) tmp2 = Eta(i,j)
    enddo
    enddo
# if defined (PARALLEL)
    call MPI_ALLREDUCE(tmp1,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
    hmax = myvar
    call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
    etamax = myvar
# endif

    ! interpolate into physical z levels
    deltz = (hmax+etamax)/float(kbbm1)
    do k = 1,kbbm1
      phy_z(k) = (float(k)-0.5)*deltz-hmax
    enddo

    rhoz = Rho0
    do i = 1,Mloc
    do j = 1,Nloc
      if(Mask(i,j)==0) cycle

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)    
        rhos(k) = Rho(i,j,k)
      enddo

      call sinter(zm,rhos,phy_z,rhoztmp,Kloc,kbbm1)

      do k = 1,kbbm1
        rhoz(i,j,k) = rhoztmp(k)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    
    NKloc = Nloc*kbbm1

    ! put the data in master processor into the global var                                                           
    if(myid==0) then
      do k = 1,kbbm1
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        rhozglob(iglob,jglob,k) = rhoz(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(rhozloc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,kbbm1
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      rhozloc(i,jk) = rhoz(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor                                                   
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = 1,kbbm1
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          jk = (k-1)*Nloc+j
          rhozglob(iglob,jglob,k) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(rhozloc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    deallocate(rhozloc)
    deallocate(xx)

    if(myid==0) then
      rhomean = zero
      do k = 1,kbbm1
        isum = 0
        do j = 1,Nglob
        do i = 1,Mglob
          if(-HCG(i,j)<=phy_z(k)) then
            isum = isum+1
            rhomean(k) = rhomean(k)+rhozglob(i,j,k)
          endif
        enddo
        enddo
        if(isum>=1) then
          rhomean(k) = rhomean(k)/float(isum)
        else
          rhomean(k) = rhomean(k-1)
        endif
      enddo
    endif

    call MPI_BCAST(rhomean,kbbm1,MPI_SP,0,MPI_COMM_WORLD,ier)

# else

    rhomean = zero
    do k = 1,kbbm1
      isum = 0
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        if(-Hc(i,j)<=phy_z(k)) then
          isum = isum+1
          rhomean(k) = rhomean(k)+rhoz(i,j,k)
        endif
      enddo
      enddo
      if(isum>=1) then
        rhomean(k) = rhomean(k)/float(isum)
      else
        rhomean(k) = rhomean(k-1)
      endif
    enddo

# endif

    ! linearly interpolate to obtain density at signa levels
    Rmean = Rho0
    do i = 1,Mloc
    do j = 1,Nloc
      if(Mask(i,j)==0) cycle

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)
      enddo

      call sinter(phy_z,rhomean,zm,rhos,kbbm1,Kloc)          

      Rmean(i,j,1:Kloc) = rhos
    enddo
    enddo
  
    return  
    end subroutine eval_dens

    subroutine sinter(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                             
!  array b.
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF (Y(I)<X(1 )) B(I) = A(1)
      IF (Y(I)>X(M1)) B(I) = A(M1)
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) * (X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter


    subroutine sinter_p(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  for baroclinic interpolation                                               
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                   
!  array b.                                                                    
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF(Y(I) < X(1 )) B(I) = A(1)-(A(2)-A(1))*(X(1)-Y(I))/(X(2)-X(1))
      IF(Y(I) > X(M1)) B(I) = A(M1)+(A(M1)-A(M1-1))*(Y(I)-X(M1))/(X(M1)-X(M1-1))                            
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) *(X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter_p


    subroutine baropg_z
!------------------------------------------------------------------------------
!   Calculate baroclinic terms in z levels
!   Called by
!      eval_duvw
!   Last Update: Gangfeng Ma, 07/05/2012
!-----------------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 201
    integer, parameter :: kbbm1 = kbb-1
    real(SP), dimension(Kloc) :: zm,rhos,pbxs,pbys
    real(SP), dimension(kbb) :: phy_z,rhoztmp,pbx,pby
    real(SP), dimension(3,3,kbb) ::pb,rhoz
    real(SP) :: Ramp1,hmax,etamax,tmp1,tmp2,myvar,deltz, &
                Rmean1,Rmean2,dz
    integer :: i,j,k,i2,j2,ic,jc

    if(TRamp==Zero) then
      Ramp1 = 1.0
    else
      Ramp1 = tanh(TIME/TRamp)
    endif

    ! subtract reference density
    Rho = Rho-Rho0

    DRhoX = Zero; DRhoY = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask9(i,j)==0) cycle

      ! find local maximum water depth 
      hmax = -large
      etamax = -large
      do j2 = j-1,j+1
      do i2 = i-1,i+1
        if(hc(i2,j2)>hmax) hmax = Hc(i2,j2)
        if(eta(i2,j2)>etamax) etamax = Eta(i2,j2)
      enddo
      enddo

      ! interpolate into physical z levels
      if(Ivgrd==1) then  
        deltz = (hmax+etamax)/float(kbbm1)
        do k = 1,kbb
          phy_z(k) = (float(k)-1.0)*deltz-hmax
        enddo
      else
        deltz = (hmax+etamax)*(Grd_R-1.0)/(Grd_R**float(kbbm1)-1.0)
        phy_z(1) = -hmax
        do k = 2,kbb
          phy_z(k) = phy_z(k-1)+deltz
          deltz = deltz*Grd_R
        enddo
      endif

      rhoz = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2  ! local index
        jc = j2-j+2  
        do k = 1,Kloc
          zm(k) = sigc(k)*D(i2,j2)-Hc(i2,j2)
          rhos(k) = Rho(i2,j2,k)
        enddo

        call sinter(zm,rhos,phy_z,rhoztmp,Kloc,kbb)

        do k = 1,kbb
          rhoz(ic,jc,k) = rhoztmp(k)
        enddo
      enddo
      enddo

      pb = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2
        jc = j2-j+2
        do k = kbbm1,1,-1
          if(phy_z(k)>=Eta(i2,j2)) then
            pb(ic,jc,k) = 0.0
          else
            if(phy_z(k)>=-Hc(i2,j2)) then
              dz = dmin1(phy_z(k+1)-phy_z(k),Eta(i2,j2)-phy_z(k))
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            elseif(phy_z(k)<-Hc(i2,j2).and.phy_z(k+1)>-Hc(i2,j2)) then
              dz = -Hc(i2,j2)-Phy_z(k)
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            else
              pb(ic,jc,k) = pb(ic,jc,k+1)
            endif
          endif
        enddo
      enddo
      enddo

      pbx = Zero; pby = Zero
      do k = 1,kbb
        if(phy_z(k)<=Eta(i,j).and.phy_z(k)>=-Hc(i,j)) then
          if(phy_z(k)<-Hc(i-1,j).and.phy_z(k)>=-Hc(i+1,j)) then
            pbx(k) = (pb(3,2,k)-pb(2,2,k))/dx
          elseif(phy_z(k)>=-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = (pb(2,2,k)-pb(1,2,k))/dx
          elseif(phy_z(k)<-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = Zero
          else
            pbx(k) = (pb(3,2,k)-pb(1,2,k))/(2.0*dx)
          endif

          if(phy_z(k)<-Hc(i,j-1).and.phy_z(k)>=-Hc(i,j+1)) then
            pby(k) = (pb(2,3,k)-pb(2,2,k))/dy
          elseif(phy_z(k)>=-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = (pb(2,2,k)-pb(2,1,k))/dy
          elseif(phy_z(k)<-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = Zero
          else
            pby(k) = (pb(2,3,k)-pb(2,1,k))/(2.0*dy)
          endif
        endif
      enddo

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)
      enddo
 
      call sinter_p(phy_z,pbx,zm,pbxs,kbb,Kloc)
      call sinter_p(phy_z,pby,zm,pbys,kbb,Kloc)

      do k = Kbeg,Kend
        DRhoX(i,j,k) = -pbxs(k)*grav*D(i,j)/Rho0*Ramp1
        DRhoY(i,j,k) = -pbys(k)*grav*D(i,j)/Rho0*Ramp1
      enddo
    enddo
    enddo

    ! Add back reference density
    Rho = Rho+Rho0

    return
    end subroutine baropg_z


    subroutine trig(alpha,beta,gama,b,x,N)
!*************************************************************!
!*                                                           *!
!*         (B1 C1                         )   (x1)     (b1)  *!
!*         (A2 B2 C2                      )   (x2)     (b2)  *!
!*         (   A3 B3 C3                   )   (x3)     (b3)  *!
!*         (      A4 B4 C4                )   (x4)====       *!
!*         (         A5 B5  C5            )   ... ====  ...  *!
!*         (            ... ... ...       )   ...       ...  *!
!*         (                An-1 Bn-1 Cn-1)   (xn-1)   (bn-1)*!
!*         (                     An   Bn  )   (xn)     (bn)  *!
!*                                                           *!
!*                                                           *!
!*************************************************************!
! where A are alpha, B are beta, C are gama
!-----------------------------------------------------------------------------
    use global, only: SP
    implicit none

    integer, intent(in) :: N
    real(SP), dimension(N), intent(in)  :: alpha,beta,gama,b
    real(SP), dimension(N), intent(out) :: x
    real(SP), dimension(N) :: betaPrime,bPrime
    real(SP) :: coeff
    integer :: II
 
    ! Perform forward elimination
    betaPrime(1) = beta(1)
    bPrime(1) = b(1)
 
    do II = 2,N
      coeff = alpha(II)/betaPrime(II-1)
      betaPrime(II) = beta(II)-coeff*gama(II-1)
      bPrime(II) = b(II)-coeff*bPrime(II-1)
    enddo

    ! Perform back substitution
    x(N) = bPrime(N) / betaPrime(N)
    do II = N-1,1,-1
      x(II) = (bPrime(II)-gama(II)*x(II+1))/betaPrime(II)
    enddo

    end subroutine trig


# if defined (VEGETATION)
    subroutine veg_height(load)
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc), intent(in) :: load
    integer :: i,j

    ! estimate height of flexible vegetation (Li and Xie, AWR, 2011)
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(load(i,j)<=100.) then
        FVegH(i,j) = -1.1093e-12*load(i,j)**6-7.9352e-12*load(i,j)**5  &
                     +7.3301e-8*load(i,j)**4-1.2141e-5*load(i,j)**3  &
                     +8.5414e-4*load(i,j)**2-3.171e-2*load(i,j)+1.0143
      else
        FVegH(i,j) = 5.4583e-18*load(i,j)**6-1.9598e-14*load(i,j)**5  &
                     +2.8577e-11*load(i,j)**4-2.1897e-8*load(i,j)**3  &
                     +9.5802e-6*load(i,j)**2-2.5010e-3*load(i,j)+0.55942
      endif
    enddo
    enddo
    FVegH = FVegH*VegH

    return
    end subroutine veg_height

    subroutine foli(load)
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc), intent(in) :: load
    integer :: i,j

    ! estimate foliage effects of flexible vegetation (Li and Xie, AWR, 2011)    
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(load(i,j)<=100.) then
        Foliage(i,j) = 8.6842e-12*load(i,j)**6-3.5010e-9*load(i,j)**5  &
                     +5.6891e-7*load(i,j)**4-4.7779e-5*load(i,j)**3  &
                     +2.2195e-3*load(i,j)**2-5.7707e-2*load(i,j)+1.0255
      else
        Foliage(i,j) = 9.0840e-13*load(i,j)**4-2.3969e-9*load(i,j)**3  &
                     +2.3272e-6*load(i,j)**2-1.056e-3*load(i,j)+0.31984 
      endif
    enddo
    enddo

    return
    end subroutine foli
# endif

# if defined (OBSTACLE)
    subroutine set_obsvel
!---------------------------------------------------------------
!   Specify or calculate obstacle velocities
!   by Gangfeng Ma, 17/08/2013
!---------------------------------------------------------------
    use global, only: Zero,obs_u,obs_v,obs_w
    implicit none

    obs_u = Zero
    obs_v = Zero
    obs_w = Zero

    return
    end subroutine set_obsvel

    subroutine set_obsflag
!----------------------------------------------------------------
!   Determine obstacle flag
!   set_flag = 1: fluid
!   set_flag = 0: obstacle
!----------------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: zc,xlin,ylin,zlin,dista,ugrad,vgrad,wgrad

    ! save flag in the previous step
    set_flag_old = set_flag

    ! default is fluid
    set_flag = 0

    ! specify obstacle 
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      zc = D(i,j)*sigc(k)

      ! obstacle 1
      if(xc(i)>=6.0.and.xc(i)<=12.0) then
        if(zc<0.05*(xc(i)-6.0)) set_flag(i,j,k) = 1
      endif

      ! obstacle 2
      if(xc(i)>=12.0.and.xc(i)<=14.0) then
        if(zc<0.3) set_flag(i,j,k) = 1
      endif

      ! obstacle 3
      if(xc(i)>=14.0.and.xc(i)<=17.0) then
        if(zc<0.1*(17.0-xc(i))) set_flag(i,j,k) = 1
      endif
    enddo
    enddo
    enddo

    ! distance from obstacle/fluid interface to neighboring fluid cell 
    set_dist_x = 1.e+20
    set_dist_y = 1.e+20
    set_dist_z = 1.e+20

    ! find the distance
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend

      ! obstacle 1
      if(xc(i)>=6.0.and.xc(i)<=12.0) then
        ! loop over all fluid cells
        if(set_flag(i,j,k)==0) then
          ! obstacle on the right
          if(set_flag(i+1,j,k)==1) then
            zc = D(i,j)*sigc(k)   
            xlin = 20.0*zc+6.0
            set_dist_x(i,j,k) = abs(xlin-xc(i))
          endif

          ! obstacle under the fluid
          if(set_flag(i,j,k-1)==1) then           
            zc = D(i,j)*sigc(k)
            zlin = 0.05*(xc(i)-6.0)
            set_dist_z(i,j,k) = abs(zc-zlin)
          endif
        endif  
      endif

      ! obstacle 2
      if(xc(i)>=12.0.and.xc(i)<=14.0) then
        ! loop over all fluid cells
        if(set_flag(i,j,k)==0) then
          ! obstacle under the fluid
          if(set_flag(i,j,k-1)==1) then
            zc = D(i,j)*sigc(k)
            set_dist_z(i,j,k) = abs(zc-0.3)
          endif
        endif
      endif

      ! obstacle 3
      if(xc(i)>=14.0.and.xc(i)<=17.0) then
        ! loop over all fluid cells
        if(set_flag(i,j,k)==0) then
          ! obstacle on the left
          if(set_flag(i-1,j,k)==1) then
            zc = D(i,j)*sigc(k)
            xlin = 17.0-10.0*zc
            set_dist_x(i,j,k) = abs(xc(i)-xlin)
          endif

          ! obstacle under the fluid
          if(set_flag(i,j,k-1)==1) then
            zc = D(i,j)*sigc(k)
            zlin = 0.1*(17.0-xc(i))
            set_dist_z(i,j,k) = abs(zc-zlin)
          endif
        endif
      endif
    enddo
    enddo
    enddo

    ! If the grid point is located inside the obstacle in the previous
    ! step and moves outside the obstacle in the current step.
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(set_flag_old(i,j,k)==1.and.set_flag(i,j,k)==0) then
        if(set_flag(i,j,k+1)==0) then
          dista = set_dist_z(i,j,k)
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          ugrad = (U(i,j,k+1)-0.0)/dista
          vgrad = (V(i,j,k+1)-0.0)/dista
          wgrad = (W(i,j,k+1)-0.0)/dista
          
          dista = set_dist_z(i,j,k)
          U(i,j,k) = ugrad*dista
          V(i,j,k) = vgrad*dista
          W(i,j,k) = wgrad*dista

          DU(i,j,k) = U(i,j,k)*D(i,j)
          DV(i,j,k) = V(i,j,k)*D(i,j)
          DW(i,j,k) = W(i,j,k)*D(i,j)
        elseif(set_flag(i,j,k-1)==0) then
          dista = set_dist_z(i,j,k)
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          ugrad= (U(i,j,k-1)-0.0)/dista
          vgrad= (V(i,j,k-1)-0.0)/dista
          wgrad = (W(i,j,k-1)-0.0)/dista
            
          dista= set_dist_z(i,j,k)
          U(i,j,k) = ugrad*dista
          V(i,j,k) = vgrad*dista
          W(i,j,k) = wgrad*dista

          DU(i,j,k) = U(i,j,k)*D(i,j)
          DV(i,j,k) = V(i,j,k)*D(i,j)
          DW(i,j,k) = W(i,j,k)*D(i,j)
        endif
      endif
    enddo
    enddo
    enddo

    return
    end subroutine set_obsflag

    subroutine imm_obs
!---------------------------------------------------------------------
!   Immersed-Boundary Method
!   Calculate forcing at the obstacle boundary
!   By Gangfeng Ma, 17/08/2013
!--------------------------------------------------------------------
    use global
    implicit none
    real(SP),dimension(:),allocatable :: target_vel,target_dist
    integer :: i,j,k,l,count
    real(SP) :: numer,denom,force_vel,dista,vel_grad

    ! calculate forcing at boundaries
    ObsForceX = Zero
    ObsForceY = Zero
    ObsForceZ = Zero

    allocate(target_vel(1:6))
    allocate(target_dist(1:6))

    ! x-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (U(i+1,j,k)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_u-U(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (U(i,j+1,k)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_u-U(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (U(i,j,k+1)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_u-U(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 1
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 1        continue
          ObsForceX(i,j,k) = (D(i,j)*force_vel-DU(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo

    ! y-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (V(i+1,j,k)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_v-V(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (V(i,j+1,k)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_v-V(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (V(i,j,k+1)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_v-V(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 2
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 2        continue
          ObsForceY(i,j,k) = (D(i,j)*force_vel-DV(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo
 
    ! z-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (W(i+1,j,k)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_w-W(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (W(i,j+1,k)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_w-W(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (W(i,j,k+1)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_w-W(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 3
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 3        continue
          ObsForceZ(i,j,k) = (D(i,j)*force_vel-DW(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo

    deallocate(target_vel)
    deallocate(target_dist)
 1000 continue

    return
    end subroutine imm_obs
# endif

# if defined (POROUSMEDIA)
    subroutine read_porosity
    use global
    implicit none
    integer  :: i,j,k
    real(SP) :: Vis_Por,zlev

    Vis_Por = 1.e-6

    Porosity = 1.0
    Ap_Por = 0.0
    Bp_Por = 0.0
    Cp_Por = 0.0

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      zlev = sigc(k)*D(i,j)
      if(xc(i)>=Por_X0.and.xc(i)<=Por_Xn.and.yc(j)>=Por_Y0.and.yc(j)<=Por_Yn  &
           .and.zlev>=Por_Z0.and.zlev<=Por_Zn) then
        Porosity(i,j,k) = Por_n
        Ap_Por(i,j,k) = alpha_por*(1-Porosity(i,j,k))**2/Porosity(i,j,k)**2*Vis_Por/D50_Por
        Bp_Por(i,j,k) = beta_por*(1-Porosity(i,j,k))/Porosity(i,j,k)**2/D50_Por
        Cp_Por(i,j,k) = 0.34*(1-Porosity(i,j,k))/Porosity(i,j,k)
      endif
    enddo
    enddo
    enddo

    return
    end subroutine read_porosity
# endif

!---------------------------------------------------------------------
!   Record the max runup of wave
!   By Fengyan Shi (Cheng Zhang), 20/11/2016
!---------------------------------------------------------------------
SUBROUTINE max_min_property
     use global
     implicit none
	integer :: i,j
 
      DO j = Jbeg,Jend
      DO i = Ibeg,Iend
       IF(OUT_M)THEN
        IF(MASK(i,j).GT.0)THEN
        IF(Eta(i,j).GT.HeightMax(i,j)) HeightMax(i,j)=Eta(i,j)
        ENDIF
       ENDIF
	   
      ENDDO
      ENDDO
 
END SUBROUTINE max_min_property

!---------------------------------------------------------------------
!   Hot start function for the simulation
!   including subroutine read_2d and subroutine read_3d
!   By Fengyan Shi (Cheng Zhang), 21/11/2016
!---------------------------------------------------------------------
SUBROUTINE hot_start
      use global
     implicit none
     integer :: i,j,k,Iglob,Jglob
     character(LEN=80) :: file=''
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: TMP_READ
     REAL(SP),DIMENSION(Mloc,Nloc) :: R1
     REAL(SP) :: TIME0,Cmiu
	 
    if (RNG) then
       cmiu = 0.085
    else
       cmiu = 0.09
    endif

     open(5,file='time0.dat')
     read(5,*) TIME0,Icount
     close(5)

     TIME=TIME0

     file=TRIM(Eta_HotStart_File)
     call read_2d(Eta,file)
     
     file=TRIM(U_HotStart_File)
     call read_3d(u,file)

     file=TRIM(V_HotStart_File)
     call read_3d(v,file)

     file=TRIM(W_HotStart_File)
     call read_3d(w,file)

     file=TRIM(P_HotStart_File)
     call read_3d(p,file)


    IF(VISCOUS_FLOW)THEN
!     file=TRIM(Rho_HotStart_File)
!     call read_3d(Rho,file)

     file=TRIM(TKE_HotStart_File)
     call read_3d(Tke,file)

     file=TRIM(EPS_HotStart_File)
     call read_3d(Eps,file)
    ENDIF

# if defined (SALINITY)
     file=TRIM(Sali_HotStart_File)
     call read_3d(Sali,file)

     file=TRIM(Temp_HotStart_File)
     call read_3d(Temp,file)

# endif

# if defined (LANDSLIDE) || defined (FLUIDSLIDE) || defined (TWOLAYERSLIDE) || defined (LANDSLIDE_COMPREHENSIVE)
     file=TRIM(Depth_HotStart_File)
     call read_2d(Hc,file)
# endif

! prepare all initials
# if defined (LANDSLIDE) || defined (FLUIDSLIDE) || defined (TWOLAYERSLIDE) || defined (LANDSLIDE_COMPREHENSIVE)
! Hc
     call phi_2D_coll(Hc)
# endif
! eta
     ! collect data into ghost cells
     call phi_2D_coll(Eta)
     Eta0 = Eta
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     D = max(Hc+Eta, MinDep)

! uvw


     call vel_bc(1)
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
# endif

# if defined (SALINITY)
! sali
     call sali_bc(1)  
# endif  

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
# if defined (SALINITY)
       DSali(I,J,K)=D(I,J)*Sali(I,J,K)
       DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif
     enddo
     enddo
     enddo

  IF(VISCOUS_FLOW)THEN
     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       DTke(I,J,K)=D(I,J)*Tke(I,J,K)
       DEps(I,J,K)=D(I,J)*Eps(I,J,K)
     ENDDO
     ENDDO
     ENDDO

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=MinDep.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    CmuHt = CmuVt   ! for IHturb>=10
   
   ELSE
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc       
       CmuVt(I,J,K) = Cmut_min
       CmuHt(I,J,K) = Cmut_min
    enddo
    enddo
    enddo

   ENDIF ! end viscous_flow


#if defined (SALINITY)
     call eval_dens
# endif

END SUBROUTINE hot_start

SUBROUTINE read_2d(phi,filename)
      use global
     implicit none
     integer :: i,j,k, Jglob,Iglob
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: TMP_READ
     REAL(SP),INTENT(OUT) :: PHI(Mloc,Nloc)
     character(len=80),INTENT(in) :: filename

     ALLOCATE (TMP_READ(Mglob,Nglob))

    
     OPEN(5,FILE=trim(filename),STATUS='OLD')
        DO J=1,Nglob
          READ(5,*)(TMP_READ(I,J),I=1,Mglob)
        ENDDO
     CLOSE(5)

# if defined (PARALLEL)
 
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = npx*(Mloc-2*Nghost)+i-Nghost
           Jglob = npy*(Nloc-2*Nghost)+j-Nghost
         phi(I,J)=TMP_READ(Iglob,Jglob)
       ENDDO
       ENDDO
# else

       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = i-Nghost
           Jglob = j-Nghost
           phi(I,J)=TMP_READ(Iglob,Jglob)
       ENDDO
       ENDDO
  
 
# endif

 
     DEALLOCATE (TMP_READ)
END SUBROUTINE read_2d

SUBROUTINE read_3d(phi,filename)
      use global
     implicit none
     integer :: i,j,k,Jglob,Iglob
     REAL(SP),INTENT(OUT) :: PHI(Mloc,Nloc,Kloc)
     character(len=80),INTENT(in) :: filename
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: TMP_READ
  
     ALLOCATE (TMP_READ(Mglob,Nglob,Kglob))
     
    
     OPEN(5,FILE=trim(filename),STATUS='OLD')
       DO K=1,Kglob
        DO J=1,Nglob
          READ(5,*)(TMP_READ(I,J,K),I=1,Mglob)
        ENDDO
       ENDDO
     CLOSE(5)

# if defined (PARALLEL)
       DO K=Kbeg,Kend
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = npx*(Mloc-2*Nghost)+i-Nghost
           Jglob = npy*(Nloc-2*Nghost)+j-Nghost
         phi(I,J,K)=TMP_READ(Iglob,Jglob,K-Nghost)
       ENDDO
       ENDDO
      
     ENDDO
# else
       DO K=Kbeg,Kend
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = i-Nghost
           Jglob = j-Nghost
           phi(I,J,K)=TMP_READ(Iglob,Jglob,K-Nghost)
       ENDDO
       ENDDO
  
     ENDDO
# endif


     DEALLOCATE (TMP_READ)
END SUBROUTINE read_3d

! added by Cheng for nesting (following two subroutines)
# if defined (COUPLING)
!---------------------------------------------------
!    This subroutine is used to read nesting file at the first time
!    Called by 
!       main
!    Last update: 05/15/2013, fyshi
!---------------------------------------------------
SUBROUTINE READ_NESTING_FILE
     use global
     use input_util
     implicit none
     INTEGER :: I,J,K
       OPEN(11,FILE=TRIM(COUPLING_FILE))
         READ(11,*)  ! title
         READ(11,*)  ! boundary info
! boundary basic info including point number of coupling, start point, etc
! east
         READ(11,*)  ! east
         READ(11,*) N_COUPLING_EAST,J_START_EAST
! west 
         READ(11,*)  ! west
         READ(11,*) N_COUPLING_WEST,J_START_WEST
! south 
         READ(11,*)  ! south
         READ(11,*) N_COUPLING_SOUTH,I_START_SOUTH
! north 
         READ(11,*)  ! north
         READ(11,*) N_COUPLING_NORTH,I_START_NORTH

! read time and variable at the first level

         READ(11,*) ! time start title
         READ(11,*) TIME_COUPLING_1 
! initialize time_2
         TIME_COUPLING_2 = TIME_COUPLING_1

! east
         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               V_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               W_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               Z_COUPLING_EAST(N_COUPLING_EAST,2), &
               P_COUPLING_EAST(N_COUPLING_EAST,Kglob,2))
# if defined (SALINITY)
           ALLOCATE(S_COUPLING_EAST(N_COUPLING_EAST,Kglob,2))
# endif	
# if defined (TEMPERATURE)	
           ALLOCATE(T_COUPLING_EAST(N_COUPLING_EAST,Kglob,2))
# endif
               
             READ(11,*)   ! east
             READ(11,119)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)((U_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((V_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((W_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((P_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# if defined (SALINITY)
             READ(11,119)((S_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
             READ(11,119)((T_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# endif	
!   initialize first step
             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)
# endif
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_east
119      FORMAT(5E16.6)

! west
         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               V_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               W_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               Z_COUPLING_WEST(N_COUPLING_WEST,2),&
               P_COUPLING_WEST(N_COUPLING_WEST,Kglob,2))
# if defined (SALINITY)
           ALLOCATE(S_COUPLING_WEST(N_COUPLING_WEST,Kglob,2))
# endif	
# if defined (TEMPERATURE)
           ALLOCATE(T_COUPLING_WEST(N_COUPLING_WEST,Kglob,2))
# endif

             READ(11,*)   ! west
             READ(11,119)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)((U_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((V_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((W_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((P_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# if defined (SALINITY)
             READ(11,119)((S_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
             READ(11,119)((T_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# endif
!   initialize first step
             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)
# endif	
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_west
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               W_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               P_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2))
# if defined (SALINITY)
           ALLOCATE(S_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2))
# endif	
# if defined (TEMPERATURE)
           ALLOCATE(T_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2))
# endif	
             READ(11,*)   ! south
             READ(11,119)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)((U_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((V_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((W_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((P_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# if defined (SALINITY)
             READ(11,119)((S_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
             READ(11,119)((T_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# endif
!   initialize first step
             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
# if defined (SALINITY)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)
# endif
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_south
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               W_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,2), &
               P_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2))
# if defined (SALINITY)
           ALLOCATE(S_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2))
# endif	
# if defined (TEMPERATURE)
           ALLOCATE(T_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2))
# endif
             READ(11,*)   ! north
             READ(11,119)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)((U_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((V_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((W_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((P_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# if defined (SALINITY)
             READ(11,119)((S_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
             READ(11,119)((T_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# endif	
!   initialize first step
             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
# if defined (SALINITY)
             S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)
# endif	
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_north


! specify boundary start points

! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
      Kstart_WEST=J_START_WEST+Nghost -npy*Nglob/py
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1 -npy*Nglob/py
      IF((Kstart_WEST>Nghost.AND.Kstart_WEST<Nloc-Nghost+1).OR.&
     (Kend_WEST<Nloc-Nghost+1.AND.Kend_WEST>Nghost))THEN
       IF(Kstart_WEST<Nghost+1)THEN
         Kshift_WEST= -Kstart_WEST+Nghost+1
         Kstart_WEST=Nghost+1
       ELSE
         Kshift_WEST=-(Kstart_WEST-Nghost)+1
       ENDIF
       IF(Kend_WEST>Nloc-Nghost)THEN
         Kend_WEST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_WEST=.TRUE.
      ELSE
       IF(Kstart_WEST<=Nghost.AND.Kend_WEST>=Nloc-Nghost+1)THEN
         Kshift_WEST=-Kstart_WEST+Nghost+1
         Kstart_WEST = Nghost+1
         Kend_WEST = Nloc-Nghost
         IN_DOMAIN_WEST=.TRUE.
       ELSE
         IN_DOMAIN_WEST=.FALSE.
       ENDIF
      ENDIF
! check print*,myid,Kshift_WEST,Kstart_WEST,Kend_WEST
     endif
# else
      Kstart_WEST=J_START_WEST+Nghost
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1
      Kshift_WEST = -(Kstart_WEST-Nghost)+1
      IN_DOMAIN_WEST = .TRUE.
# endif

   ENDIF

! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
      Kstart_EAST=J_START_EAST+Nghost -npy*Nglob/py
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1 -npy*Nglob/py
      IF((Kstart_EAST>Nghost.AND.Kstart_EAST<Nloc-Nghost+1).OR.&
     (Kend_EAST<Nloc-Nghost+1.AND.Kend_EAST>Nghost))THEN
       IF(Kstart_EAST<Nghost+1)THEN
         Kshift_EAST= -Kstart_EAST+Nghost+1
         Kstart_EAST=Nghost+1
       ELSE
         Kshift_EAST=-(Kstart_EAST-Nghost)+1
       ENDIF
       IF(Kend_EAST>Nloc-Nghost)THEN
         Kend_EAST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_EAST=.TRUE.
      ELSE

       IF(Kstart_EAST<=Nghost.AND.Kend_EAST>=Nloc-Nghost+1)THEN
         Kshift_EAST=-Kstart_EAST+Nghost+1
         Kstart_EAST = Nghost+1
         Kend_EAST = Nloc-Nghost
         IN_DOMAIN_EAST=.TRUE.
       ELSE
         IN_DOMAIN_EAST=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_EAST=J_START_EAST+Nghost
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1
      Kshift_EAST = -(Kstart_EAST-Nghost)+1
      IN_DOMAIN_EAST = .TRUE.
# endif
    ENDIF

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
      Kstart_SOUTH=I_START_SOUTH+Nghost -npx*Mglob/px
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1 -npx*Mglob/px
      IF((Kstart_SOUTH>Nghost.AND.Kstart_SOUTH<Mloc-Nghost+1).OR.&
     (Kend_SOUTH<Mloc-Nghost+1.AND.Kend_SOUTH>Nghost))THEN
       IF(Kstart_SOUTH<Nghost+1)THEN
         Kshift_SOUTH= -Kstart_SOUTH+Nghost+1
         Kstart_SOUTH=Nghost+1
       ELSE
         Kshift_SOUTH=-(Kstart_SOUTH-Nghost)+1
       ENDIF
       IF(Kend_SOUTH>Mloc-Nghost)THEN
         Kend_SOUTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_SOUTH=.TRUE.
      ELSE

       IF(Kstart_SOUTH<=Nghost.AND.Kend_SOUTH>=Mloc-Nghost+1)THEN
         Kshift_SOUTH=-Kstart_SOUTH+Nghost+1
         Kstart_SOUTH = Nghost+1
         Kend_SOUTH = Mloc-Nghost
         IN_DOMAIN_SOUTH=.TRUE.
       ELSE
         IN_DOMAIN_SOUTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_SOUTH=I_START_SOUTH+Nghost
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1
      Kshift_SOUTH = -(Kstart_SOUTH-Nghost)+1
      IN_DOMAIN_SOUTH = .TRUE.
# endif
   ENDIF

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
      Kstart_NORTH=I_START_NORTH+Nghost -npx*Mglob/px
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1 -npx*Mglob/px
      IF((Kstart_NORTH>Nghost.AND.Kstart_NORTH<Mloc-Nghost+1).OR.&
     (Kend_NORTH<Mloc-Nghost+1.AND.Kend_NORTH>Nghost))THEN
       IF(Kstart_NORTH<Nghost+1)THEN
         Kshift_NORTH= -Kstart_NORTH+Nghost+1
         Kstart_NORTH=Nghost+1
       ELSE
         Kshift_NORTH=-(Kstart_NORTH-Nghost)+1
       ENDIF
       IF(Kend_NORTH>Mloc-Nghost)THEN
         Kend_NORTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_NORTH=.TRUE.
      ELSE

       IF(Kstart_NORTH<=Nghost.AND.Kend_NORTH>=Mloc-Nghost+1)THEN
         Kshift_NORTH=-Kstart_NORTH+Nghost+1
         Kstart_NORTH = Nghost+1
         Kend_NORTH = Mloc-Nghost
         IN_DOMAIN_NORTH=.TRUE.
       ELSE
         IN_DOMAIN_NORTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_NORTH=I_START_NORTH+Nghost
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1
      Kshift_NORTH = -(Kstart_NORTH-Nghost)+1
      IN_DOMAIN_NORTH = .TRUE.
# endif
   ENDIF

END SUBROUTINE READ_NESTING_FILE
# endif

# if defined (COUPLING)
!-------------------------------------------------------------------
!   This subroutine is used to pass coupling variables into ghost cells                                                         
!   Called by
!      MAIN
!   Update: 05/15/2013 Fengyan Shi, University of Delaware                                       
!-------------------------------------------------------------------
SUBROUTINE OneWayCoupling
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: I,J,K
    REAL(SP) :: tmp1,tmp2
    LOGICAL, SAVE :: FirstReadCoupling = .TRUE.

119      FORMAT(5E16.6)  ! this is a fixed format for I/O
 
! determine time slot

    IF(TIME>TIME_COUPLING_1.AND.TIME>TIME_COUPLING_2) THEN
         TIME_COUPLING_1=TIME_COUPLING_2
         
         READ(11,*,END=120) ! time coupling
         READ(11,*,END=120) TIME_COUPLING_2 
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             READ(11,*,END=120)   ! east

             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)
# endif

             READ(11,119,END=120)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)((U_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# if defined (SALINITY)
             READ(11,119,END=120)((S_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
			 READ(11,119,END=120)((T_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
# endif
			 
            IF (FirstReadCoupling)THEN
             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)
# endif
            ENDIF

         ELSE
             READ(11,*,END=120)   ! east            
         ENDIF
! west
         IF(N_COUPLING_WEST.GT.0)THEN
             READ(11,*,END=120)   ! west

             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)
# endif

             READ(11,119,END=120)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)((U_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# if defined (SALINITY)
			 READ(11,119,END=120)((S_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
			 READ(11,119,END=120)((T_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
# endif	
			IF (FirstReadCoupling)THEN
             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
# if defined (SALINITY)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)
# endif
            ENDIF

         ELSE
             READ(11,*,END=120)   ! west            
         ENDIF
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             READ(11,*,END=120)   ! south

             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
# if defined (SALINITY)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)
# endif

             READ(11,119,END=120)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)((U_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# if defined (SALINITY)
			 READ(11,119,END=120)((S_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
			 READ(11,119,END=120)((T_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
# endif
			IF (FirstReadCoupling)THEN
             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
# if defined (SALINITY)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)
# endif
            ENDIF

         ELSE
             READ(11,*,END=120)   ! south            
         ENDIF
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             READ(11,*,END=120)   ! north

             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
# if defined (SALINITY)
             S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)
# endif	

             READ(11,119,END=120)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)((U_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# if defined (SALINITY)
			 READ(11,119,END=120)((S_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# endif	
# if defined (TEMPERATURE)
			 READ(11,119,END=120)((T_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
# endif	
		  IF (FirstReadCoupling)THEN
             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
# if defined (SALINITY)
			 S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
# endif	
# if defined (TEMPERATURE)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)
# endif	
           ENDIF

         ELSE
             READ(11,*,END=120)   ! north            
         ENDIF

         FirstReadCoupling = .FALSE.

    ENDIF  !stime>time_2 and time_1

120 CONTINUE

    tmp2=1.0
    tmp1=ZERO

    IF(TIME>TIME_COUPLING_1)THEN
      IF(TIME_COUPLING_1.EQ.TIME_COUPLING_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=1.0
      ELSE
      tmp2=(TIME_COUPLING_2-TIME) &
            /MAX(SMALL, ABS(TIME_COUPLING_2-TIME_COUPLING_1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1


! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_WEST)THEN

      DO J=Kstart_WEST,Kend_WEST 
      DO I=1,Nghost
        ETA(I,J)=Z_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +Z_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost
        U(I,J,K)=U_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +U_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +V_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2  
        W(I,J,K)=W_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +W_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2  
        P(I,J,K)=P_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +P_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2 
# if defined (SALINITY)
        Sali(I,J,K)=S_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +S_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2
# endif	
# if defined (TEMPERATURE)
        Temp(I,J,K)=T_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +T_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2
# endif	

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)   
# if defined (SALINITY)		
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K) 
# endif	
# if defined (TEMPERATURE)		
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K)   
# endif		

      ENDDO
      ENDDO
      ENDDO

!print*,'intp',Sali(1,Jbeg+1,Kend),Sali(2,Jbeg+1,Kend),Sali(3,Jbeg+1,Kend),Sali(4,Jbeg+1,Kend)
!J=Kstart_WEST
!print*,'swet',S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,2),S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,1)
!print*,'calc',S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,2)*tmp1+S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,1)*tmp2

     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_west>0


! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_EAST)THEN
      DO J=Kstart_EAST,Kend_EAST  
      DO I=Iend+1,Iend+Nghost
        ETA(I,J)=Z_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +Z_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost
        U(I,J,K)=U_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +U_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +V_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +W_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +P_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
# if defined (SALINITY)	
        Sali(I,J,K)=S_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +S_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
# endif	
# if defined (TEMPERATURE)
        Temp(I,J,K)=T_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +T_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
# endif	

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)   
# if defined (SALINITY)			
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
# endif	
# if defined (TEMPERATURE)
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
# endif	
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_east>0

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_SOUTH)THEN
      DO I=Kstart_SOUTH,Kend_SOUTH  
      DO J=1,Nghost
        ETA(I,J)=Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost    
        U(I,J,K)=U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +W_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +P_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
# if defined (SALINITY)	
        Sali(I,J,K)=S_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +S_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
# endif	
# if defined (TEMPERATURE)
        Temp(I,J,K)=T_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +T_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
# endif	

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)  
# if defined (SALINITY)		
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
# endif	
# if defined (TEMPERATURE)
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
# endif	
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_south>0

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_NORTH)THEN
      DO I=Kstart_NORTH,Kend_NORTH  
      DO J=Jend+1,Jend+Nghost
        ETA(I,J)=Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost       
        U(I,J,K)=U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +W_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +P_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
# if defined (SALINITY)
        Sali(I,J,K)=S_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +S_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
# endif	
# if defined (TEMPERATURE)
        Temp(I,J,K)=T_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +T_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
# endif	

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)  
# if defined (SALINITY)		
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
# endif	
# if defined (TEMPERATURE)
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
# endif
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_north>0

END SUBROUTINE OneWayCoupling
# endif 
! end nesting

!------------------------The End----------------------------------------------
