! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
    subroutine update_maskvs
!------------------------------------------------------  
!    This subroutine is used to update maskvs for wetting-drying
!    Called by                                                
!       update_bathymetry
!    Last update: 26/03/2017, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	integer,dimension(:,:),allocatable :: Masktmp
	 
    allocate(Masktmp(Mloc,Nloc))
    Masktmp = Maskvs
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maskvs(i,j)==0) then
        if(Maskvs(i-1,j)==1.and.Hvs(i-1,j)<Hvs(i,j)) Masktmp(i,j)=1
        if(Maskvs(i+1,j)==1.and.Hvs(i+1,j)<Hvs(i,j)) Masktmp(i,j)=1
        if(Maskvs(i,j-1)==1.and.Hvs(i,j-1)<Hvs(i,j)) Masktmp(i,j)=1
        if(Maskvs(i,j+1)==1.and.Hvs(i,j+1)<Hvs(i,j)) Masktmp(i,j)=1
      else
        ! drying (wet->dry)               
        if(Dvs(i,j)-SLIDE_MINTHICK<=1.e-8) then
          Masktmp(i,j) = 0
          Dvs(i,j) = SLIDE_MINTHICK
        endif
      endif
	  Hvs(i,j) = Hc0(i,j)-Dvs(i,j)
    enddo
    enddo
	Maskvs = Masktmp

# if defined (PARALLEL)
     ! collect mask into ghost cells  
     call phi_int_exch(Maskvs)
# endif
	 
	 deallocate(Masktmp)

    end subroutine update_maskvs
	
	
     subroutine update_vars_vs
!------------------------------------------------------ 
!    This subroutine is used to save fluid slide variables at 
!    last time step
!    Called by   
!       update_bathymetry
!    Last update: 27/03/2017, Cheng Zhang 
!----------------------------------------------------- 
     use global        
     implicit none

     Dvs0 = Dvs
     Uvs0 = Uvs
     Vvs0 = Vvs
     DUvs0 = DUvs
     DVvs0 = DVvs

     end subroutine update_vars_vs
	 
	
     subroutine fluxes_vs
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces for fluid slide
!    Called by
!       update_bathymetry
!    Last update: 26/03/2017, Cheng Zhang
!------------------------------------------------
     use global
     implicit none

     ! second order construction
     call delxyzFun_vs
     call construction_vs

     ! calculate wave speed
     call wave_speed_vs

     ! calculate fluxes at faces
     call fluxes_at_faces_HLL_vs

     ! impose boundary conditions
     call flux_bc_vs

     end subroutine fluxes_vs

     subroutine delxyzFun_vs
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes_vs 
!    Last update: 27/03/2017, Cheng Zhang
!------------------------------------------
     use global
     implicit none
     
	 call delxFun_vs_2D(Hvs,DelxHvs)
     call delxFun_vs_2D(Uvs,DelxUvs)
     call delxFun_vs_2D(Vvs,DelxVvs)
     call delxFun_vs_2D(DUvs,DelxDUvs)
     call delxFun_vs_2D(DVvs,DelxDVvs)

	 call delyFun_vs_2D(Hvs,DelyHvs)
     call delyFun_vs_2D(Uvs,DelyUvs)
     call delyFun_vs_2D(Vvs,DelyVvs)
     call delyFun_vs_2D(DUvs,DelyDUvs)
     call delyFun_vs_2D(DVvs,DelyDVvs)

     end subroutine delxyzFun_vs
	 
     subroutine delxFun_vs_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun_vs
!    Last update: 27/03/2017, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Maskvs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Maskvs(i,j)==0) then
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
     end subroutine delxFun_vs_2D
	 
     subroutine delyFun_vs_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y
!    Called by 
!       delxyzFun_vs
!    Last update: 27/03/2017, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Maskvs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Maskvs(i,j)==0) then 
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
     end subroutine delyFun_vs_2D
	 
	 subroutine construction_vs
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes_vs
!    Last update: 27/03/2017, Cheng Zhang
!-----------------------------------------
     use global
     implicit none
     integer :: i,j

	 call construct_2D_x(Hvs,DelxHvs,HvsxL,HvsxR)
     call construct_2D_x(Uvs,DelxUvs,UvsxL,UvsxR)
     call construct_2D_x(Vvs,DelxVvs,VvsxL,VvsxR)
     call construct_2D_x(DUvs,DelxDUvs,DUvsxL,DUvsxR)
     call construct_2D_x(DVvs,DelxDVvs,DVvsxL,DVvsxR)

     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DvsxL(i,j) = H0fx(i,j)-HvsxL(i,j)
       DvsxR(i,j) = H0fx(i,j)-HvsxR(i,j)
       EvsxL(i,j) = DUvsxL(i,j)
       EvsxR(i,j) = DUvsxR(i,j)
       FvsxL(i,j) = 6.0/5.0*DUvsxL(i,j)*UvsxL(i,j)+0.5*Grav*(HvsxL(i,j)*HvsxL(i,j)-2.0*HvsxL(i,j)*H0fx(i,j))
       FvsxR(i,j) = 6.0/5.0*DUvsxR(i,j)*UvsxR(i,j)+0.5*Grav*(HvsxR(i,j)*HvsxR(i,j)-2.0*HvsxR(i,j)*H0fx(i,j))
       GvsxL(i,j) = 6.0/5.0*DUvsxL(i,j)*VvsxL(i,j)
       GvsxR(i,j) = 6.0/5.0*DUvsxR(i,j)*VvsxR(i,j)
     enddo
     enddo

	 call construct_2D_y(Hvs,DelyHvs,HvsyL,HvsyR)
     call construct_2D_y(Uvs,DelyUvs,UvsyL,UvsyR)
     call construct_2D_y(Vvs,DelyVvs,VvsyL,VvsyR)
     call construct_2D_y(DUvs,DelyDUvs,DUvsyL,DUvsyR)
     call construct_2D_y(DVvs,DelyDVvs,DVvsyL,DVvsyR)

     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DvsyL(i,j) = H0fy(i,j)-HvsyL(i,j)
       DvsyR(i,j) = H0fy(i,j)-HvsyR(i,j)
       EvsyL(i,j) = DVvsyL(i,j)
       EvsyR(i,j) = DVvsyR(i,j)
       FvsyL(i,j) = 6.0/5.0*DVvsyL(i,j)*UvsyL(i,j)
       FvsyR(i,j) = 6.0/5.0*DVvsyR(i,j)*UvsyR(i,j)
       GvsyL(i,j) = 6.0/5.0*DVvsyL(i,j)*VvsyL(i,j)+0.5*Grav*(HvsyL(i,j)*HvsyL(i,j)-2.0*HvsyL(i,j)*H0fy(i,j))
       GvsyR(i,j) = 6.0/5.0*DVvsyR(i,j)*VvsyR(i,j)+0.5*Grav*(HvsyR(i,j)*HvsyR(i,j)-2.0*HvsyR(i,j)*H0fy(i,j))
     enddo
     enddo

     end subroutine construction_vs
	 
     subroutine wave_speed_vs
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    of fluid slide
!    Called by
!       fluxes_vs
!    Last update: 27/03/2017, Cheng Zhang
!-----------------------------------------------
     use global
     implicit none
     integer  :: i,j
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S,UL,UR,VL,VR
     
     ! x-faces
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Maskvs(i-1,j)==1.and.Maskvs(i,j)==1) then
         SQR_PHI_L = sqrt(6.0/25.0*UvsxL(i,j)**2+Grav*abs(DvsxL(i,j)))
         SQR_PHI_R = sqrt(6.0/25.0*UvsxR(i,j)**2+Grav*abs(DvsxR(i,j)))
		 UL = 6.0/5.0*UvsxL(i,j)
		 UR = 6.0/5.0*UvsxR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL-UR)
         U_S = 0.5*(UL+UR)+SQR_PHI_L-SQR_PHI_R
         SvsxL(i,j) = min(UL-SQR_PHI_L,U_S-SQR_PHI_S)
         SvsxR(i,j) = max(UR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskvs(i-1,j)==0.and.Maskvs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(6.0/25.0*UvsxR(i,j)**2+Grav*abs(DvsxR(i,j)))
		 UR = 6.0/5.0*UvsxR(i,j)
         SvsxL(i,j) = UR-2.0*SQR_PHI_R
         SvsxR(i,j) = UR+SQR_PHI_R
       elseif(Maskvs(i-1,j)==1.and.Maskvs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(6.0/25.0*UvsxL(i,j)**2+Grav*abs(DvsxL(i,j)))
		 UL = 6.0/5.0*UvsxL(i,j)
         SvsxL(i,j) = UL-SQR_PHI_L
         SvsxR(i,j) = UL+2.0*SQR_PHI_L
       endif
     enddo
     enddo

     ! y-faces
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Maskvs(i,j-1)==1.and.Maskvs(i,j)==1) then
         SQR_PHI_L = sqrt(6.0/25.0*VvsyL(i,j)**2+Grav*abs(DvsyL(i,j)))
         SQR_PHI_R = sqrt(6.0/25.0*VvsyR(i,j)**2+Grav*abs(DvsyR(i,j)))
		 VL = 6.0/5.0*VvsyL(i,j)
		 VR = 6.0/5.0*VvsyR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL-VR)
         U_S = 0.5*(VL+VR)+SQR_PHI_L-SQR_PHI_R
         SvsyL(i,j) = min(VL-SQR_PHI_L,U_S-SQR_PHI_S)
         SvsyR(i,j) = max(VR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskvs(i,j-1)==0.and.Maskvs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(6.0/25.0*VvsyR(i,j)**2+Grav*abs(DvsyR(i,j)))
		 VR = 6.0/5.0*VvsyR(i,j)
         SvsyL(i,j) = VR-2.0*SQR_PHI_R
         SvsyR(i,j) = VR+SQR_PHI_R
       elseif(Maskvs(i,j-1)==1.and.Maskvs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(6.0/25.0*VvsyL(i,j)**2+Grav*abs(DvsyL(i,j)))
		 VL = 6.0/5.0*VvsyL(i,j)
         SvsyL(i,j) = VL-SQR_PHI_L
         SvsyR(i,j) = VL+2.0*SQR_PHI_L
       endif
     enddo
     enddo

     end subroutine wave_speed_vs
	 
     subroutine fluxes_at_faces_HLL_vs
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes_vs
!    Last update: 27/03/2017, Cheng Zhang
!---------------------------------------------
     use global
     implicit none

     ! horizontal fluxes
     call HLL_2D(Mloc1,Nloc,SvsxL,SvsxR,EvsxL,EvsxR,DvsxL, DvsxR, Evsx)
     call HLL_2D(Mloc,Nloc1,SvsyL,SvsyR,EvsyL,EvsyR,DvsyL, DvsyR, Evsy)
     call HLL_2D(Mloc1,Nloc,SvsxL,SvsxR,FvsxL,FvsxR,DUvsxL,DUvsxR,Fvsx)
     call HLL_2D(Mloc,Nloc1,SvsyL,SvsyR,FvsyL,FvsyR,DUvsyL,DUvsyR,Fvsy)
     call HLL_2D(Mloc1,Nloc,SvsxL,SvsxR,GvsxL,GvsxR,DVvsxL,DVvsxR,Gvsx)
     call HLL_2D(Mloc,Nloc1,SvsyL,SvsyR,GvsyL,GvsyR,DVvsyL,DVvsyR,Gvsy)   

     return
     end subroutine fluxes_at_faces_HLL_vs

     subroutine HLL_2D(M,N,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLL reconstruction for 2D model
!    Called by
!       fluxes_at_faces_HLL_vs
!    Last update: 27/03/2017, Cheng Zhang
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N
     REAL(SP),INTENT(IN),DIMENSION(M,N)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N)::FOUT
	 INTEGER :: I,J

      DO J=1,N
      DO I=1,M     
      IF(SL(I,J)>=ZERO) THEN
        FOUT(I,J)=FL(I,J)
      ELSEIF(SR(I,J)<=ZERO) THEN
        FOUT(I,J)=FR(I,J)
      ELSE
        FOUT(I,J)=SR(I,J)*FL(I,J)-SL(I,J)*FR(I,J)+SL(I,J)*SR(I,J)*(UR(I,J)-UL(I,J))
        IF((ABS(SR(I,J)-SL(I,J)))<SMALL)THEN
         FOUT(I,J)=FOUT(I,J)/SMALL
        ELSE
         FOUT(I,J)=FOUT(I,J)/(SR(I,J)-SL(I,J))
        ENDIF
      ENDIF
      ENDDO
      ENDDO

     end subroutine HLL_2D

     subroutine flux_bc_vs
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    for fluid slide
!    Called by
!       fluxes_vs
!    Last update: 27/03/2017, Cheng Zhang
!--------------------------------------------
     use global
     implicit none
     integer :: i,j

     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
	 do j = Jbeg,Jend
       Evsx(Ibeg,j) = Zero
	   Fvsx(Ibeg,j) = 0.5*Grav*(HvsxR(Ibeg,j)*HvsxR(Ibeg,j)-  &
	            2.0*HvsxR(Ibeg,j)*H0fx(Ibeg,j))
       Gvsx(Ibeg,j) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
       Evsx(Iend1,j) = Zero
       Fvsx(Iend1,j) = 0.5*Grav*(HvsxL(Iend1,j)*HvsxL(Iend1,j)-  &
	            2.0*HvsxL(Iend1,j)*H0fx(Iend1,j))	
       Gvsx(Iend1,j) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
       Evsy(i,Jbeg) = Zero
       Fvsy(i,Jbeg) = Zero
       Gvsy(i,Jbeg) = 0.5*Grav*(HvsyR(i,Jbeg)*HvsyR(i,Jbeg)-  &
	          2.0*HvsyR(i,Jbeg)*H0fy(i,Jbeg)) 
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
       Evsy(i,Jend1) = Zero
       Fvsy(i,Jend1) = Zero
	   Gvsy(i,Jend1) = 0.5*Grav*(HvsyL(i,Jend1)*HvsyL(i,Jend1)-  &
	           2.0*HvsyL(i,Jend1)*H0fy(i,Jend1))
     enddo
# if defined (PARALLEL)
     endif
# endif

    do j = Jbeg-1,Jend+1
    do i = Ibeg-1,Iend+1
      if(Maskvs(i,j)==0) then
        if(Maskvs(i-1,j)==1) then !right boundary
          Evsx(i,j) = Zero
          Fvsx(i,j) = 0.5*Grav*(HvsxL(i,j)*HvsxL(i,j)-  &
	            2.0*HvsxL(i,j)*H0fx(i,j))
          Gvsx(i,j) = Zero
        elseif(Maskvs(i+1,j)==1) then !left boundary
          Evsx(i+1,j) = Zero
          Fvsx(i+1,j) = 0.5*Grav*(HvsxR(i+1,j)*HvsxR(i+1,j)-  &
	            2.0*HvsxR(i+1,j)*H0fx(i+1,j))
          Gvsx(i+1,j) = Zero
        elseif(Maskvs(i,j-1)==1) then !back boundary
          Evsy(i,j) = Zero
          Fvsy(i,j) = Zero
          Gvsy(i,j) = 0.5*Grav*(HvsyL(i,j)*HvsyL(i,j)-  &
	           2.0*HvsyL(i,j)*H0fy(i,j))
        elseif(Maskvs(i,j+1)==1) then !front boundary
          Evsy(i,j+1) = Zero
          Fvsy(i,j+1) = Zero
          Gvsy(i,j+1) = 0.5*Grav*(HvsyR(i,j+1)*HvsyR(i,j+1)-  &
	          2.0*HvsyR(i,j+1)*H0fy(i,j+1)) 
        endif
      endif
    enddo
    enddo

     end subroutine flux_bc_vs
	 
	 
     subroutine source_terms_vs
!------------------------------------------------
!    This subroutine is used to evaluate source terms
!    in momentum eqs of fluid slide
!    Called by
!       update_bathymetry
!    Last update: 27/03/2017, Cheng Zhang
!------------------------------------------------
     use global
     implicit none
     integer :: i,j
	 real(SP), dimension(:,:),allocatable :: DelxPD,DelyPD,DelxPH,DelyPH,  &
	                                         DiffvsX,DiffvsY,Dwvs,PB
	 
	! source terms for water pressure
     allocate(DelxPD(Mloc,Nloc))
     allocate(DelyPD(Mloc,Nloc))
     allocate(DelxPH(Mloc,Nloc))
     allocate(DelyPH(Mloc,Nloc))
	 allocate(Dwvs(Mloc,Nloc))
	 allocate(PB(Mloc,Nloc))
	 DelxPD = zero;DelyPD = zero;DelxPH = zero;DelyPH = zero;Dwvs = zero;PB = zero
     do j = 1,Nloc
     do i = 1,Mloc
       PB(i,j) = P(i,j,kbeg)! dynamic pressure over the slide
	   Dwvs(i,j) = Rho0*grav*(Eta(i,j)+Hvs(i,j))! hydrostatic pressure over the slide
     enddo
     enddo
	 call delxFun_2D(PB,DelxPD) 
	 call delyFun_2D(PB,DelyPD)
	 call delxFun_2D(Dwvs,DelxPH)
	 call delyFun_2D(Dwvs,DelyPH)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrcpvsX(i,j) = -Dvs(i,j)/SLIDE_DENSITY*(DelxPD(i,j)+DelxPH(i,j))*Mask(i,j)
       SrcpvsY(i,j) = -Dvs(i,j)/SLIDE_DENSITY*(DelyPD(i,j)+DelyPH(i,j))*Mask(i,j)
     enddo
     enddo
	 
     ! source terms for hydrostatic pressure of slide
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrchvsX(i,j) = -Grav*Hvs(i,j)*DelxH0(i,j)*Maskvs(i,j)
       SrchvsY(i,j) = -Grav*Hvs(i,j)*DelyH0(i,j)*Maskvs(i,j)
     enddo
     enddo
	 
     ! source terms for diffusion
	 allocate(DiffvsX(Mloc,Nloc))
     allocate(DiffvsY(Mloc,Nloc))
	 DiffvsX = zero;DiffvsY = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   DiffvsX(i,j) = (Uvs(i+1,j)-2.*Uvs(i,j)+Uvs(i-1,j))/dx**2+  &
	                  (Uvs(i,j+1)-2.*Uvs(i,j)+Uvs(i,j-1))/dy**2
	   DiffvsY(i,j) = (Vvs(i+1,j)-2.*Vvs(i,j)+Vvs(i-1,j))/dx**2+  &
	                  (Vvs(i,j+1)-2.*Vvs(i,j)+Vvs(i,j-1))/dy**2
       SrcdvsX(i,j) = SLIDE_VISCOSITY*(Dvs(i,j)*DiffvsX(i,j))*Maskvs(i,j)
       SrcdvsY(i,j) = SLIDE_VISCOSITY*(Dvs(i,j)*DiffvsY(i,j))*Maskvs(i,j)
     enddo
     enddo
	 
	 deallocate(DelxPD)
     deallocate(DelyPD)
     deallocate(DelxPH)
     deallocate(DelyPH)
	 deallocate(Dwvs)
	 deallocate(PB)
     deallocate(DiffvsX)
	 deallocate(DiffvsY)

     end subroutine source_terms_vs
	 

	 subroutine eval_duv_vs(ISTEP)
!-----------------------------------------------
!    Update all variables Dvs,Uvs,Vvs for fluid slide
!    Called by
!       update_bathymetry
!    Last update: 27/03/2017, Cheng Zhang
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP),dimension(:,:),allocatable :: R1,R2,R3
     integer :: i,j
     REAL(SP) :: FROUDECAP_VS,FroudeU,DUU,Dangle

     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc))
     allocate(R3(Mloc,Nloc))
	 
     ! solve total slide depth Dvs
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskvs(i,j)==0) cycle
       R1(i,j) = -1.0/dx*(Evsx(i+1,j)-Evsx(i,j))-1.0/dy*(Evsy(i,j+1)-Evsy(i,j)) 
       Dvs(i,j) = ALPHA(ISTEP)*Dvs0(i,j)+BETA(ISTEP)*(Dvs(i,j)+dt*R1(i,j)) 
     enddo
     enddo
     ! update Dvs and Hvs          
     Dvs = max(Dvs,SLIDE_MINTHICK)
     call wl_bc_vs
     Hvs = Hc0-Dvs
	 
	 ! solve DUvs
     R2 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskvs(i,j)==0) cycle  
       R2(i,j) = -1.0/dx*(Fvsx(i+1,j)-Fvsx(i,j))-1.0/dy*(Fvsy(i,j+1)-Fvsy(i,j)) &
                 +SrcpvsX(i,j)+SrchvsX(i,j)+SrcdvsX(i,j)
	   DUvs(i,j) = DUvs(i,j)+dt*R2(i,j)
     enddo
     enddo
	 
	 ! solve DVvs
     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskvs(i,j)==0) cycle  
       R3(i,j) = -1.0/dx*(Gvsx(i+1,j)-Gvsx(i,j))-1.0/dy*(Gvsy(i,j+1)-Gvsy(i,j)) &
                 +SrcpvsY(i,j)+SrchvsY(i,j)+SrcdvsY(i,j)
	   DVvs(i,j) = DVvs(i,j)+dt*R3(i,j)
     enddo
     enddo
	 
     ! SSP Runge-Kutta time stepping
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	 ! implicit vertical stress (without manning friction)
	   DUvs(i,j) = DUvs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY/Dvs(i,j)**2)
	   DVvs(i,j) = DVvs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY/Dvs(i,j)**2)
	 ! implicit vertical stress (with manning friction)
!       DUvs(i,j) = DUvs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY/Dvs(i,j)**2 & 
!	               +Grav*0.01**2*sqrt(DUvs(i,j)**2+DVvs(i,j)**2)/Dvs(i,j)**(7./3.))
!       DVvs(i,j) = DVvs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY/Dvs(i,j)**2 &
!	               +Grav*0.01**2*sqrt(DUvs(i,j)**2+DVvs(i,j)**2)/Dvs(i,j)**(7./3.))
       DUvs(i,j) = ALPHA(ISTEP)*DUvs0(i,j)+BETA(ISTEP)*DUvs(i,j)
       DVvs(i,j) = ALPHA(ISTEP)*DVvs0(i,j)+BETA(ISTEP)*DVvs(i,j)
     enddo
     enddo
	 
     ! limit the maximum Froude number
     FROUDECAP_VS=0.5
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Maskvs(I,J)>0)THEN
         FroudeU=SQRT(grav*Dvs(I,J))*FROUDECAP_VS*Dvs(I,J)
         DUU=SQRT(DUvs(I,J)**2+DVvs(I,J)**2)
         IF(DUU>FroudeU)THEN
          Dangle=atan2(DVvs(I,J),DUvs(I,J))
          DUvs(I,J)=FroudeU*COS(Dangle)
          DVvs(I,J)=FroudeU*SIN(Dangle)
         ENDIF
       ENDIF
     ENDDO
     ENDDO

     ! update mask for fluid slide
	 call update_maskvs

     ! boundary conditions and final velocity
     call get_UV_vs

	 deallocate(R1)
     deallocate(R2)
     deallocate(R3)
	 
	 end subroutine eval_duv_vs
	 
     subroutine wl_bc_vs
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duv_vs
!    Last update: 27/03/2017, Cheng Zhang
!-----------------------------------------------------------
     use global
     implicit none
     integer :: i,j

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do i = 1,Nghost
       Dvs(Ibeg-i,j) = Dvs(Ibeg+i-1,j)
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
       Dvs(Iend+i,j) = Dvs(Iend-i+1,j)
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
         Dvs(i,j) = Dvs(i,Jbeg+Nghost-j)
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
         Dvs(i,Jend+j) = Dvs(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(Dvs)
# endif
     
     return
     end subroutine wl_bc_vs
	 
     subroutine get_UV_vs
!------------------------------------------------
!    Obtain Uvs,Vvs for fluid slide
!    Called by
!       eval_duv_vs
!    Last update: 27/03/2017, Cheng Zhang
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j
	 
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   if(Maskvs(i,j)==1) then
         Uvs(i,j) = DUvs(i,j)/Dvs(i,j)
         Vvs(i,j) = DVvs(i,j)/Dvs(i,j)
	   else
         Uvs(i,j) = Zero
         Vvs(i,j) = Zero
         DUvs(i,j) = Zero
         DVvs(i,j) = Zero
	   endif
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc_vs
# if defined (PARALLEL)
     call phi_2D_exch(Uvs)
     call phi_2D_exch(Vvs)
     call phi_2D_exch(DUvs)
     call phi_2D_exch(DVvs)
# endif

     end subroutine get_UV_vs

     subroutine vel_bc_vs
!----------------------------------------------------
!    Boundary conditions for slide velocity
!    Called by 
!       get_UV_vs
!    Last update: 27/03/2017, Cheng Zhang
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,imask

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do i = 1,Nghost
       Uvs(Ibeg-i,j) = -Uvs(Ibeg+i-1,j)
       Vvs(Ibeg-i,j) = Vvs(Ibeg+i-1,j)
       DUvs(Ibeg-i,j) = -DUvs(Ibeg+i-1,j)
       DVvs(Ibeg-i,j) = DVvs(Ibeg+i-1,j)
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
       Uvs(Iend+i,j) = -Uvs(Iend-i+1,j)
       Vvs(Iend+i,j) = Vvs(Iend-i+1,j)
       DUvs(Iend+i,j) = -DUvs(Iend-i+1,j)
       DVvs(Iend+i,j) = DVvs(Iend-i+1,j)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do j = 1,Nghost
       Uvs(i,Jbeg-j) = Uvs(i,Jbeg+j-1)
       Vvs(i,Jbeg-j) = -Vvs(i,Jbeg+j-1)
       DUvs(i,Jbeg-j) = DUvs(i,Jbeg+j-1)
       DVvs(i,Jbeg-j) = -DVvs(i,Jbeg+j-1)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do j = 1,Nghost
       Uvs(i,Jend+j) = Uvs(i,Jend-j+1)
       Vvs(i,Jend+j) = -Vvs(i,Jend-j+1)
       DUvs(i,Jend+j) = DUvs(i,Jend-j+1)
       DVvs(i,Jend+j) = -DVvs(i,Jend-j+1)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     end subroutine vel_bc_vs
# endif