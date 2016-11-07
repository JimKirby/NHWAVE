!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!
!    two_layer_slide.F
!
!    This file contains the subroutines which compute the slope-oriented, two-layer slide
!    based on granular flow rheology.  This does NOT contain Dmitri's viscous lower layer.
!
!    fluxes_ll
!    source_terms_ll
!    eval_huv_ll
!    convert_Ha_to_Hs
!
!    James Kirby, 6/27/16
!
!    Model described in:
!
!    Ma, G., Kirby, J. T., Hsu, T.-J. and Shi, F., 2015, "A two-layer granular landslide model for
!        tsunami wave generation: Theory and computation", Ocean Modelling, 93, 40-55, 
!        doi:10.1016/j.ocemod.2015.07.012 
!
!    This file is part of NHWAVE.
!
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
!---------------------------------------------------------------------------------------
!
# if defined (TWOLAYERSLIDE)
!
!---------------------------------------------------------------------------------------
!   The following subroutines are designed for granular flow landslide.
!   Update: Gangfeng Ma, Oct.15, 2013
!---------------------------------------------------------------------------------------
    subroutine source_terms_ll
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: sinx,siny,sdens,vfluid,CosX,CosY,Dudy,Dvdx,Umag,Ubar,Vbar
    real(SP), dimension(:,:), allocatable :: Delx2P,Dely2P,Delx2D,Dely2D,Pbot,DelxHc0,DelyHc0
    real(SP) :: tanbeta,cosbeta

    allocate(Pbot(Mloc,Nloc))
    allocate(Delx2P(Mloc,Nloc))
    allocate(Dely2P(Mloc,Nloc))
    allocate(Delx2D(Mloc,Nloc))
    allocate(Dely2D(Mloc,Nloc))
    allocate(DelxHc0(Mloc,Nloc))
    allocate(DelyHc0(Mloc,Nloc))

    SorceaX = Zero
    SorceaY = Zero

! source terms for lower layer
    call delxFun_2D(Hc0,DelxHc0)
    call delyFun_2D(Hc0,DelyHc0)
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      SorceaX(i,j) = ((1.0-SlideLambda)*Kap(i,j)+SlideLambda)*Grav_LL(i,j)*Ha(i,j)*DelxHc0(i,j)
      SorceaY(i,j) = ((1.0-SlideLambda)*Kap(i,j)+SlideLambda)*Grav_LL(i,j)*Ha(i,j)*DelyHc0(i,j)
    enddo
    enddo

! shear stress at the landslide surface
    ForceX_UL = 0.0; ForceY_UL = 0.0
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0.or.Mask(i,j)==0) cycle
      Umag = sqrt((Ua(i,j)-U(i,j,Kbeg))**2+(Va(i,j)-V(i,j,Kbeg))**2)
      ForceX_UL(i,j) = 0.5*Cf_ul*Umag*(Ua(i,j)-U(i,j,Kbeg))
      ForceY_UL(i,j) = 0.5*Cf_ul*Umag*(Va(i,j)-V(i,j,Kbeg))
      SorceaX(i,j) = SorceaX(i,j)-ForceX_UL(i,j)
      SorceaY(i,j) = SorceaY(i,j)-ForceY_UL(i,j)
    enddo
    enddo

    ! Coulomb friction
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      Umag = sqrt(Ua(i,j)**2+Va(i,j)**2)

      SorceaX(i,j) = SorceaX(i,j)-(1.0-SlideLambda)*Grav_LL(i,j)*Ha(i,j)*tan(PhiBed)*Ua(i,j)/(Umag+1.e-16)
      SorceaY(i,j) = SorceaY(i,j)-(1.0-SlideLambda)*Grav_LL(i,j)*Ha(i,j)*tan(PhiBed)*Va(i,j)/(Umag+1.e-16)

!      Sxy = 0.5*(DelyUa(i,j)+DelxVa(i,j))
!      if(Sxy.ne.0.0) then
!        SorceaX(i,j) = SorceaX(i,j)-sign(1.0,Sxy)* &
!            (1.0-SlideLambda)*Grav_LL(i,j)*Ha(i,j)*Kap(i,j)*DelyHa(i,j)*sin(PhiInt)                            
!        SorceaY(i,j) = SorceaY(i,j)-sign(1.0,Sxy)*  &
!            (1.0-SlideLambda)*Grav_LL(i,j)*Ha(i,j)*Kap(i,j)*DelxHa(i,j)*sin(PhiInt)                            
!      endif
    enddo
    enddo

    ! dynamic+hydrostatic pressure at the interface
    do j = 1,Nloc
    do i = 1,Mloc
      Pbot(i,j) = P(i,j,Kbeg)
    enddo
    enddo

    call delxFun_2D(Pbot,Delx2P)
    call delyFun_2D(Pbot,Dely2P)
    call delxFun_2D(D,Delx2D)
    call delyFun_2D(D,Dely2D)
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0.or.Mask(i,j)==0) cycle
      SorceaX(i,j) = SorceaX(i,j)-Ha(i,j)/SlideDens*Delx2P(i,j)-  &
                        Rho0/SlideDens*Grav_LL(i,j)*Ha(i,j)*Delx2D(i,j)
      SorceaY(i,j) = SorceaY(i,j)-Ha(i,j)/SlideDens*Dely2P(i,j)-  &
                        Rho0/SlideDens*Grav_LL(i,j)*Ha(i,j)*Dely2D(i,j)
    enddo
    enddo

    deallocate(Pbot)
    deallocate(Delx2P)
    deallocate(Dely2P)
    deallocate(Delx2D)
    deallocate(Dely2D)
    deallocate(DelxHc0)
    deallocate(DelyHc0)

    end subroutine source_terms_ll
!
!------------------------------------------------------------------------------------
!
    subroutine update_hc_ul
!
    use global
    implicit none
    integer :: i,j

!   save old bathymetry
    Ho = Hc

!   update the thickness of the upper layer
    do j = 1,Nloc
    do i = 1,Mloc
      if(Maska(i,j)==1) then
        Hc(i,j) = Hc0(i,j)-Ha(i,j)
      else
        Hc(i,j) = Hc0(i,j)
      endif
    enddo
    enddo

!   reconstruct depth at x-y faces
!
    do j = 1,Nloc
    do i = 2,Mloc
      Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i-1,j))
    enddo
    Hfx(1,j) = Hfx(2,j)
    Hfx(Mloc1,j) = Hfx(Mloc,j)
    enddo

    do i = 1,Mloc
    do j = 2,Nloc
      Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j-1))
    enddo
    Hfy(i,1) = Hfy(i,2)
    Hfy(i,Nloc1) = Hfy(i,Nloc)
    enddo
!
!   derivatives of water depth at cell center
!
    do j = 1,Nloc
    do i = 1,Mloc
      DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
      DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
    enddo
    enddo
!
!   time derivative of water depth
!
    DeltHo = DeltH

    DeltH = zero
    do j = 1,Nloc
    do i = 1,Mloc
      DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
    enddo
    enddo
!
!   second-order time derivative
!
    if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

    end subroutine update_hc_ul
!
!--------------------------------------------------------------------------------------------------------
!
    subroutine update_maska
    use global
    implicit none
    integer :: i,j

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) then
        if(Maska(i-1,j)==1.and.Ha(i-1,j)>Ha(i,j)) Maska(i,j)=1
        if(Maska(i+1,j)==1.and.Ha(i+1,j)>Ha(i,j)) Maska(i,j)=1
        if(Maska(i,j-1)==1.and.Ha(i,j-1)>Ha(i,j)) Maska(i,j)=1
        if(Maska(i,j+1)==1.and.Ha(i,j+1)>Ha(i,j)) Maska(i,j)=1
      else
        ! drying (wet->dry)               
        if(Ha(i,j)-0.001<=1.e-6) then
          Maska(i,j) = 0
          Ha(i,j) = 0.001
        endif
      endif
    enddo
    enddo

# if defined (PARALLEL)
     ! collect mask into ghost cells  
     call phi_int_exch(Maska)
# endif

    end subroutine update_maska
!
!-------------------------------------------------------------------------------------------------------
!
    
    subroutine fluxes_ll
!
    use global
    implicit none

    ! calculate granular flow parameters
    call calc_kap

    ! second order construction
    call delxyFun_ll
    call construction_ll

    ! calculate wave speed
    call wave_speed_ll

    ! calculate fluxes at faces
    call fluxes_at_faces_HLL_ll

    ! flux bc
    call flux_bc_ll

    end subroutine fluxes_ll
!
!--------------------------------------------------------------------------------------------------------
!

    subroutine eval_huv_ll(ISTEP)
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:), allocatable :: R1,R2,R3 
    integer :: i,j
    real(SP) :: Vslide
    
    allocate(R1(Mloc,Nloc))
    allocate(R2(Mloc,Nloc))
    allocate(R3(Mloc,Nloc))
    
    R1 = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      R1(i,j) = -1.0/dx*(Eax(i+1,j)-Eax(i,j))-1.0/dy*(Eay(i,j+1)-Eay(i,j))
      Ha(i,j) = ALPHA(ISTEP)*Ha0(i,j)+BETA(ISTEP)*(Ha(i,j)+dt*R1(i,j))
    enddo
    enddo

    R2 = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      R2(i,j) = -1.0/dx*(Fax(i+1,j)-Fax(i,j))-1.0/dy*(Fay(i,j+1)-Fay(i,j))+SorceaX(i,j)
      HUa(i,j) = ALPHA(ISTEP)*HUa0(i,j)+BETA(ISTEP)*(HUa(i,j)+dt*R2(i,j))
    enddo
    enddo


    R3 = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      R3(i,j) = -1.0/dx*(Gax(i+1,j)-Gax(i,j))-1.0/dy*(Gay(i,j+1)-Gay(i,j))+SorceaY(i,j)
      HVa(i,j) = ALPHA(ISTEP)*HVa0(i,j)+BETA(ISTEP)*(HVa(i,j)+dt*R3(i,j))
    enddo
    enddo

    ! simple bc
    call phi_2D_coll(Ha)
    call phi_2D_coll(HUa)
    call phi_2D_coll(HVa)

    call update_maska

    do j = 1,Nloc
    do i = 1,Mloc
      if(Maska(i,j)==1) then
        Ua(i,j) = HUa(i,j)/Ha(i,j)
        Va(i,j) = HVa(i,j)/Ha(i,j)
      else
        Ua(i,j) = Zero
        Va(i,j) = Zero
        HUa(i,j) = Zero
        HVa(i,j) = Zero
      endif
    enddo
    enddo

    deallocate(R1)
    deallocate(R2)
    deallocate(R3)

    return
    end subroutine eval_huv_ll
!
!------------------------------------------------------------------------------------------------------------
!

    subroutine flux_bc_ll
    use global
    implicit none
    integer :: i,j

    ! left and right side
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
      Eax(Ibeg,j) = Zero
      Fax(Ibeg,j) = 0.5*GravxR(Ibeg,j)*((1.0-SlideLambda)*KapxR(Ibeg,j)+  &
                    SlideLambda)*HaxR(Ibeg,j)*HaxR(Ibeg,j)
      Gax(Ibeg,j) = Zero
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
      Eax(Iend1,j) = Zero
      Fax(Iend1,j) = 0.5*GravxL(Iend1,j)*((1.0-SlideLambda)*KapxL(Iend1,j)+  &
                     SlideLambda)*HaxL(Iend1,j)*HaxL(Iend1,j)
      Gax(Iend1,j) = Zero
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
      Eay(i,Jbeg) = Zero
      Fay(i,Jbeg) = Zero
      Gay(i,Jbeg) = 0.5*GravyR(i,Jbeg)*((1.0-SlideLambda)*KapyR(i,Jbeg)+  &
                    SlideLambda)*HayR(i,Jbeg)*HayR(i,Jbeg)
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
      Eay(i,Jend1) = Zero
      Fay(i,Jend1) = Zero
      Gay(i,Jend1) = 0.5*GravyL(i,Jend1)*((1.0-SlideLambda)*KapyL(i,Jend1)+  &
                     SlideLambda)*HayL(i,Jend1)*HayL(i,Jend1)
    enddo

# if defined (PARALLEL)
    endif
# endif

    do j = Jbeg-1,Jend+1
    do i = Ibeg-1,Iend+1
      if(Maska(i,j)==0) then
        if(Maska(i-1,j)==1) then
          Eax(i,j) = Zero
          Fax(i,j) = 0.5*GravxL(i,j)*((1.0-SlideLambda)*Kap(i-1,j)+  & 
                     SlideLambda)*HaxL(i,j)*HaxL(i,j)
          Gax(i,j) = Zero
        elseif(Maska(i+1,j)==1) then
          Eax(i+1,j) = Zero
          Fax(i+1,j) = 0.5*GravxR(i+1,j)*((1.0-SlideLambda)*Kap(i+1,j)+  & 
                       SlideLambda)*HaxR(i+1,j)*HaxR(i+1,j)
          Gax(i+1,j) = Zero
        elseif(Maska(i,j-1)==1) then
          Eay(i,j) = Zero
          Fay(i,j) = Zero
          Gay(i,j) = 0.5*GravyL(i,j)*((1.0-SlideLambda)*Kap(i,j-1)+  & 
                     SlideLambda)*HayL(i,j)*HayL(i,j)
        elseif(Maska(i,j+1)==1) then
          Eay(i,j+1) = Zero
          Fay(i,j+1) = Zero
          Gay(i,j+1) = 0.5*GravyR(i,j+1)*((1.0-SlideLambda)*Kap(i,j+1)+  & 
                       SlideLambda)*HayR(i,j+1)*HayR(i,j+1)
        endif
      endif
    enddo
    enddo

    end subroutine flux_bc_ll
!
!---------------------------------------------------------------------------------------------------------
!

    subroutine delxyFun_ll
    use global
    implicit none

    call delxFun_2D(Ha,DelxHa)
    call delxFun_2D(Ua,DelxUa)
    call delxFun_2D(Va,DelxVa)
    call delxFun_2D(HUa,DelxHUa)
    call delxFun_2D(HVa,DelxHVa)
    call delxFun_2D(Kap,DelxKap)
    call delxFun_2D(Grav_LL,DelxGrav)

    call delyFun_2D(Ha,DelyHa)
    call delyFun_2D(Ua,DelyUa)
    call delyFun_2D(Va,DelyVa)
    call delyFun_2D(HUa,DelyHUa)
    call delyFun_2D(HVa,DelyHVa)
    call delyFun_2D(Kap,DelyKap)
    call delyFun_2D(Grav_LL,DelyGrav)

    end subroutine delxyFun_ll
!
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!
    subroutine adjust_grav
!
!---------------------------------------------------------------------------------------------------------
!
!   adjust_grav
!
!   Adjust gravity to account for non-hydrostatic effects as in Denlinger and Iverson (2004)
!
!----------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j
    real(SP),dimension(:,:),allocatable :: DelxHc0,DelyHc0,wbar_n

    allocate(DelxHc0(Mloc,Nloc))
    allocate(DelyHc0(Mloc,Nloc))
    allocate(wbar_n(Mloc,Nloc))

    call delxFun_2D(Hc0,DelxHc0)
    call delyFun_2D(Hc0,DelyHc0)

    wbar_n = wbar

    wbar = zero
    do j = 1,Nloc
    do i = 1,Mloc
      if(maska(i,j)==1) then
        wbar(i,j) = 0.5*(-DeltH(i,j)-Ua(i,j)*DelxH(i,j)-Va(i,j)*DelyH(i,j))+  &
                  0.5*(-Ua(i,j)*DelxHc0(i,j)-Va(i,j)*DelyHc0(i,j))
      endif
    enddo
    enddo

    do j = 1,Nloc
    do i = 1,Mloc
      if(maska(i,j)==1) then
        write(11,*) (wbar(i,j)-wbar_n(i,j))/dt
      endif
    enddo
    enddo

    do j = 1,Nloc
    do i = 1,Mloc
      Grav_LL(i,j) = Grav
    enddo
    enddo

    deallocate(DelxHc0)
    deallocate(DelyHc0)
    deallocate(wbar_n)

    return
    end subroutine adjust_grav
!
!---------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
!
    subroutine calc_kap
!
!--------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j
    real(SP) :: Div,Umag

    call delxFun_2D(Ua,DelxUa)
    call delyFun_2D(Va,DelyVa)
!
!   The Earth pressure coefficient depends on the divergence of the flow.
!
    do j = Jbeg,Jend
    do i = Ibeg,Iend
       if(Maska(i,j)==0) cycle
      Umag = sqrt(Ua(i,j)**2+Va(i,j)**2)
      if(Umag>1.e-6) then
        ! calculate divergence
        Div = DelxUa(i,j)+DelyVa(i,j)

        if(Div>1.e-6) then
          Kap(i,j) = 2.0*(1.0-sqrt(1.0-cos(PhiInt)**2*(1.0+tan(PhiBed)**2)))/cos(PhiInt)**2-1.0              
        elseif(Div<-1.e-6) then
          Kap(i,j) = 2.0*(1.0+sqrt(1.0-cos(PhiInt)**2*(1.0+tan(PhiBed)**2)))/cos(PhiInt)**2-1.0              
        else
          Kap(i,j) = 1.0
        endif
      else
        Kap(i,j) = 1.0
      endif
    enddo
    enddo
    
    call phi_2D_coll(Kap)

    return

    end subroutine calc_kap
!
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!
    subroutine construction_ll
!
!-------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j

    call construct_2D_x(Ha,DelxHa,HaxL,HaxR)
    call construct_2D_x(Ua,DelxUa,UaxL,UaxR)
    call construct_2D_x(Va,DelxVa,VaxL,VaxR)
    call construct_2D_x(HUa,DelxHUa,HUaxL,HUaxR)
    call construct_2D_x(HVa,DelxHVa,HVaxL,HVaxR)
    call construct_2D_x(Kap,DelxKap,KapxL,KapxR)
    call construct_2D_x(Grav_LL,DelxGrav,GravxL,GravxR)

    do j = Jbeg,Jend
    do i = Ibeg,Iend1
      EaxL(i,j) = HUaxL(i,j)
      EaxR(i,j) = HUaxR(i,j)
      FaxL(i,j) = HUaxL(i,j)*UaxL(i,j)+0.5*GravxL(i,j)*((1.0-SlideLambda)*KapxL(i,j)+  &
                  SlideLambda)*HaxL(i,j)*HaxL(i,j)
      FaxR(i,j) = HUaxR(i,j)*UaxR(i,j)+0.5*GravxR(i,j)*((1.0-SlideLambda)*KapxR(i,j)+  &
                  SlideLambda)*HaxR(i,j)*HaxR(i,j)
      GaxL(i,j) = HaxL(i,j)*UaxL(i,j)*VaxL(i,j)
      GaxR(i,j) = HaxR(i,j)*UaxR(i,j)*VaxR(i,j)
    enddo
    enddo

    call construct_2D_y(Ha,DelyHa,HayL,HayR)
    call construct_2D_y(Ua,DelyUa,UayL,UayR)
    call construct_2D_y(Va,DelyVa,VayL,VayR)
    call construct_2D_y(HUa,DelyHUa,HUayL,HUayR)
    call construct_2D_y(HVa,DelyHVa,HVayL,HVayR)
    call construct_2D_y(Kap,DelyKap,KapyL,KapyR)
    call construct_2D_y(Grav_LL,DelyGrav,GravyL,GravyR)

    do j = Jbeg,Jend1
    do i = Ibeg,Iend
      EayL(i,j) = HVayL(i,j)
      EayR(i,j) = HVayR(i,j)
      FayL(i,j) = HayL(i,j)*UayL(i,j)*VayL(i,j)
      FayR(i,j) = HayR(i,j)*UayR(i,j)*VayR(i,j)
      GayL(i,j) = HVayL(i,j)*VayL(i,j)+0.5*GravyL(i,j)*((1.0-SlideLambda)*KapyL(i,j)+  &
                  SlideLambda)*HayL(i,j)*HayL(i,j)
      GayR(i,j) = HVayR(i,j)*VayR(i,j)+0.5*GravyR(i,j)*((1.0-SlideLambda)*KapyR(i,j)+  &
                  SlideLambda)*HayR(i,j)*HayR(i,j)
    enddo
    enddo

    end subroutine construction_ll    
!
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
!

    subroutine wave_speed_ll
!
!--------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j
    real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S

    do j = Jbeg,Jend
    do i = Ibeg,Iend1
      SQR_PHI_L = sqrt(GravxL(i,j)*((1.0-SlideLambda)*KapxL(i,j)+SlideLambda)*abs(HaxL(i,j)))
      SQR_PHI_R = sqrt(GravxR(i,j)*((1.0-SlideLambda)*KapxR(i,j)+SlideLambda)*abs(HaxR(i,j)))
      SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UaxL(i,j)-UaxR(i,j))
      U_S = 0.5*(UaxL(i,j)+UaxR(i,j))+SQR_PHI_L-SQR_PHI_R
      SaxL(i,j) = min(UaxL(i,j)-SQR_PHI_L,U_S-SQR_PHI_S)
      SaxR(i,j) = max(UaxR(i,j)+SQR_PHI_R,U_S+SQR_PHI_S)
      SaxS(i,j) = U_S
    enddo
    enddo

    do j = Jbeg,Jend1
    do i = Ibeg,Iend
      SQR_PHI_L = sqrt(GravyL(i,j)*((1.0-SlideLambda)*KapyL(i,j)+SlideLambda)*abs(HayL(i,j)))
      SQR_PHI_R = sqrt(GravyR(i,j)*((1.0-SlideLambda)*KapyR(i,j)+SlideLambda)*abs(HayR(i,j)))
      SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VayL(i,j)-VayR(i,j))
      U_S = 0.5*(VayL(i,j)+VayR(i,j))+SQR_PHI_L-SQR_PHI_R
      SayL(i,j) = min(VayL(i,j)-SQR_PHI_L,U_S-SQR_PHI_S)
      SayR(i,j) = max(VayR(i,j)+SQR_PHI_R,U_S+SQR_PHI_S)
      SayS(i,j) = U_S
    enddo
    enddo

    end subroutine wave_speed_ll
!
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!
    subroutine fluxes_at_faces_HLL_ll
!
!------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j

    ! horizontal fluxes
    call HLL2D(Mloc1,Nloc,SaxL,SaxR,EaxL,EaxR,HaxL,HaxR,Eax)
    call HLL2D(Mloc,Nloc1,SayL,SayR,EayL,EayR,HayL,HayR,Eay)
    call HLL2D(Mloc1,Nloc,SaxL,SaxR,FaxL,FaxR,HUaxL,HUaxR,Fax)
    call HLL2D(Mloc,Nloc1,SayL,SayR,FayL,FayR,HUayL,HUayR,Fay)
    call HLL2D(Mloc1,Nloc,SaxL,SaxR,GaxL,GaxR,HVaxL,HVaxR,Gax)
    call HLL2D(Mloc,Nloc1,SayL,SayR,GayL,GayR,HVayL,HVayR,Gay)

    end subroutine fluxes_at_faces_HLL_ll
!
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!
    subroutine HLL2D(M,N,SL,SR,FL,FR,UL,UR,FOUT)
!
!--------------------------------------------------------------------------------------------------
!
    use global, only: SP,ZERO,SMALL
    implicit none
    INTEGER,INTENT(IN)::M,N
    REAL(SP),INTENT(IN),DIMENSION(M,N)::SL,SR,FL,FR,UL,UR
    REAL(SP),INTENT(OUT),DIMENSION(M,N)::FOUT
    INTEGER :: I,J

    DO J = 1,N
    DO I = 1,M
      IF(SL(I,J)>=ZERO) THEN
        FOUT(I,J) = FL(I,J)
      ELSEIF(SR(I,J)<=ZERO) THEN
        FOUT(I,J) = FR(I,J)
      ELSE
        FOUT(I,J) = SR(I,J)*FL(I,J)-SL(I,J)*FR(I,J)+  &
              SL(I,J)*SR(I,J)*(UR(I,J)-UL(I,J))
        IF((ABS(SR(I,J)-SL(I,J)))<SMALL)THEN
          FOUT(I,J) = FOUT(I,J)/SMALL
        ELSE
          FOUT(I,J) = FOUT(I,J)/(SR(I,J)-SL(I,J))
        ENDIF
      ENDIF
    ENDDO
    ENDDO

    return

    end subroutine HLL2D
!
!-------------------------------------------------------------------------------------------------------------

# endif
