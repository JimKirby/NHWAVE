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
!
    subroutine specify_slide
!
!   Specify initial shape of slide mass
!
!   This looks like it just does a triangular shape slide at this point.  Where is the
!   Fritz-like slide done?
!
!----------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j
    real(SP) :: alpha0,zt
!
!   triangle shape landslide (Kimmoun)
!
    alpha0 = SlopeAngle*pi/180.
    do j = 1,Nloc
    do i = 1,Mloc
!
!   Test 1
!
      if(xc(i)>=0.042.and.xc(i)<=0.097) then
        zt = (xc(i)-0.042)/cos(alpha0)*tan(alpha0)
      elseif(xc(i)>0.097.and.xc(i)<=0.152) then
        zt = (0.152-xc(i))/sin(alpha0)
      else
        zt = 0.0
      endif

! Test 2
!      if(xc(i)>=0.015.and.xc(i)<=0.0825) then
!        zt = (xc(i)-0.015)/cos(alpha0)*tan(alpha0)                                                                               
!      elseif(xc(i)>0.0825.and.xc(i)<=0.15) then                                                                                  
!        zt = (0.15-xc(i))/sin(alpha0)                                                                                           
!      else                                                                                                                       
!        zt = 0.0                                                                                                                 
!      endif

      Ha(i,j) = zt
      if(Ha(i,j)-0.001<=1.e-6) then                                                                                    
        Maska(i,j) = 0                                                                                                 
        Ha(i,j) = 0.001                                                                                                
      endif                                                                                                            
    enddo                                                                                                              
    enddo                                                                                                              
  
    Ha0 = Ha

    do j = 1,Nloc
    do i = 1,Mloc
      if(Maska(i,j)==1) then
        Ua(i,j) = SlideIniU
      endif
    enddo
    enddo
!
!   update the thickness of upper layer
!
    do j = 1,Nloc
    do i = 1,Mloc
      if(Maska(i,j)==1) then
        Hc(i,j) = Hc0(i,j)-Ha(i,j)
      else
        Hc(i,j) = Hc0(i,j)
      endif
    enddo
    enddo

    call phi_2D_coll(Ua)
    Ua0 = Ua
    HUa = Ha*Ua
    HUa0 = HUa
!
!   reconstruct depth at x-y faces
!
    do j = 1,Nloc
    do i = 2,Mloc
      Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
    enddo
    Hfx(1,j) = Hc(1,j)
    Hfx(Mloc1,j) = Hc(Mloc,j)
    enddo

    do i = 1,Mloc
    do j = 2,Nloc
      Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
    enddo
    Hfy(i,1) = Hc(i,1)
    Hfy(i,Nloc1) = Hc(i,Nloc)
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

    return

    end subroutine specify_slide
!
!---------------------------------------------------------------------------------
!
    subroutine adjust_grav
!
!---------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j
    real(SP) :: DepAngx,DepAngy,CosAng,SinAng

    do j = Jbeg,Jend
    do i = Ibeg,Iend
!
!      DepAngx = (Hc0(i+1,j)-Hc0(i-1,j))/(2.0*dx)
!
      if(xc(i)<=0.3) then
        DepAngx = SlopeAngle*pi/180.
      else
        DepAngx = 0.0
      endif
      CosAng = cos(DepAngx)
      SinAng = sin(DepAngx)
      Grav_Lx(i,j) = Grav*SinAng
      Grav_Ly(i,j) = 0.0
      Grav_Lz(i,j) = Grav*CosAng
    enddo
    enddo
    call phi_2D_coll(Grav_Lx)
    call phi_2D_coll(Grav_Ly)
    call phi_2D_coll(Grav_Lz)

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then  ! account for underwater pore pressure
        Slambda(i,j) = SlideLambda
      else ! no pore pressure in subaerial setting
        Slambda(i,j) = 0.0
      endif
    enddo
    enddo
    call phi_2D_coll(Slambda)

    return
    end subroutine adjust_grav
!
!---------------------------------------------------------------------------------
!
    subroutine source_terms_ll
!
!---------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: sdens,vfluid,CosX,CosY,Dudy,Dvdx,Umag,Ubar,Vbar
    real(SP), dimension(:,:), allocatable :: Delx2P,Dely2P,Pbot
    real(SP) :: Sxy,DepAngx,CosAng

    allocate(Pbot(Mloc,Nloc))
    allocate(Delx2P(Mloc,Nloc))
    allocate(Dely2P(Mloc,Nloc))

    SorceaX = Zero
    SorceaY = Zero
!
!   source terms for lower layer
!
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      SorceaX(i,j) = Grav_Lx(i,j)*Ha(i,j)
      SorceaY(i,j) = Grav_Ly(i,j)*Ha(i,j)
    enddo
    enddo

!    ! contributions from fluid stresses (minor effects)
!    ForceX_UL = 0.0; ForceY_UL = 0.0
!    do j = Jbeg,Jend
!    do i = Ibeg,Iend
!      if(Maska(i,j)==0.or.Mask(i,j)==0) cycle
!      Umag = sqrt((Ua(i,j)-U(i,j,Kbeg))**2+(Va(i,j)-V(i,j,Kbeg))**2)
!      ForceX_UL(i,j) = 0.5*Cf_ul*Umag*(Ua(i,j)-U(i,j,Kbeg))
!      ForceY_UL(i,j) = 0.5*Cf_ul*Umag*(Va(i,j)-V(i,j,Kbeg))
!      SorceaX(i,j) = SorceaX(i,j)-ForceX_UL(i,j)
!      SorceaY(i,j) = SorceaY(i,j)-ForceY_UL(i,j)
!    enddo
!    enddo

    ! contributions from solid stresses                                                              
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle
      
      Umag = sqrt(Ua(i,j)**2+Va(i,j)**2)
      SorceaX(i,j) = SorceaX(i,j)-(1.0-Slambda(i,j))*Grav_Lz(i,j)*Ha(i,j)*  &
                       tan(PhiBed)*Ua(i,j)/(Umag+1.e-16)
      SorceaY(i,j) = SorceaY(i,j)-(1.0-Slambda(i,j))*Grav_Lz(i,j)*Ha(i,j)*  &
                       tan(PhiBed)*Va(i,j)/(Umag+1.e-16)

      Sxy = 0.5*(DelyUa(i,j)+DelxVa(i,j))
      if(Sxy.ne.0.0) then
        SorceaX(i,j) = SorceaX(i,j)-sign(1.0,Sxy)*  & 
                       (1.0-Slambda(i,j))*Grav_Lz(i,j)*Ha(i,j)*Kap(i,j)*DelyHa(i,j)*sin(PhiInt)
        SorceaY(i,j) = SorceaY(i,j)-sign(1.0,Sxy)*  &  
                       (1.0-Slambda(i,j))*Grav_Lz(i,j)*Ha(i,j)*Kap(i,j)*DelxHa(i,j)*sin(PhiInt)
      endif
    enddo
    enddo

    k = Kbeg
    do j = 1,Nloc
    do i = 1,Mloc
      Pbot(i,j) = P(i,j,k)
    enddo
    enddo

    Delx2P = Zero; Dely2P = Zero
    call DelxFun_2D(Pbot,Delx2P)
    call DelyFun_2D(Pbot,Dely2P)

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0.or.Mask(i,j)==0) cycle
      SorceaX(i,j) = SorceaX(i,j)-Rho0/SlideDens*Grav*Ha(i,j)*(DelxH(i,j)+DelxEta(i,j))  &
                          -Ha(i,j)/SlideDens*Delx2P(i,j)
      SorceaY(i,j) = SorceaY(i,j)-Rho0/SlideDens*Grav*Ha(i,j)*(DelyH(i,j)+DelyEta(i,j))  &
                          -Ha(i,j)/SlideDens*Dely2P(i,j)                     
    enddo
    enddo

    deallocate(Pbot)
    deallocate(Delx2P)
    deallocate(Dely2P)

    end subroutine source_terms_ll
!
!--------------------------------------------------------------------------------------------
!
    subroutine convert_Ha_to_Hs
!
!   Convert the height of landslide in slope-oriented coordinate (Ha) to Cartesian
!   coordinate (Hs)
!
!   waiting for improvements
!
!---------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j,i1,ic
    real(SP) :: SlopeAng,SearchAng

    ! initial values
    Hs = 0.0

    ! conversion depending on SlopeAng
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maska(i,j)==0) cycle

      if(Ha(i,j)<0.01) then
        Hs(i,j) = Ha(i,j)
      else

      if(xc(i)<=1.724) then
        SlopeAng = SlopeAngle*pi/180.
      else
        SlopeAng = 0.0
      endif
      ! no conversion if slope is too small
      SearchAng = 0.0
      if(SlopeAng>1.e-6) then  
        ! find the location of vertical height Hs
        ic = 1
        do i1 = 1,20  ! search ten grids
          if(Maska(i-i1,j)==1) then
            SearchAng = float(i1)*dx/cos(SlopeAng)/Ha(i-i1,j)
            if(SearchAng>=SlopeAng) then
              ic = i1
              exit
            endif
          endif   
        enddo
        if(Maska(i-1,j)==0) then
          Hs(i,j) = Ha(i,j)
        else
          Hs(i,j) = 0.5*(Ha(i-ic,j)+Ha(i-ic+1,j))/cos(SlopeAng)
        endif
      elseif(SlopeAng<-1.e-6) then
        ic = 1
        do i1 = 1,20
          if(Maska(i+i1,j)==1) then
            SearchAng = float(i1)*dx/cos(abs(SlopeAng))/Ha(i+i1,j)
            if(SearchAng>=SlopeAng) then
              ic = i1
              exit
            endif
          endif       
        enddo
        if(Maska(i+1,j)==0) then
          Hs(i,j) = Ha(i,j)
        else
          Hs(i,j) = 0.5*(Ha(i+ic,j)+Ha(i+ic-1,j))/cos(abs(SlopeAng))
        endif
      else
        Hs(i,j) = Ha(i,j)
      endif
      endif 
    enddo
    enddo   

    end subroutine convert_Ha_to_Hs
!
!-----------------------------------------------------------------------------------
!
    subroutine update_hc_ul
!
!-----------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j

    ! save old bathymetry
    Ho = Hc

    ! convert Ha
!    call convert_Ha_to_Hs

!   update the thickness of the upper layer
!
    do j = 1,Nloc
    do i = 1,Mloc
      if(Maska(i,j)==1) then
        Hc(i,j) = Hc0(i,j)-Ha(i,j)
      else
        Hc(i,j) = Hc0(i,j)
      endif
    enddo
    enddo

    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==0) then
        Eta(i,j) = MinDep-Hc(i,j)
        D(i,j) = Eta(i,j)+Hc(i,j)
      endif
    enddo
    enddo

!   reconstruct depth at x-y faces
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
      if(Mask(i,j)==1) then
        DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
      endif
    enddo
    enddo

    ! second-order time derivative
    if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

    end subroutine update_hc_ul

!
!--------------------------------------------------------------------------------------------------------
!
    subroutine update_maska
!
!---------------------------------------------------------------------------------------------------------
!
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
!-------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
!
!   calculate granular flow parameters
!
    call calc_kap
!
!   second order construction
!
    call delxyFun_ll
    call construction_ll
!
!   calculate wave speed
!
    call wave_speed_ll
!
!   calculate fluxes at faces
!
    call fluxes_at_faces_HLL_ll
!
!   flux bc
!
    call flux_bc_ll

    end subroutine fluxes_ll
!
!--------------------------------------------------------------------------------------------------------
!
    subroutine eval_huv_ll(ISTEP)
!
!--------------------------------------------------------------------------------------------------------
!
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
!
!------------------------------------------------------------------------------------------------------------
!
    use global
    implicit none
    integer :: i,j

    ! left and right side
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
      Eax(Ibeg,j) = Zero
      Fax(Ibeg,j) = 0.5*GravxR(Ibeg,j)*((1.0-Slambda(Ibeg,j))*KapxR(Ibeg,j)+Slambda(Ibeg,j))*HaxR(Ibeg,j)*HaxR(Ibeg,j)
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
      Fax(Iend1,j) = 0.5*GravxL(Iend1,j)*((1.0-Slambda(Iend1,j))*KapxL(Iend1,j)+Slambda(Iend1,j))*HaxL(Iend1,j)*HaxL(Iend1,j)
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
      Gay(i,Jbeg) = 0.5*GravyR(i,Jbeg)*((1.0-Slambda(i,Jbeg))*KapyR(i,Jbeg)+Slambda(i,Jbeg))*HayR(i,Jbeg)*HayR(i,Jbeg)
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
      Gay(i,Jend1) = 0.5*GravyL(i,Jend1)*((1.0-Slambda(i,Jend1))*KapyL(i,Jend1)+Slambda(i,Jend1))*HayL(i,Jend1)*HayL(i,Jend1)
    enddo
# if defined (PARALLEL)
    endif
# endif

    do j = Jbeg-1,Jend+1
    do i = Ibeg-1,Iend+1
      if(Maska(i,j)==0) then
        if(Maska(i-1,j)==1) then
          Eax(i,j) = Zero
          Fax(i,j) = 0.5*Grav_Lz(i-1,j)*((1.0-Slambda(i-1,j))*Kap(i-1,j)+  &
                     Slambda(i-1,j))*HaxL(i,j)*HaxL(i,j)
          Gax(i,j) = Zero
        elseif(Maska(i+1,j)==1) then
          Eax(i+1,j) = Zero
          Fax(i+1,j) = 0.5*Grav_Lz(i+1,j)*((1.0-Slambda(i+1,j))*Kap(i+1,j)+  &
                       Slambda(i+1,j))*HaxR(i+1,j)*HaxR(i+1,j)
          Gax(i+1,j) = Zero
        elseif(Maska(i,j-1)==1) then
          Eay(i,j) = Zero
          Fay(i,j) = Zero
          Gay(i,j) = 0.5*Grav_Lz(i,j-1)*((1.0-Slambda(i,j-1))*Kap(i,j-1)+  &
                     Slambda(i,j-1))*HayL(i,j)*HayL(i,j)
        elseif(Maska(i,j+1)==1) then
          Eay(i,j+1) = Zero
          Fay(i,j+1) = Zero
          Gay(i,j+1) = 0.5*Grav_Lz(i,j+1)*((1.0-Slambda(i,j+1))*Kap(i,j+1)+  &
                       Slambda(i,j+1))*HayR(i,j+1)*HayR(i,j+1)
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
    call delxFun_2D(Grav_Lz,DelxGravLz)

    call delyFun_2D(Ha,DelyHa)
    call delyFun_2D(Ua,DelyUa)
    call delyFun_2D(Va,DelyVa)
    call delyFun_2D(HUa,DelyHUa)
    call delyFun_2D(HVa,DelyHVa)
    call delyFun_2D(Kap,DelyKap)
    call delyFun_2D(Grav_Lz,DelyGravLz)

    end subroutine delxyFun_ll
!
!---------------------------------------------------------------------------------------------------------
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
    call construct_2D_x(Grav_Lz,DelxGravLz,GravxL,GravxR)

    do j = Jbeg,Jend
    do i = Ibeg,Iend1
      EaxL(i,j) = HUaxL(i,j)
      EaxR(i,j) = HUaxR(i,j)
      FaxL(i,j) = HUaxL(i,j)*UaxL(i,j)+0.5*GravxL(i,j)*((1.0-Slambda(i,j))*KapxL(i,j)+ &
                  Slambda(i,j))*HaxL(i,j)*HaxL(i,j)
      FaxR(i,j) = HUaxR(i,j)*UaxR(i,j)+0.5*GravxR(i,j)*((1.0-Slambda(i,j))*KapxR(i,j)+ &
                  Slambda(i,j))*HaxR(i,j)*HaxR(i,j)
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
    call construct_2D_y(Grav_Lz,DelyGravLz,GravyL,GravyR)

    do j = Jbeg,Jend1
    do i = Ibeg,Iend
      EayL(i,j) = HVayL(i,j)
      EayR(i,j) = HVayR(i,j)
      FayL(i,j) = HayL(i,j)*UayL(i,j)*VayL(i,j)
      FayR(i,j) = HayR(i,j)*UayR(i,j)*VayR(i,j)
      GayL(i,j) = HVayL(i,j)*VayL(i,j)+0.5*GravyL(i,j)*((1.0-Slambda(i,j))*KapyL(i,j)+ &
                  Slambda(i,j))*HayL(i,j)*HayL(i,j)
      GayR(i,j) = HVayR(i,j)*VayR(i,j)+0.5*GravyR(i,j)*((1.0-Slambda(i,j))*KapyR(i,j)+ &
                  Slambda(i,j))*HayR(i,j)*HayR(i,j)
    enddo
    enddo

    end subroutine construction_ll
!
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
!
    subroutine wave_speed_ll
    use global
    implicit none
    integer :: i,j
    real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S

    do j = Jbeg,Jend
    do i = Ibeg,Iend1
      SQR_PHI_L = sqrt(GravxL(i,j)*((1.0-Slambda(i,j))*KapxL(i,j)+ &
                  Slambda(i,j))*abs(HaxL(i,j)))
      SQR_PHI_R = sqrt(GravxR(i,j)*((1.0-Slambda(i,j))*KapxR(i,j)+ &
                  Slambda(i,j))*abs(HaxR(i,j)))
      SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UaxL(i,j)-UaxR(i,j))
      U_S = 0.5*(UaxL(i,j)+UaxR(i,j))+SQR_PHI_L-SQR_PHI_R
      SaxL(i,j) = min(UaxL(i,j)-SQR_PHI_L,U_S-SQR_PHI_S)
      SaxR(i,j) = max(UaxR(i,j)+SQR_PHI_R,U_S+SQR_PHI_S)
      SaxS(i,j) = U_S
    enddo
    enddo

    do j = Jbeg,Jend1
    do i = Ibeg,Iend
      SQR_PHI_L = sqrt(GravyL(i,j)*((1.0-Slambda(i,j))*KapyL(i,j)+ &
                  Slambda(i,j))*abs(HayL(i,j)))
      SQR_PHI_R = sqrt(GravyR(i,j)*((1.0-Slambda(i,j))*KapyR(i,j)+ &
                  Slambda(i,j))*abs(HayR(i,j)))
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
