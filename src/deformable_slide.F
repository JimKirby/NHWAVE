! added by Cheng for deformable slide
# if defined (DEFORMABLESLIDE)
    subroutine update_maskgs
!------------------------------------------------------  
!    This subroutine is used to update maskgs for wetting-drying
!    Called by                                                
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	integer,dimension(:,:),allocatable :: Masktmp
	 
    allocate(Masktmp(Mloc,Nloc))
    Masktmp = Maskgs
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Maskgs(i,j)==0) then
        if(Maskgs(i-1,j)==1.and.Hgs(i-1,j)<Hgs(i,j).and. &
		          Dgs(i-1,j)-SLIDE_MINTHICK>SLIDE_MINTHICK) Masktmp(i,j)=1
        if(Maskgs(i+1,j)==1.and.Hgs(i+1,j)<Hgs(i,j).and. &
		          Dgs(i+1,j)-SLIDE_MINTHICK>SLIDE_MINTHICK) Masktmp(i,j)=1
        if(Maskgs(i,j-1)==1.and.Hgs(i,j-1)<Hgs(i,j).and. &
		          Dgs(i,j-1)-SLIDE_MINTHICK>SLIDE_MINTHICK) Masktmp(i,j)=1
        if(Maskgs(i,j+1)==1.and.Hgs(i,j+1)<Hgs(i,j).and. &
		          Dgs(i,j+1)-SLIDE_MINTHICK>SLIDE_MINTHICK) Masktmp(i,j)=1
      else
        ! drying (wet->dry)               
        if(Dgs(i,j)-SLIDE_MINTHICK<=Small) then
          Masktmp(i,j) = 0
          Dgs(i,j) = SLIDE_MINTHICK
        endif
      endif
	  Hgs(i,j) = Hc0(i,j)-Dgs(i,j)
    enddo
    enddo
	Maskgs = Masktmp
	
	!boundary conditions
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         Maskgs(i,j) = Maskgs(Ibeg+Nghost-i,j)
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
         Maskgs(Iend+i,j) = Maskgs(Iend-i+1,j)
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
         Maskgs(i,j) = Maskgs(i,Jbeg+Nghost-j)
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
         Maskgs(i,Jend+j) = Maskgs(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     ! collect maskgs into ghost cells  
     call phi_int_exch(Maskgs)
# endif
	 
	 deallocate(Masktmp)

    end subroutine update_maskgs
	
	
     subroutine update_vars_gs
!------------------------------------------------------ 
!    This subroutine is used to save deformable slide variables at 
!    last time step
!    Called by   
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang 
!----------------------------------------------------- 
     use global        
     implicit none
	 integer :: i,j

     Dgs0 = Dgs
     Ugs0 = Ugs
     Vgs0 = Vgs
	 Wgs0 = Wgs
	 Wgs00 = Wgs0
     DUgs0 = DUgs
     DVgs0 = DVgs
	 DWgs0 = DWgs
	 Rhogs0=Rhogs
	 if(trim(RHEO_OPT)=='GRANULAR') then
!	    Rhogs=Zero
		PhiInt=Zero;PhiBed=Zero
	    do j = Jbeg,Jend
        do i = Ibeg,Iend
		  if(Mask(i,j)==1) then
!	        Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC + Rho0*(1.0-SLIDE_CONC)
	        PhiInt(i,j) = PhiInt_F
			PhiBed(i,j) = PhiBed_F
		  else
!		    Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC
	        PhiInt(i,j) = PhiInt_A
			PhiBed(i,j) = PhiBed_A
		  endif
		enddo
		enddo
!		call phi_2D_coll(Rhogs)
		call phi_2D_coll(PhiInt)
		call phi_2D_coll(PhiBed)
	 endif

     end subroutine update_vars_gs
	 
    subroutine flux_coeff_gs
!------------------------------------------------------  
!    This subroutine is used to update coefficients for flux
!    Called by                                                
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	real(SP) :: TMP1,TMP2,DispT
    real(SP),dimension(:,:),allocatable :: TMP3
	 
    allocate(TMP3(Mloc,Nloc))
	 
    if(trim(RHEO_OPT)=='VISCOUS') then
	  Cxgs = Grav
      Cygs = Grav
	elseif(trim(RHEO_OPT)=='GRANULAR') then
	  !Get local coordinate
	  call transform_gs
	  ! Obtain DwDt and Qsgs in non-hydrostatic model
	  if (NON_HYDRO_SLD.and.REDU_GRAV_SLD) then
	    ! Obtain dynamic pressure at slide surface
	    if (NON_HYDRO_UP) then
          do j = 1,Nloc
          do i = 1,Mloc
            Qtgs(i,j) = P(i,j,kbeg)*Mask(i,j)! dynamic pressure over the slide
          enddo
          enddo
		else
		  Qtgs=Zero
	    endif
		TMP3=Zero
        do j = Jbeg,Jend
        do i = Ibeg,Iend
		  TMP3(i,j)=(Qbgs(i,j)-Qtgs(i,j))/Rhogs(i,j)/Dgs(i,j)
		  TMP1=SLIDE_CONC*GRAIN_DENSITY/Rhogs(i,j)*Grav
  		  TMP3(i,j)=max(TMP3(i,j),-TMP1)
  		  TMP3(i,j)=min(TMP3(i,j),TMP1)
		  !Qsgs
		  Qsgs(i,j)=TMP3(i,j)*Rhogs(i,j)*Dgs(i,j)*SLIDE_LAMBDA
		  !DwDt(including fluid non-hydrostatic pressure and shear stress in z momentum eq)
  		  DwDt(i,j)=TMP3(i,j)-(TMP3(i,j)*Rhogs(i,j)*Dgs(i,j)-Qsgs(i,j))/Rhogs(i,j)/Dgs(i,j) 
!		            +(-SrctgsZ(i,j)*Rhogs(i,j)-Tezz(i,j)+Pss(i,j))/Dgs(i,j)/Rhogs(i,j)
        enddo
        enddo
		call phi_2D_coll(Qsgs)
	    call phi_2D_coll(DwDt)
	  endif
	  
	  ! Obtain coefficients of flux
      do j = Jbeg,Jend
      do i = Ibeg,Iend
  	    TMP1=SLIDE_CONC*GRAIN_DENSITY/Rhogs(i,j)*Grav
  		TMP2=(1.0-SLIDE_CONC*GRAIN_DENSITY/Rhogs(i,j))*Grav
		Cxgs(i,j)=(Cxx(i,j)/max(Czz(i,j),Small)*TMP1+TMP2)*Maskgs(i,j)
		Cygs(i,j)=(Cyy(i,j)/max(Czz(i,j),Small)*TMP1+TMP2)*Maskgs(i,j)
      enddo
      enddo
      do j = Jbeg,Jend
      do i = Ibeg,Iend
	    ! ensure the hyperbolic property in x
        TMP1=Cxx(i,j)/max(Czz(i,j),Small)*DwDt(i,j)-Qsgs(i,j)/Rhogs(i,j)/Dgs(i,j)
        if (TMP1<=-Cxgs(i,j)) then
          TMP1=-Cxgs(i,j)
        elseif (TMP1>=Cxgs(i,j)) then
          TMP1=Cxgs(i,j)
        endif
        Cxgs(i,j)=Cxgs(i,j)+TMP1
		! ensure the hyperbolic property in y
		TMP2=Cyy(i,j)/max(Czz(i,j),Small)*DwDt(i,j)-Qsgs(i,j)/Rhogs(i,j)/Dgs(i,j)
        if (TMP2<=-Cygs(i,j)) then
          TMP2=-Cygs(i,j)
        elseif (TMP2>=Cygs(i,j)) then
          TMP2=Cygs(i,j)
        endif
        Cygs(i,j)=Cygs(i,j)+TMP2
      enddo
      enddo
	endif
	
	! Dispersion correction term in non-hydrostatic model
    if (NON_HYDRO_SLD.and.DISP_CORR_SLD) then
	  Call DelxH_gs
	  call delxFun_gs_2D(Hgs,DelxHgs)
	  call delyFun_gs_2D(Hgs,DelyHgs)
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DispT=DelxH0(i,j)*DelxHgs(i,j)/3.0*Grav+DelyH0(i,j)*DelyHgs(i,j)/3.0*Grav
		! ensure the hyperbolic property in x
        if (DispT<=-Cxgs(i,j)) then
          DispT=-Cxgs(i,j)
        elseif (DispT>=Cxgs(i,j)) then
          DispT=Cxgs(i,j)
        endif
        Cxgs(i,j)=Cxgs(i,j)+DispT
		! ensure the hyperbolic property in y
        if (DispT<=-Cygs(i,j)) then
          DispT=-Cygs(i,j)
        elseif (DispT>=Cygs(i,j)) then
          DispT=Cygs(i,j)
        endif
        Cygs(i,j)=Cygs(i,j)+DispT
      enddo
      enddo
    endif

	call phi_2D_coll(Cxgs)
	call phi_2D_coll(Cygs)

    deallocate(TMP3)

    end subroutine flux_coeff_gs
	
    subroutine transform_gs
!------------------------------------------------------  
!    This subroutine is used to construct local coordinate
!    Called by                                                
!       flux_coeff_gs
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	real(SP),dimension(:,:),allocatable :: C11,C12,C13,C21,C22,C23,C31,C32,C33, &
	                                       DelxUgsP,DelyUgsP,DelxVgsP,DelyVgsP, &
										   DelxWgsP,DelyWgsP
	real(SP) :: A,B,DxhDyh,Udhhc0,crit,Int_factor,DugsDx,DvgsDy,DwgsDz,DIV,Kxa,Kxp
	real(SP) :: UvecPC(3),UvecPL(3),UvecPR(3),UvecPF(3),UvecPB(3)
	 
    allocate(C11(Mloc,Nloc))
	allocate(C12(Mloc,Nloc))
	allocate(C13(Mloc,Nloc))
    allocate(C21(Mloc,Nloc))
	allocate(C22(Mloc,Nloc))
	allocate(C23(Mloc,Nloc))
    allocate(C31(Mloc,Nloc))
	allocate(C32(Mloc,Nloc))
	allocate(C33(Mloc,Nloc))
    allocate(DelxUgsP(Mloc,Nloc))
	allocate(DelyUgsP(Mloc,Nloc))
	allocate(DelxVgsP(Mloc,Nloc))
    allocate(DelyVgsP(Mloc,Nloc))
	allocate(DelxWgsP(Mloc,Nloc))
	allocate(DelyWgsP(Mloc,Nloc))
	
    C11=Zero;C12=Zero;C13=Zero;C21=Zero;C22=Zero;C23=Zero;C31=Zero;C32=Zero;C33=Zero
	DelxUgsP=Zero;DelyUgsP=Zero;DelxVgsP=Zero;DelyVgsP=Zero;DelxWgsP=Zero;DelyWgsP=Zero
	call DelxH_gs
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ! generate transformation matrix
       A=sqrt(1.+DelxH0(i,j)**2+DelyH0(i,j)**2);
       if(sqrt(Ugs(i,j)**2+Vgs(i,j)**2+Wgs(i,j)**2)<1.e-8) then
         if(abs(DelxH0(i,j))<1.e-8.and.abs(DelyH0(i,j))<1.e-8) then
           B=sqrt(1.+DelxH0(i,j)**2)
           C11(i,j)=1./max(B,Small)
           C21(i,j)=Zero
           C31(i,j)=-DelxH0(i,j)/max(B,Small)
           C12(i,j)=-DelxH0(i,j)*DelyH0(i,j)/max(A,Small)/max(B,Small)
           C22(i,j)=(1.+DelxH0(i,j)**2)/max(A,Small)/max(B,Small)
           C32(i,j)=-DelyH0(i,j)/max(A,Small)/max(B,Small)
         else 
           DxhDyh=DelxH0(i,j)**2+DelyH0(i,j)**2
           B=sqrt(DelxH0(i,j)**2+DelyH0(i,j)**2+DxhDyh**2)
           C11(i,j)=DelxH0(i,j)/max(B,Small)
           C21(i,j)=DelyH0(i,j)/max(B,Small)
           C31(i,j)=-DxhDyh/max(B,Small)
           C12(i,j)=-DelyH0(i,j)*A/max(B,Small)
           C22(i,j)=DelxH0(i,j)*A/max(B,Small)
           C32(i,j)=Zero
         endif  
       else
         Udhhc0=Ugs(i,j)*DelxH0(i,j)+Vgs(i,j)*DelyH0(i,j)+Wgs(i,j)
         B=sqrt((Ugs(i,j)-Udhhc0*DelxH0(i,j)/max(A,Small)**2)**2 &
		       +(Vgs(i,j)-Udhhc0*DelyH0(i,j)/max(A,Small)**2)**2 &
			   +(Wgs(i,j)-Udhhc0/max(A,Small)**2)**2)
         C11(i,j)=(Ugs(i,j)-Udhhc0*DelxH0(i,j)/max(A,Small)**2)/max(B,Small)
         C21(i,j)=(Vgs(i,j)-Udhhc0*DelyH0(i,j)/max(A,Small)**2)/max(B,Small)
         C31(i,j)=(Wgs(i,j)-Udhhc0/max(A,Small)**2)/max(B,Small)
         C12(i,j)=(DelyH0(i,j)*Wgs(i,j)-Vgs(i,j))/max(A,Small)/max(B,Small)
         C22(i,j)=(-DelxH0(i,j)*Wgs(i,j)+Ugs(i,j))/max(A,Small)/max(B,Small)
         C32(i,j)=(DelxH0(i,j)*Vgs(i,j)-DelyH0(i,j)*Ugs(i,j))/max(A,Small)/max(B,Small)
       endif
	   C13(i,j)=DelxH0(i,j)/max(A,Small)
       C23(i,j)=DelyH0(i,j)/max(A,Small)
       C33(i,j)=1./max(A,Small)
       ! local velocity and derivative
	   UvecPC(1)=C11(i,j)*Ugs(i,j)+C21(i,j)*Vgs(i,j)+C31(i,j)*Wgs(i,j)
	   UvecPC(2)=C12(i,j)*Ugs(i,j)+C22(i,j)*Vgs(i,j)+C32(i,j)*Wgs(i,j)
	   UvecPC(3)=C13(i,j)*Ugs(i,j)+C23(i,j)*Vgs(i,j)+C33(i,j)*Wgs(i,j)
       UgsP(i,j)=UvecPC(1);VgsP(i,j)=UvecPC(2);WgsP(i,j)=UvecPC(3)
	   UvecPL(1)=C11(i,j)*Ugs(i-1,j)+C21(i,j)*Vgs(i-1,j)+C31(i,j)*Wgs(i-1,j)
	   UvecPL(2)=C12(i,j)*Ugs(i-1,j)+C22(i,j)*Vgs(i-1,j)+C32(i,j)*Wgs(i-1,j)
	   UvecPL(3)=C13(i,j)*Ugs(i-1,j)+C23(i,j)*Vgs(i-1,j)+C33(i,j)*Wgs(i-1,j)
	   UvecPR(1)=C11(i,j)*Ugs(i+1,j)+C21(i,j)*Vgs(i+1,j)+C31(i,j)*Wgs(i+1,j)
	   UvecPR(2)=C12(i,j)*Ugs(i+1,j)+C22(i,j)*Vgs(i+1,j)+C32(i,j)*Wgs(i+1,j)
	   UvecPR(3)=C13(i,j)*Ugs(i+1,j)+C23(i,j)*Vgs(i+1,j)+C33(i,j)*Wgs(i+1,j)
	   UvecPF(1)=C11(i,j)*Ugs(i,j-1)+C21(i,j)*Vgs(i,j-1)+C31(i,j)*Wgs(i,j-1)
	   UvecPF(2)=C12(i,j)*Ugs(i,j-1)+C22(i,j)*Vgs(i,j-1)+C32(i,j)*Wgs(i,j-1)
	   UvecPF(3)=C13(i,j)*Ugs(i,j-1)+C23(i,j)*Vgs(i,j-1)+C33(i,j)*Wgs(i,j-1)
	   UvecPB(1)=C11(i,j)*Ugs(i,j+1)+C21(i,j)*Vgs(i,j+1)+C31(i,j)*Wgs(i,j+1)
	   UvecPB(2)=C12(i,j)*Ugs(i,j+1)+C22(i,j)*Vgs(i,j+1)+C32(i,j)*Wgs(i,j+1)
	   UvecPB(3)=C13(i,j)*Ugs(i,j+1)+C23(i,j)*Vgs(i,j+1)+C33(i,j)*Wgs(i,j+1)
	   if(Maskgs(i,j)==0) then
	     DelxUgsP(i,j)=Zero;DelyUgsP(i,j)=Zero
		 DelxVgsP(i,j)=Zero;DelyVgsP(i,j)=Zero
		 DelxWgsP(i,j)=Zero;DelyWgsP(i,j)=Zero
	   else
	     call delFunP_gs_2D(UvecPC,UvecPL,UvecPR,UvecPF,UvecPB,DelxUgsP(i,j),DelyUgsP(i,j), &
	                        DelxVgsP(i,j),DelyVgsP(i,j),DelxWgsP(i,j),DelyWgsP(i,j))
	   endif
     enddo
     enddo
     call phi_2D_coll(DelxUgsP)
	 call phi_2D_coll(DelyUgsP)
     call phi_2D_coll(DelxVgsP)
	 call phi_2D_coll(DelyVgsP)
	 call phi_2D_coll(DelxWgsP)
	 call phi_2D_coll(DelyWgsP)
	 call delxFun_gs_2D(Ugs,DelxUgs)
	 call delyFun_gs_2D(Vgs,DelyVgs)
     crit=0.2 !for linear intepreting between Ka and Kp 
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==1) then
         DIV=DelxUgs(i,j)+DelyVgs(i,j)
         DugsDx = C11(i,j)*DelxUgsP(i,j)+C21(i,j)*DelyUgsP(i,j)
         DvgsDy = C12(i,j)*DelxVgsP(i,j)+C22(i,j)*DelyVgsP(i,j)
         DwgsDz = C13(i,j)*DelxWgsP(i,j)+C23(i,j)*DelyWgsP(i,j)
         if(PhiInt(i,j)>=PhiBed(i,j)) then
           Kxa=2./cos(PhiInt(i,j))**2*(1.-sqrt(1.-cos(PhiInt(i,j))**2/cos(PhiBed(i,j))**2))-1.
           Kxp=2./cos(PhiInt(i,j))**2*(1.+sqrt(1.-cos(PhiInt(i,j))**2/cos(PhiBed(i,j))**2))-1.
         else
           Kxa=(1.+sin(PhiInt(i,j))**2)/(1.-sin(PhiInt(i,j))**2)
           Kxp=Kxa
         endif
		 ! old rheology kx=ky
     !     if(abs(UgsP(i,j))>1.e-8) then
     !       if(DIV>=crit) then
     !         Kx(i,j)=Kxa
     !         Ky(i,j)=Kxa
     !       elseif(DIV<=-crit) then
     !         Kx(i,j)=Kxp
     !         Ky(i,j)=Kxp
     !       else
     !         Int_factor = (crit-abs(DIV))/(crit-1.e-8)
     !         if(DIV>1.e-8.and.DIV<crit) then
     !           Kx(i,j)=Int_factor+(1.0-Int_factor)*Kxa
     !           Ky(i,j)=Kx(i,j)
     !         elseif(DIV<-1.e-8.and.DIV>-crit) then
     !           Kx(i,j)=Int_factor+(1.0-Int_factor)*Kxp
     !           Ky(i,j)=Kx(i,j)
     !         else
     !           Kx(i,j)=1.0
     !           Ky(i,j)=1.0
     !         endif
     !       endif
     !     else
     !       Kx(i,j)=1.0
     !       Ky(i,j)=1.0
     !     endif
	     ! new rheology kx~=ky
         if(abs(UgsP(i,j))>1.e-8) then
           if(DugsDx>=crit) then
             Kx(i,j)=Kxa
           elseif(DugsDx<=-crit) then
             Kx(i,j)=Kxp
           else
             Int_factor = (crit-abs(DugsDx))/(crit-1.e-8)
             if(DugsDx>1.e-8.and.DugsDx<crit) then
               Kx(i,j)=Int_factor+(1.0-Int_factor)*Kxa
             elseif(DugsDx<-1.e-8.and.DugsDx>-crit) then
               Kx(i,j)=Int_factor+(1.0-Int_factor)*Kxp
             else
               Kx(i,j)=1.0
             endif
           endif
           if(DvgsDy>crit) then
             Ky(i,j)=0.5*(Kx(i,j)+1.0-sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed(i,j))**2))
           elseif(DvgsDy<-crit) then
             Ky(i,j)=0.5*(Kx(i,j)+1.0+sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed(i,j))**2))
           else
             Int_factor = (crit-abs(DvgsDy))/(crit-1.e-8)
             if(DvgsDy>1.e-8.and.DvgsDy<crit) then
               Ky(i,j)=Int_factor+(1.0-Int_factor) &
			           *(0.5*(Kx(i,j)+1.0-sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed(i,j))**2)))
             elseif(DvgsDy<-1.e-8.and.DvgsDy>-crit) then
               Ky(i,j)=Int_factor+(1.0-Int_factor) &
			           *(0.5*(Kx(i,j)+1.0+sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed(i,j))**2)))
             else
               Ky(i,j)=1.0
             endif
           endif 
         else
           Kx(i,j)=1.0
           Ky(i,j)=1.0
         endif
         Cxx(i,j)=C11(i,j)*C11(i,j)*Kx(i,j)+C12(i,j)*C12(i,j)*Ky(i,j)+C13(i,j)*C13(i,j)
         Cyx(i,j)=C21(i,j)*C11(i,j)*Kx(i,j)+C22(i,j)*C12(i,j)*Ky(i,j)+C23(i,j)*C13(i,j)
         Cyy(i,j)=C21(i,j)*C21(i,j)*Kx(i,j)+C22(i,j)*C22(i,j)*Ky(i,j)+C23(i,j)*C23(i,j)
         Cxz(i,j)=C11(i,j)*C31(i,j)*Kx(i,j)+C12(i,j)*C32(i,j)*Ky(i,j)+C13(i,j)*C33(i,j)
         Cyz(i,j)=C21(i,j)*C31(i,j)*Kx(i,j)+C22(i,j)*C32(i,j)*Ky(i,j)+C23(i,j)*C33(i,j)
         Czz(i,j)=C31(i,j)*C31(i,j)*Kx(i,j)+C32(i,j)*C32(i,j)*Ky(i,j)+C33(i,j)*C33(i,j)
         Czz0(i,j)=Czz(i,j)-(C33(i,j)*C31(i,j)+C31(i,j)*C33(i,j))*tan(PhiBed(i,j))
         if(abs(UgsP(i,j))>1.e-8) then
           Cxx(i,j)=Cxx(i,j)-(C13(i,j)*C11(i,j)+C11(i,j)*C13(i,j))*tan(PhiBed(i,j))
           Cyx(i,j)=Cyx(i,j)-(C23(i,j)*C11(i,j)+C21(i,j)*C13(i,j))*tan(PhiBed(i,j))
           Cyy(i,j)=Cyy(i,j)-(C23(i,j)*C21(i,j)+C21(i,j)*C23(i,j))*tan(PhiBed(i,j))
           Cxz(i,j)=Cxz(i,j)-(C13(i,j)*C31(i,j)+C11(i,j)*C33(i,j))*tan(PhiBed(i,j))
           Cyz(i,j)=Cyz(i,j)-(C23(i,j)*C31(i,j)+C21(i,j)*C33(i,j))*tan(PhiBed(i,j))
           Czz(i,j)=Czz(i,j)-(C33(i,j)*C31(i,j)+C31(i,j)*C33(i,j))*tan(PhiBed(i,j))
         endif
         Cxy(i,j)=Cyx(i,j)
         Czx(i,j)=Cxz(i,j)
         Czy(i,j)=Cyz(i,j)
       endif
     enddo
     enddo
     call phi_2D_coll(Cxx)
	 call phi_2D_coll(Cxy)
	 call phi_2D_coll(Cxz)
     call phi_2D_coll(Cyx)
	 call phi_2D_coll(Cyy)
	 call phi_2D_coll(Cyz)
     call phi_2D_coll(Czx)
	 call phi_2D_coll(Czy)
	 call phi_2D_coll(Czz)
     call phi_2D_coll(Czz0)
	 
    deallocate(C11)
	deallocate(C12)
	deallocate(C13)
    deallocate(C21)
	deallocate(C22)
	deallocate(C23)
    deallocate(C31)
	deallocate(C32)
	deallocate(C33)
    deallocate(DelxUgsP)
	deallocate(DelyUgsP)
	deallocate(DelxVgsP)
    deallocate(DelyVgsP)
	deallocate(DelxWgsP)
	deallocate(DelyWgsP)

    end subroutine transform_gs
	
     subroutine delFunP_gs_2D(UvecPC,UvecPL,UvecPR,UvecPF,UvecPB,DXUP,DYUP,DXVP,DYVP,DXWP,DYWP)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       transform_gs
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,dy
     implicit none
     real(SP),intent(in),dimension(3) :: UvecPC,UvecPL,UvecPR,UvecPF,UvecPB
     real(SP),intent(out) :: DXUP,DYUP,DXVP,DYVP,DXWP,DYWP
     real(SP) :: TMP1,TMP2,LIMITER_GS
    
     ! pu'/px'
     TMP1 = (UvecPR(1)-UvecPC(1))/dx
     TMP2 = (UvecPC(1)-UvecPL(1))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXUP = Zero
     else
       DXUP = LIMITER_GS(TMP1,TMP2)
     endif
     ! pu'/py'
     TMP1 = (UvecPB(1)-UvecPC(1))/dy
     TMP2 = (UvecPC(1)-UvecPF(1))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYUP = Zero
     else
       DYUP = LIMITER_GS(TMP1,TMP2)
     endif
     ! pv'/px'
     TMP1 = (UvecPR(2)-UvecPC(2))/dx
     TMP2 = (UvecPC(2)-UvecPL(2))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXVP = Zero
     else
       DXVP = LIMITER_GS(TMP1,TMP2)
     endif
     ! pv'/py'
     TMP1 = (UvecPB(2)-UvecPC(2))/dy
     TMP2 = (UvecPC(2)-UvecPF(2))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYVP = Zero
     else
       DYVP = LIMITER_GS(TMP1,TMP2)
     endif
     !pw'/px'
     TMP1 = (UvecPR(3)-UvecPC(3))/dx
     TMP2 = (UvecPC(3)-UvecPL(3))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXWP = Zero
     else
       DXWP = LIMITER_GS(TMP1,TMP2)
     endif
     !pw'/py'
     TMP1 = (UvecPB(3)-UvecPC(3))/dy
     TMP2 = (UvecPC(3)-UvecPF(3))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYWP = Zero
     else
       DYWP = LIMITER_GS(TMP1,TMP2)
     endif

     return
     end subroutine delFunP_gs_2D
	 
    subroutine DelxH_gs
!------------------------------------------------------  
!    This subroutine is used to obtain derivatives and
!    cell face flux of bathymetry
!    Called by                                                
!       transform_gs
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	real(SP) :: TMP1,TMP2
	 
     ! reconstruct base depth below the slide at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       H0fx(i,j) = 0.5*(Hc0(i-1,j)+Hc0(i,j))
     enddo
     H0fx(1,j) = Hc0(1,j)
     H0fx(Mloc1,j) = Hc0(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       H0fy(i,j) = 0.5*(Hc0(i,j-1)+Hc0(i,j))
     enddo
     H0fy(i,1) = Hc0(i,1)
     H0fy(i,Nloc1) = Hc0(i,Nloc)
     enddo

     ! derivatives of base depth below the slide at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH0(i,j) = (H0fx(i+1,j)-H0fx(i,j))/dx
       DelyH0(i,j) = (H0fy(i,j+1)-H0fy(i,j))/dy
     enddo
     enddo
    
    ! second derivatives of depth at cell center
    do i = 2,Mloc-1
    do j = 1,Nloc
       TMP1 = (Hc0(i+1,j)-Hc0(i,j))/dx
       TMP2 = (Hc0(i,j)-Hc0(i-1,j))/dx
       Delx2H0(i,j) = (TMP1-TMP2)/dx
    enddo
    enddo
    do j = 1,Nloc
      Delx2H0(1,j) = Delx2H0(2,j)
      Delx2H0(Mloc,j) = Delx2H0(Mloc-1,j)
    enddo
    do i = 1,Mloc
    do j = 2,Nloc-1
       TMP1 = (Hc0(i,j+1)-Hc0(i,j))/dy
       TMP2 = (Hc0(i,j)-Hc0(i,j-1))/dy
       Dely2H0(i,j) = (TMP1-TMP2)/dy
    enddo
    enddo
    do i = 1,Mloc
      Dely2H0(i,1) = Dely2H0(i,2)
      Dely2H0(i,Nloc) = Dely2H0(i,Nloc-1)
    enddo  
	
# if defined (PARALLEL)
	call phi_2D_exch(DelxH0)
	call phi_2D_exch(DelyH0)
	call phi_2D_exch(Delx2H0)
	call phi_2D_exch(Dely2H0)
# endif

    end subroutine DelxH_gs
	
     subroutine fluxes_gs
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces for deformable slide
!    Called by
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------------
     use global
     implicit none

     ! second order construction
     call delxyzFun_gs
     call construction_gs

     ! calculate wave speed
     call wave_speed_gs

     ! calculate fluxes at faces
     call fluxes_at_faces_HLL_gs

     ! impose boundary conditions
     call flux_bc_gs

     end subroutine fluxes_gs

     subroutine delxyzFun_gs
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes_gs 
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global
     implicit none
     
	 call delxFun_gs_2D(Dgs,DelxDgs)
	 call delxFun_gs_2D(Hgs,DelxHgs)
     call delxFun_gs_2D(Ugs,DelxUgs)
     call delxFun_gs_2D(Vgs,DelxVgs)
	 call delxFun_gs_2D(Wgs,DelxWgs)
     call delxFun_gs_2D(DUgs,DelxDUgs)
     call delxFun_gs_2D(DVgs,DelxDVgs)
	 call delxFun_gs_2D(DWgs,DelxDWgs)

	 call delyFun_gs_2D(Dgs,DelyDgs)
	 call delyFun_gs_2D(Hgs,DelyHgs)
     call delyFun_gs_2D(Ugs,DelyUgs)
     call delyFun_gs_2D(Vgs,DelyVgs)
	 call delyFun_gs_2D(Wgs,DelyWgs)
     call delyFun_gs_2D(DUgs,DelyDUgs)
     call delyFun_gs_2D(DVgs,DelyDVgs)
	 call delyFun_gs_2D(DWgs,DelyDWgs)

     end subroutine delxyzFun_gs
	 
     subroutine delxFun_gs_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x with limiter
!    Called by
!       delxyzFun_gs
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER_GS
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Maskgs(i,j)==0) then
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER_GS(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_gs_2D
	 
     subroutine delyFun_gs_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y with limiter
!    Called by 
!       delxyzFun_gs
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER_GS
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Maskgs(i,j)==0) then 
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER_GS(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_gs_2D
	 
!------------------------------------------
!    Slope limiter used in landslide model
!    Called by 
!       delxFun_gs_2D,delyFun_gs_2D
!    Last update: 26/06/2019, Cheng Zhang
!-----------------------------------------
     FUNCTION LIMITER_GS(A,B)
     use global, only: SP,Zero,One,Small
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: LIMITER_GS

!     ! minmod limiter
      if (abs(A)<=abs(B).and.A*B>Zero) then
        LIMITER_GS=A
      elseif (abs(A)>abs(B).and.A*B>Zero) then
        LIMITER_GS=B
      else
        LIMITER_GS=Zero
      endif

!     ! van Leer limiter
!     LIMITER_GS=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))

     RETURN
     END FUNCTION LIMITER_GS
	 
	 subroutine construction_gs
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes_gs
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------
     use global
     implicit none
     integer :: i,j

	 call construct_2D_x(Dgs,DelxDgs,DgsxL,DgsxR)
	 call construct_2D_x(Hgs,DelxHgs,HSgsxL,HSgsxR)
     call construct_2D_x(Ugs,DelxUgs,UgsxL,UgsxR)
     call construct_2D_x(Vgs,DelxVgs,VgsxL,VgsxR)
	 call construct_2D_x(Wgs,DelxWgs,WgsxL,WgsxR)
     call construct_2D_x(DUgs,DelxDUgs,DUgsxL,DUgsxR)
     call construct_2D_x(DVgs,DelxDVgs,DVgsxL,DVgsxR)
	 call construct_2D_x(DWgs,DelxDWgs,DWgsxL,DWgsxR)

     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       Cxgsx(i,j)=0.5*(Cxgs(i,j)+Cxgs(i-1,j))
     enddo
     enddo
	 
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if (Maskgs(i,j)==1) then
	     if (Maskgs(i+1,j)==0.and.Maskgs(i-1,j)==1) then
           DgsxL(i+1,j)=Dgs(i,j)
           HSgsxL(i+1,j)=Hgs(i,j)
           DUgsxL(i+1,j)=Zero
           DVgsxL(i+1,j)=Zero
           DWgsxL(i+1,j)=Zero
           UgsxL(i+1,j)=Zero
           VgsxL(i+1,j)=Zero
           WgsxL(i+1,j)=Zero
           Cxgsx(i+1,j)=Cxgs(i,j)
		 elseif (Maskgs(i+1,j)==1.and.Maskgs(i-1,j)==0) then
           DgsxR(i,j)=Dgs(i,j)
           HSgsxR(i,j)=Hgs(i,j)
           DUgsxR(i,j)=Zero
           DVgsxR(i,j)=Zero
           DWgsxR(i,j)=Zero
           UgsxR(i,j)=Zero
           VgsxR(i,j)=Zero
           WgsxR(i,j)=Zero
           Cxgsx(i,j)=Cxgs(i,j)
		 elseif (Maskgs(i+1,j)==0.and.Maskgs(i-1,j)==0) then
           DgsxL(i+1,j)=Dgs(i,j)
           HSgsxL(i+1,j)=Hgs(i,j)
           DUgsxL(i+1,j)=Zero
           DVgsxL(i+1,j)=Zero
           DWgsxL(i+1,j)=Zero
           UgsxL(i+1,j)=Zero
           VgsxL(i+1,j)=Zero
           WgsxL(i+1,j)=Zero
           Cxgsx(i+1,j)=Cxgs(i,j)
           DgsxR(i,j)=Dgs(i,j)
           HSgsxR(i,j)=Hgs(i,j)
           DUgsxR(i,j)=Zero
           DVgsxR(i,j)=Zero
           DWgsxR(i,j)=Zero
           UgsxR(i,j)=Zero
           VgsxR(i,j)=Zero
           WgsxR(i,j)=Zero
           Cxgsx(i,j)=Cxgs(i,j)
		 endif
	   elseif (Maskgs(i,j)==0) then
	     if (Maskgs(i-1,j)==0.and.Maskgs(i+1,j)==1) then
           DgsxL(i+1,j)=Dgs(i+1,j)
           HSgsxL(i+1,j)=Hgs(i+1,j)
           DUgsxL(i+1,j)=Zero
           DVgsxL(i+1,j)=Zero
           DWgsxL(i+1,j)=Zero
           UgsxL(i+1,j)=Zero
           VgsxL(i+1,j)=Zero
           WgsxL(i+1,j)=Zero
           Cxgsx(i+1,j)=Cxgs(i+1,j)
           DgsxR(i,j)=Dgs(i,j)
           HSgsxR(i,j)=Hgs(i,j)
           DUgsxR(i,j)=Zero
           DVgsxR(i,j)=Zero
           DWgsxR(i,j)=Zero
           UgsxR(i,j)=Zero
           VgsxR(i,j)=Zero
           WgsxR(i,j)=Zero
           Cxgsx(i,j)=Cxgs(i,j)
		 elseif (Maskgs(i-1,j)==1.and.Maskgs(i+1,j)==0) then
           DgsxL(i+1,j)=Dgs(i,j)
           HSgsxL(i+1,j)=Hgs(i,j)
           DUgsxL(i+1,j)=Zero
           DVgsxL(i+1,j)=Zero
           DWgsxL(i+1,j)=Zero
           UgsxL(i+1,j)=Zero
           VgsxL(i+1,j)=Zero
           WgsxL(i+1,j)=Zero
           Cxgsx(i+1,j)=Cxgs(i,j)
           DgsxR(i,j)=Dgs(i-1,j)
           HSgsxR(i,j)=Hgs(i-1,j)
           DUgsxR(i,j)=Zero
           DVgsxR(i,j)=Zero
           DWgsxR(i,j)=Zero
           UgsxR(i,j)=Zero
           VgsxR(i,j)=Zero
           WgsxR(i,j)=Zero
           Cxgsx(i,j)=Cxgs(i-1,j)
		 elseif (Maskgs(i-1,j)==0.and.Maskgs(i+1,j)==0) then
           DgsxL(i+1,j)=Dgs(i,j)
           HSgsxL(i+1,j)=Hgs(i,j)
           DUgsxL(i+1,j)=Zero
           DVgsxL(i+1,j)=Zero
           DWgsxL(i+1,j)=Zero
           UgsxL(i+1,j)=Zero
           VgsxL(i+1,j)=Zero
           WgsxL(i+1,j)=Zero
           Cxgsx(i+1,j)=Cxgs(i,j)
           DgsxR(i,j)=Dgs(i,j)
           HSgsxR(i,j)=Hgs(i,j)
           DUgsxR(i,j)=Zero
           DVgsxR(i,j)=Zero
           DWgsxR(i,j)=Zero
           UgsxR(i,j)=Zero
           VgsxR(i,j)=Zero
           WgsxR(i,j)=Zero
           Cxgsx(i,j)=Cxgs(i,j)
		 endif
	   endif
     enddo
     enddo

	 DgsxL = max(DgsxL,Small)
	 DgsxR = max(DgsxR,Small)

     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       EgsxL(i,j) = DUgsxL(i,j)
       EgsxR(i,j) = DUgsxR(i,j)
       FgsxL(i,j) = SLIDE_GAMMA*DUgsxL(i,j)*DUgsxL(i,j)/DgsxL(i,j)+0.5*Cxgsx(i,j)*DgsxL(i,j)*DgsxL(i,j)
       FgsxR(i,j) = SLIDE_GAMMA*DUgsxR(i,j)*DUgsxR(i,j)/DgsxR(i,j)+0.5*Cxgsx(i,j)*DgsxR(i,j)*DgsxR(i,j)
       GgsxL(i,j) = SLIDE_GAMMA*DUgsxL(i,j)*DVgsxL(i,j)/DgsxL(i,j)
       GgsxR(i,j) = SLIDE_GAMMA*DUgsxR(i,j)*DVgsxR(i,j)/DgsxR(i,j)
       HgsxL(i,j) = SLIDE_GAMMA*DUgsxL(i,j)*DWgsxL(i,j)/DgsxL(i,j)
       HgsxR(i,j) = SLIDE_GAMMA*DUgsxR(i,j)*DWgsxR(i,j)/DgsxR(i,j)
     enddo
     enddo

	 call construct_2D_y(Dgs,DelyDgs,DgsyL,DgsyR)
	 call construct_2D_y(Hgs,DelyHgs,HSgsyL,HSgsyR)
     call construct_2D_y(Ugs,DelyUgs,UgsyL,UgsyR)
     call construct_2D_y(Vgs,DelyVgs,VgsyL,VgsyR)
	 call construct_2D_y(Wgs,DelyWgs,WgsyL,WgsyR)
     call construct_2D_y(DUgs,DelyDUgs,DUgsyL,DUgsyR)
     call construct_2D_y(DVgs,DelyDVgs,DVgsyL,DVgsyR)
	 call construct_2D_y(DWgs,DelyDWgs,DWgsyL,DWgsyR)

     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       Cygsy(i,j)=0.5*(Cygs(i,j)+Cygs(i,j-1));
     enddo
     enddo
	 
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if (Maskgs(i,j)==1) then
	     if (Maskgs(i,j+1)==0.and.Maskgs(i,j-1)==1) then
           DgsyL(i,j+1)=Dgs(i,j)
           HSgsyL(i,j+1)=Hgs(i,j)
           DUgsyL(i,j+1)=Zero
           DVgsyL(i,j+1)=Zero
           DWgsyL(i,j+1)=Zero
           UgsyL(i,j+1)=Zero
           VgsyL(i,j+1)=Zero
           WgsyL(i,j+1)=Zero
           Cygsy(i,j+1)=Cygs(i,j)
		 elseif (Maskgs(i,j+1)==1.and.Maskgs(i,j-1)==0) then
           DgsyR(i,j)=Dgs(i,j)
           HSgsyR(i,j)=Hgs(i,j)
           DUgsyR(i,j)=Zero
           DVgsyR(i,j)=Zero
           DWgsyR(i,j)=Zero
           UgsyR(i,j)=Zero
           VgsyR(i,j)=Zero
           WgsyR(i,j)=Zero
           Cygsy(i,j)=Cygs(i,j)
		 elseif (Maskgs(i,j+1)==0.and.Maskgs(i,j-1)==0) then
           DgsyL(i,j+1)=Dgs(i,j)
           HSgsyL(i,j+1)=Hgs(i,j)
           DUgsyL(i,j+1)=Zero
           DVgsyL(i,j+1)=Zero
           DWgsyL(i,j+1)=Zero
           UgsyL(i,j+1)=Zero
           VgsyL(i,j+1)=Zero
           WgsyL(i,j+1)=Zero
           Cygsy(i,j+1)=Cygs(i,j)
           DgsyR(i,j)=Dgs(i,j)
           HSgsyR(i,j)=Hgs(i,j)
           DUgsyR(i,j)=Zero
           DVgsyR(i,j)=Zero
           DWgsyR(i,j)=Zero
           UgsyR(i,j)=Zero
           VgsyR(i,j)=Zero
           WgsyR(i,j)=Zero
           Cygsy(i,j)=Cygs(i,j)
		 endif
	   elseif (Maskgs(i,j)==0) then
	     if (Maskgs(i,j-1)==0.and.Maskgs(i,j+1)==1) then
           DgsyL(i,j+1)=Dgs(i,j+1)
           HSgsyL(i,j+1)=Hgs(i,j+1)
           DUgsyL(i,j+1)=Zero
           DVgsyL(i,j+1)=Zero
           DWgsyL(i,j+1)=Zero
           UgsyL(i,j+1)=Zero
           VgsyL(i,j+1)=Zero
           WgsyL(i,j+1)=Zero
           Cygsy(i,j+1)=Cygs(i,j+1)
           DgsyR(i,j)=Dgs(i,j)
           HSgsyR(i,j)=Hgs(i,j)
           DUgsyR(i,j)=Zero
           DVgsyR(i,j)=Zero
           DWgsyR(i,j)=Zero
           UgsyR(i,j)=Zero
           VgsyR(i,j)=Zero
           WgsyR(i,j)=Zero
           Cygsy(i,j)=Cygs(i,j)
		 elseif (Maskgs(i,j-1)==1.and.Maskgs(i,j+1)==0) then
           DgsyL(i,j+1)=Dgs(i,j)
           HSgsyL(i,j+1)=Hgs(i,j)
           DUgsyL(i,j+1)=Zero
           DVgsyL(i,j+1)=Zero
           DWgsyL(i,j+1)=Zero
           UgsyL(i,j+1)=Zero
           VgsyL(i,j+1)=Zero
           WgsyL(i,j+1)=Zero
           Cygsy(i,j+1)=Cygs(i,j)
           DgsyR(i,j)=Dgs(i,j-1)
           HSgsyR(i,j)=Hgs(i,j-1)
           DUgsyR(i,j)=Zero
           DVgsyR(i,j)=Zero
           DWgsyR(i,j)=Zero
           UgsyR(i,j)=Zero
           VgsyR(i,j)=Zero
           WgsyR(i,j)=Zero
           Cygsy(i,j)=Cygs(i,j-1)
		 elseif (Maskgs(i,j-1)==0.and.Maskgs(i,j+1)==0) then
           DgsyL(i,j+1)=Dgs(i,j)
           HSgsyL(i,j+1)=Hgs(i,j)
           DUgsyL(i,j+1)=Zero
           DVgsyL(i,j+1)=Zero
           DWgsyL(i,j+1)=Zero
           UgsyL(i,j+1)=Zero
           VgsyL(i,j+1)=Zero
           WgsyL(i,j+1)=Zero
           Cygsy(i,j+1)=Cygs(i,j)
           DgsyR(i,j)=Dgs(i,j)
           HSgsyR(i,j)=Hgs(i,j)
           DUgsyR(i,j)=Zero
           DVgsyR(i,j)=Zero
           DWgsyR(i,j)=Zero
           UgsyR(i,j)=Zero
           VgsyR(i,j)=Zero
           WgsyR(i,j)=Zero
           Cygsy(i,j)=Cygs(i,j)
		 endif
	   endif
     enddo
     enddo

	 DgsyL = max(DgsyL,Small)
	 DgsyR = max(DgsyR,Small)

     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       EgsyL(i,j) = DVgsyL(i,j)
       EgsyR(i,j) = DVgsyR(i,j)
       FgsyL(i,j) = SLIDE_GAMMA*DVgsyL(i,j)*DUgsyL(i,j)/DgsyL(i,j)
       FgsyR(i,j) = SLIDE_GAMMA*DVgsyR(i,j)*DUgsyR(i,j)/DgsyR(i,j)
       GgsyL(i,j) = SLIDE_GAMMA*DVgsyL(i,j)*DVgsyL(i,j)/DgsyL(i,j)+0.5*Cygsy(i,j)*DgsyL(i,j)*DgsyL(i,j)
       GgsyR(i,j) = SLIDE_GAMMA*DVgsyR(i,j)*DVgsyR(i,j)/DgsyR(i,j)+0.5*Cygsy(i,j)*DgsyR(i,j)*DgsyR(i,j)
       HgsyL(i,j) = SLIDE_GAMMA*DVgsyL(i,j)*DWgsyL(i,j)/DgsyL(i,j)
       HgsyR(i,j) = SLIDE_GAMMA*DVgsyR(i,j)*DWgsyR(i,j)/DgsyR(i,j)
     enddo
     enddo

     end subroutine construction_gs
	 
     subroutine wave_speed_gs
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    of deformable slide
!    Called by
!       fluxes_gs
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------
     use global
     implicit none
     integer  :: i,j
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S,UL,UR,VL,VR
     
     ! x-faces
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Maskgs(i-1,j)==1.and.Maskgs(i,j)==1) then
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxL(i,j)**2 &
		                  +Cxgsx(i,j)*abs(DgsxL(i,j)))
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxR(i,j)**2 &
		                  +Cxgsx(i,j)*abs(DgsxR(i,j)))
		 UL = SLIDE_GAMMA*UgsxL(i,j)
		 UR = SLIDE_GAMMA*UgsxR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL-UR)
         U_S = 0.5*(UL+UR)+SQR_PHI_L-SQR_PHI_R
         SgsxL(i,j) = min(UL-SQR_PHI_L,U_S-SQR_PHI_S)
         SgsxR(i,j) = max(UR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskgs(i-1,j)==0.and.Maskgs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxR(i,j)**2 &
		                  +Cxgsx(i,j)*abs(DgsxR(i,j)))
		 UR = SLIDE_GAMMA*UgsxR(i,j)
         SgsxL(i,j) = UR-2.0*SQR_PHI_R
         SgsxR(i,j) = UR+SQR_PHI_R
       elseif(Maskgs(i-1,j)==1.and.Maskgs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxL(i,j)**2 &
		                  +Cxgsx(i,j)*abs(DgsxL(i,j)))
		 UL = SLIDE_GAMMA*UgsxL(i,j)
         SgsxL(i,j) = UL-SQR_PHI_L
         SgsxR(i,j) = UL+2.0*SQR_PHI_L
       endif
     enddo
     enddo

     ! y-faces
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Maskgs(i,j-1)==1.and.Maskgs(i,j)==1) then
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyL(i,j)**2 &
		                  +Cygsy(i,j)*abs(DgsyL(i,j)))
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyR(i,j)**2 &
		                  +Cygsy(i,j)*abs(DgsyR(i,j)))
		 VL = SLIDE_GAMMA*VgsyL(i,j)
		 VR = SLIDE_GAMMA*VgsyR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL-VR)
         U_S = 0.5*(VL+VR)+SQR_PHI_L-SQR_PHI_R
         SgsyL(i,j) = min(VL-SQR_PHI_L,U_S-SQR_PHI_S)
         SgsyR(i,j) = max(VR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskgs(i,j-1)==0.and.Maskgs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyR(i,j)**2 &
		                  +Cygsy(i,j)*abs(DgsyR(i,j)))
		 VR = SLIDE_GAMMA*VgsyR(i,j)
         SgsyL(i,j) = VR-2.0*SQR_PHI_R
         SgsyR(i,j) = VR+SQR_PHI_R
       elseif(Maskgs(i,j-1)==1.and.Maskgs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyL(i,j)**2 &
		                  +Cygsy(i,j)*abs(DgsyL(i,j)))
		 VL = SLIDE_GAMMA*VgsyL(i,j)
         SgsyL(i,j) = VL-SQR_PHI_L
         SgsyR(i,j) = VL+2.0*SQR_PHI_L
       endif
     enddo
     enddo

     end subroutine wave_speed_gs
	 
     subroutine fluxes_at_faces_HLL_gs
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes_gs
!    Last update: 10/07/2018, Cheng Zhang
!---------------------------------------------
     use global
     implicit none

     ! horizontal fluxes
     call HLL_2D(Mloc1,Nloc,SgsxL,SgsxR,EgsxL,EgsxR,DgsxL, DgsxR, Egsx)
     call HLL_2D(Mloc,Nloc1,SgsyL,SgsyR,EgsyL,EgsyR,DgsyL, DgsyR, Egsy)
     call HLL_2D(Mloc1,Nloc,SgsxL,SgsxR,FgsxL,FgsxR,DUgsxL,DUgsxR,Fgsx)
     call HLL_2D(Mloc,Nloc1,SgsyL,SgsyR,FgsyL,FgsyR,DUgsyL,DUgsyR,Fgsy)
     call HLL_2D(Mloc1,Nloc,SgsxL,SgsxR,GgsxL,GgsxR,DVgsxL,DVgsxR,Ggsx)
     call HLL_2D(Mloc,Nloc1,SgsyL,SgsyR,GgsyL,GgsyR,DVgsyL,DVgsyR,Ggsy)  
     call HLL_2D(Mloc1,Nloc,SgsxL,SgsxR,HgsxL,HgsxR,DWgsxL,DWgsxR,Hgsx)
     call HLL_2D(Mloc,Nloc1,SgsyL,SgsyR,HgsyL,HgsyR,DWgsyL,DWgsyR,Hgsy) 	 

     return
     end subroutine fluxes_at_faces_HLL_gs

     subroutine HLL_2D(M,N,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLL reconstruction for 2D model
!    Called by
!       fluxes_at_faces_HLL_gs
!    Last update: 10/07/2018, Cheng Zhang
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

     subroutine flux_bc_gs
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    for deformable slide
!    Called by
!       fluxes_gs
!    Last update: 10/07/2018, Cheng Zhang
!--------------------------------------------
     use global
     implicit none
     integer :: i,j

     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
	 do j = Jbeg,Jend
       Egsx(Ibeg,j) = Zero
	   Fgsx(Ibeg,j) = 0.5*Cxgsx(Ibeg,j)*DgsxR(Ibeg,j)*DgsxR(Ibeg,j)
       Ggsx(Ibeg,j) = Zero
	   Hgsx(Ibeg,j) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
       Egsx(Iend1,j) = Zero
       Fgsx(Iend1,j) = 0.5*Cxgsx(Iend1,j)*DgsxL(Iend1,j)*DgsxL(Iend1,j)
       Ggsx(Iend1,j) = Zero
	   Hgsx(Iend1,j) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
       Egsy(i,Jbeg) = Zero
       Fgsy(i,Jbeg) = Zero
       Ggsy(i,Jbeg) = 0.5*Cygsy(i,Jbeg)*DgsyR(i,Jbeg)*DgsyR(i,Jbeg)
	   Hgsy(i,Jbeg) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
       Egsy(i,Jend1) = Zero
       Fgsy(i,Jend1) = Zero
	   Ggsy(i,Jend1) = 0.5*Cygsy(i,Jend1)*DgsyL(i,Jend1)*DgsyL(i,Jend1)
	   Hgsy(i,Jend1) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

    do j = Jbeg-1,Jend+1
    do i = Ibeg-1,Iend+1
      if(Maskgs(i,j)==0) then
        if(Maskgs(i-1,j)==1) then !right boundary
          Egsx(i,j) = Zero
          Fgsx(i,j) = 0.5*Cxgs(i-1,j)*DgsxL(i,j)*DgsxL(i,j)
          Ggsx(i,j) = Zero
		  Hgsx(i,j) = Zero
        elseif(Maskgs(i+1,j)==1) then !left boundary
          Egsx(i+1,j) = Zero
          Fgsx(i+1,j) = 0.5*Cxgs(i+1,j)*DgsxR(i+1,j)*DgsxR(i+1,j)
          Ggsx(i+1,j) = Zero
		  Hgsx(i+1,j) = Zero
        elseif(Maskgs(i,j-1)==1) then !back boundary
          Egsy(i,j) = Zero
          Fgsy(i,j) = Zero
          Ggsy(i,j) = 0.5*Cygs(i,j-1)*DgsyL(i,j)*DgsyL(i,j)
		  Hgsx(i,j) = Zero
        elseif(Maskgs(i,j+1)==1) then !front boundary
          Egsy(i,j+1) = Zero
          Fgsy(i,j+1) = Zero
          Ggsy(i,j+1) = 0.5*Cygs(i,j+1)*DgsyR(i,j+1)*DgsyR(i,j+1)
		  Hgsx(i,j+1) = Zero
        endif
      endif
    enddo
    enddo

     end subroutine flux_bc_gs
	 
	 
     subroutine source_terms_gs
!------------------------------------------------
!    This subroutine is used to evaluate source terms
!    in momentum eqs of deformable slide
!    Called by
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------------
     use global
     implicit none
     integer :: i,j
	 real(SP), dimension(:,:),allocatable :: PH,DelxPH,DelyPH,  &
	                                         Taxxd,Tayxd,Taxyd,Tayyd,Taxzd,Tayzd, &
											 DelxTaxxd,DelyTayxd,DelxTaxyd,DelyTayyd,DelxTaxzd,DelyTayzd
     real(SP) :: alpha0,alpha1,Hc0L,Hc0R
	 
	 ! source terms for hydrostatic pressure
	 ! well-balanced upwind scheme 
	   !x direction
       do j = Jbeg,Jend
       do i = Ibeg,Iend1
         if ((abs(SgsxR(i,j)-SgsxL(i,j)))<Small) then
           alpha0=(SgsxR(i,j)*abs(SgsxL(i,j))-SgsxL(i,j)*abs(SgsxR(i,j)))/Small
           alpha1=(abs(SgsxR(i,j))-abs(SgsxL(i,j)))/Small
         else
           alpha0=(SgsxR(i,j)*abs(SgsxL(i,j))-SgsxL(i,j)*abs(SgsxR(i,j)))/(SgsxR(i,j)-SgsxL(i,j))
           alpha1=(abs(SgsxR(i,j))-abs(SgsxL(i,j)))/(SgsxR(i,j)-SgsxL(i,j))
         endif
         Hc0L=DgsxL(i,j)+HSgsxL(i,j)
         Hc0R=DgsxR(i,j)+HSgsxR(i,j)
         SrcmgsxL(i,j)= 0.5*alpha0*(Hc0R-Hc0L)
         SrcmgsxR(i,j)=-0.5*alpha0*(Hc0R-Hc0L)
         SrcpgsxL(i,j)= 0.5*(alpha1-1.)*Grav*(DgsxR(i,j)+DgsxL(i,j))/2.0*(Hc0R-Hc0L) &
		               -Grav*(DgsxL(i,j)+DgsxR(i-1,j))/2.*Hc0L
         SrcpgsxR(i,j)=-0.5*(alpha1+1.)*Grav*(DgsxR(i,j)+DgsxL(i,j))/2.0*(Hc0R-Hc0L) &
		               +Grav*(DgsxL(i+1,j)+DgsxR(i,j))/2.*Hc0R
       enddo
       enddo
	   !y direction
       do j = Jbeg,Jend1
       do i = Ibeg,Iend
         if ((abs(SgsyR(i,j)-SgsyL(i,j)))<Small) then
           alpha0=(SgsyR(i,j)*abs(SgsyL(i,j))-SgsyL(i,j)*abs(SgsyR(i,j)))/Small
           alpha1=(abs(SgsyR(i,j))-abs(SgsyL(i,j)))/Small
         else
           alpha0=(SgsyR(i,j)*abs(SgsyL(i,j))-SgsyL(i,j)*abs(SgsyR(i,j)))/(SgsyR(i,j)-SgsyL(i,j))
           alpha1=(abs(SgsyR(i,j))-abs(SgsyL(i,j)))/(SgsyR(i,j)-SgsyL(i,j))
         endif
         Hc0L=DgsyL(i,j)+HSgsyL(i,j)
         Hc0R=DgsyR(i,j)+HSgsyR(i,j)
         SrcmgsyL(i,j)= 0.5*alpha0*(Hc0R-Hc0L)
         SrcmgsyR(i,j)=-0.5*alpha0*(Hc0R-Hc0L)
         SrcpgsyL(i,j)= 0.5*(alpha1-1.)*Grav*(DgsyR(i,j)+DgsyL(i,j))/2.0*(Hc0R-Hc0L) &
		               -Grav*(DgsyL(i,j)+DgsyR(i,j-1))/2.*Hc0L
         SrcpgsyR(i,j)=-0.5*(alpha1+1.)*Grav*(DgsyR(i,j)+DgsyL(i,j))/2.0*(Hc0R-Hc0L) &
		               +Grav*(DgsyL(i,j+1)+DgsyR(i,j))/2.*Hc0R
       enddo
       enddo

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrcmgsX(i,j) = -(SrcmgsxL(i+1,j)+SrcmgsxR(i,j))/dx*Maskgs(i,j)
       SrcmgsY(i,j) = -(SrcmgsyL(i,j+1)+SrcmgsyR(i,j))/dy*Maskgs(i,j)
       SrcpgsX(i,j) = -(SrcpgsxL(i+1,j)+SrcpgsxR(i,j))/dx*Maskgs(i,j)
       SrcpgsY(i,j) = -(SrcpgsyL(i,j+1)+SrcpgsyR(i,j))/dy*Maskgs(i,j)
     enddo
     enddo
	 
	 ! source terms from upper layer pressure
     allocate(DelxPH(Mloc,Nloc))
     allocate(DelyPH(Mloc,Nloc))
	 allocate(PH(Mloc,Nloc))
	 DelxPH = zero;DelyPH = zero;PH = zero
     do j = 1,Nloc
     do i = 1,Mloc
	   PH(i,j) = Rho0*Grav*(Eta(i,j)+Hgs(i,j))! hydrostatic pressure over the slide
     enddo
     enddo
	 call delxFun_2D(PH,DelxPH)
	 call delyFun_2D(PH,DelyPH)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrcpgsX(i,j) = SrcpgsX(i,j)-Dgs(i,j)/Rhogs(i,j)*DelxPH(i,j)*Maskgs(i,j)*Mask(i,j)
       SrcpgsY(i,j) = SrcpgsY(i,j)-Dgs(i,j)/Rhogs(i,j)*DelyPH(i,j)*Maskgs(i,j)*Mask(i,j)
     enddo
     enddo
	 
     ! source terms for rhelogy stress
     allocate(Taxxd(Mloc,Nloc))
     allocate(Tayxd(Mloc,Nloc))
	 allocate(Taxyd(Mloc,Nloc))
	 allocate(Tayyd(Mloc,Nloc))
	 allocate(Taxzd(Mloc,Nloc))
	 allocate(Tayzd(Mloc,Nloc))
     allocate(DelxTaxxd(Mloc,Nloc))
     allocate(DelyTayxd(Mloc,Nloc))
	 allocate(DelxTaxyd(Mloc,Nloc))
	 allocate(DelyTayyd(Mloc,Nloc))
	 allocate(DelxTaxzd(Mloc,Nloc))
	 allocate(DelyTayzd(Mloc,Nloc))
     if(trim(RHEO_OPT)=='VISCOUS') then
	   call laminar(SLIDE_VISCOSITY,Rhogs,Tbxx,Tbxy,Tbxz,Tbyx,Tbyy,Tbyz,Tbzx,Tbzy,Tbzz)
	   Taxx=0.5*Tbxx;Tayx=0.5*Tbyx;Taxy=0.5*Tbxy;Tayy=0.5*Tbyy;Taxz=0.5*Tbxz;Tayz=0.5*Tbyz
	   Tbxx=zero;Tbyx=zero;Tbzx=zero;Tbxy=zero;Tbyy=zero;Tbzy=zero;Tbxz=zero;Tbyz=zero;Tbzz=zero;
	 elseif(trim(RHEO_OPT)=='GRANULAR') then
	   call coulomb
	 endif
	 Taxxd=Taxx*Dgs;Tayxd=Tayx*Dgs;Taxyd=Taxy*Dgs;Tayyd=Tayy*Dgs;Taxzd=Taxz*Dgs;Tayzd=Tayz*Dgs
	 DelxTaxxd=zero;DelyTayxd=zero;DelxTaxyd=zero;DelyTayyd=zero;DelxTaxzd=zero;DelyTayzd=zero
	 call delxFun_gs_2D(Taxxd,DelxTaxxd)
	 call delyFun_gs_2D(Tayxd,DelyTayxd)
	 call delxFun_gs_2D(Taxyd,DelxTaxyd)
	 call delyFun_gs_2D(Tayyd,DelyTayyd)
	 call delxFun_gs_2D(Taxzd,DelxTaxzd)
	 call delyFun_gs_2D(Tayzd,DelyTayzd)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrctgsX(i,j) = 1.0/Rhogs(i,j)*(DelxTaxxd(i,j)+DelyTayxd(i,j) &
	                 -Tbxx(i,j)*DelxH0(i,j)-Tbyx(i,j)*DelyH0(i,j)-Tbzx(i,j))*Maskgs(i,j)
       SrctgsY(i,j) = 1.0/Rhogs(i,j)*(DelxTaxyd(i,j)+DelyTayyd(i,j) &
	                 -Tbxy(i,j)*DelxH0(i,j)-Tbyy(i,j)*DelyH0(i,j)-Tbzy(i,j))*Maskgs(i,j)
       SrctgsZ(i,j) = 1.0/Rhogs(i,j)*(DelxTaxzd(i,j)+DelyTayzd(i,j) &
	                 -Tbxz(i,j)*DelxH0(i,j)-Tbyz(i,j)*DelyH0(i,j)-Tbzz(i,j))*Maskgs(i,j)
     enddo
     enddo
	 
     deallocate(DelxPH)
     deallocate(DelyPH)
	 deallocate(PH)
     deallocate(Taxxd)
     deallocate(Tayxd)
	 deallocate(Taxyd)
	 deallocate(Tayyd)
	 deallocate(Taxzd)
	 deallocate(Tayzd)
     deallocate(DelxTaxxd)
     deallocate(DelyTayxd)
	 deallocate(DelxTaxyd)
	 deallocate(DelyTayyd)
	 deallocate(DelxTaxzd)
	 deallocate(DelyTayzd)

     end subroutine source_terms_gs
	 
     subroutine laminar(Nu,Rho,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz)
!-----------------------------------------
!    Rheology of laminar flow
!    Called by 
!       source_terms_gs, coulomb
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Ugs,Vgs,Dgs,DelxH0,DelyH0,Zero
     implicit none
     real(SP),intent(in) :: NU
	 real(SP),intent(in),dimension(Mloc,Nloc) :: Rho
     real(SP),intent(out),dimension(Mloc,Nloc) :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
     integer :: i,j

     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   if (abs(Ugs(i,j))>=1e-8.and.abs(Vgs(i,j))>=1e-8) then
         Txx(i,j)=6.0*Nu*Rho(i,j)*Ugs(i,j)/Dgs(i,j)*DelxH0(i,j)
         Tyy(i,j)=6.0*Nu*Rho(i,j)*Vgs(i,j)/Dgs(i,j)*DelyH0(i,j)
         Tzz(i,j)=-6.0*Nu*Rho(i,j)/Dgs(i,j)*(Ugs(i,j)*DelxH0(i,j)+Vgs(i,j)*DelyH0(i,j))
         Txy(i,j)=3.0*Nu*Rho(i,j)/Dgs(i,j)*(Ugs(i,j)*DelyH0(i,j)+Vgs(i,j)*DelxH0(i,j))
         Txz(i,j)=3.0*Nu*Rho(i,j)*Ugs(i,j)/Dgs(i,j) &
	             -3.0*Nu*Rho(i,j)/Dgs(i,j)*(Ugs(i,j)*DelxH0(i,j) &
	                                       +Vgs(i,j)*DelyH0(i,j))*DelxH0(i,j)
         Tyz(i,j)=3.0*Nu*Rho(i,j)*Vgs(i,j)/Dgs(i,j) &
	             -3.0*Nu*Rho(i,j)/Dgs(i,j)*(Ugs(i,j)*DelxH0(i,j) &
	                                       +Vgs(i,j)*DelyH0(i,j))*DelyH0(i,j)
	    else
         Txx(i,j)=Zero
         Tyy(i,j)=Zero
         Tzz(i,j)=Zero
         Txy(i,j)=Zero
         Txz(i,j)=Zero
         Tyz(i,j)=Zero
        endif
     enddo
     enddo
     call phi_2D_coll(Txx)
	 call phi_2D_coll(Tyy)
	 call phi_2D_coll(Tzz)
     call phi_2D_coll(Txy)
	 call phi_2D_coll(Txz)
	 call phi_2D_coll(Tyz)
     Tyx=Txy;Tzx=Txz;Tzy=Tyz

     return
     end subroutine laminar
	 
     subroutine coulomb
!------------------------------------------------
!    This subroutine is used to evaluate source terms
!    in momentum eqs of deformable slide
!    Called by
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------------
     use global
     implicit none
     integer :: i,j
	 real(SP), dimension(:,:),allocatable :: Tbfxx,Tbfxy,Tbfxz,Tbfyx,Tbfyy, &
	                                         Tbfyz,Tbfzx,Tbfzy,Tbfzz, &
											 Rhof
	 real(SP) :: Nugse
	 
	 allocate(Tbfxx(Mloc,Nloc))
     allocate(Tbfxy(Mloc,Nloc))
	 allocate(Tbfxz(Mloc,Nloc))
	 allocate(Tbfyx(Mloc,Nloc))
     allocate(Tbfyy(Mloc,Nloc))
	 allocate(Tbfyz(Mloc,Nloc))
	 allocate(Tbfzx(Mloc,Nloc))
     allocate(Tbfzy(Mloc,Nloc))
	 allocate(Tbfzz(Mloc,Nloc))
	 allocate(Rhof(Mloc,Nloc))
	 Tbfxx = Zero;Tbfxy = Zero;Tbfxz = Zero;Tbfyx = Zero;Tbfyy = Zero;Tbfyz = Zero
	 Tbfzx = Zero;Tbfzy = Zero;Tbfzz = Zero;Tezz = Zero
	 
	 ! fluid stress
	 Nugse=(1.-SLIDE_CONC)*Visc
	 Rhof=Rho0
     call laminar(Nugse,Rhof,Tbfxx,Tbfxy,Tbfxz,Tbfyx,Tbfyy,Tbfyz,Tbfzx,Tbfzy,Tbfzz)
	 
	 ! solid stress+fluid stress
     do j = 1,Nloc
     do i = 1,Mloc
	   Tezz(i,j) = Rhogs(i,j)*DwDt(i,j)*Dgs(i,j)+SLIDE_CONC*GRAIN_DENSITY*Grav*Dgs(i,j)
	   Pss(i,j) = SLIDE_CONC*GRAIN_DENSITY*Grav*Dgs(i,j)+Qsgs(i,j)
	   !x momentum
	   Taxx(i,j) = 0.5*Tbfxx(i,j)*Maskgs(i,j)!solid stress is in Cx_gs
	   Tayx(i,j) = 0.5*(Tbfyx(i,j)-Cyx(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbxx(i,j) = (Tbfxx(i,j)-Cxx(i,j)/max(Czz(i,j),Small)*Tezz(i,j)+Pss(i,j))*Maskgs(i,j)
	   Tbyx(i,j) = (Tbfyx(i,j)-Cyx(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbzx(i,j) = Tbfzx(i,j)*Maskgs(i,j)
	   !y momentum
	   Taxy(i,j) = 0.5*(Tbfxy(i,j)-Cxy(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tayy(i,j) = 0.5*Tbfyy(i,j)*Maskgs(i,j)!solid stress is in Cy_gs
	   Tbxy(i,j) = (Tbfxy(i,j)-Cxy(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbyy(i,j) = (Tbfyy(i,j)-Cyy(i,j)/max(Czz(i,j),Small)*Tezz(i,j)+Pss(i,j))*Maskgs(i,j)
	   Tbzy(i,j) = Tbfzy(i,j)*Maskgs(i,j)!solid stress determined in eval_duvw_gs
	   !z momentum
	   Taxz(i,j) = 0.5*(Tbfxz(i,j)-Cxz(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tayz(i,j) = 0.5*(Tbfyz(i,j)-Cyz(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbxz(i,j) = (Tbfxz(i,j)-Cxz(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbyz(i,j) = (Tbfyz(i,j)-Cyz(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbzz(i,j) = (Tbfzz(i,j)-Tezz(i,j)+Pss(i,j))*Maskgs(i,j) !Ts=-Te+Ps
     enddo
     enddo
	 
	 call phi_2D_coll(Tezz)
	 call phi_2D_coll(Pss)
     call phi_2D_coll(Taxx)
	 call phi_2D_coll(Tayx)
	 call phi_2D_coll(Tbxx)
	 call phi_2D_coll(Tbyx)
	 call phi_2D_coll(Tbzx)
	 call phi_2D_coll(Taxy)
     call phi_2D_coll(Tayy)
	 call phi_2D_coll(Tbxy)
	 call phi_2D_coll(Tbyy)
     call phi_2D_coll(Tbzy)
	 call phi_2D_coll(Taxz)
     call phi_2D_coll(Tayz)
	 call phi_2D_coll(Tbxz)
	 call phi_2D_coll(Tbyz)
     call phi_2D_coll(Tbzz)
	 
	 deallocate(Tbfxx)
     deallocate(Tbfxy)
	 deallocate(Tbfxz)
	 deallocate(Tbfyx)
     deallocate(Tbfyy)
	 deallocate(Tbfyz)
	 deallocate(Tbfzx)
     deallocate(Tbfzy)
	 deallocate(Tbfzz)
	 deallocate(Rhof)

     end subroutine coulomb
	 

	 subroutine eval_duvw_gs(ISTEP)
!-----------------------------------------------
!    Update all variables Dgs,Ugs,Vgs for deformable slide
!    Called by
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP),dimension(:,:),allocatable :: R1,R2,R3,R4,Thetax,Thetay,DelxTheta,DelyTheta
     integer :: i,j
     REAL(SP) :: FROUDECAP_gs,FroudeU,DUU,Dangle,TauC,Xib

     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc))
     allocate(R3(Mloc,Nloc))
     allocate(R4(Mloc,Nloc))
     allocate(Thetax(Mloc,Nloc))
	 allocate(Thetay(Mloc,Nloc))
     allocate(DelxTheta(Mloc,Nloc))
     allocate(DelyTheta(Mloc,Nloc))
	 
     ! solve total slide depth Dgs
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==0) cycle
       R1(i,j) = -1.0/dx*(Egsx(i+1,j)-Egsx(i,j))-1.0/dy*(Egsy(i,j+1)-Egsy(i,j)) &
	             +SrcmgsX(i,j)+SrcmgsY(i,j)
# if defined (OBSTACLE)
       if(set_flag_gs(i,j)==0) then
         Dgs(i,j) = Dgs(i,j)+dt*R1(i,j)
       endif
# else
	   Dgs(i,j) = Dgs(i,j)+dt*R1(i,j)
# endif
     enddo
     enddo
     ! update Dgs and Hgs          
     Dgs = max(Dgs,SLIDE_MINTHICK)
     call wl_bc_gs
     Hgs = Hc0-Dgs
	 
	 ! solve DUgs
     R2 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskgs(i,j)==0) cycle  
       R2(i,j) = -1.0/dx*(Fgsx(i+1,j)-Fgsx(i,j))-1.0/dy*(Fgsy(i,j+1)-Fgsy(i,j)) &
                 +SrcpgsX(i,j)+SrctgsX(i,j)
# if defined (OBSTACLE)
       if(set_flag_gs(i,j)==1) then
         DUgs(i,j) = Dgs(i,j)*obs_u
       else
         DUgs(i,j) = DUgs(i,j)+dt*R2(i,j)
       endif
# else
	   DUgs(i,j) = DUgs(i,j)+dt*R2(i,j)
# endif
     enddo
     enddo
	 
	 ! solve DVgs
     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskgs(i,j)==0) cycle  
       R3(i,j) = -1.0/dx*(Ggsx(i+1,j)-Ggsx(i,j))-1.0/dy*(Ggsy(i,j+1)-Ggsy(i,j)) &
                 +SrcpgsY(i,j)+SrctgsY(i,j)
# if defined (OBSTACLE)
       if(set_flag_gs(i,j)==1) then
         DVgs(i,j) = Dgs(i,j)*obs_v
       else
         DVgs(i,j) = DVgs(i,j)+dt*R3(i,j)
       endif
# else
	   DVgs(i,j) = DVgs(i,j)+dt*R3(i,j)
# endif
     enddo
     enddo
	 
	 ! solve DWgs
     R4 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskgs(i,j)==0) cycle  
       R4(i,j) = -1.0/dx*(Hgsx(i+1,j)-Hgsx(i,j))-1.0/dy*(Hgsy(i,j+1)-Hgsy(i,j)) &
                 +SrctgsZ(i,j)
# if defined (OBSTACLE)
       if(set_flag_gs(i,j)==1) then
         DWgs(i,j) = Dgs(i,j)*obs_w
       else
         DWgs(i,j) = DWgs(i,j)+dt*R4(i,j)
       endif
# else
	   DWgs(i,j) = DWgs(i,j)+dt*R4(i,j)
# endif
     enddo
     enddo
	 
     ! bottom shear stress
     if(trim(RHEO_OPT)=='GRANULAR') then
        call delxFun_gs_2D(Dgs,DelxDgs)
        call delyFun_gs_2D(Dgs,DelyDgs)
        call construct_2D_x(Dgs,DelxDgs,DgsxL,DgsxR)
        call construct_2D_y(Dgs,DelyDgs,DgsyL,DgsyR)
        call DelxH_gs
		Thetax=Zero;Thetay=Zero;DelxTheta=Zero;DelyTheta=Zero
!		do j = 1,Nloc
!        do i = 1,Mloc
!          Thetax(i,j) = atan(DelxH0(i,j))
!	      Thetay(i,j) = atan(DelyH0(i,j))
!        enddo
!        enddo
!        call delxFun_gs_2D(Thetax,DelxTheta)
!        call delyFun_gs_2D(Thetay,DelyTheta)
		do j = Jbeg,Jend
        do i = Ibeg,Iend
          if(Maskgs(i,j)==1) then
		    TauC=tan(PhiBed(i,j))/max(Czz0(i,j),Small) &
			     *((Rhogs(i,j)*DwDt(i,j)+SLIDE_CONC*GRAIN_DENSITY*Grav) &
				   *((DgsxL(i,j)+DgsxR(i+1,j))/2.+(DgsyL(i,j)+DgsyR(i,j+1))/2.)/2.)
            if(sqrt(DUgs(i,j)**2+DVgs(i,j)**2)>dt/Rhogs(i,j)*TauC) then
              Tbzx(i,j) = -Czx(i,j)/max(Czz(i,j),Small) &
			             *((Rhogs(i,j)*DwDt(i,j)+SLIDE_CONC*GRAIN_DENSITY*Grav) &
						   *(DgsxL(i,j)+DgsxR(i+1,j))/2.) &
						 *(1.-(DUgs(i,j)/Dgs(i,j))**2*DelxTheta(i,j)/Grav)
              Tbzy(i,j) = -Czy(i,j)/max(Czz(i,j),Small) &
			             *((Rhogs(i,j)*DwDt(i,j)+SLIDE_CONC*GRAIN_DENSITY*Grav) &
						   *(DgsyL(i,j)+DgsyR(i,j+1))/2.) &
						 *(1.-(DVgs(i,j)/Dgs(i,j))**2*DelyTheta(i,j)/Grav)
            else
              Dgs(i,j) = Dgs(i,j)-dt*R1(i,j)
              Tbzx(i,j) = Rhogs(i,j)/dt*DUgs(i,j)
              Tbzy(i,j) = Rhogs(i,j)/dt*DVgs(i,j)
			  DWgs(i,j) = DWgs(i,j)-dt*R4(i,j)
            endif
          endif
        enddo
        enddo
       ! update Dgs and Hgs          
       Dgs = max(Dgs,SLIDE_MINTHICK)
       call wl_bc_gs
       Hgs = Hc0-Dgs
	   ! update DUgs and DVgs
	   do j = Jbeg,Jend
       do i = Ibeg,Iend
         DUgs(i,j) = DUgs(i,j)-dt/Rhogs(i,j)*Tbzx(i,j)
         DVgs(i,j) = DVgs(i,j)-dt/Rhogs(i,j)*Tbzy(i,j)
       enddo
       enddo
	 elseif(trim(RHEO_OPT)=='VISCOUS') then
	   call DelxH_gs
       do j = Jbeg,Jend
       do i = Ibeg,Iend
	     Xib = sqrt(1.0+DelxH0(i,j)**2+DelyH0(i,j)**2)
	     ! implicit basal stress 
	     DUgs(i,j) = DUgs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY*Xib**2/Dgs(i,j)**2)
	     DVgs(i,j) = DVgs(i,j)/dt/(1./dt+3.*SLIDE_VISCOSITY*Xib**2/Dgs(i,j)**2)
		 DWgs(i,j) = DWgs(i,j)+3.*dt*SLIDE_VISCOSITY/Dgs(i,j)**2 &
		             *(DelxH0(i,j)*Xib**2*DUgs(i,j)+DelyH0(i,j)*Xib**2*DVgs(i,j))
       enddo
       enddo
     endif
	 
     ! SSP Runge-Kutta time stepping
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   Dgs(i,j) = ALPHA(ISTEP)*Dgs0(i,j)+BETA(ISTEP)*Dgs(i,j)
       DUgs(i,j) = ALPHA(ISTEP)*DUgs0(i,j)+BETA(ISTEP)*DUgs(i,j)
       DVgs(i,j) = ALPHA(ISTEP)*DVgs0(i,j)+BETA(ISTEP)*DVgs(i,j)
       DWgs(i,j) = ALPHA(ISTEP)*DWgs0(i,j)+BETA(ISTEP)*DWgs(i,j)
     enddo
     enddo
	 
     ! update Dgs and Hgs          
     Dgs = max(Dgs,SLIDE_MINTHICK)
     call wl_bc_gs
     Hgs = Hc0-Dgs

	 ! obtain hydrostatic velocity
	 call get_UVW_gs
	 
	 ! vertical velocity in hydrostatic model
	 if(.not.NON_HYDRO_SLD) then
	   ! solve Wgs
	   call Kbc_Surface_gs
       call Kbc_Bottom_gs
	   do j = Jbeg,Jend
       do i = Ibeg,Iend
	     if (Maskgs(i,j)==1) then 
           Wgs(i,j) = 0.5*(Wtgs(i,j)+Wbgs(i,j))
           DWgs(i,j) = Wgs(i,j)*Dgs(i,j)
		 endif
	   enddo
	   enddo
	   
	 ! avoid non-zero DUgs/DVgs/DWgs in grid with Mask(i,j)==0
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==0) then
         DUgs(i,j) = Zero
         DVgs(i,j) = Zero
         DWgs(i,j) = Zero
       endif
     enddo
     enddo
	 
	 ! run non-hydrostatic simulation
	 elseif(NON_HYDRO_SLD.and.ISTEP==It_Order) then
	   
       ! solve dynamic pressure 
       call poisson_solver_gs

       ! correct velocity field  
       call projection_corrector_gs
	   
	 endif
	 
# if defined (OBSTACLE)
     ! update velocities for calculating IB forces
     call get_UVW_gs     

     ! calculate forcing at obstacle boundary                                                                                    
     call imm_obs_gs

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==0) cycle

       if(set_flag_gs(i,j)==1) then
         DUgs(i,j) = Dgs(i,j)*obs_u
         DVgs(i,j) = Dgs(i,j)*obs_v
         DWgs(i,j) = Dgs(i,j)*obs_w
       else
         DUgs(i,j) = DUgs(i,j)+dt*ObsForceX_gs(i,j)
         DVgs(i,j) = DVgs(i,j)+dt*ObsForceY_gs(i,j)
         DWgs(i,j) = DWgs(i,j)+dt*ObsForceZ_gs(i,j)
       endif
     enddo
     enddo
# endif
	 
     ! limit the maximum Froude number
!     FROUDECAP_gs=5.0
!     DO J=Jbeg,Jend
!     DO I=Ibeg,Iend
!       IF(Maskgs(I,J)>0)THEN
!         FroudeU=SQRT(grav*Dgs(I,J))*FROUDECAP_gs*Dgs(I,J)
!         DUU=SQRT(DUgs(I,J)**2+DVgs(I,J)**2+DWgs(I,J)**2)
!         IF(DUU>FroudeU)THEN
!          Dangle=atan2(DVgs(I,J),DUgs(I,J))
!          DUgs(I,J)=FroudeU*COS(Dangle)
!          DVgs(I,J)=FroudeU*SIN(Dangle)
!		  DWgs(i,j)=Zero
!         ENDIF
!       ENDIF
!     ENDDO
!     ENDDO

     ! update maskgs for deformable slide
	 ! call update_maskgs

     ! boundary conditions and final velocity
     call get_UVW_gs

	 deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(R4)
     deallocate(Thetax)
	 deallocate(Thetay)
     deallocate(DelxTheta)
     deallocate(DelyTheta)
	 
	 end subroutine eval_duvw_gs
	 
     subroutine wl_bc_gs
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duvw_gs
!    Last update: 10/07/2018, Cheng Zhang
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
       Dgs(Ibeg-i,j) = Dgs(Ibeg+i-1,j)
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
       Dgs(Iend+i,j) = Dgs(Iend-i+1,j)
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
         Dgs(i,j) = Dgs(i,Jbeg+Nghost-j)
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
         Dgs(i,Jend+j) = Dgs(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(Dgs)
# endif
     
     return
     end subroutine wl_bc_gs
	 
     subroutine KBC_surface_gs
!-------------------------------------------------
!    Applied kinematic boundary conditions at free surface                                                                        
!    of deformable slide                                                                                   
!    Called by
!       eval_duvw_gs                                                                                                                  
!    Last update: 10/07/2018, Cheng Zhang                                                                                    
!------------------------------------------------                                                                                
     use global
     implicit none
     integer :: i,j
     real(SP),dimension(:,:),allocatable :: Ugstmp,Vgstmp

     allocate(Ugstmp(Mloc,Nloc))
	 allocate(Vgstmp(Mloc,Nloc))

	 Wtgs=Zero
!	 call delxFun_gs_2D(Hgs,DelxHgs)
!	 call delyFun_gs_2D(Hgs,DelyHgs)
!	 do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       Wtgs(i,j) = -(Hgs(i,j)-(Hc0(i,j)-Dgs0(i,j)))/dt &
!	               -DUgs(i,j)/Dgs(i,j)*DelxHgs(i,j) &
!				   -DVgs(i,j)/Dgs(i,j)*DelyHgs(i,j)                                                     
!     enddo
!     enddo

     ! give better results when bed is not moving (no time derivative involved)
	 Ugstmp=Zero;Vgstmp=Zero
	 do j = 1,Nloc
     do i = 1,Mloc
       Ugstmp(i,j) = DUgs(i,j)/Dgs(i,j)
	   Vgstmp(i,j) = DVgs(i,j)/Dgs(i,j)
     enddo
     enddo
	 call delxFun_gs_2D(Ugstmp,DelxUgs)
	 call delyFun_gs_2D(Vgstmp,DelyVgs)
	 do j = Jbeg,Jend
     do i = Ibeg,Iend
       Wtgs(i,j) = -Dgs(i,j)*(DelxUgs(i,j)+DelyVgs(i,j)) &
	               -(Ugstmp(i,j)*DelxH0(i,j)+Vgstmp(i,j)*DelyH0(i,j))                                                     
     enddo
     enddo
     deallocate(Ugstmp)
	 deallocate(Vgstmp)

     end subroutine KBC_surface_gs

     subroutine KBC_bottom_gs
!-------------------------------------------------  
!    Applied kinematic boundary conditions at bottom                                                                        
!    of deformable slide                                                                                   
!    Called by
!       eval_duvw_gs                                                                                                                  
!    Last update: 10/07/2018, Cheng Zhang                                                                                    
!------------------------------------------------                                                                                 
     use global
     implicit none
     integer :: i,j

	 Wbgs=Zero
	 call DelxH_gs
	 do j = Jbeg,Jend
     do i = Ibeg,Iend
       Wbgs(i,j) = -DUgs(i,j)/Dgs(i,j)*DelxH0(i,j) &
	               -DVgs(i,j)/Dgs(i,j)*DelyH0(i,j)                                                     
     enddo
     enddo

     end subroutine KBC_bottom_gs
	 
     subroutine get_UVW_gs
!------------------------------------------------
!    Obtain Ugs,Vgs for deformable slide
!    Called by
!       eval_duvw_gs
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j
	 
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   if(Maskgs(i,j)==1) then
         Ugs(i,j) = DUgs(i,j)/Dgs(i,j)
         Vgs(i,j) = DVgs(i,j)/Dgs(i,j)
		 Wgs(i,j) = DWgs(i,j)/Dgs(i,j)
	   else
         Ugs(i,j) = Zero
         Vgs(i,j) = Zero
		 Wgs(i,j) = Zero
         DUgs(i,j) = Zero
         DVgs(i,j) = Zero
		 DWgs(i,j) = Zero
	   endif
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc_gs
# if defined (PARALLEL)
     call phi_2D_exch(Ugs)
     call phi_2D_exch(Vgs)
	 call phi_2D_exch(Wgs)
     call phi_2D_exch(DUgs)
     call phi_2D_exch(DVgs)
	 call phi_2D_exch(DWgs)
# endif

     end subroutine get_UVW_gs

     subroutine vel_bc_gs
!----------------------------------------------------
!    Boundary conditions for slide velocity
!    Called by 
!       get_UVW_gs
!    Last update: 10/07/2018, Cheng Zhang
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do i = 1,Nghost
       Ugs(Ibeg-i,j) = -Ugs(Ibeg+i-1,j)
       Vgs(Ibeg-i,j) = Vgs(Ibeg+i-1,j)
	   Wgs(Ibeg-i,j) = Wgs(Ibeg+i-1,j)
       DUgs(Ibeg-i,j) = -DUgs(Ibeg+i-1,j)
       DVgs(Ibeg-i,j) = DVgs(Ibeg+i-1,j)
	   DWgs(Ibeg-i,j) = DWgs(Ibeg+i-1,j)
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
       Ugs(Iend+i,j) = -Ugs(Iend-i+1,j)
       Vgs(Iend+i,j) = Vgs(Iend-i+1,j)
	   Wgs(Iend+i,j) = Wgs(Iend-i+1,j)
       DUgs(Iend+i,j) = -DUgs(Iend-i+1,j)
       DVgs(Iend+i,j) = DVgs(Iend-i+1,j)
	   DWgs(Iend+i,j) = DWgs(Iend-i+1,j)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do j = 1,Nghost
       Ugs(i,Jbeg-j) = Ugs(i,Jbeg+j-1)
       Vgs(i,Jbeg-j) = -Vgs(i,Jbeg+j-1)
	   Wgs(i,Jbeg-j) = Wgs(i,Jbeg+j-1)
       DUgs(i,Jbeg-j) = DUgs(i,Jbeg+j-1)
       DVgs(i,Jbeg-j) = -DVgs(i,Jbeg+j-1)
	   DWgs(i,Jbeg-j) = DWgs(i,Jbeg+j-1)
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
       Ugs(i,Jend+j) = Ugs(i,Jend-j+1)
       Vgs(i,Jend+j) = -Vgs(i,Jend-j+1)
	   Wgs(i,Jend+j) = Wgs(i,Jend-j+1)
       DUgs(i,Jend+j) = DUgs(i,Jend-j+1)
       DVgs(i,Jend+j) = -DVgs(i,Jend-j+1)
	   DWgs(i,Jend+j) = DWgs(i,Jend-j+1)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     end subroutine vel_bc_gs
	 
     subroutine poisson_solver_gs
!--------------------------------------------
!    Solve poisson equation for dynamic pressure
!    of deformable slide
!    Called by
!       eval_duvw_gs
!    Last update: 10/07/2018, Cheng Zhang
!----------------------------------------------
     use global
     implicit none
     integer :: i,j,imask
# if !defined (PARALLEL)
!     ! variables for serial computation
!     real(SP), dimension(:), allocatable :: Wksp
!     integer,  dimension(:), allocatable :: IWksp
!     real(SP), dimension(neqns) :: Phi
!     real(SP) :: RPARM(30),Pbar(1)
!     integer :: IPARM(30),S(1),IS(1),nwksp,inwksp,Ndim,Mdim,N,Maxnz,ierr,neq
!     external :: MIC3,IC3,SOR3,GMRES,CG,BCGS
# endif

     ! Obtain dynamic pressure at slide surface
	 if (NON_HYDRO_UP) then
       do j = 1,Nloc
       do i = 1,Mloc
         Qtgs(i,j) = P(i,j,kbeg)*Mask(i,j)! dynamic pressure over the slide
       enddo
       enddo
	 else
	   Qtgs=Zero
	 endif
   
     ! generate coefficient matrix and rhs
     call generate_coef_rhs_gs

# if defined (PARALLEL)
     ! use HYPRE package for parallel computation
     call hypre_pres_solver_gs
# else
     ! use NSPCG package for serial computation (waiting for developement)
!     call dfault(IPARM,RPARM)
!
!     ! reset default values
!     IPARM(2) = itmax
!     IPARM(3) = 3
!     IPARM(4) = 33
!     RPARM(1) = tol
!
!     Ndim = 5*neqns
!     Mdim = 5*15
!     N = neqns
!     Maxnz = 15
!     nwksp = 30*neqns
!     inwksp = 10*neqns
!
!     allocate(Wksp(nwksp))
!     allocate(Iwksp(inwksp))
!   
!     ! initial guess
!     neq = 0
!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       neq = neq+1
!       Phi(neq) = P(i,j,k)
!     enddo
!     enddo
!     enddo
!     
!     ! solve poisson equation using preconditioned CG or GMRES
!     if(isolver==1) then
!       call nspcg(MIC3,CG,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
!           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
!     elseif(isolver==2) then
!       call nspcg(IC3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
!           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)  
!     elseif(isolver==3) then
!       call nspcg(SOR3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
!           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
!     endif   
!
!     neq = 0
!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       neq = neq+1
!       P(i,j,k) = Phi(neq)
!     enddo
!     enddo
!     enddo
!
!     deallocate(Wksp)
!     deallocate(Iwksp)
# endif

     ! collect into ghost cell corners
     call phiC_2D_coll(QbgsC)

     ! interpolate cell corner results to cell center
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Qbgs(i,j)=0.25*(QbgsC(i+1,j+1)+QbgsC(i,j+1)+QbgsC(i+1,j)+QbgsC(i,j))
     enddo
     enddo 

     ! collect into ghost cells
     call phi_2D_coll(Qbgs)

     end subroutine poisson_solver_gs
	 
     subroutine generate_coef_rhs_gs
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    (deformable slide)
!    Called by 
!       poisson_solver_gs
!    Last update: 10/07/2018, Cheng Zhang
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,neq,n,ic
     real(SP), dimension(:,:), allocatable :: DelxQtgs,DelyQtgs
     integer,  dimension(:,:), allocatable :: indx
	 REAL(SP) :: DgsC,DelxH0C,DelyH0C,QtgsC,NHC,&
	             DelxQtgsR,QtgsR,DelxHgsR,DelxDgsR,DgsR,NHR, &
	             DelxQtgsL,QtgsL,DelxHgsL,DelxDgsL,DgsL,NHL, &
				 DelyQtgsB,QtgsB,DelyHgsB,DelyDgsB,DgsB,NHB, &
				 DelyQtgsF,QtgsF,DelyHgsF,DelyDgsF,DgsF,NHF, &
				 QC,DelxH0R,QR,DelxH0L,QL,DelyH0B,QB,DelyH0F,QF

     allocate(DelxQtgs(Mloc,Nloc))
     allocate(DelyQtgs(Mloc,Nloc))
	 allocate(indx(MlocC,NlocC))
	 
	 DelxQtgs=Zero;DelyQtgs=Zero
	 call DelxH_gs
     call delxFun_nl_gs_2D(Dgs,DelxDgs)
	 call delyFun_nl_gs_2D(Dgs,DelyDgs)
     call delxFun_nl_gs_2D(Qtgs,DelxQtgs)
	 call delyFun_nl_gs_2D(Qtgs,DelyQtgs)
     call delxFun_nl_gs_2D(Hgs,DelxHgs)
	 call delyFun_nl_gs_2D(Hgs,DelyHgs)
	 
     ! generate coefficient matrix
     neq = 0; indx=0
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       neq = neq+1
       indx(i,j) = neq
     enddo
     enddo 

     ! generate source term 
     Rhsgs = Zero; WgsC = Zero
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       DgsC=0.25*(Dgs(i,j)+Dgs(i-1,j)+Dgs(i,j-1)+Dgs(i-1,j-1))
       DelxH0C=1.0/(2.0*dx)*(Hc0(i,j)-Hc0(i-1,j)+Hc0(i,j-1)-Hc0(i-1,j-1))
       DelyH0C=1.0/(2.0*dy)*(Hc0(i,j)+Hc0(i-1,j)-Hc0(i,j-1)-Hc0(i-1,j-1))
	   !NHR
       DelxQtgsR=0.5*(DelxQtgs(i,j)+DelxQtgs(i,j-1))
       QtgsR=0.5*(Qtgs(i,j)+Qtgs(i,j-1))
       DelxHgsR=0.5*(DelxHgs(i,j)+DelxHgs(i,j-1))
       DelxDgsR=0.5*(DelxDgs(i,j)+DelxDgs(i,j-1))
       DgsR=0.5*(Dgs(i,j)+Dgs(i,j-1))
       NHR=-(dt/DgsR)*((1.0/3.0)*DgsR*DelxQtgsR+(1.0/3.0)*QtgsR*DelxDgsR+QtgsR*DelxHgsR)
	   !NHL
       DelxQtgsL=0.5*(DelxQtgs(i-1,j)+DelxQtgs(i-1,j-1))
       QtgsL=0.5*(Qtgs(i-1,j)+Qtgs(i-1,j-1))
       DelxHgsL=0.5*(DelxHgs(i-1,j)+DelxHgs(i-1,j-1))
       DelxDgsL=0.5*(DelxDgs(i-1,j)+DelxDgs(i-1,j-1))
       DgsL=0.5*(Dgs(i-1,j)+Dgs(i-1,j-1))
       NHL=-(dt/DgsL)*((1.0/3.0)*DgsL*DelxQtgsL+(1.0/3.0)*QtgsL*DelxDgsL+QtgsL*DelxHgsL)
	   !NHB
       DelyQtgsB=0.5*(DelyQtgs(i,j)+DelyQtgs(i-1,j))
       QtgsB=0.5*(Qtgs(i,j)+Qtgs(i-1,j))
       DelyHgsB=0.5*(DelyHgs(i,j)+DelyHgs(i-1,j))
       DelyDgsB=0.5*(DelyDgs(i,j)+DelyDgs(i-1,j))
       DgsB=0.5*(Dgs(i,j)+Dgs(i-1,j))
       NHB=-(dt/DgsB)*((1.0/3.0)*DgsB*DelyQtgsB+(1.0/3.0)*QtgsB*DelyDgsB+QtgsB*DelyHgsB)
	   !NHF
       DelyQtgsF=0.5*(DelyQtgs(i,j-1)+DelyQtgs(i-1,j-1))
       QtgsF=0.5*(Qtgs(i,j-1)+Qtgs(i-1,j-1))
       DelyHgsF=0.5*(DelyHgs(i,j-1)+DelyHgs(i-1,j-1))
       DelyDgsF=0.5*(DelyDgs(i,j-1)+DelyDgs(i-1,j-1))
       DgsF=0.5*(Dgs(i,j-1)+Dgs(i-1,j-1))
       NHF=-(dt/DgsF)*((1.0/3.0)*DgsF*DelyQtgsF+(1.0/3.0)*QtgsF*DelyDgsF+QtgsF*DelyHgsF)
	   !NHC
       QtgsC=0.25*(Qtgs(i,j)+Qtgs(i-1,j)+Qtgs(i,j-1)+Qtgs(i-1,j-1))
       NHC=-(2.0*dt/DgsC)*QtgsC
	   WgsC(i,j)=0.25*(Wgs(i,j)+Wgs(i-1,j)+Wgs(i,j-1)+Wgs(i-1,j-1))
       Rhsgs(indx(i,j))=(-0.5*(1./dx+DelxH0C/DgsC)*(Rhogs(i,j)*(Ugs(i,j)  +Ugs(i,j-1)  )+2.*NHR) &
                         +0.5*(1./dx-DelxH0C/DgsC)*(Rhogs(i,j)*(Ugs(i-1,j)+Ugs(i-1,j-1))+2.*NHL) &
                         -0.5*(1./dy+DelyH0C/DgsC)*(Rhogs(i,j)*(Vgs(i,j)  +Vgs(i-1,j)  )+2.*NHB) &
                         +0.5*(1./dy-DelyH0C/DgsC)*(Rhogs(i,j)*(Vgs(i,j-1)+Vgs(i-1,j-1))+2.*NHF) &
                         -1./DgsC*(2.*Rhogs(i,j)*WgsC(i,j)+NHC))*DgsC!multiply DgsC to have better matrix property			 
     enddo
     enddo

	 ! generate coefficient matrix
     Coefgs = Zero
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       DgsC=0.25*(Dgs(i,j)+Dgs(i-1,j)+Dgs(i,j-1)+Dgs(i-1,j-1))
       DelxH0C=1.0/(2.0*dx)*(Hc0(i,j)-Hc0(i-1,j)+Hc0(i,j-1)-Hc0(i-1,j-1))
       DelyH0C=1.0/(2.0*dy)*(Hc0(i,j)+Hc0(i-1,j)-Hc0(i,j-1)-Hc0(i-1,j-1))
       !QR
       DelxDgsR=0.5*(DelxDgs(i,j)+DelxDgs(i,j-1))
       DelxH0R=0.5*(DelxH0(i,j)+DelxH0(i,j-1))
       DgsR=0.5*(Dgs(i,j)+Dgs(i,j-1))
       QR=(-dt*(1./dx+DelxH0C/DgsC)*(2./(3.*dx)+(1./3.)*(DelxDgsR/DgsR)-0.5*(DelxH0R/DgsR)))*DgsC!multiply DgsC to have better matrix property
       !QL
       DelxDgsL=0.5*(DelxDgs(i-1,j)+DelxDgs(i-1,j-1))
       DelxH0L=0.5*(DelxH0(i-1,j)+DelxH0(i-1,j-1))
       DgsL=0.5*(Dgs(i-1,j)+Dgs(i-1,j-1))
       QL=(-dt*(1./dx-DelxH0C/DgsC)*(2./(3.*dx)-(1./3.)*(DelxDgsL/DgsL)+0.5*(DelxH0L/DgsL)))*DgsC
       !QB
       DelyDgsB=0.5*(DelyDgs(i,j)+DelyDgs(i-1,j))
       DelyH0B=0.5*(DelyH0(i,j)+DelyH0(i-1,j))
       DgsB=0.5*(Dgs(i,j)+Dgs(i-1,j))
       QB=(-dt*(1./dy+DelyH0C/DgsC)*(2./(3.*dy)+(1./3.)*(DelyDgsB/DgsB)-0.5*(DelyH0B/DgsB)))*DgsC
       !QF
       DelyDgsF=0.5*(DelyDgs(i,j-1)+DelyDgs(i-1,j-1))
       DelyH0F=0.5*(DelyH0(i,j-1)+DelyH0(i-1,j-1))
       DgsF=0.5*(Dgs(i,j-1)+Dgs(i-1,j-1))
       QF=(-dt*(1./dy-DelyH0C/DgsC)*(2./(3.*dy)-(1./3.)*(DelyDgsF/DgsF)+0.5*(DelyH0F/DgsF)))*DgsC
	   !QC
       QC=(dt*(1./dx+DelxH0C/DgsC)*(2./(3.*dx)-(1./3.)*(DelxDgsR/DgsR)+0.5*(DelxH0R/DgsR)) &
          +dt*(1./dx-DelxH0C/DgsC)*(2./(3.*dx)+(1./3.)*(DelxDgsL/DgsL)-0.5*(DelxH0L/DgsL)) &
          +dt*(1./dy+DelyH0C/DgsC)*(2./(3.*dy)-(1./3.)*(DelyDgsB/DgsB)+0.5*(DelyH0B/DgsB)) &
          +dt*(1./dy-DelyH0C/DgsC)*(2./(3.*dy)+(1./3.)*(DelyDgsF/DgsF)-0.5*(DelyH0F/DgsF)) &
          +2.*dt/DgsC**2)*DgsC
       Coefgs(indx(i,j),1) = QC    ! (i,j)
       Coefgs(indx(i,j),2) = QL    ! (i-1,j)
       Coefgs(indx(i,j),3) = QR    ! (i+1,j)
       Coefgs(indx(i,j),4) = QF    ! (i,j-1)
       Coefgs(indx(i,j),5) = QB    ! (i,j+1)
     enddo
     enddo

!     ! boundary conditions at masks face
!     do i = Ibeg+1,Iend-1
!     do j = Jbeg+1,Jend-1
!       if(Maskgs(i,j)==0) then
!         ! left 
!         if(Maskgs(i+1,j)==1) then
!           ic = indx(i+1,j)
!           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,2)
!           Coefgs(ic,2) = Zero
!         ! right 
!         elseif(Maskgs(i-1,j)==1) then
!           ic = indx(i-1,j)
!           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,3)
!           Coefgs(ic,3) = Zero
!         ! south
!         elseif(Maskgs(i,j+1)==1) then
!           ic = indx(i,j+1)
!           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,4)
!           Coefgs(ic,4) = Zero
!         ! north
!         elseif(Maskgs(i,j-1)==1) then
!           ic = indx(i,j-1)
!           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,5)
!           Coefgs(ic,5) = Zero
!         endif ! end Maskgs+1=1 
!       endif ! end Maskgs=0
!     enddo
!     enddo

     ! boundary conditions
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     i = Ibeg
     do j = Jbeg,JendC
       ic = indx(i,j)
       Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,2)
       Coefgs(ic,2) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! right side
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     i = IendC
     do j = Jbeg,JendC
       ic = indx(i,j)
       Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,3)
       Coefgs(ic,3) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     j = Jbeg
     do i = Ibeg,IendC
       ic = indx(i,j)         
       Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,4)
       Coefgs(ic,4) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! back side
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     j = JendC
     do i = Ibeg,IendC
       ic = indx(i,j)
       Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,5)
       Coefgs(ic,5) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! left front corner
# if defined (PARALLEL)
     if((n_west.eq.MPI_PROC_NULL).and.(n_suth.eq.MPI_PROC_NULL)) then
# endif
     i = Ibeg
     j = Jbeg
     ic = indx(i,j)
     Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,2)+Coefgs(ic,4)
     Coefgs(ic,2) = Zero
	 Coefgs(ic,4) = Zero
# if defined (PARALLEL)
     endif
# endif

     ! right front corner
# if defined (PARALLEL)
     if((n_east.eq.MPI_PROC_NULL).and.(n_suth.eq.MPI_PROC_NULL)) then
# endif
     i = IendC
     j = Jbeg
     ic = indx(i,j)
     Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,3)+Coefgs(ic,4)
     Coefgs(ic,3) = Zero
	 Coefgs(ic,4) = Zero
# if defined (PARALLEL)
     endif
# endif

     ! left back corner
# if defined (PARALLEL)
     if((n_west.eq.MPI_PROC_NULL).and.(n_nrth.eq.MPI_PROC_NULL)) then
# endif
     i = Ibeg
     j = JendC
     ic = indx(i,j)         
     Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,2)+Coefgs(ic,5)
     Coefgs(ic,2) = Zero
	 Coefgs(ic,5) = Zero
# if defined (PARALLEL)
     endif
# endif

     ! right back side
# if defined (PARALLEL)
     if((n_east.eq.MPI_PROC_NULL).and.(n_nrth.eq.MPI_PROC_NULL)) then
# endif
     i = IendC
     j = JendC
     ic = indx(i,j)
     Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,3)+Coefgs(ic,5)
     Coefgs(ic,3) = Zero
     Coefgs(ic,5) = Zero
# if defined (PARALLEL)
     endif
# endif

     ! take (i=2,j=2) to obtain the diagonal information (Serial Computing used)
     JCoefgs(1) = indx(Ibeg+1,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i,j)
     JCoefgs(2) = indx(Ibeg  ,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i-1,j)
     JCoefgs(3) = indx(Ibeg+2,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i+1,j) 
	 JCoefgs(4) = indx(Ibeg+1,Jbeg  )-indx(Ibeg+1,Jbeg+1)  ! (i,j-1)
     JCoefgs(5) = indx(Ibeg+1,Jbeg+2)-indx(Ibeg+1,Jbeg+1)  ! (i,j+1)

     deallocate(DelxQtgs)
     deallocate(DelyQtgs)
	 deallocate(indx)

     return
     end subroutine generate_Coef_rhs_gs

     subroutine hypre_pres_solver_gs
!---------------------------------------------
!    solve pressure using hypre package 
!    (deformable slide)
!    called by
!       poisson_solver_gs
!    Last update: 10/07/2018, Cheng Zhang
!---------------------------------------------
# if defined (PARALLEL)
     use global
     implicit none
     integer, parameter :: ndim=2
     integer, parameter :: nentries=5
     integer :: i,j,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(MlocC),j_glob(NlocC)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,-1,0,1,0,0,-1,0,1/

     ! set up a two dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     iupper(1) = i_glob(IendC)
     iupper(2) = j_glob(JendC)

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)

     call HYPRE_StructGridAssemble(grid,ierr)

     ! define the discretization stencil
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)

     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

     ! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)

     call HYPRE_StructMatrixInitialize(matrix,ierr)

     ! set the matrix Coefficient
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(MlocC,NlocC))
 
     neq = 0
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       neq = neq+1
       indx(i,j) = neq
     enddo
     enddo
    
     nvalues = (IendC-Ibeg+1)*(JendC-Jbeg+1)*nentries !(grid points)*(stencil points)
     allocate(values(nvalues))

     ivalues = 0
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coefgs(indx(i,j),n)
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

     ! set the vector Coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhsgs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqnsgs))
     do j = Jbeg,JendC
     do i = Ibeg,IendC
       Phi(indx(i,j)) = QbgsC(i,j)
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

     do j = Jbeg,JendC
     do i = Ibeg,IendC
       QbgsC(i,j) = Phi(indx(i,j))
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
     end subroutine hypre_pres_solver_gs
	 
     subroutine projection_corrector_gs
!-------------------------------------------
!    Correct the velocity field for slide
!    Called by
!       eval_duvw_gs
!    Last update: 10/07/2018, Cheng Zhang
!-------------------------------------------
     use global
     implicit none
     integer :: i,j
     real(SP), dimension(:,:),allocatable :: DelxQtgs,DelyQtgs,WtgsC,WbgsC,DIV, &
	                                         UgsR,UgsL,VgsB,VgsF
	 REAL(SP) :: UgsC,VgsC,DelxH0C,DelyH0C,DgsC,QtgsC,MAXDIV, &
	             DelxQtgsR,DelxQbgsR,QtgsR,QbgsR,DelxHgsR,DelxH0R,DelxDgsR,DgsR,NHTR,NHBR, &
				 DelxQtgsL,DelxQbgsL,QtgsL,QbgsL,DelxHgsL,DelxH0L,DelxDgsL,DgsL,NHTL,NHBL, &
				 DelyQtgsB,DelyQbgsB,QtgsB,QbgsB,DelyHgsB,DelyH0B,DelyDgsB,DgsB,NHTB,NHBB, &
				 DelyQtgsF,DelyQbgsF,QtgsF,QbgsF,DelyHgsF,DelyH0F,DelyDgsF,DgsF,NHTF,NHBF
# if defined (PARALLEL)
     real(SP) :: myvar
# endif
 
     allocate(DelxQtgs(Mloc,Nloc))
     allocate(DelyQtgs(Mloc,Nloc))
     allocate(WtgsC(MlocC,NlocC))
     allocate(WbgsC(MlocC,NlocC))
	 allocate(DIV(MlocC,NlocC))
	 allocate(UgsR(MlocC,NlocC))
	 allocate(UgsL(MlocC,NlocC))
	 allocate(VgsB(MlocC,NlocC))
	 allocate(VgsF(MlocC,NlocC))
	 
	 DelxQtgs=Zero;DelyQtgs=Zero
	 call DelxH_gs
     call delxFun_nl_gs_2D(Dgs,DelxDgs)
	 call delyFun_nl_gs_2D(Dgs,DelyDgs)
     call delxFun_nl_gs_2D(Qtgs,DelxQtgs)
	 call delyFun_nl_gs_2D(Qtgs,DelyQtgs)
     call delxFun_nl_gs_2D(Hgs,DelxHgs)
	 call delyFun_nl_gs_2D(Hgs,DelyHgs)
	 
	 ! update Ugs and Vgs at cell face, Wgs at cell corner
	 WtgsC=Zero;WbgsC=Zero;WgsC=Zero
	 UgsR=Zero;UgsL=Zero;VgsB=Zero;VgsF=Zero
     do j = Jbeg,JendC
	 do i = Ibeg,IendC
	     !UgsR
         DelxQtgsR=0.5*(DelxQtgs(i,j)+DelxQtgs(i,j-1))
		 DelxQbgsR=(QbgsC(i+1,j)-QbgsC(i,j))/dx
         QtgsR=0.5*(Qtgs(i,j)+Qtgs(i,j-1))
		 QbgsR=0.5*(QbgsC(i+1,j)+QbgsC(i,j))
         DelxHgsR=0.5*(DelxHgs(i,j)+DelxHgs(i,j-1))
		 DelxH0R=0.5*(DelxH0(i,j)+DelxH0(i,j-1))
         DelxDgsR=0.5*(DelxDgs(i,j)+DelxDgs(i,j-1))
         DgsR=0.5*(Dgs(i,j)+Dgs(i,j-1))
         NHTR=-(dt/DgsR/Rhogs(i,j))*((1.0/3.0)*DgsR*DelxQtgsR+(1.0/3.0)*QtgsR*DelxDgsR+QtgsR*DelxHgsR)
		 NHBR=-(dt/DgsR/Rhogs(i,j))*((2.0/3.0)*DgsR*DelxQbgsR+(2.0/3.0)*QbgsR*DelxDgsR-QbgsR*DelxH0R)
		 UgsR(i,j)=0.5*(Ugs(i,j)+Ugs(i,j-1))+NHTR+NHBR
	     !UgsL
         DelxQtgsL=0.5*(DelxQtgs(i-1,j)+DelxQtgs(i-1,j-1))
		 DelxQbgsL=(QbgsC(i,j)-QbgsC(i-1,j))/dx
         QtgsL=0.5*(Qtgs(i-1,j)+Qtgs(i-1,j-1))
		 QbgsL=0.5*(QbgsC(i,j)+QbgsC(i-1,j))
         DelxHgsL=0.5*(DelxHgs(i-1,j)+DelxHgs(i-1,j-1))
		 DelxH0L=0.5*(DelxH0(i-1,j)+DelxH0(i-1,j-1))
         DelxDgsL=0.5*(DelxDgs(i-1,j)+DelxDgs(i-1,j-1))
         DgsL=0.5*(Dgs(i-1,j)+Dgs(i-1,j-1))
         NHTL=-(dt/DgsL/Rhogs(i,j))*((1.0/3.0)*DgsL*DelxQtgsL+(1.0/3.0)*QtgsL*DelxDgsL+QtgsL*DelxHgsL)
		 NHBL=-(dt/DgsL/Rhogs(i,j))*((2.0/3.0)*DgsL*DelxQbgsL+(2.0/3.0)*QbgsL*DelxDgsL-QbgsL*DelxH0L)
		 UgsL(i,j)=0.5*(Ugs(i-1,j)+Ugs(i-1,j-1))+NHTL+NHBL
	     !VgsB
         DelyQtgsB=0.5*(DelyQtgs(i,j)+DelyQtgs(i-1,j))
		 DelyQbgsB=(QbgsC(i,j+1)-QbgsC(i,j))/dy
         QtgsB=0.5*(Qtgs(i,j)+Qtgs(i-1,j))
		 QbgsB=0.5*(QbgsC(i,j+1)+QbgsC(i,j))
         DelyHgsB=0.5*(DelyHgs(i,j)+DelyHgs(i-1,j))
		 DelyH0B=0.5*(DelyH0(i,j)+DelyH0(i-1,j))
         DelyDgsB=0.5*(DelyDgs(i,j)+DelyDgs(i-1,j))
         DgsB=0.5*(Dgs(i,j)+Dgs(i-1,j))
         NHTB=-(dt/DgsB/Rhogs(i,j))*((1.0/3.0)*DgsB*DelyQtgsB+(1.0/3.0)*QtgsB*DelyDgsB+QtgsB*DelyHgsB)
		 NHBB=-(dt/DgsB/Rhogs(i,j))*((2.0/3.0)*DgsB*DelyQbgsB+(2.0/3.0)*QbgsB*DelyDgsB-QbgsB*DelyH0B)
		 VgsB(i,j)=0.5*(Vgs(i,j)+Vgs(i-1,j))+NHTB+NHBB
	     !VgsF
         DelyQtgsF=0.5*(DelyQtgs(i,j-1)+DelyQtgs(i-1,j-1))
		 DelyQbgsF=(QbgsC(i,j)-QbgsC(i,j-1))/dy
         QtgsF=0.5*(Qtgs(i,j-1)+Qtgs(i-1,j-1))
		 QbgsF=0.5*(QbgsC(i,j)+QbgsC(i,j-1))
         DelyHgsF=0.5*(DelyHgs(i,j-1)+DelyHgs(i-1,j-1))
		 DelyH0F=0.5*(DelyH0(i,j-1)+DelyH0(i-1,j-1))
         DelyDgsF=0.5*(DelyDgs(i,j-1)+DelyDgs(i-1,j-1))
         DgsF=0.5*(Dgs(i,j-1)+Dgs(i-1,j-1))
         NHTF=-(dt/DgsF/Rhogs(i,j))*((1.0/3.0)*DgsF*DelyQtgsF+(1.0/3.0)*QtgsF*DelyDgsF+QtgsF*DelyHgsF)
		 NHBF=-(dt/DgsF/Rhogs(i,j))*((2.0/3.0)*DgsF*DelyQbgsF+(2.0/3.0)*QbgsF*DelyDgsF-QbgsF*DelyH0F)
		 VgsF(i,j)=0.5*(Vgs(i,j-1)+Vgs(i-1,j-1))+NHTF+NHBF
	     !WbgsC&WtgsC
         UgsC=0.5*(UgsR(i,j)+UgsL(i,j))
         VgsC=0.5*(VgsB(i,j)+VgsF(i,j))
         WgsC(i,j)=0.25*(Wgs(i,j)+Wgs(i-1,j)+Wgs(i,j-1)+Wgs(i-1,j-1))
         DelxH0C=1./(2.*dx)*(Hc0(i,j)-Hc0(i-1,j)+Hc0(i,j-1)-Hc0(i-1,j-1))
         DelyH0C=1./(2.*dy)*(Hc0(i,j)+Hc0(i-1,j)-Hc0(i,j-1)-Hc0(i-1,j-1))
         DgsC=0.25*(Dgs(i,j)+Dgs(i-1,j)+Dgs(i,j-1)+Dgs(i-1,j-1))
         QtgsC=0.25*(Qtgs(i,j)+Qtgs(i-1,j)+Qtgs(i,j-1)+Qtgs(i-1,j-1))   
         WbgsC(i,j)=-UgsC*DelxH0C-VgsC*DelyH0C
         WtgsC(i,j)=2.*WgsC(i,j)-WbgsC(i,j)-2.*dt/SLIDE_DENSITY/DgsC*(QtgsC-QbgsC(i,j))
         WgsC(i,j)=0.5*(WtgsC(i,j)+WbgsC(i,j))
     enddo
     enddo
	 call phiC_2D_coll(UgsR)
	 call phiC_2D_coll(UgsL)
	 call phiC_2D_coll(VgsB)
	 call phiC_2D_coll(VgsF)
	 call phiC_2D_coll(WbgsC)
	 call phiC_2D_coll(WtgsC)
	 call phiC_2D_coll(WgsC)
	 
	 ! Obtain Ugs, Vgs and Wgs at cell center by interpolations
     do j = Jbeg,Jend
	 do i = Ibeg,Iend
         Ugs(i,j)=0.5*(UgsR(i,j+1)+UgsR(i,j))
         DUgs(i,j)=Ugs(i,j)*Dgs(i,j)
         Vgs(i,j)=0.5*(VgsB(i+1,j)+VgsB(i,j))
         DVgs(i,j)=Vgs(i,j)*Dgs(i,j)
	     Wtgs(i,j)=0.25*(WtgsC(i+1,j+1)+WtgsC(i,j+1)+WtgsC(i+1,j)+WtgsC(i,j))
		 Wbgs(i,j)=0.25*(WbgsC(i+1,j+1)+WbgsC(i,j+1)+WbgsC(i+1,j)+WbgsC(i,j))
         Wgs(i,j)=0.25*(WgsC(i+1,j+1)+WgsC(i,j+1)+WgsC(i+1,j)+WgsC(i,j))
         DWgs(i,j) = Wgs(i,j)*Dgs(i,j)
     enddo
     enddo
	 call phi_2D_coll(Wtgs)
	 call phi_2D_coll(Wbgs)
	 
	 ! collect data into ghost cells
     call vel_bc_gs
# if defined (PARALLEL)
     call phi_2D_exch(Ugs)
     call phi_2D_exch(Vgs)
	 call phi_2D_exch(Wgs)
     call phi_2D_exch(DUgs)
     call phi_2D_exch(DVgs)
	 call phi_2D_exch(DWgs)
# endif

     ! check divergence free in cell corner
	 DIV=Zero
     do j = Jbeg,JendC
	 do i = Ibeg,IendC
            DgsC=0.25*(Dgs(i,j)+Dgs(i-1,j)+Dgs(i,j-1)+Dgs(i-1,j-1))
            DIV(i,j)=(UgsR(i,j)-UgsL(i,j))/dx+(VgsB(i,j)-VgsF(i,j))/dy+(WtgsC(i,j)-WbgsC(i,j))/DgsC
     enddo
     enddo
	 MAXDIV = MAXVAL(abs(DIV(Ibeg:IendC,Jbeg:JendC)))
# if defined (PARALLEL)
     call MPI_ALLREDUCE(MAXDIV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MAXDIV = myvar
# endif
     !print
# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     WRITE(*,*),'MAXDIV:  ',MAXDIV
# if defined (PARALLEL)
     endif
# endif
 
     deallocate(DelxQtgs)
     deallocate(DelyQtgs)
     deallocate(WtgsC)
	 deallocate(WbgsC)
	 deallocate(DIV)
	 deallocate(UgsR)
	 deallocate(UgsL)
	 deallocate(VgsB)
	 deallocate(VgsF)

     return
     end subroutine projection_corrector_gs

# if defined (OBSTACLE)
    subroutine imm_obs_gs
!---------------------------------------------------------------------
!   Immersed-Boundary Method
!   Calculate forcing at the obstacle boundary
!   by Cheng Zhang, 21/07/2020
!--------------------------------------------------------------------
    use global
    implicit none
    real(SP),dimension(:),allocatable :: target_vel,target_dist
    integer :: i,j,l,count
    real(SP) :: numer,denom,force_vel,dista,vel_grad

    ! calculate forcing at boundaries
    ObsForceX_gs = Zero
    ObsForceY_gs = Zero
    ObsForceZ_gs = Zero

    allocate(target_vel(1:4))
    allocate(target_dist(1:4))

    !no-slip
    ! x-direction
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for slide cells
      if(set_flag_gs(i,j)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j) is slide, (i-1,j) is obstacle
        if(set_flag_gs(i-1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i+1,j) to obstacle/slide interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (Ugs(i+1,j)-obs_u)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i+1,j) is obstacle
        if(set_flag_gs(i+1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j) 
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i-1,j) to obstacle/slide interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_u-Ugs(i-1,j))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j-1) is obstacle
        if(Nloc>=2.and.set_flag_gs(i,j-1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j+1) to obstacle/slide interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (Ugs(i,j+1)-obs_u)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j+1) is obstacle 
        if((Nloc>=2).and.set_flag_gs(i,j+1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j-1) to obstacle/slide interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_u-Ugs(i,j-1))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(yc(j)-yc(j-1))
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
          ObsForceX_gs(i,j) = (Dgs(i,j)*force_vel-DUgs(i,j))/dt
        endif
      endif
    enddo
    enddo

    ! y-direction
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for slide cells
      if(set_flag_gs(i,j)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j) is slide, (i-1,j) is obstacle
        if(set_flag_gs(i-1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i+1,j) to obstacle/slide interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (Vgs(i+1,j)-obs_v)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i+1,j) is obstacle
        if(set_flag_gs(i+1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j) 
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i-1,j) to obstacle/slide interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_v-Vgs(i-1,j))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j-1) is obstacle
        if(Nloc>=2.and.set_flag_gs(i,j-1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j+1) to obstacle/slide interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (Vgs(i,j+1)-obs_v)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j+1) is obstacle 
        if((Nloc>=2).and.set_flag_gs(i,j+1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j-1) to obstacle/slide interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_v-Vgs(i,j-1))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(yc(j)-yc(j-1))
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
          ObsForceY_gs(i,j) = (Dgs(i,j)*force_vel-DVgs(i,j))/dt
        endif
      endif
    enddo
    enddo
 
    ! z-direction
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for slide cells
      if(set_flag_gs(i,j)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j) is slide, (i-1,j) is obstacle
        if(set_flag_gs(i-1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i+1,j) to obstacle/slide interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (Wgs(i+1,j)-obs_w)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i+1,j) is obstacle
        if(set_flag_gs(i+1,j)==1.and.set_dist_x_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j) 
          dista = set_dist_x_gs(i,j)
          ! compute velocity gradient from (i-1,j) to obstacle/slide interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_w-Wgs(i-1,j))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j-1) is obstacle
        if(Nloc>=2.and.set_flag_gs(i,j-1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j+1) to obstacle/slide interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (Wgs(i,j+1)-obs_w)/dista
          ! interpolate velocity from obstacle/slide interface to (i,j) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j) is slide, (i,j+1) is obstacle 
        if((Nloc>=2).and.set_flag_gs(i,j+1)==1.and.set_dist_y_gs(i,j)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j)
          dista = set_dist_y_gs(i,j)
          ! compute velocity gradient from (i,j-1) to obstacle/slide interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_w-Wgs(i,j-1))/dista
          ! interpolate velocity from obstacle/slide interface to (i,j)
          dista = dista-(yc(j)-yc(j-1))
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
          ObsForceZ_gs(i,j) = (Dgs(i,j)*force_vel-DWgs(i,j))/dt
        endif
      endif
    enddo
    enddo
	
    deallocate(target_vel)
    deallocate(target_dist)
 1000 continue

    return
    end subroutine imm_obs_gs
# endif
	 
     subroutine delxFun_nl_gs_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x without limiter
!    Called by
!       poisson_solver_gs,projection_corrector_gs
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
!       if(Maskgs(i,j)==0) then
!         Dout(i,j) = Zero
!       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1+TMP2)/2.0
!         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_nl_gs_2D
	 
     subroutine delyFun_nl_gs_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y without limiter
!    Called by 
!       poisson_solver_gs,projection_corrector_gs
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
!       if(Maskgs(i,j)==0) then 
!         Dout(i,j) = Zero
!       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1+TMP2)/2.0
         endif
!       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_nl_gs_2D
	 
     subroutine delx2Fun_gs_2D(Din,Dout)
!-------------------------------------------
!    Second derivative in x
!    Called by
!       generate_Coefgs_rhs_gs
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Maskgs(i,j)==0) then
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1-TMP2)/dx
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = Dout(2,j)
       Dout(Mloc,j) = Dout(Mloc-1,j)
     enddo  

     return
     end subroutine delx2Fun_gs_2D
	 
     subroutine dely2Fun_gs_2D(Din,Dout)
!-----------------------------------------
!    Second derivative in y
!    Called by 
!       generate_Coefgs_rhs_gs
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Maskgs
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Maskgs(i,j)==0) then 
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1-TMP2)/dy
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = Dout(i,2)
       Dout(i,Nloc) = Dout(i,Nloc)
     enddo

     return
     end subroutine dely2Fun_gs_2D
	 
     subroutine delx2Fun_2D(Din,Dout)
!-------------------------------------------
!    Second derivative in x
!    Called by
!       generate_Coefgs_rhs_gs
!    Last update: 10/07/2018, Cheng Zhang
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
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
           Dout(i,j) = (TMP1-TMP2)/dx
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = Dout(2,j)
       Dout(Mloc,j) = Dout(Mloc-1,j)
     enddo  

     return
     end subroutine delx2Fun_2D
	 
     subroutine dely2Fun_2D(Din,Dout)
!-----------------------------------------
!    Second derivative in y
!    Called by 
!       generate_Coefgs_rhs_gs
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2
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
           Dout(i,j) = (TMP1-TMP2)/dy
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = Dout(i,2)
       Dout(i,Nloc) = Dout(i,Nloc)
     enddo

     return
     end subroutine dely2Fun_2D
	 
     subroutine phiC_2D_coll(phi)
!-----------------------------------------------------
!    This subroutine is used to collect data into ghost cell corners
!    Called by
!       poisson_solver_gs
!    Last update: 28/06/2019, Cheng Zhang
!-----------------------------------------------------
     use global
     implicit none
     real(SP), intent(inout) :: phi(MlocC,NlocC)
     integer :: i,j

     ! x-direction
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,JendC
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
       do j = Jbeg,JendC
       do i = 1,Nghost
         phi(IendC+i,j) = phi(IendC-i+1,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif
 
     ! y-direction and corners
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif     
       do i = 1,MlocC
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
       do i = 1,MlocC
       do j = 1,Nghost
         phi(i,JendC+j) = phi(i,JendC-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phiC_2D_exch(phi)
# endif    

     end subroutine phiC_2D_coll
	 
# if defined (PARALLEL)
    subroutine phiC_2D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(MlocC,NlocC)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    REAL(SP),DIMENSION(MlocC,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(NlocC,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = NlocC * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, NlocC
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
       do j = 1, NlocC
       do i = 1, Nghost
          sEmsg(j,i) = PHI(IendC-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, NlocC
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, NlocC
       do i = 1, Nghost
          PHI(IendC+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = MlocC * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, MlocC
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
       do i = 1, MlocC
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,JendC-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, MlocC
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, MlocC
       do j = 1, Nghost
          PHI(i,JendC+j) = rNmsg(i,j)
       enddo
       enddo
    endif

    return
    END SUBROUTINE phiC_2D_exch
# endif
# endif