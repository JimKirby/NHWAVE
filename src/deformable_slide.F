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
        if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
          Masktmp(i,j) = 0
          Dgs(i,j) = SLIDE_MINTHICK
        endif
      endif
	  Hgs(i,j) = Hc0(i,j)-Dgs(i,j)
    enddo
    enddo
	Maskgs = Masktmp

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
!	 if(trim(RHEO_OPT)=='GRANULAR') then
!	    Rhogs=Zero
!	    do j = Jbeg,Jend
!        do i = Ibeg,Iend
!		  if(Mask(i,j)==1) then
!	        Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC + Rho0*(1.0-SLIDE_CONC)
!		  else
!		    Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC
!		  endif
!		enddo
!		enddo
!		call phi_2D_coll(Rhogs)
!	 endif

     end subroutine update_vars_gs
	 
    subroutine flux_coff_gs
!------------------------------------------------------  
!    This subroutine is used to update coefficients for flux
!    Called by                                                
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!-----------------------------------------------------
    use global
    implicit none
    integer :: i,j
	real(SP),dimension(:,:),allocatable :: DwDt
	real(SP) :: TMP1,TMP2
	 
    allocate(DwDt(Mloc,Nloc))
	DwDt=Zero
    if(trim(RHEO_OPT)=='VISCOUS') then
	  Cxgs = Grav
      Cygs = Grav
	elseif(trim(RHEO_OPT)=='GRANULAR') then
	  call transform_gs
      do j = Jbeg,Jend
      do i = Ibeg,Iend
	    TMP1=SLIDE_CONC*GRAIN_DENSITY/Rhogs(i,j)*Grav
		TMP2=(1.0-SLIDE_CONC*GRAIN_DENSITY/Rhogs(i,j))*Grav
	    ! DwDt(i,j)=(Wgs0(i,j)-Wgs00(i,j))/dt
		! DwDt(i,j)=max(DwDt(i,j),-TMP1)
		! DwDt(i,j)=min(DwDt(i,j),TMP1)
		Cxgs(i,j)=(Cxx(i,j)/max(Czz(i,j),Small)*(   TMP1)+TMP2)*Maskgs(i,j)
		Cygs(i,j)=(Cyy(i,j)/max(Czz(i,j),Small)*(   TMP1)+TMP2)*Maskgs(i,j)
      enddo
      enddo
	endif
# if defined (PARALLEL)
	call phi_2D_exch(Cxgs)
	call phi_2D_exch(Cygs)
# endif
	
	deallocate(DwDt)

    end subroutine flux_coff_gs
	
    subroutine transform_gs
!------------------------------------------------------  
!    This subroutine is used to construct local coordinate
!    Called by                                                
!       flux_coff_gs
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
         if(PhiInt>=PhiBed) then
           Kxa=2./cos(PhiInt)**2*(1.-sqrt(1.-cos(PhiInt)**2/cos(PhiBed)**2))-1.
           Kxp=2./cos(PhiInt)**2*(1.+sqrt(1.-cos(PhiInt)**2/cos(PhiBed)**2))-1.
         else
           Kxa=(1.+sin(PhiInt)**2)/(1.-sin(PhiInt)**2)
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
           if(DvgsDy>1.e-8) then
             Ky(i,j)=0.5*(Kx(i,j)+1.0-sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed)**2))
           elseif(DvgsDy<-1.e-8) then
             Ky(i,j)=0.5*(Kx(i,j)+1.0+sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed)**2))
           else
             Int_factor = (crit-abs(DvgsDy))/(crit-1.e-8)
             if(DvgsDy>1.e-8.and.DvgsDy<crit) then
               Ky(i,j)=Int_factor+(1.0-Int_factor) &
			           *(0.5*(Kx(i,j)+1.0-sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed)**2)))
             elseif(DvgsDy<-1.e-8.and.DvgsDy>-crit) then
               Ky(i,j)=Int_factor+(1.0-Int_factor) &
			           *(0.5*(Kx(i,j)+1.0+sqrt((Kx(i,j)-1.0)**2+4.0*tan(PhiBed)**2)))
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
         Czz0(i,j)=Czz(i,j)-(C33(i,j)*C31(i,j)+C31(i,j)*C33(i,j))*tan(PhiBed)
         if(abs(UgsP(i,j))>1.e-8) then
           Cxx(i,j)=Cxx(i,j)-(C13(i,j)*C11(i,j)+C11(i,j)*C13(i,j))*tan(PhiBed)
           Cyx(i,j)=Cyx(i,j)-(C23(i,j)*C11(i,j)+C21(i,j)*C13(i,j))*tan(PhiBed)
           Cyy(i,j)=Cyy(i,j)-(C23(i,j)*C21(i,j)+C21(i,j)*C23(i,j))*tan(PhiBed)
           Cxz(i,j)=Cxz(i,j)-(C13(i,j)*C31(i,j)+C11(i,j)*C33(i,j))*tan(PhiBed)
           Cyz(i,j)=Cyz(i,j)-(C23(i,j)*C31(i,j)+C21(i,j)*C33(i,j))*tan(PhiBed)
           Czz(i,j)=Czz(i,j)-(C33(i,j)*C31(i,j)+C31(i,j)*C33(i,j))*tan(PhiBed)
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
     real(SP) :: TMP1,TMP2,LIMITER
    
     ! pu'/px'
     TMP1 = (UvecPR(1)-UvecPC(1))/dx
     TMP2 = (UvecPC(1)-UvecPL(1))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXUP = Zero
     else
       DXUP = LIMITER(TMP1,TMP2)
     endif
     ! pu'/py'
     TMP1 = (UvecPB(1)-UvecPC(1))/dy
     TMP2 = (UvecPC(1)-UvecPF(1))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYUP = Zero
     else
       DYUP = LIMITER(TMP1,TMP2)
     endif
     ! pv'/px'
     TMP1 = (UvecPR(2)-UvecPC(2))/dx
     TMP2 = (UvecPC(2)-UvecPL(2))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXVP = Zero
     else
       DXVP = LIMITER(TMP1,TMP2)
     endif
     ! pv'/py'
     TMP1 = (UvecPB(2)-UvecPC(2))/dy
     TMP2 = (UvecPC(2)-UvecPF(2))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYVP = Zero
     else
       DYVP = LIMITER(TMP1,TMP2)
     endif
     !pw'/px'
     TMP1 = (UvecPR(3)-UvecPC(3))/dx
     TMP2 = (UvecPC(3)-UvecPL(3))/dx
     if((abs(TMP1)+abs(TMP2))<Small) then
       DXWP = Zero
     else
       DXWP = LIMITER(TMP1,TMP2)
     endif
     !pw'/py'
     TMP1 = (UvecPB(3)-UvecPC(3))/dy
     TMP2 = (UvecPC(3)-UvecPF(3))/dy
     if((abs(TMP1)+abs(TMP2))<Small) then
       DYWP = Zero
     else
       DYWP = LIMITER(TMP1,TMP2)
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
     call delxFun_gs_2D(Ugs,DelxUgs)
     call delxFun_gs_2D(Vgs,DelxVgs)
     call delxFun_gs_2D(DUgs,DelxDUgs)
     call delxFun_gs_2D(DVgs,DelxDVgs)
	 call delxFun_gs_2D(Cxgs,DelxCxgs)

	 call delyFun_gs_2D(Dgs,DelyDgs)
     call delyFun_gs_2D(Ugs,DelyUgs)
     call delyFun_gs_2D(Vgs,DelyVgs)
     call delyFun_gs_2D(DUgs,DelyDUgs)
     call delyFun_gs_2D(DVgs,DelyDVgs)
	 call delyFun_gs_2D(Cygs,DelyCygs)

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
     real(SP) :: TMP1,TMP2,LIMITER
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
     real(SP) :: TMP1,TMP2,LIMITER
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
     end subroutine delyFun_gs_2D
	 
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
     call construct_2D_x(Ugs,DelxUgs,UgsxL,UgsxR)
     call construct_2D_x(Vgs,DelxVgs,VgsxL,VgsxR)
     call construct_2D_x(DUgs,DelxDUgs,DUgsxL,DUgsxR)
     call construct_2D_x(DVgs,DelxDVgs,DVgsxL,DVgsxR)
	 call construct_2D_x(Cxgs,DelxCxgs,CxgsxL,CxgsxR)

     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       EgsxL(i,j) = DUgsxL(i,j)
       EgsxR(i,j) = DUgsxR(i,j)
       FgsxL(i,j) = SLIDE_GAMMA*DUgsxL(i,j)*UgsxL(i,j)+0.5*CxgsxL(i,j)*DgsxL(i,j)*DgsxL(i,j)
       FgsxR(i,j) = SLIDE_GAMMA*DUgsxR(i,j)*UgsxR(i,j)+0.5*CxgsxR(i,j)*DgsxR(i,j)*DgsxR(i,j)
       GgsxL(i,j) = SLIDE_GAMMA*DUgsxL(i,j)*VgsxL(i,j)
       GgsxR(i,j) = SLIDE_GAMMA*DUgsxR(i,j)*VgsxR(i,j)
     enddo
     enddo

	 call construct_2D_y(Dgs,DelyDgs,DgsyL,DgsyR)
     call construct_2D_y(Ugs,DelyUgs,UgsyL,UgsyR)
     call construct_2D_y(Vgs,DelyVgs,VgsyL,VgsyR)
     call construct_2D_y(DUgs,DelyDUgs,DUgsyL,DUgsyR)
     call construct_2D_y(DVgs,DelyDVgs,DVgsyL,DVgsyR)
	 call construct_2D_y(Cygs,DelyCygs,CygsyL,CygsyR)

     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       EgsyL(i,j) = DVgsyL(i,j)
       EgsyR(i,j) = DVgsyR(i,j)
       FgsyL(i,j) = SLIDE_GAMMA*DVgsyL(i,j)*UgsyL(i,j)
       FgsyR(i,j) = SLIDE_GAMMA*DVgsyR(i,j)*UgsyR(i,j)
       GgsyL(i,j) = SLIDE_GAMMA*DVgsyL(i,j)*VgsyL(i,j)+0.5*CygsyL(i,j)*DgsyL(i,j)*DgsyL(i,j)
       GgsyR(i,j) = SLIDE_GAMMA*DVgsyR(i,j)*VgsyR(i,j)+0.5*CygsyR(i,j)*DgsyR(i,j)*DgsyR(i,j)
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
		                  +CxgsxL(i,j)*abs(DgsxL(i,j)))
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxR(i,j)**2 &
		                  +CxgsxR(i,j)*abs(DgsxR(i,j)))
		 UL = SLIDE_GAMMA*UgsxL(i,j)
		 UR = SLIDE_GAMMA*UgsxR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL-UR)
         U_S = 0.5*(UL+UR)+SQR_PHI_L-SQR_PHI_R
         SgsxL(i,j) = min(UL-SQR_PHI_L,U_S-SQR_PHI_S)
         SgsxR(i,j) = max(UR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskgs(i-1,j)==0.and.Maskgs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxR(i,j)**2 &
		                  +CxgsxR(i,j)*abs(DgsxR(i,j)))
		 UR = SLIDE_GAMMA*UgsxR(i,j)
         SgsxL(i,j) = UR-2.0*SQR_PHI_R
         SgsxR(i,j) = UR+SQR_PHI_R
       elseif(Maskgs(i-1,j)==1.and.Maskgs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*UgsxL(i,j)**2 &
		                  +CxgsxL(i,j)*abs(DgsxL(i,j)))
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
		                  +CygsyL(i,j)*abs(DgsyL(i,j)))
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyR(i,j)**2 &
		                  +CygsyR(i,j)*abs(DgsyR(i,j)))
		 VL = SLIDE_GAMMA*VgsyL(i,j)
		 VR = SLIDE_GAMMA*VgsyR(i,j)
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL-VR)
         U_S = 0.5*(VL+VR)+SQR_PHI_L-SQR_PHI_R
         SgsyL(i,j) = min(VL-SQR_PHI_L,U_S-SQR_PHI_S)
         SgsyR(i,j) = max(VR+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Maskgs(i,j-1)==0.and.Maskgs(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyR(i,j)**2 &
		                  +CygsyR(i,j)*abs(DgsyR(i,j)))
		 VR = SLIDE_GAMMA*VgsyR(i,j)
         SgsyL(i,j) = VR-2.0*SQR_PHI_R
         SgsyR(i,j) = VR+SQR_PHI_R
       elseif(Maskgs(i,j-1)==1.and.Maskgs(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(SLIDE_GAMMA*(SLIDE_GAMMA-1.0)*VgsyL(i,j)**2 &
		                  +CygsyL(i,j)*abs(DgsyL(i,j)))
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
	   Fgsx(Ibeg,j) = 0.5*CxgsxR(Ibeg,j)*DgsxR(Ibeg,j)*DgsxR(Ibeg,j)
       Ggsx(Ibeg,j) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
       Egsx(Iend1,j) = Zero
       Fgsx(Iend1,j) = 0.5*CxgsxL(Iend1,j)*DgsxL(Iend1,j)*DgsxL(Iend1,j)
       Ggsx(Iend1,j) = Zero
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
       Ggsy(i,Jbeg) = 0.5*CygsyR(i,Jbeg)*DgsyR(i,Jbeg)*DgsyR(i,Jbeg)
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
	   Ggsy(i,Jend1) = 0.5*CygsyL(i,Jend1)*DgsyL(i,Jend1)*DgsyL(i,Jend1)
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
        elseif(Maskgs(i+1,j)==1) then !left boundary
          Egsx(i+1,j) = Zero
          Fgsx(i+1,j) = 0.5*Cxgs(i+1,j)*DgsxR(i+1,j)*DgsxR(i+1,j)
          Ggsx(i+1,j) = Zero
        elseif(Maskgs(i,j-1)==1) then !back boundary
          Egsy(i,j) = Zero
          Fgsy(i,j) = Zero
          Ggsy(i,j) = 0.5*Cygs(i,j-1)*DgsyL(i,j)*DgsyL(i,j)
        elseif(Maskgs(i,j+1)==1) then !front boundary
          Egsy(i,j+1) = Zero
          Fgsy(i,j+1) = Zero
          Ggsy(i,j+1) = 0.5*Cygs(i,j+1)*DgsyR(i,j+1)*DgsyR(i,j+1)
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
	                                         Taxxd,Tayxd,Taxyd,Tayyd, &
											 DelxTaxxd,DelyTayxd,DelxTaxyd,DelyTayyd
	 
	 call DelxH_gs
	 ! source terms for hydrostatic pressure
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrcpgsX(i,j) = (Grav*Dgs(i,j)*DelxH0(i,j))*Maskgs(i,j);
       SrcpgsY(i,j) = (Grav*Dgs(i,j)*DelyH0(i,j))*Maskgs(i,j);
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
     allocate(DelxTaxxd(Mloc,Nloc))
     allocate(DelyTayxd(Mloc,Nloc))
	 allocate(DelxTaxyd(Mloc,Nloc))
	 allocate(DelyTayyd(Mloc,Nloc))
     if(trim(RHEO_OPT)=='VISCOUS') then
	   call laminar(SLIDE_VISCOSITY,Rhogs,Tbxx,Tbxy,Tbxz,Tbyx,Tbyy,Tbyz,Tbzx,Tbzy,Tbzz)
	   Taxx=0.5*Tbxx;Tayx=0.5*Tbyx;Taxy=0.5*Tbxy;Tayy=0.5*Tbyy
	 elseif(trim(RHEO_OPT)=='GRANULAR') then
	   call coulomb
	 endif
	 Taxxd=Taxx*Dgs;Tayxd=Tayx*Dgs;Taxyd=Taxy*Dgs;Tayyd=Tayy*Dgs
	 DelxTaxxd=zero;DelyTayxd=zero;DelxTaxyd=zero;DelyTayyd=zero
	 call delxFun_gs_2D(Taxxd,DelxTaxxd)
	 call delyFun_gs_2D(Tayxd,DelyTayxd)
	 call delxFun_gs_2D(Taxyd,DelxTaxxd)
	 call delyFun_gs_2D(Tayyd,DelyTayxd)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SrctgsX(i,j) = 1.0/Rhogs(i,j)*(DelxTaxxd(i,j)+DelyTayxd(i,j) &
	                 -Tbxx(i,j)*DelxH0(i,j)-Tbyx(i,j)*DelyH0(i,j)-Tbzx(i,j))*Maskgs(i,j)
       SrctgsY(i,j) = 1.0/Rhogs(i,j)*(DelxTaxyd(i,j)+DelyTayyd(i,j) &
	                 -Tbxy(i,j)*DelxH0(i,j)-Tbyy(i,j)*DelyH0(i,j)-Tbzy(i,j))*Maskgs(i,j)
     enddo
     enddo
	 
     deallocate(DelxPH)
     deallocate(DelyPH)
	 deallocate(PH)
     deallocate(Taxxd)
     deallocate(Tayxd)
	 deallocate(Taxyd)
	 deallocate(Tayyd)
     deallocate(DelxTaxxd)
     deallocate(DelyTayxd)
	 deallocate(DelxTaxyd)
	 deallocate(DelyTayyd)

     end subroutine source_terms_gs
	 
     subroutine laminar(Nu,Rho,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz)
!-----------------------------------------
!    Rheology of laminar flow
!    Called by 
!       source_terms_gs, coulomb
!    Last update: 10/07/2018, Cheng Zhang 
!------------------------------------------ 
     use global, only: SP,Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Ugs,Vgs,Dgs,DelxH0,DelyH0
     implicit none
     real(SP),intent(in) :: NU
	 real(SP),intent(in),dimension(Mloc,Nloc) :: Rho
     real(SP),intent(out),dimension(Mloc,Nloc) :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
     integer :: i,j

     do j = Jbeg,Jend
     do i = Ibeg,Iend
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
	                                         Tbfyz,Tbfzx,Tbfzy,Tbfzz,Tezz, &
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
	 allocate(Tezz(Mloc,Nloc))
	 allocate(Rhof(Mloc,Nloc))
	 Tbfxx = Zero;Tbfxy = Zero;Tbfxz = Zero;Tbfyx = Zero;Tbfyy = Zero;Tbfyz = Zero
	 Tbfzx = Zero;Tbfzy = Zero;Tbfzz = Zero;Tezz = Zero
	 
	 ! fluid stress
	 Nugse=(1.-SLIDE_CONC)*Visc
	 Rhof=Rho0
     call laminar(Nugse,Rhof,Tbfxx,Tbfxy,Tbfxz,Tbfyx,Tbfyy,Tbfyz,Tbfzx,Tbfzy,Tbfzz)
	 
	 ! solid stress
     do j = 1,Nloc
     do i = 1,Mloc
	   Tezz(i,j) = (    SLIDE_CONC*GRAIN_DENSITY*Grav*Dgs(i,j))
	   Taxx(i,j) = 0.5*Tbfxx(i,j)*Maskgs(i,j)
	   Tbxx(i,j) = (Tbfxx(i,j)-Cxx(i,j)/max(Czz(i,j),Small)*Tezz(i,j)+Tezz(i,j))*Maskgs(i,j)
	   Tayy(i,j) = 0.5*Tbfyy(i,j)*Maskgs(i,j)
	   Tbyy(i,j) = (Tbfyy(i,j)-Cyy(i,j)/max(Czz(i,j),Small)*Tezz(i,j)+Tezz(i,j))*Maskgs(i,j)
	   Tbxy(i,j) = (Tbfxy(i,j)-Cxy(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbyx(i,j) = (Tbfyx(i,j)-Cyx(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Taxy(i,j) = 0.5*(Tbfxy(i,j)-Cxy(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tayx(i,j) = 0.5*(Tbfyx(i,j)-Cyx(i,j)/max(Czz(i,j),Small)*Tezz(i,j))*Maskgs(i,j)
	   Tbzz(i,j) = Tbfzz(i,j)-Tezz(i,j)+SLIDE_CONC*GRAIN_DENSITY*Grav*Dgs(i,j) !Ts=-Te+Ps
     enddo
     enddo
     Tbzx=Tbfzx
     Tbzy=Tbfzy
     call phi_2D_coll(Taxx)
	 call phi_2D_coll(Tbxx)
     call phi_2D_coll(Tayy)
	 call phi_2D_coll(Tbyy)
     call phi_2D_coll(Tbxy)
	 call phi_2D_coll(Tbyx)
     call phi_2D_coll(Taxy)
	 call phi_2D_coll(Tayx)
     call phi_2D_coll(Tbzx)
	 call phi_2D_coll(Tbzy)
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
	 deallocate(Tezz)
	 deallocate(Rhof)

     end subroutine coulomb
	 

	 subroutine eval_duv_gs(ISTEP)
!-----------------------------------------------
!    Update all variables Dgs,Ugs,Vgs for deformable slide
!    Called by
!       update_bathymetry
!    Last update: 10/07/2018, Cheng Zhang
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP),dimension(:,:),allocatable :: R1,R2,R3,Thetax,Thetay,DelxTheta,DelyTheta
     integer :: i,j
     REAL(SP) :: FROUDECAP_gs,FroudeU,DUU,Dangle,TauC

     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc))
     allocate(R3(Mloc,Nloc))
     allocate(Thetax(Mloc,Nloc))
	 allocate(Thetay(Mloc,Nloc))
     allocate(DelxTheta(Mloc,Nloc))
     allocate(DelyTheta(Mloc,Nloc))
	 
     ! solve total slide depth Dgs
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==0) cycle
       R1(i,j) = -1.0/dx*(Egsx(i+1,j)-Egsx(i,j))-1.0/dy*(Egsy(i,j+1)-Egsy(i,j)) 
       Dgs(i,j) = Dgs(i,j)+dt*R1(i,j) 
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
	   DUgs(i,j) = DUgs(i,j)+dt*R2(i,j)
     enddo
     enddo
	 
	 ! solve DVgs
     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Maskgs(i,j)==0) cycle  
       R3(i,j) = -1.0/dx*(Ggsx(i+1,j)-Ggsx(i,j))-1.0/dy*(Ggsy(i,j+1)-Ggsy(i,j)) &
                 +SrcpgsY(i,j)+SrctgsY(i,j)
	   DVgs(i,j) = DVgs(i,j)+dt*R3(i,j)
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
		    TauC=tan(PhiBed)/max(Czz0(i,j),Small) &
			     *(   SLIDE_CONC*GRAIN_DENSITY*Grav*((DgsxL(i,j)+DgsxR(i+1,j))/2. &
                                                    +(DgsyL(i,j)+DgsyR(i,j+1))/2.)/2.)
            if(sqrt(DUgs(i,j)**2+DVgs(i,j)**2)>dt/Rhogs(i,j)*TauC) then
              Tbzx(i,j) = -Czx(i,j)/max(Czz(i,j),Small) &
			             *(   SLIDE_CONC*GRAIN_DENSITY*Grav*(DgsxL(i,j)+DgsxR(i+1,j))/2.) &
						 *(1.-(DUgs(i,j)/Dgs(i,j))**2*DelxTheta(i,j)/Grav)
              Tbzy(i,j) = -Czy(i,j)/max(Czz(i,j),Small) &
			             *(   SLIDE_CONC*GRAIN_DENSITY*Grav*(DgsyL(i,j)+DgsyR(i,j+1))/2.) &
						 *(1.-(DVgs(i,j)/Dgs(i,j))**2*DelyTheta(i,j)/Grav)
            else
              Dgs(i,j) = Dgs(i,j)-dt*R1(i,j)
              Tbzx(i,j) = Rhogs(i,j)/dt*DUgs(i,j)
              Tbzy(i,j) = Rhogs(i,j)/dt*DVgs(i,j)
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
     endif
	 
	 ! run non-hydrostatic simulation
	 if(NON_HYDRO_SLD) then
	   ! obtain hydrostatic velocity
	   call get_UV_gs
	   
       ! solve dynamic pressure 
       call poisson_solver_gs

       ! correct velocity field  
       call projection_corrector_gs
	 endif
	 
     ! SSP Runge-Kutta time stepping
     do j = Jbeg,Jend
     do i = Ibeg,Iend
	   Dgs(i,j) = ALPHA(ISTEP)*Dgs0(i,j)+BETA(ISTEP)*Dgs(i,j)
       DUgs(i,j) = ALPHA(ISTEP)*DUgs0(i,j)+BETA(ISTEP)*DUgs(i,j)
       DVgs(i,j) = ALPHA(ISTEP)*DVgs0(i,j)+BETA(ISTEP)*DVgs(i,j)
	   if(NON_HYDRO_SLD) then
         Wgs(i,j) = ALPHA(ISTEP)*Wgs0(i,j)+BETA(ISTEP)*Wgs(i,j);
         DWgs(i,j) = Wgs(i,j)*Dgs(i,j);
	   endif
     enddo
     enddo
	 
     ! update Dgs and Hgs          
     Dgs = max(Dgs,SLIDE_MINTHICK)
     call wl_bc_gs
     Hgs = Hc0-Dgs
	 
     ! limit the maximum Froude number
     FROUDECAP_gs=50.0
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Maskgs(I,J)>0)THEN
         FroudeU=SQRT(grav*Dgs(I,J))*FROUDECAP_gs*Dgs(I,J)
         DUU=SQRT(DUgs(I,J)**2+DVgs(I,J)**2)
         IF(DUU>FroudeU)THEN
          Dangle=atan2(DVgs(I,J),DUgs(I,J))
          DUgs(I,J)=FroudeU*COS(Dangle)
          DVgs(I,J)=FroudeU*SIN(Dangle)
         ENDIF
       ENDIF
     ENDDO
     ENDDO
	 
	 if(.not.NON_HYDRO_SLD) then
	   ! solve Wgs
	   call Kbc_Surface_gs
       call Kbc_Bottom_gs
	   do j = Jbeg,Jend
       do i = Ibeg,Iend
         Wgs(i,j) = 0.5*(Wtgs(i,j)+Wbgs(i,j))
         DWgs(i,j) = Wgs(i,j)*Dgs(i,j)
	   enddo
	   enddo
	 endif

     ! update maskgs for deformable slide
	 ! call update_maskgs

     ! boundary conditions and final velocity
     call get_UVW_gs

	 deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(Thetax)
	 deallocate(Thetay)
     deallocate(DelxTheta)
     deallocate(DelyTheta)
	 
	 end subroutine eval_duv_gs
	 
     subroutine wl_bc_gs
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duv_gs
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
!       eval_duv_gs                                                                                                                  
!    Last update: 10/07/2018, Cheng Zhang                                                                                    
!------------------------------------------------                                                                                
     use global
     implicit none
     integer :: i,j

	 Wtgs=Zero
	 call delxFun_gs_2D(Hgs,DelxHgs)
	 call delyFun_gs_2D(Hgs,DelyHgs)
	 do j = Jbeg,Jend
     do i = Ibeg,Iend
       Wtgs(i,j) = -(Hgs(i,j)-(Hc0(i,j)-Dgs0(i,j)))/dt &
	               -DUgs(i,j)/Dgs(i,j)*DelxHgs(i,j) &
				   -DVgs(i,j)/Dgs(i,j)*DelyHgs(i,j)                                                     
     enddo
     enddo

     end subroutine KBC_surface_gs

     subroutine KBC_bottom_gs
!-------------------------------------------------  
!    Applied kinematic boundary conditions at bottom                                                                        
!    of deformable slide                                                                                   
!    Called by
!       eval_duv_gs                                                                                                                  
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
	 
     subroutine get_UV_gs
!------------------------------------------------
!    Obtain Ugs,Vgs for deformable slide
!    Called by
!       eval_duv_gs
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
	   else
         Ugs(i,j) = Zero
         Vgs(i,j) = Zero
         DUgs(i,j) = Zero
         DVgs(i,j) = Zero
	   endif
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc_gs
# if defined (PARALLEL)
     call phi_2D_exch(Ugs)
     call phi_2D_exch(Vgs)
     call phi_2D_exch(DUgs)
     call phi_2D_exch(DVgs)
# endif

     end subroutine get_UV_gs
	 
     subroutine get_UVW_gs
!------------------------------------------------
!    Obtain Ugs,Vgs for deformable slide
!    Called by
!       eval_duv_gs
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
!       get_UV_gs get_UVW_gs
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
       Ugs(Ibeg-i,j) = Ugs(Ibeg+i-1,j)
       Vgs(Ibeg-i,j) = Vgs(Ibeg+i-1,j)
	   Wgs(Ibeg-i,j) = Wgs(Ibeg+i-1,j)
       DUgs(Ibeg-i,j) = DUgs(Ibeg+i-1,j)
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
       Ugs(Iend+i,j) = Ugs(Iend-i+1,j)
       Vgs(Iend+i,j) = Vgs(Iend-i+1,j)
	   Wgs(Iend+i,j) = Wgs(Iend-i+1,j)
       DUgs(Iend+i,j) = DUgs(Iend-i+1,j)
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
     do i = Ibeg,Iend
     do j = 1,Nghost
       Ugs(i,Jbeg-j) = Ugs(i,Jbeg+j-1)
       Vgs(i,Jbeg-j) = Vgs(i,Jbeg+j-1)
	   Wgs(i,Jbeg-j) = Wgs(i,Jbeg+j-1)
       DUgs(i,Jbeg-j) = DUgs(i,Jbeg+j-1)
       DVgs(i,Jbeg-j) = DVgs(i,Jbeg+j-1)
	   DWgs(i,Jbeg-j) = DWgs(i,Jbeg+j-1)
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
       Ugs(i,Jend+j) = Ugs(i,Jend-j+1)
       Vgs(i,Jend+j) = Vgs(i,Jend-j+1)
	   Wgs(i,Jend+j) = Wgs(i,Jend-j+1)
       DUgs(i,Jend+j) = DUgs(i,Jend-j+1)
       DVgs(i,Jend+j) = DVgs(i,Jend-j+1)
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
     do j = 1,Nloc
     do i = 1,Mloc
       Qtgs(i,j) = P(i,j,kbeg)! dynamic pressure over the slide
     enddo
     enddo
   
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
  
     ! boundary condition for dry cells
     ! set zero for dry set is inaccurate
     ! dry cells
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Maskgs(i,j)==0) then
         Qbgs(i,j) = Zero
         
         ! south boundary 
         if(Maskgs(i,j+1)==1)then
           do imask=1,Nghost
             Qbgs(i,j-imask+1)=Qbgs(i,j+imask)
           enddo
         ! north boundary
         elseif(Maskgs(i,j-1)==1)then
           do imask=1,Nghost
             Qbgs(i,j+imask-1)=Qbgs(i,j-imask)
           enddo
         ! west boundary
         elseif(Maskgs(i+1,j)==1)then
           do imask=1,Nghost
             Qbgs(i-imask+1,j)=Qbgs(i+imask,j)
           enddo
         ! east boundary
         elseif(Maskgs(i-1,j)==1)then
           do imask=1,Nghost
             Qbgs(i+imask-1,j)=Qbgs(i-imask,j)
           enddo
         endif
       endif 
     enddo
     enddo

     ! collect into ghost cells
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do i = 1,Nghost
       Qbgs(Ibeg-i,j) = Qbgs(Ibeg+i-1,j)
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
       Qbgs(Iend+i,j) = Qbgs(Iend-i+1,j)
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
       Qbgs(i,Jbeg-j) = Qbgs(i,Jbeg+j-1)
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
       Qbgs(i,Jend+j) = Qbgs(i,Jend-j+1)
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(Qbgs)
# endif

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
     real(SP), dimension(:,:), allocatable :: Delx2Dgs,Dely2Dgs,Delx2Hgs,Dely2Hgs, &
	                                          RhoUgs,RhoVgs,DelxRhoUgs,DelyRhoVgs, &
										      DelxQtgs,DelyQtgs,Delx2Qtgs,Dely2Qtgs
     integer,  dimension(:,:), allocatable :: indx
	 REAL(SP) :: R1,R2,R3,R4,R5,R6,R7,R8,R9,L1,L2,L3,L4,L5,L6,L7

     allocate(Delx2Dgs(Mloc,Nloc))
     allocate(Dely2Dgs(Mloc,Nloc))
     allocate(Delx2Hgs(Mloc,Nloc))
     allocate(Dely2Hgs(Mloc,Nloc))
     allocate(RhoUgs(Mloc,Nloc))
	 allocate(RhoVgs(Mloc,Nloc))
     allocate(DelxRhoUgs(Mloc,Nloc))
	 allocate(DelyRhoVgs(Mloc,Nloc))
     allocate(DelxQtgs(Mloc,Nloc))
     allocate(DelyQtgs(Mloc,Nloc))
     allocate(Delx2Qtgs(Mloc,Nloc))
     allocate(Dely2Qtgs(Mloc,Nloc))
	 allocate(indx(Mloc,Nloc))
	 
	 Delx2Dgs=Zero;Dely2Dgs=Zero;Delx2Hgs=Zero;Dely2Hgs=Zero
	 RhoUgs=Zero;RhoVgs=Zero;DelxRhoUgs=Zero;DelyRhoVgs=Zero
	 DelxQtgs=Zero;DelyQtgs=Zero;Delx2Qtgs=Zero;Dely2Qtgs=Zero
     call delxFun_nl_gs_2D(Dgs,DelxDgs)
	 call delyFun_nl_gs_2D(Dgs,DelyDgs)
	 call delx2Fun_gs_2D(Dgs,Delx2Dgs)
	 call dely2Fun_gs_2D(Dgs,Dely2Dgs)
     call delxFun_nl_gs_2D(Hgs,DelxHgs)
	 call delyFun_nl_gs_2D(Hgs,DelyHgs)
	 call delx2Fun_gs_2D(Hgs,Delx2Hgs)
	 call dely2Fun_gs_2D(Hgs,Dely2Hgs)
	 call DelxH_gs
	 RhoUgs = Rhogs*Ugs
	 RhoVgs = Rhogs*Vgs
	 call delxFun_nl_gs_2D(RhoUgs,DelxRhoUgs)
	 call delyFun_nl_gs_2D(RhoVgs,DelyRhoVgs)
     call delxFun_nl_gs_2D(Qtgs,DelxQtgs)
	 call delyFun_nl_gs_2D(Qtgs,DelyQtgs)
	 call delx2Fun_gs_2D(Qtgs,Delx2Qtgs)
	 call dely2Fun_gs_2D(Qtgs,Dely2Qtgs)
	 call KBC_bottom_gs
	 
     ! generate coefficient matrix
     neq = 0; indx=0
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j) = neq
     enddo
     enddo 

     ! generate source term 
     Rhsgs = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       R1=2.*dt/3.
       R2=dt/Dgs(i,j)*(2./3.*DelxDgs(i,j)+DelxHgs(i,j))
       R3=dt/Dgs(i,j)*(2./3.*Delx2Dgs(i,j)+Delx2Hgs(i,j) &
	                  -1./Dgs(i,j)*DelxDgs(i,j)*(2./3.*DelxDgs(i,j)+DelxHgs(i,j)))
	   R4=2.*dt/3.
       R5=dt/Dgs(i,j)*(2./3.*DelyDgs(i,j)+DelyHgs(i,j))
       R6=dt/Dgs(i,j)*(2./3.*Dely2Dgs(i,j)+Dely2Hgs(i,j) &
	                  -1./Dgs(i,j)*DelyDgs(i,j)*(2./3.*DelyDgs(i,j)+DelyHgs(i,j)))
	   R7=2.*dt/Dgs(i,j)**2
       R8=2.*dt/Dgs(i,j)**2*(   Tbzz(i,j))
       R9=-(DelxRhoUgs(i,j)+DelyRhoVgs(i,j) &
	       +2.*Rhogs(i,j)*Wgs(i,j)/Dgs(i,j)-2.*Rhogs(i,j)*Wbgs(i,j)/Dgs(i,j))
       Rhsgs(indx(i,j)) = R1*Delx2Qtgs(i,j)+R4*Dely2Qtgs(i,j) &
	                     +R2*DelxQtgs(i,j)+R5*DelyQtgs(i,j)+(R3+R6+R7)*Qtgs(i,j)+R8+R9			 
     enddo
     enddo

	 ! generate coefficient matrix
     Coefgs = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       L1=-2.*dt/3.
       L2=-dt/Dgs(i,j)*(2./3.*DelxDgs(i,j)-DelxH0(i,j))
       L3=-dt/Dgs(i,j)*(2./3.*Delx2Dgs(i,j)-Delx2H0(i,j) &
	                   -1./Dgs(i,j)*DelxDgs(i,j)*(2./3.*DelxDgs(i,j)-DelxH0(i,j)))
	   L4=-2.*dt/3.
       L5=-dt/Dgs(i,j)*(2./3.*DelyDgs(i,j)-DelyH0(i,j))
       L6=-dt/Dgs(i,j)*(2./3.*Dely2Dgs(i,j)-Dely2H0(i,j) &
	                   -1./Dgs(i,j)*DelyDgs(i,j)*(2./3.*DelyDgs(i,j)-DelyH0(i,j)));
	   L7=2.*dt/Dgs(i,j)**2;
       Coefgs(indx(i,j),1) = (L3+L6+L7)-2.*L1/dx**2-2.*L4/dy**2
       Coefgs(indx(i,j),2) = L1/dx**2-L2/2./dx
       Coefgs(indx(i,j),3) = L1/dx**2+L2/2./dx
       Coefgs(indx(i,j),4) = L4/dy**2-L5/2./dy
       Coefgs(indx(i,j),5) = L4/dy**2+L5/2./dy
     enddo
     enddo
	 
     ! boundary conditions at masks face
     do i = Ibeg+1,Iend-1
     do j = Jbeg+1,Jend-1
       if(Maskgs(i,j)==0) then
         ! left 
         if(Maskgs(i+1,j)==1) then
           ic = indx(i+1,j)
           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,2)
           Coefgs(ic,2) = Zero
         ! right 
         elseif(Maskgs(i-1,j)==1) then
           ic = indx(i-1,j)
           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,3)
           Coefgs(ic,3) = Zero
         ! south
         elseif(Maskgs(i,j+1)==1) then
           ic = indx(i,j+1)
           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,4)
           Coefgs(ic,4) = Zero
         ! north
         elseif(Maskgs(i,j-1)==1) then
           ic = indx(i,j-1)
           Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,5)
           Coefgs(ic,5) = Zero
         endif ! end Maskgs+1=1 
       endif ! end Maskgs=0
     enddo
     enddo

     ! boundary conditions
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     i = Ibeg
     do j = Jbeg,Jend
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
     i = Iend
     do j = Jbeg,Jend
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
     do i = Ibeg,Iend
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
     j = Jend
     do i = Ibeg,Iend
       ic = indx(i,j)
       Coefgs(ic,1) = Coefgs(ic,1)+Coefgs(ic,5)
       Coefgs(ic,5) = Zero
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! take (i=2,j=2) to obtain the diagonal information
     JCoefgs(1) = indx(Ibeg+1,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i,j)
     JCoefgs(2) = indx(Ibeg  ,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i-1,j)
     JCoefgs(3) = indx(Ibeg+2,Jbeg+1)-indx(Ibeg+1,Jbeg+1)  ! (i+1,j) 
	 JCoefgs(4) = indx(Ibeg+1,Jbeg  )-indx(Ibeg+1,Jbeg+1)  ! (i,j-1)
     JCoefgs(5) = indx(Ibeg+1,Jbeg+2)-indx(Ibeg+1,Jbeg+1)  ! (i,j+1)

     deallocate(Delx2Dgs)
     deallocate(Dely2Dgs)
     deallocate(Delx2Hgs)
     deallocate(Dely2Hgs)
     deallocate(RhoUgs)
	 deallocate(RhoVgs)
     deallocate(DelxRhoUgs)
	 deallocate(DelyRhoVgs)
     deallocate(DelxQtgs)
     deallocate(DelyQtgs)
     deallocate(Delx2Qtgs)
     deallocate(Dely2Qtgs)
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
     integer :: i_glob(Mloc),j_glob(Nloc)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,-1,0,1,0,0,-1,0,1/

     ! set up a two dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)

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

     allocate(indx(Mloc,Nloc))
 
     neq = 0
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j) = neq
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*nentries
     allocate(values(nvalues))

     ivalues = 0
     do j = Jbeg,Jend
     do i = Ibeg,Iend
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

     ! set the vector Coefgsficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhsgs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqnsgs))
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j)) = Qbgs(i,j)
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

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Qbgs(i,j) = Phi(indx(i,j))
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
     real(SP), dimension(:,:),allocatable :: DelxQtgs,DelyQtgs,DelxQbgs,DelyQbgs
 
     allocate(DelxQtgs(Mloc,Nloc))
     allocate(DelyQtgs(Mloc,Nloc))
     allocate(DelxQbgs(Mloc,Nloc))
     allocate(DelyQbgs(Mloc,Nloc))

     DelxQtgs = Zero;DelyQtgs = Zero
     DelxQbgs = Zero;DelyQbgs = Zero
	 call delxFun_nl_gs_2D(Dgs,DelxDgs)
	 call delyFun_nl_gs_2D(Dgs,DelyDgs)
	 call delxFun_nl_gs_2D(Qtgs,DelxQtgs)
	 call delyFun_nl_gs_2D(Qtgs,DelyQtgs)
	 call delxFun_nl_gs_2D(Qbgs,DelxQbgs)
	 call delyFun_nl_gs_2D(Qbgs,DelyQbgs)

     do j = Jbeg,Jend
	 do i = Ibeg,Iend
       if(Maskgs(i,j)==0) cycle
         Ugs(i,j)=Ugs(i,j)-dt/Rhogs(i,j)/Dgs(i,j) &
		 *(2./3.*(Dgs(i,j)*DelxQtgs(i,j)+Qtgs(i,j)*DelxDgs(i,j))+Qtgs(i,j)*DelxHgs(i,j) &
		  +2./3.*(Dgs(i,j)*DelxQbgs(i,j)+Qbgs(i,j)*DelxDgs(i,j))-Qbgs(i,j)*DelxH0(i,j))
         Vgs(i,j)=Vgs(i,j)-dt/Rhogs(i,j)/Dgs(i,j) &
		 *(2./3.*(Dgs(i,j)*DelyQtgs(i,j)+Qtgs(i,j)*DelyDgs(i,j))+Qtgs(i,j)*DelyHgs(i,j) &
		  +2./3.*(Dgs(i,j)*DelyQbgs(i,j)+Qbgs(i,j)*DelyDgs(i,j))-Qbgs(i,j)*DelyH0(i,j))
         Wtgs(i,j)=2.*Wgs(i,j)-Wbgs(i,j)-2.*dt/Rhogs(i,j)/Dgs(i,j)*(Tbzz(i,j)+Qtgs(i,j)-Qbgs(i,j));
         Wgs(i,j)=0.5*(Wtgs(i,j)+Wbgs(i,j));
         DUgs(i,j)=Ugs(i,j)*Dgs(i,j);
         DVgs(i,j)=Vgs(i,j)*Dgs(i,j);
     enddo
     enddo
 
     deallocate(DelxQtgs)
     deallocate(DelyQtgs)
     deallocate(DelxQbgs)
	 deallocate(DelyQbgs)

     return
     end subroutine projection_corrector_gs
	 
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
       if(Maskgs(i,j)==0) then
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1+TMP2)/2.0
         endif
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
       if(Maskgs(i,j)==0) then 
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = (TMP1+TMP2)/2.0
         endif
       endif
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
# endif