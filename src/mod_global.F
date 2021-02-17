     Module GLOBAL
# if defined (PARALLEL)
     use mpi
# endif
     implicit none

     ! define precision
# if defined (DOUBLE_PRECISION)
     integer, parameter :: SP=8
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_DOUBLE_PRECISION
# endif
# else
     integer, parameter :: SP=SELECTED_REAL_KIND(6,30)
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_REAL
# endif
# endif

     ! define parameters
     real(SP), parameter :: pi=3.141592653
     real(SP), parameter :: Small=1.0e-16
     real(SP), parameter :: Large=10000000.0
     real(SP), parameter :: Grav=9.81
     real(SP), parameter :: Zero=0.0
     real(SP), parameter :: One=1.0
     real(SP), parameter :: Rho0=1000.0
     real(SP), parameter :: RhoA=1.20
     real(SP), parameter :: Srho = 2650.
     real(SP), parameter :: Kappa=0.41
     integer,  parameter :: MaxNumFreq=500
     integer,  parameter :: MaxNumDir=100

     ! ghost cells (>=1)
     integer, parameter :: Nghost=2

     ! define characters
     character(len=80) :: TITLE
     character(len=80) :: RESULT_FOLDER
     character(len=80) :: HIGH_ORDER
     character(len=80) :: TIME_ORDER
     character(len=80) :: WaveMaker
     character(len=80) :: DEPTH_TYPE
     character(len=80) :: dt_constraint
     character(len=80) :: CONVECTION
     character(len=80) :: DEPTH_FILE

# if defined (PARALLEL)
     integer :: myid,ier
     integer :: comm2d
     integer :: n_west,n_east,n_suth,n_nrth
     integer :: npx,npy
     integer :: ndims=2
     integer :: NumP
     integer, dimension(2) :: dims,coords
     logical, dimension(2) :: periods
     logical :: reorder=.true.
# endif


# if defined (BUBBLE)
     integer, parameter :: Mg = 20
     real(SP), parameter :: Surface_Tension = 71.97e-3
     real(SP), parameter :: Con_b = 0.5
     real(SP), parameter :: Eps_Cric = 0.01
     real(SP), parameter :: TIM_B = 10.0
     real(SP), dimension(:), allocatable :: Rbg,DRbg,Wbg,Entrain
     real(SP), dimension(:,:,:,:), allocatable :: Nbg,DNbg,DNbg0
     real(SP), dimension(:,:,:), allocatable :: Vbg
# endif

# if defined (CORALREEF)
     real(SP) :: Creef
# endif

# if defined (SEDIMENT)
     real(SP), parameter :: Spor = 0.4
     logical :: COUPLE_FS,BED_LOAD
     character(len=80) :: Sed_Type
     real(SP), dimension(:,:),   allocatable :: SSour,Taub,Bed,Bed0,Qbedx,Qbedy, &
                                                QbxAvg,QbyAvg
     real(SP), dimension(:,:),   allocatable :: Sed_X0,Sed_Xn,Smean,TotDep
     real(SP), dimension(:,:,:), allocatable :: Conc,DConc,DConc0,Euler_Smean
     real(SP) :: Update_Bed_T,MorDt,SD50,Af,Shields_c,Tau_ce,Tau_cd,Erate,  &
                 Tim_Sedi,Wset,Sedi_Ws
     integer :: ntyws,nbed
# endif

# if defined (SALINITY)
     real(SP), dimension(:,:,:), allocatable :: Sali,DSali,DSali0
     real(SP), dimension(:,:), allocatable :: Sin_X0,Sin_Xn 
# endif

# if defined (TEMPERATURE)
     real(SP), dimension(:,:,:), allocatable :: Temp,DTemp,DTemp0
     real(SP), dimension(:,:), allocatable :: Tin_X0,Tin_Xn
# endif

# if defined (VEGETATION)
     character(len=80) :: Veg_Type
     real(SP) :: VegH,VegDens,StemD,VegDrag,VegVol,VegVM,Veg_X0,Veg_Xn,Veg_Y0,Veg_Yn,EI,Cfe,Cfk
     real(SP), dimension(:,:,:), allocatable :: Prod_v,Tke_w,Eps_w,Prod_w,Wke_p,DWke,DWke0
     real(SP), dimension(:,:), allocatable :: FVegH,Foliage
# endif

# if defined (LANDSLIDE)
     character(len=80) :: SlideType
     real(SP) :: SlideT,SlideL,SlideW,SlideAngle,SlopeAngle,SlideD,SlideX0,SlideY0
     real(SP) :: SlideUt,SlideA0,SlideDens
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
! more comprehensive slide setting including multi-slides, rectangular shape
! rotational etc
     real(SP) :: kb_slide,kw_slide, &
                 acceleration_lab

! for multiple slides
     CHARACTER(LEN=80) :: MultiSlidesFile
     CHARACTER(LEN=80) :: VoidString
     INTEGER :: NumSlides
     REAL(SP),DIMENSION(:),ALLOCATABLE :: T_slide,b_slide,w_slide,e_slide,alpha_slide, &
                 x0_slide,y0_slide,term_v,slope_slide
     REAL(SP),DIMENSION(:),ALLOCATABLE :: uslide,ac_slide
     REAL(SP),DIMENSION(:),ALLOCATABLE  :: subaer_center2water, &
                                   subaer_center2water0,Coulomb_phi
     REAL(SP),DIMENSION(:),ALLOCATABLE  :: Smax_slide,SlideStartTime

! add base bathymetry for landslide fyshi
     real(SP), dimension(:,:), allocatable :: Dep0,DepC0    
     LOGICAL :: SLIDE_SHAPE_BOX 
     LOGICAL :: SMFCut,Slide,SlideStop
     real(SP) :: sf_slide,tf_slide,ut_slide,ACC_lab
     real(SP) :: XX0_Slide,YY0_Slide ! used for grilli ellipse case

     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: EleMax,EleMin
     REAL(SP) :: SlideErrorRange = 0.01  ! avoid sharp edge in box cases 

# endif 
! end comprehensive

# if defined (BALANCE2D)
     real(SP), dimension(:,:), allocatable :: DUDT2D,DVDT2D,DEDX2D,DEDY2D,DPDX2D,DPDY2D,DIFFX2D,DIFFY2D,  &
                TAUBX2D,TAUBY2D,FVEGX2D,FVEGY2D,ADVX2D,ADVY2D
# endif

     real(SP), dimension(:,:,:), allocatable :: Porosity,Ap_Por,Bp_Por,Cp_Por
# if defined (POROUSMEDIA)
     real(SP) :: alpha_por,beta_por,D50_Por,Por_n
     real(SP) :: Por_X0,Por_Xn,Por_Y0,Por_Yn,Por_Z0,Por_Zn
# endif

# if defined (OBSTACLE)
     real(SP) :: obs_u,obs_v,obs_w
     real(SP), dimension(:,:,:), allocatable :: ObsForceX,ObsForceY,ObsForceZ
     real(SP), dimension(:,:,:), allocatable :: set_dist_x,set_dist_y,set_dist_z
     integer, dimension(:,:,:), allocatable :: set_flag,set_flag_old
# if defined (DEFORMABLESLIDE)
     real(SP), dimension(:,:), allocatable :: ObsForceX_gs,ObsForceY_gs,ObsForceZ_gs
     real(SP), dimension(:,:), allocatable :: set_dist_x_gs,set_dist_y_gs,set_dist_z_gs
     integer, dimension(:,:), allocatable :: set_flag_gs,set_flag_gs_old
# endif
# endif

# if defined (TWOLAYERSLIDE)
     real(SP) :: SlideT,SlideW,SlideL,SlideAngle,SlopeAngle,SlideX0,SlideY0,SlideDens,SlideVisc
     real(SP) :: SlideLambda,SlideIniU,PhiInt,PhiBed,Cf_ul,Cs_ul
     real(SP), dimension(:,:), allocatable :: Ha0,Ha,Hs,Ua0,Va0,Ua,Va,HUa0,HVa0,HUa,HVa, &
                                              HUaxL,HUaxR,HVaxL,HVaxR,HUayL,HUayR,HVayL,HVayR, &
                                              Eax,Eay,EaxL,EaxR,EayL,EayR,Fax,Fay,FaxL,FaxR,FayL, &
                                              FayR,Gax,Gay,GaxL,GaxR,GayL,GayR,SorceaX,SorceaY, &
                                              DelxHa,DelxUa,DelxVa,DelxHUa,DelxHVa,DelyHa,DelyUa,DelyVa, &
                                              DelyHUa,DelyHVa,HaxL,HaxR,UaxL,UaxR,VaxL,VaxR,HayL,HayR, &
                                              UayL,UayR,VayL,VayR,Kap,KapxL,KapxR,KapyL,KapyR,SaxL,SaxR, &
                                              SaxS,SayL,SayR,SayS,DelxKap,DelyKap,Slambda, &
                                              ForceX_UL,ForceY_UL,Grav_LL,DelxGravLz,DelyGravLz,GravxL,GravxR, &
                                              GravyL,GravyR,wbar,Grav_Lx,Grav_Ly,Grav_Lz
     integer, dimension(:,:), allocatable :: Maska
# endif

! define output logical parameters
     logical :: ANA_BATHY,NON_HYDRO,VISCOUS_FLOW,SPONGE_ON,OUT_H,OUT_E,OUT_U,OUT_V,OUT_W,OUT_P, &
                OUT_K,OUT_D,OUT_S,OUT_C,OUT_B,OUT_A,OUT_F,OUT_T,OUT_G,OUT_I,OUT_Z,OUT_M,PERIODIC_X,PERIODIC_Y, &
                WAVE_AVERAGE_ON,ADV_HLLC,BAROTROPIC,RIGID_LID,BED_CHANGE,EXTERNAL_FORCING,STATIONARY, &
                INITIAL_EUVW,INITIAL_SALI,RHEOLOGY_ON,RNG

! variables
     integer :: It_Order,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend,Kend1,PX,PY,IVturb,IHturb,Iws,  &
                Mglob,Nglob,Kglob,Mloc,Nloc,Kloc,Mloc1,Nloc1,Kloc1,Icount,RUN_STEP,Ivgrd,SIM_STEPS,Ibot, &
                NumFreq,NumDir,NSTAT,WaveheightID,ProdType
     integer :: Bc_X0,Bc_Xn,Bc_Y0,Bc_Yn,Bc_Z0,Bc_Zn
     real(SP) :: dt,dt_old,dt_min,dt_max,dt_ini,dx,dy,Theta,CFL,VISCOUS_NUMBER,MinDep,TIME,TOTAL_TIME,Plot_Intv,  &
                 Screen_Intv,Screen_Count,Plot_Count,Visc,Cvs,Chs,Zob,Tke_min,  &
                 Eps_min,Cmut_min,Cd0,Plot_Start,Plot_Intv_Stat, &
                 Plot_Count_Stat,xstat(200),ystat(200),zstat(200),Wave_Ave_Start,  & 
                 Wave_Ave_End,Schmidt,TRamp,Grd_R,Yield_Stress,Plastic_Visc, &
                 Mud_Visc,Water_Depth,WindU,WindV,slat,fcor,Pgrad0,DepConst
     real(SP) :: Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave,Freq(MaxNumFreq),  &
                 Dire(MaxNumDir),Wave_Spc2d(MaxNumDir,MaxNumFreq), &
                 Random_Phs(MaxNumDir,MaxNumFreq),Hm0,Tp,Freq_Min,Freq_Max,  & 
                 Jon_Spc(MaxNumFreq),RanPhs(MaxNumFreq),Cur_Wave
     real(SP) :: Sponge_West_Width,Sponge_East_Width,Sponge_South_Width,Sponge_North_Width, &
                 Xsource_West,Xsource_East,Ysource_Suth,Ysource_Nrth
     real(SP), dimension(3) :: ALPHA,BETA

     ! M.Derakhti added this for FOCUSED wavemaker 
     integer :: nwave,Component_Amp_Type
     real(SP) :: k_center,f_center,x_breaking,t_breaking,Slope_group,normalized_delta_f,  &
                 depth_comp,eta_mean,sd_return     
     real(SP), dimension(12000) :: Amp_i,f_i,k_i,Phase_i

     ! real arrays
     real(SP), dimension(:), allocatable :: x,xc,y,yc,sig,dsig,sigc,Ein_X0,Din_X0,Ein_Xn,Din_Xn
     real(SP), dimension(:,:), allocatable :: Ho,H,Hc,HCG,Hc0,Hfx,Hfy,Hfx0,Hfy0,DeltH,DeltHo,Delt2H,  &
                                              DelxH,DelyH,D,D0,Eta,Eta0,Eta00, &
                                              SourceX,SourceY,SourceC,DxL,DxR,DyL,DyR,EtaxL,EtaxR,EtayL,EtayR, &
                                              DelxEta,DelyEta,DelxD,DelyD,Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn, &
                                              Win_Xn,Bc_Prs,Sponge,Setup,WaveHeight,Umean,Vmean,Emax,Emin,WdU,WdV,Wsx,Wsy, &
											  HeightMax,Uin_Xni,Uin_Xni0
     real(SP), dimension(:,:,:), allocatable :: U,V,W,U0,V0,W0,U00,V00,W00,Omega,P,DU,DV,DW,DU0,DV0,DW0, &
                                                UxL,UxR,VxL,VxR,WxL,WxR,DUxL,DUxR,DVxL,DVxR,DWxL, &
                                                DWxR,UyL,UyR,VyL,VyR,WyL,WyR,DUyL,DUyR,DVyL,DVyR,DWyL,DWyR, &
                                                UzL,UzR,VzL,VzR,WzL,WzR,OzL,OzR,SxL,SxR,SxS,SyL,SyR,SyS,ExL,ExR,FxL, &
                                                FxR,GxL,GxR,HxL,HxR,EyL,EyR,FyL,FyR,GyL,GyR,HyL,HyR,Ex,Ey,Fx, &
                                                Fy,Fz,Gx,Gy,Gz,Hx,Hy,Hz,DelxU,DelyU,DelzU,DelxV,DelyV,DelzV, &
                                                DelxW,DelyW,DelzW,DelxDU,DelyDU,DelxDV,DelyDV,DelxDW,DelyDW, &
                                                DelzO,Uf,Vf,Wf,Cmu,CmuHt,CmuVt,CmuR,Diffxx,Diffxy,Diffxz,Diffyx,  &
                                                Diffyy,Diffyz,Diffzx,Diffzy,Diffzz,DelxSc,DelySc,Rho,Rmean,Tke,Eps,Skl, &
                                                DTke,DEps,DTke0,DEps0,Prod_s,Prod_b,Richf,Lag_Umean,Lag_Vmean,Lag_Wmean, &
                                                Euler_Umean,Euler_Vmean,Euler_Wmean,DRhoX,DRhoY,ExtForceX,ExtForceY, &
                                                UpWp,PresForceX,PresForceY,PresForceZ,Pdiff,DelxSl,DelySl
!     real(SP), dimension(:,:,:,:), allocatable :: UGrad,VGrad,WGrad  ! deleted by Cheng, they are arrays never used.

     ! integer arrays
     integer, dimension(:,:), allocatable :: Mask,Mask_Struct,Mask9,Brks,Num_Zero_Up
     integer, dimension(:,:,:), allocatable :: IsMove
     
     ! poisson solvers
     integer  :: itmax,isolver,neqns
     real(SP) :: tol
     real(SP), dimension(:),   allocatable :: Rhs
     integer,  dimension(:),   allocatable :: JCoef
     real(SP), dimension(:,:), allocatable :: Coef
	 
     ! added by Cheng for nesting
# if defined(COUPLING)
! coupling
       CHARACTER(LEN=80) COUPLING_FILE
       INTEGER :: icount_coupling = 1
       INTEGER :: N_COUPLING_EAST,&
                  N_COUPLING_WEST,&
                  N_COUPLING_SOUTH,&
                  N_COUPLING_NORTH,&
                  N_COUPLING_DATA,&
                  J_START_EAST,&
                  J_START_WEST,&
                  I_START_SOUTH,&
                  I_START_NORTH
       REAL(SP),DIMENSION(:,:),ALLOCATABLE::Z_COUPLING_EAST,&
                  Z_COUPLING_WEST,&
                  Z_COUPLING_SOUTH,&
                  Z_COUPLING_NORTH
       REAL(SP),DIMENSION(:,:,:),ALLOCATABLE::U_COUPLING_EAST,&
                  V_COUPLING_EAST,&
                  U_COUPLING_WEST,V_COUPLING_WEST,&
                  U_COUPLING_SOUTH,V_COUPLING_SOUTH,&
                  U_COUPLING_NORTH,V_COUPLING_NORTH,&
                  P_COUPLING_WEST,P_COUPLING_EAST,P_COUPLING_SOUTH,&
                  P_COUPLING_NORTH,&
                  W_COUPLING_WEST,W_COUPLING_EAST,W_COUPLING_SOUTH,&
                  W_COUPLING_NORTH
# if defined (SALINITY)
       REAL(SP),DIMENSION(:,:,:),ALLOCATABLE:: &
                  S_COUPLING_WEST,S_COUPLING_EAST,S_COUPLING_SOUTH,&
                  S_COUPLING_NORTH
# endif		  
# if defined (TEMPERATURE)	
       REAL(SP),DIMENSION(:,:,:),ALLOCATABLE:: &			
                  T_COUPLING_WEST,T_COUPLING_EAST,T_COUPLING_SOUTH,&
                  T_COUPLING_NORTH
# endif		

       REAL(SP) :: TIME_COUPLING_1,TIME_COUPLING_2
       INTEGER :: Kstart_EAST,Kend_EAST,Kshift_EAST,&
                  Kstart_WEST,Kend_WEST,Kshift_WEST,&
                  Kstart_SOUTH,Kend_SOUTH,Kshift_SOUTH,&
                  Kstart_NORTH,Kend_NORTH,Kshift_NORTH
       LOGICAL :: IN_DOMAIN_EAST,IN_DOMAIN_WEST,&
                  IN_DOMAIN_SOUTH,IN_DOMAIN_NORTH
# endif
	 
     ! hot start added by Cheng for hot start
     LOGICAL :: HOTSTART
     CHARACTER(LEN=80) :: Eta_HotStart_File,U_HotStart_File,V_HotStart_File,&
                          W_HotStart_File,P_HotStart_File,&
                          Rho_HotStart_File,TKE_HotStart_File,&
                          EPS_HotStart_File
# if defined (SALINITY)
     CHARACTER(LEN=80) :: Sali_HotStart_File,Temp_HotStart_File
# endif
# if defined (LANDSLIDE) || defined (FLUIDSLIDE) || defined (TWOLAYERSLIDE) || defined (LANDSLIDE_COMPREHENSIVE) || defined (DEFORMABLESLIDE)
     CHARACTER(LEN=80) :: Depth_HotStart_File
# endif
# if defined (FLUIDSLIDE)
     CHARACTER(LEN=80) :: Us_HotStart_File,Vs_HotStart_File
# endif

	 ! added by Cheng for limiting the maximum Froude number
# if defined (FROUDE_CAP)
     REAL(SP) :: FROUDECAP
# endif
	 
	 ! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     CHARACTER(LEN=80) :: Slide_File	
     real(SP), dimension(:,:), allocatable :: Uvs ,Vvs ,DUvs ,DVvs ,Dvs ,DelxH0,Hvs, &
	                                          Uvs0,Vvs0,DUvs0,DVvs0,Dvs0,DelyH0, &
                                              DelxUvs,DelxVvs,DelxDUvs,DelxDVvs, &
	                                          DelyUvs,DelyVvs,DelyDUvs,DelyDVvs, &
											  DelxHvs,SrcpvsX,SrchvsX, SrcdvsX,  &
											  DelyHvs,SrcpvsY,SrchvsY, SrcdvsY
     real(SP), dimension(:,:), allocatable :: DvsxL,DvsxR,UvsxL,UvsxR,VvsxL,VvsxR, &
											  H0fx,DUvsxL,DUvsxR,DVvsxL,DVvsxR, &
											  HvsxL,HvsxR,SvsxL,SvsxR,Evsx,Fvsx,Gvsx, &
											  EvsxL,EvsxR,FvsxL,FvsxR,GvsxL,GvsxR
     real(SP), dimension(:,:), allocatable :: DvsyL,DvsyR,UvsyL,UvsyR,VvsyL,VvsyR, &
	                                          H0fy,DUvsyL,DUvsyR,DVvsyL,DVvsyR, &
											  HvsyL,HvsyR,SvsyL,SvsyR,Evsy,Fvsy,Gvsy, &
											  EvsyL,EvsyR,FvsyL,FvsyR,GvsyL,GvsyR
     integer,  dimension(:,:), allocatable :: Maskvs
     real(SP) :: SLIDE_DENSITY,SLIDE_VISCOSITY,SLIDE_MINTHICK
# endif

# if defined (DEFORMABLESLIDE)
     CHARACTER(LEN=80) :: Slide_File
     CHARACTER(LEN=80) :: Us_HotStart_File,Vs_HotStart_File,Ws_HotStart_File
	 CHARACTER(LEN=80) :: RHEO_OPT
     real(SP), dimension(:,:), allocatable :: Ugs ,Vgs ,Wgs ,DUgs ,DVgs ,DWgs ,Dgs ,Rhogs ,Hgs  ,Wtgs,Qtgs, &
	                                          Ugs0,Vgs0,Wgs0,DUgs0,DVgs0,DWgs0,Dgs0,Rhogs0,Wgs00,Wbgs,Qbgs, &
                                              DelxUgs,DelxVgs,DelxWgs,DelxDUgs,DelxDVgs,DelxDWgs,DelxDgs,DelxHgs, &
	                                          DelyUgs,DelyVgs,DelyWgs,DelyDUgs,DelyDVgs,DelyDWgs,DelyDgs,DelyHgs, &
											  DelxH0,Delx2H0,Kx,Cxgs,SrcmgsX,SrcpgsX,SrctgsX,UgsP,VgsP,WgsP, &
											  DelyH0,Dely2H0,Ky,Cygs,SrcmgsY,SrcpgsY,SrctgsY,SrctgsZ,DwDt, &
											  Cxx,Cxy,Cxz,Cyx,Cyy,Cyz,Czx,Czy,Czz,Czz0, &
											  Taxx,Tayx,Taxy,Tayy,Taxz,Tayz,Tezz,Pss,Qsgs, &
											  Tbxx,Tbxy,Tbxz,Tbyx,Tbyy,Tbyz,Tbzx,Tbzy,Tbzz, &
											  PhiInt,PhiBed
     real(SP), dimension(:,:), allocatable :: DgsxL,DgsxR,HSgsxL,HSgsxR,UgsxL,UgsxR,VgsxL,VgsxR,WgsxL,WgsxR, &
											  H0fx,DUgsxL,DUgsxR,DVgsxL,DVgsxR,DWgsxL,DWgsxR, &
											  Cxgsx,SgsxL,SgsxR,Egsx,Fgsx,Ggsx,Hgsx, &
											  EgsxL,EgsxR,FgsxL,FgsxR,GgsxL,GgsxR,HgsxL,HgsxR, &
											  SrcmgsxL,SrcmgsxR,SrcpgsxL,SrcpgsxR
     real(SP), dimension(:,:), allocatable :: DgsyL,DgsyR,HSgsyL,HSgsyR,UgsyL,UgsyR,VgsyL,VgsyR,WgsyL,WgsyR, &
	                                          H0fy,DUgsyL,DUgsyR,DVgsyL,DVgsyR,DWgsyL,DWgsyR, &
											  Cygsy,SgsyL,SgsyR,Egsy,Fgsy,Ggsy,Hgsy, &
											  EgsyL,EgsyR,FgsyL,FgsyR,GgsyL,GgsyR,HgsyL,HgsyR, &
											  SrcmgsyL,SrcmgsyR,SrcpgsyL,SrcpgsyR
     real(SP), dimension(:,:), allocatable :: QbgsC,WgsC
     integer,  dimension(:,:), allocatable :: Maskgs
	 logical :: NON_HYDRO_SLD,DISP_CORR_SLD,REDU_GRAV_SLD,NON_HYDRO_UP
     real(SP) :: SLIDE_DENSITY,SLIDE_VISCOSITY,SLIDE_MINTHICK,SLIDE_GAMMA, &
	             SLIDE_CONC,GRAIN_DENSITY,PhiInt_A,PhiBed_A,SLIDE_LAMBDA, &
				 SLIDE_INIU,SLIDE_INIV,SLIDE_INIW,PhiInt_F,PhiBed_F
     ! poisson solvers
     integer  :: neqnsgs,IendC,JendC,MlocC,NlocC
     real(SP), dimension(:),   allocatable :: Rhsgs
     integer,  dimension(:),   allocatable :: JCoefgs
     real(SP), dimension(:,:), allocatable :: Coefgs
# endif
     End Module GLOBAL
