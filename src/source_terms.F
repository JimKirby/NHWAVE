    subroutine source_terms
!-----------------------------------------------------------------------------------------------------------
!    This subroutine is used to evaluate sources
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!
!    Note that the wavemaker code needs to be split between here and initialize.F
!    Wavemaker parmameters are getting recomputed every time step.  This is what they
!    call BAD.
!
!    Change log:
!           2/7/19 Kirby modifying wavenumber solver using Eckart's initial guess.
!                  (Should iterate to 10-6 accuracy in two iterations.
!                  Change Segma to Omega so I can use Sigma when there is a current
!                  entering through the boundary (not yet - Segma also appears elsewhere)
!----------------------------------------------------------------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,Iter,nn,ndir,nfreq,nk
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Source_Area,myvar, &
                 WnumX,WnumY,Phs_lag,dfreq,ddir,Angle,tmp1,tmp2,tmp3,tmp4,Umag,Ubar,Vbar
     real(SP) :: Ytrough,Mod1,Zup,Zlow,Zmid,Xstart,Zero1,cnoidal_cn,cnoidal_ck,Atmp
!
!  All wave info here needs to be moved to a routine called by initialize except for final calculation
!  of the source

!  internal wavemaker for linear wave
!
     if(WaveMaker(1:7)=='INT_LIN') then
!
!  Find wave number for linear wave (Newton-Ralphson Method)
!
       Segma = 2.0*pi/Per_Wave
       call wave_dispersion(Segma, Wave_Number)
!
!       Celerity = sqrt(Grav*Dep_Wave)
!       Wave_Length = Celerity*Per_Wave
!       Wave_Number = 2.0*pi/Wave_Length
!
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
           Phs_lag = yc(j)*WnumY
           SourceC(i,j) = Celerity*Amp_Wave/Source_Area*cos(pi/2-Segma*time+Phs_lag)
         endif
       enddo
       enddo
     endif
!
! internal wavemaker for random waves
!
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
             call wave_dispersion(Segma, Wave_Number)
!            Celerity = sqrt(Grav*Dep_Wave)
!            Wave_Length = Celerity*Per_Wave
!            Wave_Number = 2.0*pi/Wave_Length
       
!            Iter = 0
!  75        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
!            if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
!            Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
!              Grav*tanh(Wave_Number*Dep_Wave)
!             Wave_Number = Wave_Number-Fk/Fkdif
!             Iter = Iter+1
!             goto 75
!   85        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

!   adjust wave direction for periodic bc
!   change as in funwave - to do (Kirby, 2/7/19)
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

!   judge between nk-1 and nk which is closer
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

!   judge between nk-1 and nk which is closer
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

!    calculate root-mean-square wave height for each component  (there is no such thing - Kirby)
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

             Phs_lag = (dy/2.0+yc(j))*WnumY
             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
           enddo
           enddo
         endif
       enddo
       enddo
     endif
!
! Added by Cheng Zhang, internal wavemaker for irregular waves
!  Why is this different from INT_SPC? (Kirby)
!
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
             call wave_dispersion(Segma, Wave_Number)
!             Celerity = sqrt(Grav*Dep_Wave)
!             Wave_Length = Celerity*Per_Wave
!             Wave_Number = 2.0*pi/Wave_Length
       
!             Iter = 0
!   95        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
!             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 105
!             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
!               Grav*tanh(Wave_Number*Dep_Wave)
!             Wave_Number = Wave_Number-Fk/Fkdif
!             Iter = Iter+1
!             goto 95
!   105        continue
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
             Phs_lag = (dy/2.0+yc(j))*WnumY
			 
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
             call wave_dispersion(Segma, Wave_Number)
!             Celerity = sqrt(Grav*Dep_Wave)
!             Wave_Length = Celerity*Per_Wave
!             Wave_Number = 2.0*pi/Wave_Length
       
!             Iter = 0
 !76          Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
!             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
!             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
!               Grav*tanh(Wave_Number*Dep_Wave)
!             Wave_Number = Wave_Number-Fk/Fkdif
!             Iter = Iter+1
!             goto 76
! 86          continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

!   root-mean-square wave height (not defined)
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
	   !hfx0 is flux of bed, which is different from hfx at wet-dry front (Cheng)
       SourceX(i,j) = Grav*0.5*(EtaxL(i+1,j)+EtaxR(i,j))*(hfx0(i+1,j)-hfx0(i,j))/dx*Mask(i,j)
       SourceY(i,j) = Grav*0.5*(EtayL(i,j+1)+EtayR(i,j))*(hfy0(i,j+1)-hfy0(i,j))/dy*Mask(i,j)
     enddo
     enddo

     end subroutine source_terms
 
