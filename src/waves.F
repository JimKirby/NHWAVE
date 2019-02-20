!------------------------------------------------------------------------------------------------
!
!   Waves.F
!
!   Subroutines related to specification of wavemaker conditions or processing
!   of wave information.
!
!   James Kirby, 2/7/2019
!
!-----------------------------------------------------------------------------------------------
!
!   This file is part of NHWAVE, and is copyrighted under the terms of the BSD 2-Clause License
!
!-----------------------------------------------------------------------------------------------
!
!   BSD 2-Clause License
!
!   Copyright (c) 2019, NHWAVE Development Group
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!----------------------------------------------------------------------------------------------
!
!   Subroutines included in this file:
!
!   (1) wave_dispersion      Calculate wavenumber given depth anf frequency

!   cnoidal_wave_left_boundary
!
!
!----------------------------------------------------------------------------------------------
!
!   (1) Subroutine wave_dispersion
!
!       Solve the linear wave dispersion relation \omega^2 - gk\tanh kh = 0
!       by Newton-Raphson.  Initial guess due to Eckart (1951).
!
!       J. Kirby, 2/7/19, 2/11/19
!
!----------------------------------------------------------------------------------------------
!
      subroutine wave_dispersion(Segma, Wave_Number)
!
      use global, only: SP,Dep_Wave, Grav
      implicit none
      real(SP) :: Segma, Wave_Number,Wn_new, error,Fk,Fkdif
!
!  Initial guess (Eckart, 1951)
!
      Wave_Number = (Segma*Segma)/(Grav*(tanh(Segma*Segma*Dep_Wave/grav))**0.5)
      error = 100.
      do while ( error >= 1.0e-8  )
         Fk = Segma**2 - Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)
         Fkdif = -Grav*( tanh(Wave_Number*Dep_Wave) +  (Wave_Number*Dep_wave)/  &
                  cosh(Wave_Number*Dep_Wave)**2 )
         Wn_new = Wave_Number - Fk/Fkdif
         error = abs( (Wave_Number - Wn_new)/Wave_Number )
         Wave_Number = Wn_new
      end do

      end subroutine wave_dispersion
!
!-------------------------------------------------------------------------------------------
!
!    Subroutine cnoidal_wave_left_boundary
!
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 06/02/2011, Gangfeng Ma
!------------------------------------------------------------------------------------------
!
     subroutine cnoidal_wave_left_boundary
!
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer :: j,k,nn
     real(SP) :: Wave_Length,Celerity,Ytrough,Mod1,Wave_Number,Zup,Zlow,Zmid,Xstart, &
                 Zero1,Zlev,Xcn,Xsn,Xdn,cnoidal_cn,cnoidal_ck,Stokes_Drift,Fact

      call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)
!
! wave number
!
      Wave_Number = 2.0*pi/Wave_Length
!
! find zero start
!
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
          goto 200 !goto's???
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
!
! mean mass transport
!
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
!
!--------------------------------------------------------------------------------------------
!
!    Subroutine cnoidal
!
!    Compute elliptic integrals needed for cnoidal wave solution
!
!    Called by: cnoidal_wave_left_boundary
!
!    Last update: 06/02/2011, Gangfeng Ma
!
!--------------------------------------------------------------------------------------------
!
     subroutine cnoidal(height,depth,period,l,c,x2,mod)
!
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
!
!-------------------------------------------------------------------------------------------------------
!

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
!
!------------------------------------------------------------------------
!
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
!
!----------------------------------------------------------------
!

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
!
!----------------------------------------------------------------------
!
!    Subroutine random_wave_left_boundary
!
!    Specify random waves on left boundary
!
!    Called by: update_wave_bc
!
!    Last update: 04/11/2011, Gangfeng Ma
!                 02/11/2019, J. Kirby remove dispersion calculation to routine wave_dispersion
!
!--------------------------------------------------------------------------------------------------
!
     subroutine random_wave_left_boundary
!
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy,yc, &
                       MaxNumFreq,MaxNumDir,Nglob,Dep_Wave,Random_Phs
# if defined (PARALLEL)
     use global, only: myid
# endif
     implicit none
     integer :: nfreq,ndir,Iter,nk,i,j,k
     real(SP) :: Per_Wave,Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Angle,tmp1,tmp2,tmp3,tmp4,Fact, &
                 Wtheta(MaxNumDir,MaxNumFreq),Wh(MaxNumDir,MaxNumFreq),Stokes_Drift(MaxNumDir,MaxNumFreq), &
                 Wnum(MaxNumDir,MaxNumFreq),Phs_Lag,dfreq,ddir,Zlev

!    find the right angle for periodic bc
     do nfreq = 1,NumFreq
     do ndir = 1,NumDir
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       call wave_dispersion(Segma, Wave_Number)
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
         Phs_Lag = yc(j)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
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
         Phs_Lag = yc(j)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
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
!
!------------------------------------------------------------------------------------------------------------------
!

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
     call wave_dispersion(Segma, Wave_Number)
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
!
!----------------------------------------------------------------------------------------------------------------
!


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
       call wave_dispersion(dep_wave, Segma, Wave_Number)
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave
!
!   root-mean-square wave height   ! not defined for a single component
!
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
!
!------------------------------------------------------------------------------------------------------------
!
!    Subroutine focused_wave_left_boundary
!
!    Used to specify a focusing wave packet on the left boundary
!    left boundary
!
!    Called by: update_wave_bc
!
!    Last update: 03/15/2013, Derakhti Morteza (derakhti@udel.edu)
!
!------------------------------------------------------------------------------------------------------------
!
     subroutine focused_wave_left_boundary
!
     use global, only: Amp_i,k_i,f_i,nwave,x_breaking,t_breaking,depth_comp,sd_return, &
          SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
          Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,myid
!
# if defined (INTEL)
     use IFPORT, only: rand
# endif
!
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
!
!------------------------------------------------------------------------------------------------------------
!
!    Subroutine wave_current_left_boundary
!
!    Used to specify (what?) on left boundary
!
!    Called by:  update_wave_bc
!
!    Last update: 05/02/2011, Gangfeng Ma
!
!-----------------------------------------------------
!
     subroutine wave_current_left_boundary
!
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
     call wave_dispersion(Segma, Wave_Number)
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = yc(j)*sin(Theta_Wave*pi/180.)*Wave_Number
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

!
!-----------------------------------------------------------------------------------------------------------------
!

     subroutine linear_wave_right_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Iend1,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy,yc
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave

     call wave_dispersion(Segma, Wave_Number)

     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = yc(j)*sin(Theta_Wave*pi/180.)*Wave_Number
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
!
!----------------------------------------------------------------------------------------------------------------
!

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
     call wave_dispersion(Segma, Wave_Number)

     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = yc(j)*sin(Theta_Wave*pi/180.)*Wave_Number
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
!
!--------------------------------------------------------------------------------------------------------------
!
!    Subroutine solitary wave left boundary
!
!    Used to specify solitary wave on left boundary
!
!    Called by:  update_bc
!
!    Last update: 22/12/2010, Gangfeng Ma
!
!-------------------------------------------------------------------------------------------------------------
!
     subroutine solitary_wave_left_boundary
!
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



