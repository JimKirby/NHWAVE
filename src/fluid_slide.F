! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
      subroutine nlmassF90()
      use global
      implicit none
      real(SP) :: Rx, Ry, GXX, XM, XN
      integer ilo,i,iup, jlo,j,jup
      integer dilo,diup,djlo,djup
	  
      GXX=SLIDE_MINTHICK
      Rx=SLIDE_RKDT/dx
      Ry=SLIDE_RKDT/dy


      !ilo=max(Ibeg,3)
      !iup=min(Iend,Mglob-1)
      !jlo=max(Jbeg,3)
      !jup=min(Jend,Nglob-1)

      dilo=0; djlo=0; diup=0; djup=0

!# if defined (PARALLEL)
!      if(npx==0) dilo=2
!      if(npy==0) djlo=2
!      if((npx+1)*(Mloc-2*Nghost)==Mglob) diup=-1
!      if((npy+1)*(Nloc-2*Nghost)==Nglob) djup=-1
!# endif

      ilo=Ibeg+dilo; iup=Iend+diup;
      jlo=Jbeg+djlo; jup=Jend+djup;

      do j=jlo,jup
      do i=ilo,iup



        SldDZ1(I,J)=SLIDE_RKC1*SldDZ0(I,J)+SLIDE_RKC2*SldDZs(I,J) &
	    -Rx*(SldMs(i,j)-SldMs(i-1,j))-Ry*(SldNs(i,j)-SldNs(i,j-1)) &
              + (SldDZs(I-1,J)-2.0d0*SldDZs(I,J)+SldDZs(I+1,J))*Rx/dx*SLIDE_DIFF &
              + (SldDZs(I,J-1)-2.0d0*SldDZs(I,J)+SldDZs(I,J+1))*Ry/dy*SLIDE_DIFF

        SldZ1(i,j)=SldDZ1(i,j)-Hc0(i,j)

        if(SldDZ1(i,j).lt.GXX) then
          SldDZ1(i,j)=0.0d0
          SldZ1(i,j)=-Hc0(i,j)
        endif
      enddo
      enddo

      ! ghost cells
      call phi_2D_coll(SldDZ1)
      call phi_2D_coll(SldZ1)

      end subroutine nlmassF90


subroutine nlmmntF90()
      use global
      implicit none
      real(SP) :: R, GXX,CFL_Cap,sCFL
      integer ilo,i,iup, jlo,j,jup
      integer dilo,diup,djlo,djup
	  
      real(SP) SldDM1(Mloc,Nloc), SldDN1(Mloc,Nloc)

      logical LBANK(Mloc,Nloc)
      logical RBANK(Mloc,Nloc)
      logical TBANK(Mloc,Nloc)
      logical BBANK(Mloc,Nloc)


      real(SP) Rx, Ry, RHO1, RHO2
      real(SP) FN, FM, MU, FF, XM, XN, XNE, XNW, XNN, XVE, XVW, XVN, XVV
      real(SP) XMM, XMS, XMN, XME, XUM, XUS, XUN, XUU
      real(SP) MaxVel, tmpVel, tmpN, tmpM

      real(SP) SLDalpha, SLDbeta, SW_REO, CC
	  
      SW_REO=6.0d0/5.0d0
      CFL_Cap = 0.5d0

      Rx=SLIDE_RKDT/dx
      Ry=SLIDE_RKDT/dy
      FM=SLIDE_ROUGHNESS
      Mu=SLIDE_VISCOSITY
      GXX=SLIDE_MINTHICK 
      RHO1=Rho0 !WATER_DENSITY
      RHO2=SLIDE_DENSITY


      MAXVEL=0.0d0
      !ilo=max(Ibeg-1,1)
      !iup=min(Iend+1,Mglob-1)
      !jlo=max(Jbeg-1,1)
      !jup=min(Jend+1,Nglob)

!      dilo=-1; djlo=-1; diup=1; djup=1
      dilo=0; djlo=0; diup=0; djup=0

!# if defined (PARALLEL)
!      if(npx==0) dilo=0
!      if(npy==0) djlo=0
!      if((npx+1)*(Mloc-2*Nghost)==Mglob) diup=-1
!      if((npy+1)*(Nloc-2*Nghost)==Nglob) djup=0
!# endif

      ilo=Ibeg+dilo; iup=Iend+diup;
      jlo=Jbeg+djlo; jup=Jend+djup;

      DO 10 J=jlo,jup
         DO 10 I=ilo,iup
            SldDM1(I,J)=0.5*(SldDZ1(I,J)+SldDZ1(I+1,J))
           
            RBANK(I,J)=(Hc0(I,J).LE.Hc0(I+1,J)+1.0D-5)
            LBANK(I,J)=(Hc0(I+1,J).LE.Hc0(I,J)+1.0D-5)

            IF(SldDM1(I,J).GT.GXX) THEN
               SldUs(I,J)=SldMs(I,J)/SldDM1(I,J)
! FroudCap in x
               sCFL=DT*ABS(SldUs(I,J))/DX
               IF(sCFL>CFL_Cap)THEN
                 SldUs(I,J)=SIGN(CFL_Cap*DX/DT,SldUs(I,J))
               ENDIF
            ELSE
               SldUs(I,J)=0.0d0
            ENDIF
            MAXVEL=MAX1(ABS(SldUs(i,j)),MAXVEL)
      10 CONTINUE

      !ilo=max(Ibeg-1,1)
      !iup=min(Iend+1,Mglob)
      !jlo=max(Jbeg-1,1) 
      !jup=min(Jend+1,Nglob-1)

      ! dilo=-1; djlo=-1; diup=1; djup=1
	  dilo=0; djlo=0; diup=0; djup=0

!# if defined (PARALLEL)
!      if(npx==0) dilo=0
!      if(npy==0) djlo=0
!      if((npx+1)*(Mloc-2*Nghost)==Mglob) diup=0
!      if((npy+1)*(Nloc-2*Nghost)==Nglob) djup=-1
!# endif

      ilo=Ibeg+dilo; iup=Iend+diup;
      jlo=Jbeg+djlo; jup=Jend+djup;

      DO 20 J=jlo,jup
         DO 20 I=ilo,iup
            SldDN1(I,J)=0.5D0*(SldDZ1(I,J)+SldDZ1(I,J+1))
           
            TBANK(I,J)=(Hc0(I,J).LE.Hc0(I,J+1)+1.0D-5)
            BBANK(I,J)=(Hc0(I,J+1).LE.Hc0(I,J)+1.0D-5)

            IF(SldDN1(I,J).GT.GXX) THEN
               SldVs(I,J)=SldNs(I,J)/SldDN1(I,J)
! FroudCap in y
               sCFL=DT*ABS(SldVs(I,J))/DY
               IF(sCFL>CFL_Cap)THEN
                 SldVs(I,J)=SIGN(CFL_Cap*DY/DT,SldVs(I,J))
               ENDIF
            ELSE
               SldVs(I,J)=0.0d0
            ENDIF
            MAXVEL=MAX1(ABS(SldVs(i,j)),MAXVEL)
      20 CONTINUE


      IF(DT>dx/MAXVEL) THEN
          PRINT *,DT, '>', dx/MAXVEL
      ENDIF


      call phi_2D_coll(SldVs)
      call phi_2D_coll(SldUs)

      ilo=Ibeg; iup=Iend;
      jlo=Jbeg; jup=Jend;

      DO 21 J=jlo,jup
         DO 21 I=ilo,iup
               SldUs(I,J)=SUM(SldUs(I-1:I+1,J-1:J+1))/9.0d0
               SldVs(I,J)=SUM(SldVs(I-1:I+1,J-1:J+1))/9.0d0
      21 CONTINUE
      call phi_2D_coll(SldVs)
      call phi_2D_coll(SldUs)


      FN=Grav*FM**2.0

      !   ------- CAL. OF LINEAR TERMS (Y-DIRECTION) -------      
!      ilo=max(Ibeg,3       )
!      iup=min(Iend,Mglob-2 )
!      jlo=max(Ibeg,3       )
!      jup=min(Jend,Nglob-1 )

      dilo=0; djlo=0; diup=0; djup=0

!# if defined (PARALLEL)
!      if(npx==0) dilo=2
!      if(npy==0) djlo=2
!      if((npx+1)*(Mloc-2*Nghost)==Mglob) diup=-2
!      if((npy+1)*(Nloc-2*Nghost)==Nglob) djup=-1
!# endif

      ilo=Ibeg+dilo; iup=Iend+diup;
      jlo=Jbeg+djlo; jup=Jend+djup;

      DO 200 J=jlo,jup
         DO 200 I=ilo,iup

            XNN=0.25D0*(SldNs(I,J)+SldNs(I+1,J)+SldNs(I,J-1)+SldNs(I+1,J-1))
            FF=DT*(FN*SQRT(SldMs(I,J)**2+XNN**2)/MAX(SldDM1(I,J),GXX)**(7.0/3.0) + &
                                           3.0d0*MU/MAX(SldDM1(I,J),GXX)**2.0d0 )

            IF(((Hc0(I,J)+SldZ1(I+1,J).GE.0.0D0).AND.RBANK(I,J)).OR.   &
               ((SldDZ1(I,J).GT.0.0D0).AND.RBANK(I,J)).OR.            &
               ((Hc0(I+1,J)+SldZ1(I,J).GE.0.0D0).AND.LBANK(I,J)).OR.   &
               ((SldDZ1(I+1,J).GT.0.0D0).AND.LBANK(I,J))) THEN

               IF(MASK(I+1,J)*MASK(I,J).GT.0) THEN
                   SLDalpha=(RHO2-RHO1)/RHO2
                   SLDbeta=RHO1/RHO2
               ELSE
                   SLDalpha=1.0d0
                   SLDbeta=0.0d0
               ENDIF

               XM=(1.0D0-FF)*( SLIDE_RKC1*SldM0(I,J)+SLIDE_RKC2*SldMs(I,J) )- &
                  Grav*Rx*SldDM1(I,J)*( SLDalpha*(SldZ1(I+1,J)-SldZ1(I,J))  + &
                                        SLDbeta*(SldWZ(I+1,J)-SldWZ(I,J)) ) + &
                       SldDM1(I,J)*( SldUs(I-1,J)-2.0d0*SldUs(I,J)+SldUs(I+1,J))*Rx/dx*MU + &
                       SldDM1(I,J)*( SldUs(I,J-1)-2.0d0*SldUs(I,J)+SldUs(I,J+1))*Ry/dy*MU

               IF(SldMs(I,J).LT.0.0D0) THEN
                  XM=XM-Rx*(    -SldMs(I+2,J)*SldUs(I+2,J) &
                            +4.0*SldMs(I+1,J)*SldUs(I+1,J) &
                            -3.0*SldMs(I,J)*SldUs(I,J))/2.0*SW_REO
                  !XM=XM-Rx*(SldMs(I+1,J)*SldUs(I+1,J)-SldMs(I,J)*SldUs(I,J))*SW_REO
               ELSEIF(SldMs(I,J).GT.0.0D0) THEN
		  XM=XM-Rx*( 3.0*SldMs(I,J)*SldUs(I,J)     &
                            -4.0*SldMs(I-1,J)*SldUs(I-1,J) &
                                +SldMs(I-2,J)*SldUs(I-2,J))/2.0*SW_REO
                  !XM=XM-Rx*(SldMs(I,J)*SldUs(I,J)-SldMs(I-1,J)*SldUs(I-1,J))*SW_REO
               ELSE
                  tmpM=SldMs(I-1,J)+SldMs(I+1,J)
                  IF(tmpM.LT.0.0D0) THEN
                     XM=XM-Rx*(SldMs(I+1,J)*SldUs(I+1,J)-0.0D0)*SW_REO
                  ELSE
                     XM=XM-Rx*(0.0D0-SldMs(I-1,J)*SldUs(I-1,J))*SW_REO
                  ENDIF
               ENDIF

               XNE=0.25D0*(SldVs(I,J+1)+SldVs(I+1,J+1)+SldVs(I,J)+SldVs(I+1,J))
               XNN=0.25D0*(SldVs(I,J-1)+SldVs(I+1,J-1)+SldVs(I,J-2)+SldVs(I+1,J-2))
               XM=XM-Ry*(SldMs(I,J+1)*XNE-SldMs(I,J-1)*XNN)*SW_REO/2.0

               SldM1(I,J)=XM/(1.0D0+FF)
            ELSE
               SldM1(I,J)=0.0D0
            ENDIF
200      CONTINUE

      !   ------- CAL. OF LINEAR TERMS (Y-DIRECTION) -------
!      ilo=max(Ibeg,3      )
!      iup=min(Iend,Mglob-1)
!      jlo=max(Jbeg,3      )
!      jup=min(Jend,Nglob-2)

      dilo=0; djlo=0; diup=0; djup=0

!# if defined (PARALLEL)
!      if(npx==0) dilo=2
!      if(npy==0) djlo=2
!      if((npx+1)*(Mloc-2*Nghost)==Mglob) diup=-1
!      if((npy+1)*(Nloc-2*Nghost)==Nglob) djup=-2
!# endif

      ilo=Ibeg+dilo; iup=Iend+diup;
      jlo=Jbeg+djlo; jup=Jend+djup;


      DO 220 J=jlo,jup
         DO 220 I=ilo,iup

            XMM=0.25D0*(SldMs(I,J)+SldMs(I,J+1)+SldMs(I-1,J)+SldMs(I-1,J+1))
            FF=DT*(FN*DSQRT(SldNs(I,J)**2+XMM**2)/MAX(SldDN1(I,J),GXX)**(7.0/3.0) + & 
                                                3.0d0*MU/MAX(SldDN1(I,J),GXX)**2.0d0 )

            IF( ((Hc0(I,J)+SldZ1(I,J+1).GE.0.0D0).AND.TBANK(I,J)).OR. 	&
                ((SldDZ1(I,J).GT.0.0D0).AND.TBANK(I,J)).OR.             & 
                ((Hc0(I,J+1)+SldZ1(I,J).GE.0.0D0).AND.BBANK(I,J)).OR.    &
                ((SldDZ1(I,J+1).GT.0.0D0).AND.BBANK(I,J))) THEN	

               IF(MASK(I,J)*MASK(I,J+1).GT.0) THEN
                   SLDalpha=(RHO2-RHO1)/RHO2
                   SLDbeta=RHO1/RHO2
               ELSE
                   SLDalpha=1.0d0
                   SLDbeta=0.0d0
               ENDIF

               XN=(1.0D0-FF)*( SLIDE_RKC1*SldN0(I,J)+SLIDE_RKC2*SldNs(I,J) ) - &
                         Grav*Ry*SldDN1(I,J)*( SLDalpha*(SldZ1(I,J+1)-SldZ1(I,J))   + &
                                               SLDbeta*(SldWZ(I,J+1)-SldWZ(I,J)) ) + &
                              SldDN1(I,J)*( SldVs(I-1,J)-2.0d0*SldVs(I,J)+SldVs(I+1,J))*Rx/dx*MU + &
                              SldDN1(I,J)*( SldVs(I,J-1)-2.0d0*SldVs(I,J)+SldVs(I,J+1))*Ry/dy*MU

               IF(SldNs(I,J).LT.0.0D0) THEN
                  XN=XN-Ry*(    -SldNs(I,J+2)*SldVs(I,J+2) &
                            +4.0*SldNs(I,J+1)*SldVs(I,J+1) &
                            -3.0*SldNs(I,J)*SldVs(I,J))/2.0*SW_REO
                  !XN=XN-Ry*(SldNs(I,J+1)*SldVs(I,J+1)-SldNs(I,J)*SldVs(I,J))*SW_REO
               ELSEIF(SldNs(I,J).GT.0.0D0) THEN
                  XN=XN-Ry*( 3.0*SldNs(I,J)*SldVs(I,J)     &
                            -4.0*SldNs(I,J-1)*SldVs(I,J-1) &
                                +SldNs(I,J-2)*SldVs(I,J-2))/2.0*SW_REO
                  !XN=XN-Ry*(SldNs(I,J)*SldVs(I,J)-SldNs(I,J-1)*SldVs(I,J-1))*SW_REO
               ELSE
                  tmpN=SldNs(I,J-1)+SldNs(I,J+1)
                  IF(tmpN.LT.0.0D0) THEN
                     XN=XN-Ry*(SldNs(I,J+1)*SldVs(I,J+1)-0.0D0)*SW_REO
                  ELSE
                     XN=XN-Ry*(0.0D0-SldNs(I,J-1)*SldVs(I,J-1))*SW_REO
                  ENDIF
               ENDIF

               XME=0.25D0*(SldUs(I+1,J)+SldUs(I+1,J+1)+SldUs(I,J)+SldUs(I,J+1))
               XMM=0.25D0*(SldUs(I-1,J)+SldUs(I-1,J+1)+SldUs(I-2,J)+SldUs(I-2,J+1))
               XN=XN-Rx*(SldNs(I+1,J)*XME-SldNs(I-1,J)*XMM)*SW_REO/2.0

               SldN1(I,J)=XN/(1.0D0+FF)
            ELSE
               SldN1(I,J)=0.0D0
            ENDIF
220      CONTINUE

      ! ghost cells
      call phi_2D_coll(SldN1)
      call phi_2D_coll(SldM1)


      end subroutine nlmmntF90
# endif
