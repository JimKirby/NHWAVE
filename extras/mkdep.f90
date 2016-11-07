    program generate_bathymetry
    implicit none

    integer, parameter :: im = 500
    integer, parameter :: jm = 90
    integer :: i,j
    real(8) :: dx,dy,x,y,h0,h02,tht,cs,sn,pi,  &
               ytmp,xp,yp,ctmp1,ctmp2,depth(im,jm)    
    character(len=80) :: test_case

    pi = 2.*acos(0.)

    test_case = 'landslide'

    if(trim(test_case)=='submerged_bar') then
      dx = 0.02
      do j = 1,jm
      do i = 1,im
        x = (i-1)*dx
        if(x<=6.0) then
          depth(i,j) = 0.4
        elseif(x>6.0.and.x<=12.0) then
          depth(i,j) = 0.4-(x-6.0)/20.0
        elseif(x>12.0.and.x<=14.0) then
          depth(i,j) = 0.1
        elseif(x>14.0.and.x<17.0) then
          depth(i,j) = 0.4-(17.0-x)/10.0
        else
          depth(i,j) = 0.4
        endif
      enddo
      enddo
    elseif(trim(test_case)=='landslide') then
      dx = 0.02
      do j = 1,jm
      do i = 1,im
        x = (i-1)*dx+0.5*dx
        if(x<=1.5/tan(15.*pi/180.)+0.1) then
          depth(i,j) = (x-0.1)*tan(15.*pi/180.)
        else
          depth(i,j) = 1.5
        endif
      enddo
      enddo
    elseif(trim(test_case)=='breaking_solitary_wave') then
      dx = 0.02
      do j = 1,jm
      do i = 1,im
        x = (i-1)*dx
        if(x<=4.0) then
          depth(i,j) = 0.21
        else
          depth(i,j) = 0.21-(x-4.0)/20.
        endif
      enddo
      enddo
    elseif(trim(test_case)=='dam_break') then
      dx = 0.05
      do j = 1,jm
      do i = 1,im
        x = (i-1)*dx
        if(x<=25.5) then
          depth(i,j) = 0.01
        elseif(x>25.5.and.x<=28.5) then
          depth(i,j) = -(x-25.5)*0.4/3.0
        elseif(x>28.5.and.x<=31.5) then
          depth(i,j) = -(31.5-x)*0.4/3.0
        elseif(x>31.5) then
          depth(i,j) = 0.01
        endif
      enddo
      enddo
    elseif(trim(test_case)=='tidal_wave') then
      dx = 7.5
      do j = 1,jm
      do i = 1,im
        x = (i-1)*dx
        if(abs(x-750)<=187.5) then
          depth(i,j) = 8.0
        else
          depth(i,j) = 16.0
        endif
      enddo
      enddo
    elseif(trim(test_case)=='Berkhoff_Shoal') then
      dx = 0.025
      dy = 0.05
      h0 = .45
      h02 = 0.07
      tht = 70.*3.14159/180
      cs    = cos(tht)
      sn    = sin(tht)
      ytmp  = 5.84-(h0-h02)/.02
      do i = 1, im
      do j = 1, jm
        x  = float(i-1)*dx- 12.0
        y  = float(j-1)*dy- 10.0
        xp =  x*cs+y*sn
        yp = -x*sn+y*cs
        ctmp1 = (xp/4.)**2+(yp/3.)**2
        ctmp2 = (xp/5.)**2+(yp/3.75)**2
        if(yp.ge.5.84) then
          depth(i,j) = h0
        elseif (yp.le.ytmp) then
          depth(i,j) = h02
        elseif (ctmp1.gt.1.) then
          depth(i,j) = h0-0.02*(5.84-yp)
        else
          depth(i,j) = .75-.02*(5.84-yp)-0.5*sqrt(1.-ctmp2)
        endif
      enddo
      enddo
    endif

    open(11,file='depth.txt')
    do j = 1,jm
      write(11,'(5000f8.4)') (depth(i,j),i=1,im)
    enddo
    close(11)

    end
