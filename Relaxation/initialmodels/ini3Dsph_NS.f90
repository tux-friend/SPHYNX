
   PROGRAM ini3DSPH

     IMPLICIT NONE

     INTEGER, PARAMETER :: np = 100000       !*** number of desired SPH particles
     INTEGER, PARAMETER :: nm = 18400        !*** rows in your 1D data file

     DOUBLE PRECISION, PARAMETER :: nv = 100

     DOUBLE PRECISION,PARAMETER :: pi = acos(-1.d0),msun=1.989d33,rsun=6.96d10
     DOUBLE PRECISION,PARAMETER :: K2MeV=8.617342857d-11


     INTEGER i,count
     DOUBLE PRECISION pmass,massacum,ir,idens,ivrad,itemp,iye,ient
     DOUBLE PRECISION h,ran,a1,a2
     DOUBLE PRECISION xcm,ycm,zcm
     DOUBLE PRECISION, DIMENSION(nm) :: m,r,vrad,dens,temp,s,ye,u,p
     DOUBLE PRECISION, DIMENSION(np) :: x,y,z,vx,vy,vz,radius
     CHARACTER admy*7


     open(1,file='1D_star_profile')              	!*** 1D input data file
     open(2,file='neutron_star',form='unformatted')  !*** 3D output data file
     open(3,file='neutron_ascii')

! Read HEADERS
     !do i=1,23
     !   read(1,'(a1)') admy
     !enddo

     do i=1,nm
        read(1,*) r(i),dens(i),m(i),p(i),temp(i),ye(i),u(i)
     enddo
     print *,'Last shell: ',r(nm)/1.d5,m(nm)/msun
     close(1)
!*** END: adapt to your 1D data file structure


     !Mass of each particle
     pmass=m(nm)/dble(np)

     !Accumulated mass
     massacum=0.d0

     !3D particle position randomized
     do i=1,np
        massacum=massacum+pmass
        call profile_r(m,r,nm,massacum,ir) !find interpolated radius

        if(mod(i,2).ne.0)then
           call random_number(ran)
           a1=ran*2.d0*pi
           call random_number(ran)
           a2=acos(2.d0*ran-1.d0)
        else               !Symmetrizes odd particles with respect
           a1=pi+a1        !randomly distributed even particles.
           a2=pi-a2        !In this way the CM is almost at 0 always.
        endif

        x(i)=ir*sin(a2)*cos(a1)
        y(i)=ir*sin(a2)*sin(a1)
        z(i)=ir*cos(a2)

        radius(i)=ir
     enddo

     massacum=0.d0
     count=0

     xcm=0.d0
     ycm=0.d0
     zcm=0.d0

     !Interpolate particle physical properties at its radius
     do i=1,np
        massacum=massacum+pmass

        if(mod(i,100000).eq.0)print*, i,massacum/msun

        !*** adapt to interpolate the magnitudes needed
        call profile(m,dens,vrad,ye,temp,u,nm,massacum,&
             &idens,ivrad,iye,itemp,ient)

        !*** if hydrostatic velocities are 0 and not needed
        vx(i)=ivrad*x(i)/radius(i)
        vy(i)=ivrad*y(i)/radius(i)
        vz(i)=ivrad*z(i)/radius(i)

        h=(3.d0*nv*pmass/32.d0/pi/idens)**(1.d0/3.d0)

        count=count+1
        write(2) x(i),y(i),z(i),pmass,h,idens,itemp,&
			&iye,ient
        write(3,'(10(1x,es14.7))') x(i),y(i),z(i),radius(i),pmass,h,&
			&idens,itemp,iye,ient
        xcm=xcm+x(i)
        ycm=ycm+y(i)
        zcm=zcm+z(i)
     enddo
     print *,count,dble(count)/dble(np)
     print *,xcm/dble(np),ycm/dble(np),zcm/dble(np)

     close(2)
     close(3)
   end PROGRAM ini3DSPH
! **********************************************************************

   SUBROUTINE profile_r(pm,r,nm,pmasa,rrr)
     implicit none
     double precision pm,r,pmasa,value1,rrr,value2,value3
     integer i,j,nm
     dimension pm(nm),r(nm)

     if(pmasa.lt.pm(1))  then

        value2=log10(pm(2))-log10(pm(1))
        value3=log10(pmasa)-log10(pm(1))

        value1=(log10(r(2))-log10(r(1)))/value2
        rrr=10.d0**(log10(r(1))+value1*value3)

     else if(pmasa.gt.pm(nm)) then
        value2=log10(pm(nm))-log10(pm(nm-1))
        value3=log10(pmasa)-log10(pm(nm))

        value1=(log10(r(nm))-log10(r(nm-1)))/value2
        rrr=10.d0**(log10(r(nm))+value1*value3)

     else
        j=1

        do while (pmasa.ge.pm(j+1))
           j=j+1
        enddo

        value1=(log10(pmasa)-log10(pm(j)))/(log10(pm(j+1))-log10(pm(j)))
        rrr=10.d0**(log10(r(j))+(log10(r(j+1))-log10(r(j)))*value1)

     endif

     return
   end SUBROUTINE profile_r


! **********************************************************************

   !*** adapt this subroutine to interpolate only the magnitudes needed
   SUBROUTINE profile(pm,dens,vrad,ye,t,u,nm,pmass,uro,uvrad,uye,ut,uenr)
     implicit none
     double precision pm,u,t,s,dens,pmass,value1,ut,uro,uye,ye,us,&
          &     vrad,uvrad,value2,value3,pi,uenr
     integer i,j,nm,nmP,nisoP
     dimension pm(nm),dens(nm),t(nm),ye(nm),vrad(nm),u(nm)

     pi=dacos(-1.d0)
     uro=0.d0
     ut=0.d0
     uye=0.d0
     uvrad=0.d0

     if(pmass.lt.pm(1))  then
        write(*,*) 'warning low'

        !value2=log10(pm(2))-log10(pm(1))
        !value3=log10(pmass)-log10(pm(1))

        !value1=(log10(dens(2))-log10(dens(1)))/value2
        !uro=10.d0**(log10(dens(1))+value1*value3)

        value2=pm(2)-pm(1)
        value3=pmass-pm(1)

        value1=(dens(2)-dens(1))/value2
        uro=dens(1)+value1*value3

        value1=(t(2)-t(1))/value2
        ut=(t(1)+value1*value3)

        value1=(ye(2)-ye(1))/value2
        uye=(ye(1)+value1*value3)

        value1=(vrad(2)-vrad(1))/value2
        uvrad=(vrad(1)+value1*value3)
        
        value1=(u(2)-u(1))/value2
        uenr=(vrad(1)+value1*value3)
        
        !value1=(s(2)-s(1))/value2
        !us=(s(1)+value1*value3)

     else if(pmass.gt.pm(nm)) then
        write(*,*) 'warning high', pmass, pm(nm)


        value2=pm(nm)-pm(nm-1)
        value3=pmass-pm(nm)

        value1=(dens(nm)-dens(nm-1))/value2
        uro=(dens(nm)+value1*value3)

        value1=(t(nm)-t(nm-1))/value2
        ut=t(nm)+value1*value3

        value1=(ye(nm)-ye(nm-1))/value2
        uye=(ye(nm)+value1*value3)

        value1=(vrad(nm)-vrad(nm-1))/value2
        uvrad=(vrad(nm)+value1*value3)
        
        value1=(u(nm)-u(nm-1))/value2
        uenr=(u(nm)+value1*value3)

        !value1=(s(nm)-s(nm-1))/value2
        !us=(s(nm)+value1*value3)

     else
        j=1

        do while (pmass.ge.pm(j+1))
           j=j+1
        enddo

        value1=(pmass-pm(j))/(pm(j+1)-pm(j))

        uro=dens(j)+(dens(j+1)-dens(j))*value1
        ut=t(j)+(t(j+1)-t(j))*value1
        uye=ye(j)+(ye(j+1)-ye(j))*value1
        uvrad=vrad(j)+(vrad(j+1)-vrad(j))*value1
        uenr=u(j)+(u(j+1)-u(j))*value1
        !us=s(j)+(s(j+1)-s(j))*value1
     endif

     return
   end SUBROUTINE profile
