!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2018)                       !
!                                                     !
!=====================================================!
!                                                     !
!               SPHYNX: conservemod.f90               !
!                                                     !
! Checks conservation laws.                           !
!=====================================================!

subroutine conservation

  USE parameters

  IMPLICIT NONE

  INTEGER i,ii

  DOUBLE PRECISION dgrav,ax,ay,az
  DOUBLE PRECISION vmod2,ecin,eint,egrav,etot
  DOUBLE PRECISION acumlmx,acumlmy,acumlmz,acumamx,acumamy,acumamz
  DOUBLE PRECISION linmo,angmo,maxvykin,vykin
  DOUBLE PRECISION aux1,aux2,si,ci,di,sumsi,sumci,sumdi,growthrate


  call profile(0)

  if(flags)write(*,*) 'Conservation Laws'
  if(id.eq.0)open(3,file='conservelaws.d',position='append')
3 format(i10,10(1x,es17.10))

  ecin=0.d0
  eint=0.d0
  egrav=0.d0
  acumlmx=0.d0
  acumlmy=0.d0
  acumlmz=0.d0
  acumamx=0.d0
  acumamy=0.d0
  acumamz=0.d0
  sumsi=0.d0
  sumci=0.d0
  sumdi=0.d0
  maxvykin=0.d0


!$omp parallel private(i,ii,vmod2,ax,ay,az,dgrav,aux1,aux2,si,ci,di,vykin) &
!$omp          reduction(+:acumlmx,acumlmy,acumlmz)&
!$omp          reduction(+:acumamx,acumamy,acumamz)&
!$omp          reduction(+:ecin,egrav,eint)&
!$omp          reduction(+:sumsi,sumci,sumdi)&
!$omp          reduction(max:maxvykin)
!$omp do schedule(static)
  do i=1,n
     ii=1+dim*(i-1)
#ifdef EQMASS
     acumlmx=acumlmx+v(ii)
     acumlmy=acumlmy+v(ii+1)
     acumlmz=acumlmz+v(ii+2)
     ax=a(ii)-despl
     ay=a(ii+1)-despl
     az=a(ii+2)-despl
     acumamx=acumamx+(ay*v(ii+2)-az*v(ii+1))
     acumamy=acumamy+(az*v(ii)-ax*v(ii+2))
     acumamz=acumamz+(ax*v(ii+1)-ay*v(ii))
     vmod2=v(ii)*v(ii)+v(ii+1)*v(ii+1)+v(ii+2)*v(ii+2)
     ecin=ecin+.5d0*vmod2
     dgrav=-g*.5d0*ugrav(i)
     egrav=egrav+dgrav
     eint=eint+u(i)
#else
     acumlmx=acumlmx+v(ii)*mass(i)
     acumlmy=acumlmy+v(ii+1)*mass(i)
     acumlmz=acumlmz+v(ii+2)*mass(i)
     ax=a(ii)-despl
     ay=a(ii+1)-despl
     az=a(ii+2)-despl
     acumamx=acumamx+(ay*v(ii+2)-az*v(ii+1))*mass(i)
     acumamy=acumamy+(az*v(ii)-ax*v(ii+2))*mass(i)
     acumamz=acumamz+(ax*v(ii+1)-ay*v(ii))*mass(i)
     vmod2=v(ii)*v(ii)+v(ii+1)*v(ii+1)+v(ii+2)*v(ii+2)
     ecin=ecin+.5d0*mass(i)*vmod2
     dgrav=-mass(i)*g*.5d0*ugrav(i)
     egrav=egrav+dgrav
     eint=eint+u(i)*mass(i)
#endif
     aux1=exp(-4.d0*pi*abs(ay-0.25d0))
     aux2=exp(-4.d0*pi*abs((ybox-ay)-0.25d0))
     vykin=.5d0*promro(i)*v(ii+1)**2
     maxvykin=max(maxvykin,vykin)
     if(ay.lt.ybox*.5d0) then
       si=v(ii+1)*vol(i)*sin(4.d0*pi*ax)*aux1
       ci=v(ii+1)*vol(i)*cos(4.d0*pi*ax)*aux1
       di=vol(i)*aux1
     else
       si=v(ii+1)*vol(i)*sin(4.d0*pi*ax)*aux2
       ci=v(ii+1)*vol(i)*cos(4.d0*pi*ax)*aux2
       di=vol(i)*aux2
     endif
     sumsi=sumsi+si
     sumci=sumci+ci
     sumdi=sumdi+di
  enddo
!$omp end do
!$omp end parallel
  linmo=dsqrt(acumlmx**2+acumlmy**2+acumlmz**2)
  angmo=dsqrt(acumamx**2+acumamy**2+acumamz**2)

#ifdef EQMASS
  ecin=ecin*masspart
  eint=eint*masspart
  egrav=egrav*masspart
  linmo=linmo*masspart
  angmo=angmo*masspart
#endif

  etot=ecin+eint+egrav
  growthrate=2.d0*sqrt((sumsi/sumdi)**2+(sumci/sumdi)**2)

  if(id.eq.0) then
     write(3,3) l,tt,ecin,eint,egrav,etot,linmo,angmo,growthrate
     close(3)
  endif

  call profile(1)

  return
end subroutine conservation
