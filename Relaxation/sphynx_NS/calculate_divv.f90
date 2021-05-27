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
!               SPHYNX: calculate_divv.f90            !
!                                                     !
! Calculates divergence(velocity) with IAD.           !
!=====================================================!

  SUBROUTINE calculate_divv

    USE parameters

    IMPLICIT NONE

    INTEGER i,ii,k,j,jj,iii

    DOUBLE PRECISION d1,d2,d3,d05,d02,jumpx,jumpy,jumpz
    DOUBLE PRECISION v1,w1d,norm,dji1,dji2,dji3
    DOUBLE PRECISION vji1,vji2,vji3,dvx,dvy,dvz
    DOUBLE PRECISION kern11i,kern12i,kern13i
    DOUBLE PRECISION kern21i,kern22i,kern23i
    DOUBLE PRECISION kern31i,kern32i,kern33i
    DOUBLE PRECISION termA1,termA2,termA3

    DOUBLE PRECISION,DIMENSION(n)     :: curlvx,curlvy,curlvz
    
    call profile(0)
    if(flags)write(*,*) 'Div v + Curl v calculation (IAD0)'


    !Initialization
    divv=0.d0
    curlv=0.d0
    curlvx=0.d0
    curlvy=0.d0
    curlvz=0.d0

    !$omp parallel private(i,ii,iii,k,j,jj,d1,d2,d3,d02,d05,v1,w1d,&
    !$omp                  norm,kern11i,kern12i,kern13i,kern21i,kern22i,&
    !$omp                  kern23i,kern31i,kern32i,kern33i,dji1,dji2,dji3,&
    !$omp                  dvx,dvy,dvz,jumpx,jumpy,jumpz,vji1,vji2,vji3,&
    !$omp                  termA1,termA2,termA3)

    !Reshaping


    !$omp do schedule(static)
    do i=npini,npend
       ii=1+dim*(i-1)
       iii=i-npini+1
       norm=pk(i)*hm3(i)

       do k=1,nvi(i)-1
          j=neighbors(iii,k)
          jj=1+dim*(j-1)

          call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

          d1=a(ii)-a(jj)-jumpx
          d2=a(ii+1)-a(jj+1)-jumpy
          d3=a(ii+2)-a(jj+2)-jumpz
          d02=d1*d1+d2*d2+d3*d3
          d05=sqrt(d02)

          v1=d05/h(i)

          call Wkernel_noderiv(v1,indice(i),w1d)

         !  dji1=aloc(jj)+jumpx-aloc(ii)
         !  dji2=aloc(jj+1)+jumpy-aloc(ii+1)
         !  dji3=aloc(jj+2)+jumpz-aloc(ii+2)
          dji1=-d1
          dji2=-d2
          dji3=-d3

          vji1=v(jj)-v(ii)
          vji2=v(jj+1)-v(ii+1)
          vji3=v(jj+2)-v(ii+2)

          kern11i=c11(i)*dji1
          kern12i=c12(i)*dji2
          kern13i=c13(i)*dji3

          kern21i=c12(i)*dji1
          kern22i=c22(i)*dji2
          kern23i=c23(i)*dji3

          kern31i=c13(i)*dji1
          kern32i=c23(i)*dji2
          kern33i=c33(i)*dji3

          termA1=(kern11i+kern12i+kern13i)*w1d
          termA2=(kern21i+kern22i+kern23i)*w1d
          termA3=(kern31i+kern32i+kern33i)*w1d

          dvx=vji1*termA1
          dvy=vji2*termA2
          dvz=vji3*termA3

          divv(i)=divv(i)+(dvx+dvy+dvz)*xmass(j)

          curlvx(i)=curlvx(i)+(vji3*termA2-vji2*termA3)*xmass(j)
          curlvy(i)=curlvy(i)-(vji3*termA1-vji1*termA3)*xmass(j)
          curlvz(i)=curlvz(i)+(vji2*termA1-vji1*termA2)*xmass(j)
       enddo

       divv(i)=divv(i)/sumkx(i)*norm
       curlv(i)=sqrt(curlvx(i)**2+curlvy(i)**2+curlvz(i)**2)
       curlv(i)=abs(curlv(i))/sumkx(i)*norm

    enddo
    !$omp end do
    !$omp end parallel

    call profile(1)

    RETURN
  END SUBROUTINE calculate_divv
