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
!             SPHYNX: calculate_density.f90           !
!                                                     !
! Calculates density and omega                        !
!=====================================================!

  SUBROUTINE calculate_density

    USE parameters

    IMPLICIT NONE
  
    INTEGER i,j,k,ii,jj,iii,ik1

    DOUBLE PRECISION w1d,dw1,w2d,dw2,d05,v1,v2,d1,d2,d3
    DOUBLE PRECISION dter,dterh,dlw1d,dtern
    DOUBLE PRECISION dmy,jumpx,jumpy,jumpz

    DOUBLE PRECISION,DIMENSION(n3):: aloc,vloc

    if(flags)write(*,*) 'Density calculation'


    !Inicialization
    f=0.d0           !Gravitational force
    do i=npini,npend
       if(ready(i))cycle
       sumkx(i)=0.d0       !Volumen partition
       sumwh(i)=0.d0       !grad-h term
    enddo

    !Reshaping
!$omp parallel private(i,ii,k,j,jj,iii,d1,d2,d3,d05,v1,w1d,dw1,&
!$omp                  dter,dterh,dlw1d,dtern,jumpx,jumpy,jumpz)

!$omp do schedule(static)
    do i=1,n
       ii=1+dim*(i-1)
       aloc(ii)=a(i)
       aloc(ii+1)=a(i+n)
       aloc(ii+2)=a(i+n2)
       vloc(ii)=v(i)
       vloc(ii+1)=v(i+n)
       vloc(ii+2)=v(i+n2)
    enddo
!omp end do

!$omp do schedule(static)
    do i = npini,npend
       if(ready(i))cycle
       ii=1+dim*(i-1)
       iii=i-npini+1
       
       do k=1,nvi(i)
          j=neighbors(iii,k)
          jj=1+dim*(j-1)

          call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

          if(i.ne.j)then
             d1=aloc(ii)-aloc(jj)-jumpx
             d2=aloc(ii+1)-aloc(jj+1)-jumpy
             d3=aloc(ii+2)-aloc(jj+2)-jumpz
             d05=sqrt(d1*d1+d2*d2+d3*d3)
             v1=d05/h(i)

             call Wkernel(v1,indice(i),w1d,dw1)

             dter=pk(i)*hm3(i)*w1d
             dterh=-pk(i)*hm4(i)*(dim*w1d+v1*v1*dw1)   !additional v1 for dw1

             sumkx(i)=sumkx(i)+xmass(j)*dter
             sumwh(i)=sumwh(i)+xmass(j)*dterh

          endif
       enddo
          
       !adds self-contribution to density and omega   
       sumkx(i)=sumkx(i)+pk(i)*xmass(i)*hm3(i)
       sumwh(i)=sumwh(i)-dim*pk(i)*hm4(i)*xmass(i)
       
    enddo
!$omp end do

!$omp do schedule(static)
    do i=npini,npend
       if(ready(i))cycle
       promro(i)=sumkx(i)*masa(i)/xmass(i)
       sumwh(i)=sumwh(i)*masa(i)/xmass(i)
    enddo
!$omp end do

!$omp end parallel   

    RETURN

  END SUBROUTINE calculate_density
       
