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
!              SPHYNX: findneighbors.f90              !
!                                                     !
! Creates a list of neighbours for each particle.     !
!=====================================================!

  SUBROUTINE findneighbors

    USE mpi
    USE parameters

    IMPLICIT NONE

    INTEGER nodo,count,i,j,ii,iii,jj,ik1,ierr
    INTEGER,DIMENSION(npend-npini+1,int(nvmax))::auxneighbors,auxcube
    DOUBLE PRECISION axce,ayce,azce,d1,d2,d3,dc,d05
    DOUBLE PRECISION haux
    DOUBLE PRECISION v1,v2,nij,dmy,ah,jumpx,jumpy,jumpz
    DOUBLE PRECISION, DIMENSION(nt)::xcet,ycet,zcet
    DOUBLE PRECISION, DIMENSION(n)::aux,aux1
    LOGICAL recalc,jump
    LOGICAL, DIMENSION(nmax)::find_neighbors

    call profile(0)

!neighbors(i,j) has a disordered list of neighbors of part. i, including i itself!
    mindist(:)=1.d80
    recalc=.true.
    find_neighbors=.true.
    count=0
    auxneighbors(:,:)=0
    auxcube(:,:)=0
    nvi(:)=0
    nvi(npini:npend)=1     !Self-particle

    do while(recalc)
       recalc=.false.
       count=count+1
       if(flags)write(*,*) 'Finding neighbors',count,id

       do ik1=1,ncubes

          call apply_PBC(0,0,ik1,jumpx,jumpy,jumpz)

          xcet(1:ntotal)=xce(1:ntotal)+jumpx
          ycet(1:ntotal)=yce(1:ntotal)+jumpy
          zcet(1:ntotal)=zce(1:ntotal)+jumpz

          !$omp parallel private(i,ii,iii,nodo,axce,ayce,azce,d1,d2,d3,dc,j,&
          !$omp                  jj,d05,v1,v2)
          !$omp do schedule(static)
global:   do i = npini,npend
             iii=i-npini+1
             if(find_neighbors(i)) then
                ii=1+dim*(i-1)
                nodo=0
51              nodo=nodo+1
                axce=a(ii)-xcet(nodo)
                ayce=a(ii+1)-ycet(nodo)
                azce=a(ii+2)-zcet(nodo)
                d1=dabs(axce)
                d2=dabs(ayce)
                d3=dabs(azce)
                dc=4.d0*h(i)+dxx(nodo)/2.d0
                if((d1.le.dc).and.(d2.le.dc).and.(d3.le.dc)) goto 55
52              if (nodo .eq. nh2(npa(nodo))) then
                   nodo=npa(nodo)
                   if(nodo .eq. 0) then
                      if(nvi(i).le.int(nvmax).and.nvi(i).ge.int(nvmin).and.ik1.eq.ncubes)find_neighbors(i)=.false.
                      cycle global
                   endif
                   goto 52
                else
                   goto 51
                endif

!     If this condition is fulfilled means that it has found a particle.
!     If not, there is no need to continue with the operations.


55              if (dxx(nodo) .eq. 0.d0) then
                   j=nnod(nodo)
                   if(j.eq.i) goto 52  !Exclude self-particle for list of neighbors
                   d05=sqrt(d1*d1+d2*d2+d3*d3)
                   v1=d05/h(i)
                   v2=d05/h(j)

!                   if(v1.gt.2.d0) then
                   if(v1.gt.2.d0 .and. v2.gt.2.d0) then
                      goto 52
                   else
                      if(nvi(i).lt.int(nvmax)) then
!It is mandatory that the order of update is as follows
!auxneighbors and auxcube update for nvi(i)-1 because nvi(i) is
!initialized to 1 as it includes i, but auxneighbors and auxcube should .not
!Therefore, their lenght is nvi(i)-1
!To avoid that substraction every iteration, nvi(i) is updated only after
!auxneighbors and auxcube.
                         auxneighbors(iii,nvi(i))=j
                         auxcube(iii,nvi(i))=ik1
                         nvi(i)=nvi(i)+1
                      else
                         nvi(i)=nvi(i)+1
                      endif
                   endif
                   goto 52
                endif
                nodo=nh1(nodo)-1
                goto 51
             endif
          enddo global
!$omp end do
!$omp end parallel
       enddo

       recalc=any(find_neighbors(npini:npend))

       if(count.ge.10.and.recalc)recalc=.false.     !Do not iterate forever
       write(*,*) 'recalc',recalc
       if(.not.recalc) then
  check: do i=npini,npend     !We just need one particle out of limits to recalc
             if(nvi(i).gt.int(nvmax).or.nvi(i).lt.int(nvmin)) then
                recalc=.true.
                exit check
             endif
          enddo check
       endif

       if(recalc) then
!$omp parallel private(i,iii,haux)
!$omp do schedule(static)
          do i=npini,npend
             iii=i-npini+1
             if(find_neighbors(i))then
                print '(i7,1x,es12.5,1x,i4,1x,i3)',i,h(i),nvi(i),count
                haux=h(i)*.5d0*(1.d0+hfac*ni0/dble(nvi(i)))**hexp
                if(haux.gt.h(i)) then
                   h(i)=min(1.1d0*h(i),haux)
                else
                   h(i)=max(0.9d0*h(i),haux)
                endif
                nvi(i)=1
                auxneighbors(iii,:)=0
                auxcube(iii,:)=0
             endif
          enddo
!$omp end do
!$omp end parallel
       endif
    enddo

    aux1=0.d0;aux1(npini:npend)=h(npini:npend);aux=0.d0
    call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);h=aux
    call calculate_hpowers

    neighbors=0
    cube=0

    if(sortneigh) then
      call sort_neighbors(auxneighbors,auxcube)
    else
      !$omp parallel private(i,iii,j)
      !$omp do schedule(static)
      do i=npini,npend
        iii=i-npini+1
        do j=1,nvi(i)
            neighbors(iii,j)=auxneighbors(iii,j)
            cube(iii,j)=auxcube(iii,j)
        enddo
      enddo
      !$omp end do
      !$omp end parallel
    endif

    if(l.le.liniNR) then
       ballmass=0.d0
!$omp parallel private(i)
!$omp do schedule(static)
       do i=npini,npend
          ballmass(i)=promro(i)*h3(i)
       enddo
!$omp end do
!$omp end parallel
    endif
    call profile(1)
    RETURN
  END SUBROUTINE findneighbors
