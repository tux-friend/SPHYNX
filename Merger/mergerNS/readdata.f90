!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2017)                       !
!                                                     !
!=====================================================!
!                                                     !
!                 SPHYNX: readdata.f90                !
!                                                     !
! Input of data.                                      !
!=====================================================!

    SUBROUTINE readdata

      USE parameters

      IMPLICIT NONE

      DOUBLE PRECISION dmy,dnvi                  !dummy real variables
      INTEGER idmy                               !dummy integer variable
      INTEGER i,ii,j,k,ierr
      LOGICAL endhere
      DOUBLE PRECISION, DIMENSION(n3):: aloc,vloc

      if(flags) write(*,*) 'Reading data...'

!--------------------  User can change this  --------------------
      if(iterini.eq.1) then
         open(1,file='../initialmodels/'//inputfile)
         masscloudinic=0.d0
         mark(:)=.false.
         n=nmax
         do i=1,n1
#ifdef EQMASS 

            read(1,formatout)  idmy,a(i),a(i+n),a(i+n2),h(i),dmy,&
                &promro(i),dmy,dmy,dmy,radius(i),dmy,dmy,&
                &dmy,dmy,dmy,dmy,dmy,dmy,dmy,&
                &dmy,temp(i),dmy,dmy,&
                &ballmass(i),xmass(i),dmy,dnvi,dmy,&
                &dmy,dmy,dmy,dmy,dmy,&
                &dmy,dmy
            nvi(i)=int(dnvi)
#else
            read(1,'(2x,9(2x,es17.10))') a(i),a(i+n),a(i+n2),h(i),&
                &promro(i),mass(i),v(i),v(i+n),v(i+n2)
            print *,a(i)
#endif
            !print *,i,promro(i),id
            temp(i)=1.d0
            u(i)=1.d50
            ballmass(i)=promro(i)*h(i)**3
         enddo
#ifdef EQMASS
         masspart=xmass(1)
         print *,masspart
         if(id.eq.0) then
           open(22,file='REPORT',position='append')
           write(22,'(1x,a20,5x,es23.16)') 'masspart:',masspart
           close(22)
         endif
#endif

         !populate arrays with same values as star 1 (doubling the star)
         do i=1,n1
            j=i+n1
            a(j)=a(i)
            a(j+n)=a(i+n)
            a(j+n2)=a(i+n2)
            h(j)=h(i)
            promro(j)=promro(i)
            ballmass(j)=ballmass(i)
            xmass(j)=xmass(i)
            u(j)=u(i)
            nvi(j)=nvi(i)
            temp(j)=temp(i)
         enddo
         !Reshaping
         !$omp parallel private(i,ii)
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
         !$omp end do
         !$omp end parallel
         a(:)=aloc(:)
         v(:)=vloc(:)

         if(ncubes.gt.1) then   !If PBC do not center in CM
            cm(1)=xbox*.5d0     !Use fake central point of the bo
            cm(2)=ybox*.5d0
            cm(3)=zbox*.5d0
            do i=1,n
               ii=1+dim*(i-1)
               radius(i)=sqrt((a(ii)-cm(1))**2+(a(ii+1)-cm(2))**2+&
                    &(a(ii+2)-cm(3))**2)
            enddo
         else
            call masscenter(1,n)
            print *,'masscenter:', cm(:)
            do i=1,n
               ii=1+dim*(i-1)
               a(ii)=a(ii)-cm(1)
               a(ii+1)=a(ii+1)-cm(2)
               a(ii+2)=a(ii+2)-cm(3)
               radius(i)=sqrt(a(ii)*a(ii)+a(ii+1)*a(ii+1)+a(ii+2)*a(ii+2))
            enddo
         endif
         print *,'radius:',maxval(radius)

! RESTART -------------------------------------------------------
      else

        open(1,file='../initialmodels/'//inputfile)
        do i=1,n1

          read(1,formatout)  idmy,a(i),a(i+n),a(i+n2),h(i),u(i),promro(i),&
              &    v(i),v(i+n),v(i+n2),radius(i),dmy,dmy,&
              &    dmy,dmy,dmy,dmy,dmy,dmy,dmy,&
              &    dmy,temp(i),dmy,dmy,&
              &    ballmass(i),xmass(i),dmy,dnvi,dmy,&
              &    dmy,dmy,dmy,dmy,dmy,&
              &    dmy,dmy
          nvi(i)=int(dnvi)


#ifdef EQMASS
            masspart=xmass(1)
#else
            !mass(i)=1.3682937d27
#endif
         enddo
         
        open(2,file='../initialmodels/'//inputfile2)
        do i=n1+1,nmax

          read(2,formatout)  idmy,a(i),a(i+n),a(i+n2),h(i),u(i),promro(i),&
              &    v(i),v(i+n),v(i+n2),radius(i),dmy,dmy,&
              &    dmy,dmy,dmy,dmy,dmy,dmy,dmy,&
              &    dmy,temp(i),dmy,dmy,&
              &    ballmass(i),xmass(i),dmy,dnvi,dmy,&
              &    dmy,dmy,dmy,dmy,dmy,&
              &    dmy,dmy
          nvi(i)=int(dnvi)


#ifdef EQMASS
            masspart=xmass(1)
#else
            !mass(i)=1.3682937d27
#endif
         enddo

         !Reshaping
         !$omp parallel private(i,ii)
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
         !$omp end do
         !$omp end parallel
         a(:)=aloc(:)
         v(:)=vloc(:)

         call masscenter(1,n)
         masscloudinic=6.85189d-4/9.99805d-1
      endif


!----------------------------------------------------------------

      close(1)
      close(2)

      !Tests validity of initial data
      endhere=.false.
      do i=1,n*dim
         if(isnan(a(i))) then
            write(*,*) 'NAN position at: (',i,',',i/n,')',a(i)
            endhere=.true.
         else if(isnan(v(i))) then
            write(*,*) 'NAN velocity at: (',i,',',i/n,')',v(i)
            endhere=.true.
         endif
         if(endhere)stop
      enddo

#ifdef EQMASS
      if(isnan(masspart).or.masspart.le.0.d0) then
          write(*,*) 'Incorrect masspart',i,masspart
          endhere=.true.
      endif
#endif

      do i=1,n
         if(isnan(u(i))) then
            write(*,*) 'Incorrect U at i=',i,u(i)
            endhere=.true.
         endif
#ifdef EQMASS
#else
         if(isnan(mass(i)).or.mass(i).le.0.d0) then
            write(*,*) 'Incorrect mass at i=',i,mass(i)
            endhere=.true.
         endif
#endif
         if(isnan(h(i)).or.h(i).le.0.d0) then
            write(*,*) 'Incorrect smoothing-length at i=',i,h(i)
            endhere=.true.
         endif
         if(endhere)stop
      enddo

      if(flags) write(*,*) 'Reading data... Done!'


      RETURN
    END SUBROUTINE readdata
