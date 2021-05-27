SUBROUTINE sort_neighbors(auxneighbors,auxcube)

  USE parameters,only: neighbors,cube,npini,npend,nvmax,nvi,id

  INTEGER i,iii,k,j,l
  INTEGER,PARAMETER::invmax=int(nvmax)
  INTEGER,DIMENSION(npend-npini+1,invmax),intent(in)   :: auxneighbors,auxcube
  INTEGER,DIMENSION(invmax) :: sorted_index

  !$omp parallel private(i,iii,k,sorted_index,j,l)
  !$omp do schedule(static)
  do i=npini,npend
     iii=i-npini+1
     call indexxi(invmax,auxneighbors(iii,:),sorted_index)
     k=1
     do while (auxneighbors(iii,sorted_index(k)).eq.0)
        k=k+1
     enddo
     if(invmax-k+1.ne.nvi(i)-1)print *,i,invmax-k+1,nvi(i),k,invmax
!     if(iii.eq.2584.and.id.eq.12) then
!       print *, invmax-k+1,nvi(i)-1
!            print *, auxneighbors(iii,:)
!            print *, '****************'
!            print *,sorted_index(:)
!            do j=k,invmax
!              print *, auxneighbors(iii,sorted_index(j))
!            enddo
!            stop
!      endif
     l=0
     do j=k,invmax
       l=l+1
       neighbors(iii,l)=auxneighbors(iii,sorted_index(j))
       cube(iii,l)=auxcube(iii,sorted_index(j))
     enddo
  enddo
  !$omp end do
  !$omp end parallel

END SUBROUTINE sort_neighbors
