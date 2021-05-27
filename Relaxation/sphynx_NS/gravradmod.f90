    SUBROUTINE gravrad(theta,phi)

      USE parameters,only: l,tt,g,cvel,TD,masspart,n1,nmax

      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION,INTENT(in)::theta,phi
      DOUBLE PRECISION,parameter::gwunits=g/cvel**4/3.08568025d22   !GW units at 10kpc
      DOUBLE PRECISION sqsint,sqsinp,sqcost,sqcosp,sin2t,sin2p,cos2t,cos2p,&
           &  sint,sinp,cost,cosp,dot2ibartt,dot2ibarpp,dot2ibartp,ixx,iyy,&
           &  izz,ixy,ixz,iyz,trazai,httplus,httcross,units,ixx3,iyy3,&
           &  izz3,ixy3,ixz3,iyz3

      call profile(0)

      call d2quadmom(1,1,ixx)
      call d2quadmom(2,2,iyy)
      call d2quadmom(3,3,izz)
      call d2quadmom(1,2,ixy)
      call d2quadmom(1,3,ixz)
      call d2quadmom(2,3,iyz)


      if(theta.eq.0.d0.and.phi.eq.0.d0) then
         dot2ibartt=ixx
         dot2ibarpp=iyy
         dot2ibartp=ixy
      else
         sin2t=sin(2.d0*theta)
         sin2p=sin(2.d0*phi)
         cos2t=cos(2.d0*theta)
         cos2p=cos(2.d0*phi)
         sint=sin(theta)
         sinp=sin(phi)
         cost=cos(theta)
         cosp=cos(phi)
         sqsint=sint*sint
         sqsinp=sinp*sinp
         sqcost=cost*cost
         sqcosp=cosp*cosp

         dot2ibartt=(ixx*sqcosp+iyy*sqsinp+ixy*sin2p)*sqcost+izz*sqsint&
              &-(ixz*cosp+iyz*sinp)*sin2t
         dot2ibarpp=ixx*sqsinp+iyy*sqcosp-ixy*sin2p
         dot2ibartp=.5d0*(iyy-ixx)*cost*sin2p+ixy*cost*cos2p+&
              &(ixz*sinp-iyz*cosp)*sint
      endif

      httplus=(dot2ibartt-dot2ibarpp)*gwunits
      httcross=2.d0*dot2ibartp*gwunits
      
      open(1,file='radgrav.d',access='append')
      write(1,18) l,tt,httplus,httcross,ixx,iyy,ixy,ixz,iyz,izz
      close(1)
 18   format(1x,i6,16(1x,es14.7))

      call profile(1)

      return
    end SUBROUTINE gravrad

! ***********************************************************************

      SUBROUTINE d2quadmom(i1,i2,dot2i)

      USE parameters,only:nmax,a,v,ac,third,masspart,cm,despl,dim

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(out)::dot2i
      DOUBLE PRECISION vl,vm,xl,xm,al,am,factor1,factor2
      DOUBLE PRECISION,DIMENSION(3*nmax)::pos
      INTEGER,INTENT(in)::i1,i2
      INTEGER i,k1,k2,n,n2,ii

      n=nmax
      n2=2*n
      dot2i=0.d0

      k1=i1-1
      k2=i2-1

      do i=1,n
         ii=1+dim*(i-1)
         pos(ii)=a(ii)-despl!-cm(1)
         pos(ii+1)=a(ii+1)-despl!-cm(2)
         pos(ii+2)=a(ii+2)-despl!-cm(3)
      enddo

      if(i1.ne.i2) then
         do i=1,n
            ii=1+dim*(i-1)
            vl=v(ii+k1)
            vm=v(ii+k2)
            xl=pos(ii+k1)
            xm=pos(ii+k2)
            al=ac(ii+k1)
            am=ac(ii+k2)
            dot2i=dot2i+(2.d0*vl*vm+al*xm+xl*am)
         enddo
         dot2i=dot2i*masspart
      else
         do i=1,n
            ii=1+dim*(i-1)
            factor1=v(ii)*v(ii)+v(ii+1)*v(ii+1)+v(ii+2)*v(ii+2)
            factor2=pos(ii)*ac(ii)+pos(ii+1)*ac(ii+1)+pos(ii+2)*ac(ii+2)
            vl=v(ii+k1)
            xl=pos(ii+k1)
            al=ac(ii+k1)
            dot2i=dot2i+3.d0*(vl*vl+xl*al)-factor1-factor2
         enddo
         dot2i=dot2i*2.d0*third*masspart
      endif

      return
    end SUBROUTINE d2quadmom
