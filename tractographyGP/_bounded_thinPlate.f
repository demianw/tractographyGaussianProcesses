      SUBROUTINE thinPlate2D_mat(C,nx,ny,R,cmin,cmax,symm)

cf2py intent(inplace) C
cf2py intent(hide) nx,ny
cf2py double precision intent(in), optional :: R=1
cf2py logical intent(in), optional:: symm=0
cf2py integer intent(in), optional :: cmin=0
cf2py integer intent(in), optional :: cmax=-1
cf2py threadsafe

      INTEGER nx,ny,i,j,cmin,cmax
      DOUBLE PRECISION C(nx,ny)
      LOGICAL symm
      DOUBLE PRECISION R
      
      if (cmax.EQ.-1) then
          cmax = ny
      end if

      if(symm) then

        do j=cmin+1,cmax

          C(j,j)=1.0D0
        
          do i=1,j-1
            C(i,j) = thinPlate2D(C(i,j),R)
          enddo
        enddo

      else

        do j=cmin+1,cmax
          do i=1,nx
            C(i,j) = thinPlate2D(C(i,j),R)
          enddo
        enddo

      endif


      return
      END


      SUBROUTINE thinPlate3D_mat(C,nx,ny,R,cmin,cmax,symm)

cf2py intent(inplace) C
cf2py intent(hide) nx,ny
cf2py double precision intent(in), optional :: R=1
cf2py logical intent(in), optional:: symm=0
cf2py integer intent(in), optional :: cmin=0
cf2py integer intent(in), optional :: cmax=-1
cf2py threadsafe

      INTEGER nx,ny,i,j,cmin,cmax
      DOUBLE PRECISION C(nx,ny)
      LOGICAL symm
      DOUBLE PRECISION R
      
      if (cmax.EQ.-1) then
          cmax = ny
      end if

      if(symm) then

        do j=cmin+1,cmax

          C(j,j)=1.0D0
        
          do i=1,j-1
            C(i,j) = thinPlate3D(C(i,j),R)
          enddo
        enddo

      else

        do j=cmin+1,cmax
          do i=1,nx
            C(i,j) = thinPlate3D(C(i,j),R)
          enddo
        enddo

      endif


      return
      END


c
      DOUBLE PRECISION FUNCTION thinPlate2D(d,R)

      DOUBLE PRECISION d,R
      DOUBLE PRECISION t
cf2py threadsafe

c       Aux variables
      DOUBLE PRECISION d2,R2

      if (d .EQ. 0.0D0) then
        t = 1.0D0
      else
        if (d .LT. R) then
          d2 = d*d
          R2 = R*R
          t = (2*d2*dlog(d)-(1+2*dlog(R))*d2+R2)/R2
        else
          t = 0.0D0 
        endif
      endif

      thinPlate2D=t
      RETURN

      END
c
      FUNCTION thinPlate3D(d,R)

      DOUBLE PRECISION d,R
      DOUBLE PRECISION t
cf2py threadsafe

c       Aux variables
      DOUBLE PRECISION d2,d3,R2,R3

      if (d .EQ. 0.0D0) then
        t = 1.0D0
      else
        if (d .LT. R) then
          d2 = d*d
          d3 = dabs(d)*d2
          R2 = R*R
          R3 = dabs(R)*R2
          t = (2*d3-3*R*d2+R3)/R3
        else
          t = 0.0D0 
        endif
      endif

      thinPlate3D=t
      RETURN

      END





      SUBROUTINE innerProduct_thinPlate3D_normalized(C,nx,ny,Q,R,symm)

cf2py intent(inplace) C
cf2py intent(hide) nx,ny
cf2py double precision intent(in), optional :: R=1
cf2py double precision intent(in), optional :: Q=1
cf2py logical intent(in), optional:: symm=0
cf2py threadsafe

      INTEGER nx,ny,i,j
      DOUBLE PRECISION C(nx,ny)
      DOUBLE PRECISION R,Q
      DOUBLE PRECISION a
      LOGICAL symm
c      DOUBLE PRECISION Rs(10)
c      DOUBLE PRECISION Qs(10)

      if (R .LT. Q) then
        a = Q
        Q = R
        R = a
      end if

c      Rs(1) = R
c      Qs(1) = Q
c      do i=1,9
c        Rs(i+1) = Rs(i) * R
c        Qs(i+1) = Qs(i) * Q
c      enddo

      if (symm) then
        do j=1,ny
          do i=1,j
             C(i,j) = covIntegralThinPlateR3Normalized( C(i,j), Q, R)
             C(j,i) = C(i,j)
          enddo
        enddo
      else 
        do j=1,ny
          do i=1,nx
             C(i,j) = covIntegralThinPlateR3Normalized( C(i,j), Q, R)
          enddo
        enddo
      endif

      return
      END
 
      DOUBLE PRECISION FUNCTION covIntegralThinPlateR3Normalized(
     &                                                    w,Q,R)

      double precision w,Q,R


      double precision aux

      if ( Q .GT. R ) then
        aux = R
        R = Q
        Q = aux
      endif

      if ( w .EQ. 0 ) then
       covIntegralThinPlateR3Normalized = 
     & dacos( -1.0D0 ) * dble(Q ** 3) * dble(84 * R ** 3 - 81 *
     & R * Q ** 2 + 35 * Q ** 3) / dble(R ** 3) / 0.315D3 
       return
      elseif ( w .LE. (R-Q) ) then
        covIntegralThinPlateR3Normalized = 
     &          covIntTPR3NwLTRminusQ(w,Q,R)
      elseif ( ( (R-Q) .LT. w ) .and. ( w .LE. R) ) then
        covIntegralThinPlateR3Normalized =
     &          covIntTP3NRminusQLTw(w,Q,R)
      elseif ( w .GT. R ) then
        covIntegralThinPlateR3Normalized =
     &          covIntTP3NwGTR(w,Q,R)
      endif

      return

      end


      DOUBLE PRECISION function covIntTPR3NwLTRminusQ (w, Q, R)
        double precision w
        double precision Q
        double precision R

        if ( w .LE. Q ) then
          covIntTPR3NwLTRminusQ = dacos( -1.0D0 ) * dble(-
     &405 * Q ** 8 * R - 1260 * Q ** 6 * R * w ** 2 + 420 * Q ** 6 * R *
     &* 3 + 15 * Q * w ** 8 + 175 * Q ** 9 + 900 * Q ** 7 * w ** 2 + 378
     & * Q ** 5 * w ** 4 - 60 * Q ** 3 * w ** 6 - 4 * w ** 9) / dble(Q *
     &* 3) / dble(R ** 3) / 0.1575D4
        else 
          covIntTPR3NwLTRminusQ =  dacos( -1.0D0 ) * dble(Q 
     &** 3) * dble(-135 * Q ** 2 * w * R - 420 * w ** 3 * R + 140 * w * 
     &R ** 3 + 180 * Q ** 2 * w ** 2 + 280 * w ** 4 + 8 * Q ** 4) / dble
     &(w) / dble(R ** 3) / 0.525D3
        endif
        return 
      end

      double precision function covIntTP3NRminusQLTw (w, Q, R)
        double precision w
        double precision Q
        double precision R
        double precision pi
        pi = dacos( -1.0D0 )

        if ( w .LE. Q ) then 
         covIntTP3NRminusQLTw = -0.40D1 / 0.525D3 * pi *
     & (w ** 10 / 0.2D1 + (-0.5D1 / 0.8D1 * R - 0.5D1 / 0.8D1 * Q) * w *
     &* 9 - 0.135D3 / 0.64D2 * w ** 8 * R * Q + (0.5D1 / 0.2D1 * R ** 3 
     &+ 0.5D1 / 0.2D1 * Q ** 3) * w ** 7 + 0.105D3 / 0.16D2 * Q * R * (R
     & ** 2 + Q ** 2) * w ** 6 + (-0.63D2 / 0.4D1 * Q ** 5 - 0.63D2 / 0.
     &4D1 * R ** 5) * w ** 5 + (-0.945D3 / 0.32D2 * Q * R ** 5 - 0.175D3
     & / 0.16D2 * Q ** 3 * R ** 3 - 0.945D3 / 0.32D2 * Q ** 5 * R + 0.35
     &D2 * R ** 6 + 0.35D2 * Q ** 6) * w ** 4 + (0.105D3 / 0.2D1 * R ** 
     &6 * Q - 0.75D2 / 0.2D1 * R ** 7 + 0.105D3 / 0.2D1 * Q ** 6 * R - 0
     &.75D2 / 0.2D1 * Q ** 7) * w ** 3 + 0.45D2 / 0.2D1 * ( Q - R )** 4
     &* (Q ** 4 + 0.17D2 / 0.8D1 * Q ** 3 * R + 0.5D1 / 0.2D1 * Q ** 2 *
     & R ** 2 + 0.17D2 / 0.8D1 * Q * R ** 3 + R ** 4) * w ** 2 + (0.135D
     &3 / 0.8D1 * R ** 8 * Q - 0.175D3 / 0.24D2 * Q ** 9 - 0.35D2 /0.2D1
     & * Q ** 6 * R ** 3 - 0.175D3 / 0.24D2 * R ** 9 - 0.35D2 / 0.2D1 *
     &R ** 6 * Q ** 3 + 0.135D3 / 0.8D1 * Q ** 8 * R)*w + ( Q - R ) **
     &6 * (Q ** 4 + 0.209D3 / 0.64D2 * Q ** 3 * R + 0.147D3 / 0.32D2 * 
     &Q ** 2 * R ** 2 + 0.209D3 / 0.64D2 * Q * R ** 3 + R ** 4)) / Q ** 
     &3 / w / R ** 3 
      return

      else

      covIntTP3NRminusQLTw =  pi * (0.192D3 * Q ** 10 - 0.
     &192D3 * R ** 10 - 0.32D2 * w ** 10 + 0.8100D4 * w ** 2 * Q ** 7 * 
     &R - 0.10080D5 * w ** 3 * Q ** 6 * R + 0.8100D4 * Q * w ** 2 * R **
     & 7 - 0.3780D4 * w ** 2 * Q ** 5 * R ** 3 + 0.3360D4 * R ** 6 * Q *
     &* 3 * w - 0.3240D4 * w * Q ** 8 * R - 0.900D3 * Q ** 7 * R ** 3 + 
     &0.480D3 * w ** 7 * Q ** 3 + 0.405D3 * w ** 8 * R * Q - 0.3780D4 * 
     &Q ** 3 * w ** 2 * R ** 5 - 0.1260D4 * w ** 6 * Q ** 3 * R + 0.5670
     &D4 * Q * w ** 4 * R ** 5 + 0.5670D4 * w ** 4 * Q ** 5 * R - 0.3240
     &D4 * R ** 8 * Q * w + 0.3360D4 * w * Q ** 6 * R ** 3 + 0.2100D4 * 
     &w ** 4 * Q ** 3 * R ** 3 - 0.1260D4 * w ** 6 * Q * R ** 3 - 0.1008
     &0D5 * w ** 3 * R ** 6 * Q + 0.525D3 * Q ** 9 * R - 0.480D3 * w ** 
     &7 * R ** 3 + 0.120D3 * w ** 9 * R - 0.120D3 * w ** 9 * Q - 0.1400D
     &4 * w * Q ** 9 + 0.4320D4 * w ** 2 * Q ** 8 - 0.7200D4 * w ** 3 * 
     &Q ** 7 + 0.6720D4 * w ** 4 * Q ** 6 + 0.525D3 * Q * R ** 9 + 0.113
     &4D4 * Q ** 5 * R ** 5 - 0.900D3 * Q ** 3 * R ** 7 + 0.7200D4 * w *
     &* 3 * R ** 7 + 0.1400D4 * w * R ** 9 - 0.6720D4 * w ** 4 * R ** 6 
     &+ 0.3024D4 * w ** 5 * R ** 5 - 0.4320D4 * w ** 2 * R ** 8 - 0.3024
     &D4 * w ** 5 * Q ** 5) / Q ** 3 / w / R ** 3 / 0.25200D5
      return
      endif

      end

      double precision function covIntTP3NwGTR (w, Q, R)
        double precision w
        double precision Q
        double precision R

        if ( w .LE. R+Q ) then
          covIntTP3NwGTR = 0.4D1 / 0.525D3 * dacos( -1.0D0 ) 
     & * (R + Q - w) ** 6 * (w ** 4 / 0.6D1 + (0.3D1 / 0.8D1 * Q + 0.3D
     &1 / 0.8D1 * R) * w ** 3 + (0.103D3 / 0.64D2 * Q * R - Q ** 2 / 0.4
     &D1 - R ** 2 / 0.4D1) * w ** 2 - 0.31D2 / 0.24D2 * (R ** 2 - 0.247D
     &3 / 0.124D3 * Q * R + Q ** 2) * (R + Q) * w + R ** 4 - 0.209D3 / 0
     &.64D2 * Q ** 3 * R + 0.147D3 / 0.32D2 * Q ** 2 * R ** 2 + Q ** 4 -
     & 0.209D3 / 0.64D2 * Q * R ** 3) / Q ** 3 / w / R ** 3
        else
          covIntTP3NwGTR = 0
        endif
        return
      end

