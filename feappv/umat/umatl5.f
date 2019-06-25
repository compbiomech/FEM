!$Id:$
      subroutine umatl5(f,theta,td,d,ud,hn,h1,nh,ii,istrt, sig,dd,isw)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved
!....  Prepared by Dr. Kewei Li, TU Graz, Jun 25, 2019
!
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Neo Hookean Hyperelastic Model

!     Input:
!          f(*)    -  Current strains at point      (small deformation)
!                  -  Deformation gradient at point (finite deformation)
!          theta   -  Trace of strain at point
!                  -  Determinant of deforamtion gradient
!          td      -  Temperature change
!          d(*)    -  Program material parameters (ndd)
!          ud(*)   -  User material parameters (nud)
!          hn(nh)  -  History terms at point: t_n
!          h1(nh)  -  History terms at point: t_n+1
!          nh      -  Number of history terms
!          ii      -  Current point number
!          istrt   -  Start state: 0 = elastic; 1 = last solution
!          isw     -  Solution option from element

!     Output:
!          sig(*)  -  Stresses at point.
!                     N.B. 1-d models use only sig(1)
!          dd(6,*) -  Current material tangent moduli
!                     N.B. 1-d models use only dd(1,1) and dd(2,1)
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none
      
      include 'iofile.h'
      
      
      integer  nh,istrt,isw, ii, ntm, i, j
      real*8   td
      real*8   f(3,*),theta(*),d(*),ud(*),hn(nh),h1(nh), sig(*),dd(6,*)

      real*8    detf, detfi, j23, trbe3, bei, mub1, mub2, mub3
      real*8    u, up, upp, ha, hp, hpp, press
      real*8    one3, two3
      real*8    be(6)
      real*8    mu, bk

      data      one3 /0.3333333333333333d0/
      data      two3 /0.6666666666666667d0/


      sig(6) = 0.0d0
      
      dd(6,6) = 0.0d0
      
!     Compute deviatoric be
      detf = theta(1)
      detfi = 1.d0/detf
      j23   = detfi**two3      
      
      ! mixed element is 2d, but the storage is 3D 
      ntm = 6 
      
      mu = ud(1) 
      bk = ud(2)
      
!     Compute Left Cauchy-Green deformation tensor
      
      ! if (isw .eq. 4)   write(iow, *), "f(3, 1)=", f(1:3,1:3)

      be(1) = f(1,1)*f(1,1) + f(1,2)*f(1,2) + f(1,3)*f(1,3)
      be(2) = f(2,1)*f(2,1) + f(2,2)*f(2,2) + f(2,3)*f(2,3)
      be(3) = f(3,1)*f(3,1) + f(3,2)*f(3,2) + f(3,3)*f(3,3)
      be(4) = f(1,1)*f(2,1) + f(1,2)*f(2,2) + f(1,3)*f(2,3)
      be(5) = f(2,1)*f(3,1) + f(2,2)*f(3,2) + f(2,3)*f(3,3)
      be(6) = f(1,1)*f(3,1) + f(1,2)*f(3,2) + f(1,3)*f(3,3)
 
      do i = 1,ntm
        be(i) = be(i) * j23
      end do

      trbe3  = (be(1) + be(2) + be(3)) * one3
      be(1) = be(1) - trbe3
      be(2) = be(2) - trbe3
      be(3) = be(3) - trbe3

!     Compute pressure and volumetric moduli

      !u    = bk*( detf**2 - 1.d0 ) - bk*log(abs(detf) )
      !up   = bk*( detf - 1.d0/detf    )
      !upp  = bk*( 1.d0 + 1.d0/detf**2 )

      u    = bk*log(abs(detf))**2*0.5d0
      up   = bk*log(abs(detf)) / detf
      upp  = ( bk/detf  - up )/detf
        
!     Pressure and tangent (not mixed pressure)

      press =  up  
      upp   = upp  * detf

!     Compute Kirchhoff stress tensor.

      mub1 = mu
      do i = 1,ntm
        sig(i) = mub1 * be(i)
      end do

!     Compute tangent tensor
!                                  __             __     _
!     Rank one update: -2/3 mu * ( be x g +  g x  be ) / J

      mub3 = two3 * mub1
      do i = 1,ntm
        bei = mub3 * be(i)
        do j = 1,3
          dd(i,j) =  dd(i,j) - bei
          dd(j,i) =  dd(j,i) - bei
        end do
      end do
!                       __                     _
!     Deviatoric term 2 mu [ I - 1/3 g x g ] / J

      mub1 = mub1 * trbe3
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        dd(i  ,i  ) = dd(i  ,i  ) + mub2
        dd(i+3,i+3) = dd(i+3,i+3) + mub1
        do j = 1,3
          dd(i ,j ) = dd(i ,j )   - mub3
        end do
      end do

!     Compute deviatoric material moduli

      do i = 1,ntm
        do j = 1,i
          dd(i,j) = dd(i,j) * detfi
        end do
      end do

!     Add volumetric correction to aa

      dd(1,1) = dd(1,1) - press + upp
      dd(1,2) = dd(1,2) + press + upp
      dd(1,3) = dd(1,3) + press + upp
      dd(2,2) = dd(2,2) - press + upp
      dd(2,3) = dd(2,3) + press + upp
      dd(3,3) = dd(3,3) - press + upp
      dd(4,4) = dd(4,4) - press
      dd(5,5) = dd(5,5) - press
      dd(6,6) = dd(6,6) - press

!     Compute lower part of aa by symmetry

      do i = 2,ntm
        do j = 1,i-1
          dd(i,j) = dd(j,i)
        end do
      end do

!     Compute Cauchy stress and add pressure terms

      do i = 1,3
        sig(i) = sig(i) * detfi + press
        sig(i+3) = sig(i+3) * detfi   
      end do

!     Compute stored energy density

      !engy = u + (d(22)*trbe3 - d(22))*1.5d0

      
      
      end
