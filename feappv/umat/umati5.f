!$Id:$
      subroutine umati5(vtype,vv, d, ud, n1,n3)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Dummy user material model routine

!      Inputs:
!         vtype  - Name of material model
!         vv(5)  - Command line real data
!         d(*)   - Program material parameter data

!      Outputs:
!         ud(*)  - Material parameter data for model
!         n1     - Number of history items/point (time   dependent)
!         n3     - Number of history items/point (time independent)
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      
      include 'iofile.h' 

      logical   pcomp
      character vtype*15
      integer   n1,n3, i
      real*8    vv(5),d(*),ud(*)

!     Set command name

      if(pcomp(vtype,'mat5',4)) then     ! Default  form DO NOT CHANGE
       vtype = 'hook'                     ! Specify new 'name'

c     Input user data and save in ud(*) array

      else                              ! Perform input for user data

        write(iow,100) 'M e c h a n i c a l   P r o p e r t i e s'
        write(iow,*)
        write(iow,100) 'Neo Hookean Hyperelastic Model loaded'
        write(iow,*)

c       Read and write parameters
        do i = 1,2
          ud(i) = vv(i)
        end do

c       Print user parameter information        
        write(iow,200) 'Parameter mu           :', ud(1)
        write(iow,200) 'Parameter k            :', ud(2)

        write(iow,*)
        
        
      endif

100   format(1x,a)
200   format(1x,a,1x,F16.5)


      end
