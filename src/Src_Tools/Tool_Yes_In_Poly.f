!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      subroutine Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_INOUT)
C     Determine whether a point is inside a polygon. Note that the polygon must be closed.
c     npol equals the number of edges plus one or the number of vertices plus one
      use Global_Float_Type
      implicit none
      integer i,j,npol
      real(kind=FT) x,y,xpol(npol),ypol(npol)
      logical Yes_INOUT

      Yes_INOUT = .False.
      if((x.gt.maxval(xpol)) .or.
     &   (x.lt.minval(xpol)) .or.
     &   (y.gt.maxval(ypol)) .or.
     &   (y.lt.minval(ypol))) then
             Yes_INOUT = .False.
         goto 100
         
      else
          j = npol-1 
          do i=1,npol-1
              if ( ((ypol(i).gt.y).neqv. (ypol(j).gt.y)) .and.
     &               (x .lt. (xpol(j)-xpol(i)) * (y-ypol(i)) 
     &                / (ypol(j)-ypol(i)) + xpol(i)) ) then
                  Yes_INOUT = .NOT. Yes_INOUT
              end if
              j = i
          end do
      end if
      
  100 continue
  
      return 
      end SUBROUTINE Tool_Yes_In_Poly                          
