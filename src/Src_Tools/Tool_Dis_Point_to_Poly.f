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
 
      subroutine Tool_Dis_Point_to_Poly(x,y,
     &                                  xpol,ypol,npol,
     &                                  Dis)
C     Calculate the distance from point P to the polygon; it is positive if outside the polygon and negative if inside the polygon.
c     Note: The polygons here are in closed form, npol = number of sides 1

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type     
      implicit none
      integer,intent(in)::npol
      real(kind=FT),intent(in)::x,y,xpol(npol),ypol(npol)
      real(kind=FT),intent(out)::Dis
      logical Yes_in_Poly,Yes_on_Poly
      real(kind=FT) Dis_to_Edges(npol-1),Line_Edge(2,2)
      integer i_Edge
      call Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_in_Poly)
      
      call Tool_Yes_On_Poly(x,y,xpol,ypol,npol,Yes_ON_Poly)
      if(Yes_ON_Poly)then
          Dis = ZR
          return
      endif
      
      Dis_to_Edges(1:npol-1) = ZR
      do i_Edge = 1,npol-1
          Line_Edge(1,1) = xpol(i_Edge)
          Line_Edge(1,2) = ypol(i_Edge)
          Line_Edge(2,1) = xpol(i_Edge+1)
          Line_Edge(2,2) = ypol(i_Edge+1)
          call Tool_Dis_Point_to_Seg(Line_Edge(1,1:2),
     &                               Line_Edge(2,1:2),
     &                               [x,y],
     &                               Dis_to_Edges(i_Edge))
     
      enddo
      if(Yes_in_Poly)then
          Dis = -minval(Dis_to_Edges(1:npol-1))
      else
          Dis =  minval(Dis_to_Edges(1:npol-1))
      endif
       
      return 
      end SUBROUTINE Tool_Dis_Point_to_Poly                       
