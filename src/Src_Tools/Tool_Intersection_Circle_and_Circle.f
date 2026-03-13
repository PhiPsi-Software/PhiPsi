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
 
      subroutine Tool_Intersection_Circle_and_Circle(
     &                      x1,y1,r1,
     &                      x2,y2,r2,Circles_Status,
     &                      num_Inters,Inters)
C     Calculate the intersection points of two circles
c     Circles_Status
c     = 0, no intersection, far apart;
c     = -999, no intersection, one contains the other
c     = 1, an intersection point, inscribed;
c     =-1, an intersection point, externally tangent;
c      = 2, two points of intersection.
c     https://stackoverflow.com/questions/3349125/circle-circle-intersection-points

      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::x1,y1,r1,x2,y2,r2
      real(kind=FT),intent(out)::Inters(2,2)
      integer,intent(out)::Circles_Status
      integer,intent(out)::num_Inters
      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) c_Dis
      real(kind=FT) Line_AB(2,2), new_Line_AB(2,2)
      real(kind=FT) a,h
      real(kind=FT) P0(2),P1(2),P2(2)
      
      Inters(1:2,1:2) =ZR     
      
      c_Dis = Tool_Function_2Point_Dis([x1,y1],[x2,y2])
      
      if(c_Dis < abs(r1-r2))then
          num_Inters     = 0
          Circles_Status = -999
          return
          
      elseif(c_Dis == abs(r1-r2))then
          num_Inters     = 1
          Circles_Status = 1
          if(r1 > r2)then
              Line_AB(1,1:2) = [x1,y1]
              Line_AB(2,1:2) = [x2,y2]
              call Tool_Shorten_or_Extend_Line(Line_AB,r2,'B',
     &                                     new_Line_AB,Inters(1,1:2))
          elseif(r2 > r1)then
              Line_AB(1,1:2) = [x2,y2]
              Line_AB(2,1:2) = [x1,y1]
              call Tool_Shorten_or_Extend_Line(Line_AB,r1,'B',
     &                                     new_Line_AB,Inters(1,1:2))
              
          endif
          return
      endif
      
      if(c_Dis==(r1+r2))then
          num_Inters =1
          Circles_Status = -1
          Line_AB(1,1:2) = [x1,y1]
          Line_AB(2,1:2) = [x2,y2]
          call Tool_Shorten_or_Extend_Line(Line_AB,-r1,'A',
     &                                     new_Line_AB,Inters(1,1:2))
      elseif(c_Dis>(r1+r2))then
          num_Inters =0
          Circles_Status = 0
      elseif(c_Dis<(r1+r2))then
          num_Inters = 2
          Circles_Status = 2
          if (r1 >= r2)then
              a       = (r1**2 - r2**2 + c_Dis**2)/(2*c_Dis)
              h       = sqrt(r1**2 - a**2)
              P0      = [x1,y1]
              P1      = [x2,y2]
              P2(1:2) = P0 + a*(P1-P0)/c_Dis
              Inters(1,1) = P2(1) + h* ( y2 - y1 ) / c_Dis
              Inters(1,2) = P2(2) - h* ( x2 - x1 ) / c_Dis
              Inters(2,1) = P2(1) - h* ( y2 - y1 ) / c_Dis
              Inters(2,2) = P2(2) + h* ( x2 - x1 ) / c_Dis
          elseif (r1 < r2) then
              a       = (r2**2 - r1**2 + c_Dis**2)/(2*c_Dis)
              h       = sqrt(r2**2 - a**2)
              P0      = [x2,y2]
              P1      = [x1,y1]
              P2(1:2) = P0 + a*(P1-P0)/c_Dis
              Inters(1,1) = P2(1) + h* ( y2 - y1 ) / c_Dis
              Inters(1,2) = P2(2) - h* ( x2 - x1 ) / c_Dis
              Inters(2,1) = P2(1) - h* ( y2 - y1 ) / c_Dis
              Inters(2,2) = P2(2) + h* ( x2 - x1 ) / c_Dis
          endif
      endif
      
      
      
      return 
      end SUBROUTINE Tool_Intersection_Circle_and_Circle                       
