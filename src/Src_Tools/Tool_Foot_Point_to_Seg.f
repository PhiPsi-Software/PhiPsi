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
 
      subroutine Tool_Foot_Point_to_Seg(a,b,c,Foot)
C     Calculate the foot of the perpendicular from point c to the line segment ab.
c     http://stackoverflow.com/questions/10301001/perpendicular-on-a-line-segment-from-a-given-point
C      function getSpPoint(A,B,C){
C      var x1=A.x, y1=A.y, x2=B.x, y2=B.y, x3=C.x, y3=C.y;
C      var px = x2-x1, py = y2-y1, dAB = px*px + py*py;
C      var u = ((x3 - x1) * px + (y3 - y1) * py) / dAB;
C      var x = x1 + u * px, y = y1 + u * py;
C      return {x:x, y:y}; //this is D
C      }

      use Global_Float_Type 
      implicit none
      real(kind=FT), intent(in) :: a(2), b(2), c(2)
      real(kind=FT), intent(out) :: Foot(2)
      real(kind=FT):: x1,y1,x2,y2,x3,y3,x,y
      real(kind=FT):: px,py,dAB,u
      
      x1=a(1)
      y1=a(2)
      x2=b(1)
      y2=b(2)
      x3=c(1)
      y3=c(2)
      px = x2-x1
      py = y2-y1
      dAB = px*px + py*py
      u = ((x3 - x1) * px + (y3 - y1) * py) / dAB
      x = x1 + u * px
      y = y1 + u * py
      Foot(1:2) = [x,y]
      
      return 
      end SUBROUTINE Tool_Foot_Point_to_Seg                       
