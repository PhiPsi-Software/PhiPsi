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
 
      subroutine Tool_Yes_On_Arc(x,y,c_Arc_Crack_Coor,Yes_ON)
C	  Determine whether a point P(x,y) is on the arc segment (Added on 2017-07-16)
c     c_Arc_Crack_Coor(10):o_x,o_y,r,Radian_Start,Radian_End,Radian,Point_Start_x,Point_Start_y,Point_End_x,Point_End_y

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in):: x,y,c_Arc_Crack_Coor(11)
      logical,intent(out):: Yes_ON
      

      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) o_x,o_y
      real(kind=FT) r,c_r
      real(kind=FT) c_Arc_Center(2)
      real(kind=FT) c_Angle
      real(kind=FT) c_Arc_Radian_S,c_Arc_Radian_E
      real(kind=FT) c_Arc_Direction
      
      Yes_ON = .False.
      
      o_x = c_Arc_Crack_Coor(1)
      o_y = c_Arc_Crack_Coor(2)
      c_Arc_Center(1:2) = [o_x,o_Y]
      c_Arc_Direction  = c_Arc_Crack_Coor(3)
      r   = c_Arc_Crack_Coor(4)
      
      c_r = Tool_Function_2Point_Dis(c_Arc_Center,[x,y])
      
      c_Arc_Radian_S = c_Arc_Crack_Coor(5)
      c_Arc_Radian_E = c_Arc_Crack_Coor(6)
      
      if(abs(c_r - r)<=Tol_11)then
          call Tool_Angle_of_Inclination_of_VectorAB(
     &                          c_Arc_Center,[x,y],c_Angle)
          
          if(c_Arc_Direction > HLF)then
              if (c_Arc_Radian_E >=c_Arc_Radian_S) then
                  if ((c_Angle >=c_Arc_Radian_S) .and.
     &                (c_Angle <=c_Arc_Radian_E)) then
                      Yes_ON = .True.
                  endif
              elseif(c_Arc_Radian_E < c_Arc_Radian_S)then
                  if ((c_Angle >=c_Arc_Radian_S) .and.
     &                (c_Angle <=Con_360)) then
                      Yes_ON = .True.
                  endif
                  if ((c_Angle >=ZR) .and.
     &                (c_Angle <=c_Arc_Radian_E)) then
                      Yes_ON = .True.
                  endif
              endif
          elseif(c_Arc_Direction < (-HLF))then
              if (c_Arc_Radian_E <=c_Arc_Radian_S) then
                  if ((c_Angle >=c_Arc_Radian_E) .and.
     &                (c_Angle <=c_Arc_Radian_S)) then
                      Yes_ON = .True.
                  endif
              elseif(c_Arc_Radian_E > c_Arc_Radian_S)then
                  if ((c_Angle >=c_Arc_Radian_E) .and.
     &                (c_Angle <=Con_360)) then
                      Yes_ON = .True.
                  endif
                  if ((c_Angle >=ZR) .and.
     &                (c_Angle <=c_Arc_Radian_S)) then
                      Yes_ON = .True.
                  endif
              endif  
          endif
      else
          Yes_ON = .False.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_On_Arc                          
