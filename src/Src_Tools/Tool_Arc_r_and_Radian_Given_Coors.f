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
 
      subroutine Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_End_1,c_Arc_End_2,c_Arc_Center,
     &                           c_Yes_feasible,
     &                           c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,
     &                           c_Arc_Radian)
c     Related to arc calculations:
c     Given the endpoint coordinates of an arc and the coordinates of the arc's center, as well as the arc direction (c_Arc_Direction,
c     Counterclockwise is 1, clockwise is -1), first determine whether the data is feasible (whether the radius is consistent),
c     If feasible, calculate the radius, as well as the starting arc angle and the arc length (c_Arc_Radian)

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in):: c_Arc_Direction,
     &                           c_Arc_End_1(2),c_Arc_End_2(2),
     &                           c_Arc_Center(2)
      real(kind=FT),intent(out)::c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,
     &                           c_Arc_Radian
      logical,intent(out)::c_Yes_feasible
      
      real(kind=FT) r1,r2
      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) c_Angle_1,c_Angle_2
      
      
      c_Yes_feasible =.False.
      r1 = Tool_Function_2Point_Dis(c_Arc_Center,c_Arc_End_1)
      r2 = Tool_Function_2Point_Dis(c_Arc_Center,c_Arc_End_2)
      
      if(abs(r1-r2)<=Tol_11)then
          c_Yes_feasible =.True.
          c_Arc_r = r1
          call Tool_Angle_of_Inclination_of_VectorAB(
     &                          c_Arc_Center,c_Arc_End_1,c_Angle_1)
          call Tool_Angle_of_Inclination_of_VectorAB(
     &                          c_Arc_Center,c_Arc_End_2,c_Angle_2)
     
          if(abs(c_Angle_1)<=Tol_11)then   
              c_Angle_1 = Con_360
          endif
          if(abs(c_Angle_2)<=Tol_11)then
              c_Angle_2 = Con_360
          endif     
          
          c_Arc_Radian_S = c_Angle_1
          c_Arc_Radian_E = c_Angle_2

          if ((c_Arc_Radian_E-c_Arc_Radian_S) >= ZR) then
              c_Arc_Radian = c_Arc_Radian_E-c_Arc_Radian_S
          else
              c_Arc_Radian = c_Arc_Radian_E +(Con_360-c_Arc_Radian_S)
          endif
          if (c_Arc_Radian>Con_360)then
              c_Arc_Radian =c_Arc_Radian - Con_360
          endif
          
          if(c_Arc_Direction < (-HLF))then
              c_Arc_Radian = Con_360 - c_Arc_Radian
          endif
      endif    
      
      return 
      end SUBROUTINE Tool_Arc_r_and_Radian_Given_Coors                    
