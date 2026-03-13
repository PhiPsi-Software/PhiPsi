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
 
      subroutine Tool_Intersection_Line_and_Circle(Circle_x0,Circle_y0,
     &                                             R,Seg_A,Seg_B,
     &                                             num_Inter,State,
     &                                             Inter)
C     Calculate the coordinates of the intersection points between line segment AB and the circle
C     A and B are the two endpoints of the line segment, and num_Inter represents the number of intersections.
c     State = -1 - Disjoint (inside the circle);
c     0 - Non-intersecting (outside the circle);
c     1-Tangent
c     2-Intersection


      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: Circle_x0,Circle_y0,Seg_A(2),
     &                              Seg_B(2),R
      integer,intent(out)::num_Inter,State
      real(kind=FT),intent(out)::Inter(2,2)
      
      real(kind=FT) Tool_Function_2Point_Dis
      logical Yes_Seg_A_in,Yes_Seg_B_in
      real(kind=FT) Line_AB(2,2),S_Distance
      real(kind=FT) c_Foot(2),c_Dis_Foot_Inter
      real(kind=FT) Line_A_foot(2,2),Line_B_foot(2,2)
      real(kind=FT) tem_new_Line(2,2),c_tem_new_Point(2)
      real(kind=FT) Length_AB,Length_BT,Length_OB
      real(kind=FT) Length_OA,Length_AT
      
      Inter(1:2,1:2) = ZR 
      
      Line_AB(1,1:2) = Seg_A
      Line_AB(2,1:2) = Seg_B
      
      Length_AB = Tool_Function_2Point_Dis(Seg_A,Seg_B)
    
      
      
      
      Yes_Seg_A_in = .False.
      if(Tool_Function_2Point_Dis([Circle_x0,Circle_y0],Seg_A)<=R)then
          Yes_Seg_A_in  = .True.
      endif
      Yes_Seg_B_in = .False.
      if(Tool_Function_2Point_Dis([Circle_x0,Circle_y0],Seg_B)<=R)then
          Yes_Seg_B_in  = .True.
      endif
      
      if(Yes_Seg_A_in .and. Yes_Seg_B_in)then
          State = -1
          num_Inter = 0
          return
      elseif((Yes_Seg_A_in .EQV. .FALSE.) .and. 
     &       (Yes_Seg_B_in .EQV. .FALSE.)) then
          Length_OB = Tool_Function_2Point_Dis(
     &                [Circle_x0,Circle_y0],Seg_B)
          Length_OA = Tool_Function_2Point_Dis(
     &                [Circle_x0,Circle_y0],Seg_A)
          Length_BT = sqrt(Length_OB**2 - R**2)
          Length_AT = sqrt(Length_OA**2 - R**2)
          if(Length_OB>=Length_OA)then
              if(Length_AB<=Length_BT)then
                  State = 0
                  num_Inter = 0
                  return
              endif
          else
              if(Length_AB<=Length_AT)then
                  State = 0
                  num_Inter = 0
                  return
              endif
          endif
          call Cal_Signed_Distance(
     &         Line_AB,[Circle_x0,Circle_y0],S_Distance)
          if(abs(S_Distance)<R) then
              State = 2
              num_Inter = 2
              call Tool_Foot_Point_to_Seg(Seg_A,Seg_B,
     &                         [Circle_x0,Circle_y0],c_Foot)
              c_Dis_Foot_Inter = sqrt(R**2 - S_Distance**2)
              Line_A_foot(1,1:2) =  Seg_A
              Line_A_foot(2,1:2) =  c_Foot
              call Tool_Shorten_or_Extend_Line(
     &             Line_A_foot,-c_Dis_Foot_Inter,
     &                           'B',tem_new_Line,c_tem_new_Point)
              Inter(1,1:2) = c_tem_new_Point
              Line_B_foot(1,1:2) =  Seg_B
              Line_B_foot(2,1:2) =  c_Foot
              call Tool_Shorten_or_Extend_Line(
     &             Line_B_foot,-c_Dis_Foot_Inter,
     &                           'B',tem_new_Line,c_tem_new_Point)
              Inter(2,1:2) = c_tem_new_Point
          elseif(abs(S_Distance)==R) then
              
              State = 1
              num_Inter = 1
              call Tool_Foot_Point_to_Seg(Seg_A,Seg_B,
     &                         [Circle_x0,Circle_y0],c_Foot)
              Inter(1,1:2) = c_Foot
          elseif(abs(S_Distance)>R) then
              State = 0
              num_Inter = 0
              return
          endif
      elseif(Yes_Seg_A_in .or. Yes_Seg_B_in) then
          State = 2
          num_Inter = 1
          call Cal_Signed_Distance(
     &         Line_AB,[Circle_x0,Circle_y0],S_Distance)
          call Tool_Foot_Point_to_Seg(Seg_A,Seg_B,
     &                             [Circle_x0,Circle_y0],c_Foot)
          c_Dis_Foot_Inter = sqrt(R**2 - S_Distance**2)
          if(Yes_Seg_A_in)then
              Line_B_foot(1,1:2) =  Seg_B
              Line_B_foot(2,1:2) =  c_Foot
              call Tool_Shorten_or_Extend_Line(
     &             Line_B_foot,-c_Dis_Foot_Inter,
     &                           'B',tem_new_Line,c_tem_new_Point)
              Inter(1,1:2) = c_tem_new_Point
          elseif(Yes_Seg_B_in)then
              Line_A_foot(1,1:2) =  Seg_A
              Line_A_foot(2,1:2) =  c_Foot
              call Tool_Shorten_or_Extend_Line(
     &             Line_A_foot,-c_Dis_Foot_Inter,
     &                           'B',tem_new_Line,c_tem_new_Point)
              Inter(1,1:2) = c_tem_new_Point
          endif
      endif
      return 
      end SUBROUTINE Tool_Intersection_Line_and_Circle   
