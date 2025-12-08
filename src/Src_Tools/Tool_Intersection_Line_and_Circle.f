 
      subroutine Tool_Intersection_Line_and_Circle(Circle_x0,Circle_y0,
     &                                             R,Seg_A,Seg_B,
     &                                             num_Inter,State,
     &                                             Inter)


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
