 
      subroutine Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_End_1,c_Arc_End_2,c_Arc_Center,
     &                           c_Yes_feasible,
     &                           c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,
     &                           c_Arc_Radian)

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
