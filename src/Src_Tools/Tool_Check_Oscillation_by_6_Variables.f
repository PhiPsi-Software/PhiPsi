 
      subroutine Tool_Check_Oscillation_by_6_Variables(Input_Var,
     &                                                 Yes_Oscill)

      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: Input_Var(6)
      logical,intent(out):: Yes_Oscill
      
      real(kind=FT) Error_1_3_5,Error_2_4_6
      real(kind=FT) Yita1,Yita2
      real(kind=FT) sum_1_3_5,sum_2_4_6,Error_sum
      real(kind=FT) Error_1_4,Error_2_5,Error_3_6
      real(kind=FT) Error_1_2,Error_2_3,Error_1_3
      
      Yes_Oscill = .False.
      Yita1 = 5.0D-2;
      Yita2 = 10.0D-2;
      
      call Tool_Error_Between_3_Variables(Input_Var(1),
     &                                        Input_Var(3),
     &                                        Input_Var(5),
     &                                        Error_1_3_5)
      call Tool_Error_Between_3_Variables(Input_Var(2),
     &                                        Input_Var(4),
     &                                        Input_Var(6),
     &                                        Error_2_4_6)
      if(Error_1_3_5 <= Yita1  .and. Error_2_4_6 <= Yita1 )then
          sum_1_3_5 = Input_Var(1) + Input_Var(3) + Input_Var(5)
          sum_2_4_6 = Input_Var(2) + Input_Var(4) + Input_Var(6)
          Error_sum = abs(sum_1_3_5- sum_2_4_6 ) /
     &                max(sum_1_3_5, sum_2_4_6 )
          if(Error_sum >=Yita2 )then
              Yes_Oscill = .True.
              return
          endif
      endif
      
      Error_1_4 = abs(Input_Var(1)- Input_Var(4) ) /
     &            max(Input_Var(1), Input_Var(4) )
      Error_2_5 = abs(Input_Var(2)- Input_Var(5) ) /
     &            max(Input_Var(2), Input_Var(5) )
      Error_3_6 = abs(Input_Var(3)- Input_Var(6) ) /
     &            max(Input_Var(3), Input_Var(6) )
      if(Error_1_4 <= Yita1  .and. 
     &   Error_2_5 <= Yita1  .and.
     &   Error_3_6 <= Yita1 )    then
          Error_1_2 = abs(Input_Var(1)- Input_Var(2) ) /
     &                max(Input_Var(1), Input_Var(2) )
          Error_2_3 = abs(Input_Var(2)- Input_Var(3) ) /
     &                max(Input_Var(2), Input_Var(3) )
          Error_1_3 = abs(Input_Var(1)- Input_Var(3) ) /
     &                max(Input_Var(1), Input_Var(3) )
          if(Error_1_2 >=Yita2  .or.
     &       Error_2_3 >=Yita2  .or.
     &       Error_1_3 >=Yita2 )then
              Yes_Oscill = .True.
              return
          endif
      endif
      
      return 
      end SUBROUTINE Tool_Check_Oscillation_by_6_Variables                         
