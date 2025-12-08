 
      subroutine Cal_HF_S(Counter,S,total_time)
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      integer,intent(in)::Counter
      real(kind=FT),intent(in)::total_time
      real(kind=FT),intent(out)::S(num_Tol_CalP_Water)
      integer i_C,i_CalP,Num_Div_Points,num_Count,i_Pair,CalP_1,CalP_2 
      real(kind=FT) Inject_Q,Length_Out,CalP_x,CalP_y
      real(kind=FT) Theor_Time,Theor_Pres,delta_S,delta_Time,
     &                 Theor_Aper
      
      
      S(1:num_Tol_CalP_Water)=ZR
      
      if(total_time > maxval(Inject_Q_Time))then
           print *,'    ~~~~~~!~~~~~~!~~~~~~~!~~~~~~!~~~~~~!~~~~~~!~~~~'
           print *,'    Error:: time beyonds the range of '
     &                  //'Inject_Q_Time!'
           print *,'    ~~~~~~!~~~~~~!~~~~~~~!~~~~~~!~~~~~~!~~~~~~!~~~~'
           call Warning_Message('S',Keywords_Blank)  
      endif
      call Tool_Get_Value_from_x_y_Curve(Inject_Q_Time,Inject_Q_Val,
     &                                   200,total_time,Inject_Q)
      num_Count = 0
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then  
              Num_Div_Points = Cracks_CalP_Num(i_C)
              do i_CalP = 1,Num_Div_Points
                  num_Count = num_Count + 1
                  CalP_x=Cracks_CalP_Coors(i_C,i_CalP,1)
                  CalP_y=Cracks_CalP_Coors(i_C,i_CalP,2)
                  if(Key_Symm_HF==1)then
                      if (i_CalP == 1) then
                          if (i_C == Inject_Crack_Num) then
                              S(num_Count) = Inject_Q*HLF
                          end if
                      else
                          S(num_Count) = ZR
                      end if
                  elseif(Key_Symm_HF==0)then
                      if (i_CalP == CalP_num_InjP_Local) then
                          if (i_C == Inject_Crack_Num) then
                              S(num_Count) = Inject_Q
                          end if
                      else
                          S(num_Count) = ZR
                      end if
                  end if
                  if(Key_Leakoff==1)then
                      call Tool_Length_CalP_to_InjP(i_C,
     &                             Length_Out,CalP_x,CalP_y,i_CalP)
                      call Cal_HF_Theoretic_Pressure(1,Counter,
     &                    Length_Out,Theor_Time,Theor_Pres,Theor_Aper)
                      delta_Time = total_time-Theor_Time
                      if(delta_Time<=ZR)then
                          delta_Time = 1.0D-1
                      endif
                      delta_S =TWO*Coeff_Leak/sqrt(delta_Time)
                      S(num_Count)=S(num_Count)- delta_S   
                  end if
              end do
          end if
      end do
      
      
      return 
      end SUBROUTINE Cal_HF_S               
