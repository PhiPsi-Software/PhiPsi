 
      subroutine Cal_HF_Flow_Quan(ifra,iter,Counter,total_time) 
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      integer,intent(in)::iter,ifra,Counter
      real(kind=FT),intent(in)::total_time
      integer i_C,i_CalP,Num_Div_Points
      real(kind=FT) c_CalP_w,c_CalP_x,c_CalP_y,n_CalP_x,n_CalP_y
      real(kind=FT) n_CalP_P,c_CalP_P,u,Q
      real(kind=FT) l_CalP_x,l_CalP_y,l_CalP_P,Q_1
      real(kind=FT) Point_1(2),Point_2(2),Point_3(2),a,b,c
      real(kind=FT) CalP_2_x,CalP_3_x,CalP_4_x
      real(kind=FT) CalP_2_y,CalP_3_y,CalP_4_y,CalP_1_x,CalP_1_y
      real(kind=FT) delta_1_2,delta_2_3,delta_3_4,Pa_2_x
      real(kind=FT) Pa_3_x,Pa_4_x,Pa_1_x,delta_x_n,delta_x_l
      real(kind=FT) Max_Aperture,Pre_Gradi,c_CalP_c,Real_Visco
      real(kind=FT) Last_Cr_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP)
      real(kind=FT) u_tem,Inject_Q
      
      Cracks_CalP_Velo(1:Max_Num_Cr,1:Max_Num_Cr_CalP)= ZR 
      Cracks_CalP_Quan(1:Max_Num_Cr,1:Max_Num_Cr_CalP)= ZR 
      
      Last_Cr_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP) =ZR
      if(ifra==1)then
      elseif(ifra > 1)then
          Last_Cr_CalP_Conc = Map_L_Cracks_CalP_Conc
      endif
      call Tool_Get_Value_from_x_y_Curve(Inject_Q_Time,Inject_Q_Val,
     &                                   200,total_time,Inject_Q)
      if(Key_Symm_HF==1)then
          Inject_Q = Inject_Q *HLF
      endif
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then
              Num_Div_Points = Cracks_CalP_Num(i_C)
              do i_CalP=2,Num_Div_Points-1
                  c_CalP_w = Cracks_CalP_Aper(i_C,i_CalP) 
                  c_CalP_P = Cracks_CalP_Pres(i_C,i_CalP)
                  l_CalP_P = Cracks_CalP_Pres(i_C,i_CalP-1)
                  c_CalP_x = Cracks_CalP_Coors(i_C,i_CalP,1)
                  c_CalP_y = Cracks_CalP_Coors(i_C,i_CalP,2)
                  l_CalP_x = Cracks_CalP_Coors(i_C,i_CalP-1,1)
                  l_CalP_y = Cracks_CalP_Coors(i_C,i_CalP-1,2)   
                  n_CalP_P = Cracks_CalP_Pres(i_C,i_CalP+1)
                  n_CalP_x = Cracks_CalP_Coors(i_C,i_CalP+1,1)
                  n_CalP_y = Cracks_CalP_Coors(i_C,i_CalP+1,2) 
                  c_CalP_c = Last_Cr_CalP_Conc(i_C,i_CalP)
                  if(Key_Visco_Type==1)then
                      Real_Visco = Viscosity
                  elseif(Key_Visco_Type==2)then
                      Real_Visco = Viscosity*
     &                   (ONE-c_CalP_c/Max_c)**(-Viscosity_Par_m)
                  endif
                  if(Real_Visco>Visco_Zoom_Factor*Viscosity)then
                      Real_Visco=Visco_Zoom_Factor*Viscosity
                  endif
                  delta_x_l = sqrt((l_CalP_x-c_CalP_x)**2 +
     &                             (l_CalP_y-c_CalP_y)**2)
                  delta_x_n = sqrt((n_CalP_x-c_CalP_x)**2 +
     &                                 (n_CalP_y-c_CalP_y)**2)
                  if(c_CalP_P>=n_CalP_P) then
                      Pre_Gradi = (c_CalP_P-l_CalP_P)/delta_x_l
                  elseif(c_CalP_P < n_CalP_P) then
                      Pre_Gradi = (n_CalP_P-c_CalP_P)/delta_x_n
                  endif
                  if((c_CalP_P>n_CalP_P).and.
     &                  (c_CalP_P>l_CalP_P))then 
                      Pre_Gradi = ZR
                  end if      
                  Cracks_CalP_Pgra(i_C,i_CalP)=Pre_Gradi
                  if(Key_Visco_Type==1)then
                      u = -c_CalP_w**2/12.0D0/Real_Visco*Pre_Gradi
                  elseif(Key_Visco_Type==2)then
                      u_tem = TWO/THR/Real_Visco
                      u = -c_CalP_w**2/12.0D0*Pre_Gradi*u_tem 
                  endif
                  Cracks_CalP_Velo(i_C,i_CalP) = u
                  Q = c_CalP_w*u
                  Cracks_CalP_Quan(i_C,i_CalP) = Q          
              end do
              Max_Aperture = 
     &                    maxval(Cracks_CalP_Aper(i_C,1:Num_Div_Points))
              if (Cracks_CalP_Aper(i_C,1)<=ZPZ1*max_aperture)then
                  Cracks_CalP_Pgra(i_C,1)=ZR
                  Cracks_CalP_Velo(i_C,1)=ZR
                  Cracks_CalP_Quan(i_C,1)=ZR
              else
                  c_CalP_c = Last_Cr_CalP_Conc(i_C,1)
                  if(Key_Visco_Type==1)then
                      Real_Visco = Viscosity
                  elseif(Key_Visco_Type==2)then
                      Real_Visco = Viscosity*
     &                   (ONE-c_CalP_c/Max_c)**(-Viscosity_Par_m)
                  endif
                  if(Real_Visco>Visco_Zoom_Factor*Viscosity)then
                      Real_Visco=Visco_Zoom_Factor*Viscosity
                  endif
                  c_CalP_w = Cracks_CalP_Aper(i_C,1)
                  if(i_C == Inject_Crack_Num.and.
     &                          CalP_num_InjP_Local==1)then
                          Cracks_CalP_Quan(i_C,1)=Inject_Q
                          Cracks_CalP_Velo(i_C,1)=Inject_Q/c_CalP_w
                          Cracks_CalP_Pgra(i_C,1) = 
     &                      -12.0D0*Real_Visco*Inject_Q/c_CalP_w**3
                  else
                      if(Num_Div_Points>=4)then
                          CalP_1_x = Cracks_CalP_Coors(i_C,1,1)
                          CalP_1_y = Cracks_CalP_Coors(i_C,1,2)
                          CalP_2_x = Cracks_CalP_Coors(i_C,2,1)
                          CalP_2_y = Cracks_CalP_Coors(i_C,2,2)
                          CalP_3_x = Cracks_CalP_Coors(i_C,3,1)
                          CalP_3_y = Cracks_CalP_Coors(i_C,3,2)
                          CalP_4_x = Cracks_CalP_Coors(i_C,4,1)
                          CalP_4_y = Cracks_CalP_Coors(i_C,4,2)
                          delta_1_2 = sqrt((CalP_2_x-CalP_1_x)**2 +
     &                                     (CalP_2_y-CalP_1_y)**2)
                          delta_2_3 = sqrt((CalP_3_x-CalP_2_x)**2 +
     &                                     (CalP_3_y-CalP_2_y)**2)
                          delta_3_4 = sqrt((CalP_4_x-CalP_3_x)**2 +
     &                                     (CalP_4_y-CalP_3_y)**2)
                          Pa_1_x = 0   
                          Pa_2_x = delta_1_2
                          Pa_3_x = delta_1_2 + delta_2_3
                          Pa_4_x = delta_1_2 + delta_2_3 + delta_3_4
                          Point_1(1) = Pa_2_x
                          Point_1(2) = Cracks_CalP_Quan(i_C,2)
                          Point_2(1) = Pa_3_x
                          Point_2(2) = Cracks_CalP_Quan(i_C,3)
                          Point_3(1) = Pa_4_x
                          Point_3(2) = Cracks_CalP_Quan(i_C,4)
                          call Tool_abc_of_Parabola(Point_1,
     &                                           Point_2,Point_3,a,b,c)
                          Q_1 = a*Pa_1_x**2 + b*Pa_1_x + c
                          Cracks_CalP_Quan(i_C,1) = Q_1
                          Cracks_CalP_Velo(i_C,1) = Q_1/c_CalP_w
                              Cracks_CalP_Pgra(i_C,1) = 
     &                          -12.0D0*Real_Visco*Q_1/c_CalP_w**3
                      else  
                        Cracks_CalP_Pgra(i_C,1)=Cracks_CalP_Pgra(i_C,2)
                        Cracks_CalP_Velo(i_C,1)=Cracks_CalP_Velo(i_C,2)
                        Cracks_CalP_Quan(i_C,1)=Cracks_CalP_Quan(i_C,2)
                        print *,'    WARNING::less than 4 fluids nodes!'
                      endif
                  endif
              end if
              if (Cracks_CalP_Aper(i_C,Num_Div_Points)
     &                                <= ZPZ1*max_aperture) then
                  Cracks_CalP_Pgra(i_C,Num_Div_Points)=ZR
                  Cracks_CalP_Velo(i_C,Num_Div_Points)=ZR
                  Cracks_CalP_Quan(i_C,Num_Div_Points)=ZR
              else
                  c_CalP_c = Last_Cr_CalP_Conc(i_C,Num_Div_Points)
                  if(Key_Visco_Type==1)then
                      Real_Visco = Viscosity
                  elseif(Key_Visco_Type==2)then
                      Real_Visco = Viscosity*
     &                   (ONE-c_CalP_c/Max_c)**(-Viscosity_Par_m)
                  endif
                  if(Real_Visco>Visco_Zoom_Factor*Viscosity)then
                      Real_Visco=Visco_Zoom_Factor*Viscosity
                  endif
                  c_CalP_w = Cracks_CalP_Aper(i_C,Num_Div_Points)
                  if(i_C == Inject_Crack_Num.and.
     &                        CalP_num_InjP_Local==Num_Div_Points)then
                      Cracks_CalP_Quan(i_C,Num_Div_Points)=Inject_Q
                      Cracks_CalP_Velo(i_C,Num_Div_Points)=
     &                                 Inject_Q/c_CalP_w
                      Cracks_CalP_Pgra(i_C,Num_Div_Points) = 
     &                      -12.0D0*Real_Visco*Inject_Q/c_CalP_w**3
                  else
                      if(Num_Div_Points>=4)then
                          CalP_1_x = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-3,1)
                          CalP_1_y = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-3,2)
                          CalP_2_x = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-2,1)
                          CalP_2_y = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-2,2)
                          CalP_3_x = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-1,1)
                          CalP_3_y = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points-1,2)
                          CalP_4_x = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points,1)
                          CalP_4_y = 
     &                        Cracks_CalP_Coors(i_C,Num_Div_Points,2)
                          delta_1_2 = sqrt((CalP_2_x-CalP_1_x)**2 +
     &                                     (CalP_2_y-CalP_1_y)**2)
                          delta_2_3 = sqrt((CalP_3_x-CalP_2_x)**2 +
     &                                     (CalP_3_y-CalP_2_y)**2)
                          delta_3_4 = sqrt((CalP_4_x-CalP_3_x)**2 +
     &                                     (CalP_4_y-CalP_3_y)**2)
                          Pa_1_x = 0   
                          Pa_2_x = delta_1_2
                          Pa_3_x = delta_1_2 + delta_2_3
                          Pa_4_x = delta_1_2 + delta_2_3 + delta_3_4
                          Point_1(1) = Pa_1_x
                          Point_1(2) = 
     &                           Cracks_CalP_Quan(i_C,Num_Div_Points-3)
                          Point_2(1) = Pa_2_x
                          Point_2(2) = 
     &                           Cracks_CalP_Quan(i_C,Num_Div_Points-2)
                          Point_3(1) = Pa_3_x
                          Point_3(2) = 
     &                           Cracks_CalP_Quan(i_C,Num_Div_Points-1)
                          call Tool_abc_of_Parabola(Point_1,
     &                                           Point_2,Point_3,a,b,c)
                          Q_1 = a*Pa_4_x**2 + b*Pa_4_x + c
                          Cracks_CalP_Quan(i_C,Num_Div_Points) = Q_1
                          Cracks_CalP_Velo(i_C,Num_Div_Points) = 
     &                                     Q_1/c_CalP_w
                              Cracks_CalP_Pgra(i_C,Num_Div_Points) = 
     &                          -12.0D0*Real_Visco*Q_1/c_CalP_w**3
                      else   
                        Cracks_CalP_Pgra(i_C,Num_Div_Points)=
     &                            Cracks_CalP_Pgra(i_C,Num_Div_Points-1)
                        Cracks_CalP_Velo(i_C,Num_Div_Points)=
     &                            Cracks_CalP_Velo(i_C,Num_Div_Points-1)
                        Cracks_CalP_Quan(i_C,Num_Div_Points)=
     &                            Cracks_CalP_Quan(i_C,Num_Div_Points-1)
                        print *,'    WARNING::less than 4 fluids nodes!'
                      endif
                  end if
              end if
          end if
      end do
      
      return 
      end SUBROUTINE Cal_HF_Flow_Quan    
