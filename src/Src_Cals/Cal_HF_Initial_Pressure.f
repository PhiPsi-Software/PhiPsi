 
      subroutine Cal_HF_Initial_Pressure(iter,ifra,Counter,
     &                                   Yes_Growth,
     &                                   Last_Cracks_CalP_Num_ifra,
     &                                   Last_Cr_CalP_Pres_ifra,
     &                                   Initial_Cr_CalP_Pres) 
     
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      
      integer,intent(in)::iter,ifra,Counter
      real(kind=FT),intent(in)::Last_Cr_CalP_Pres_ifra(Max_Num_Cr,
     &                                   Max_Num_Cr_CalP)
      integer,intent(in)::Last_Cracks_CalP_Num_ifra(Max_Num_Cr)
      real(kind=FT),intent(out)::Initial_Cr_CalP_Pres(Max_Num_Cr,
     &                                   Max_Num_Cr_CalP) 
      logical,intent(in)::Yes_Growth(Max_Num_Cr,2)
      integer i_C,i_CalP,Num_Div_Points,num_Count,j_C
      real(kind=FT) Theor_Press,Inject_Crack_L
      real(kind=FT) Inj_x,Inj_y,c_iCalp_Dis,c_iCalp_x,c_iCalp_y,
     &                 l_Crack,Kesi,tem,tem1,tem2,tem3,cc_time,
     &                 A0,A1,B1,E_1,u_1
      real(kind=FT) c_Pressure,factor
      integer Positive_Cr,Num_Div_Points_P_Cr,Num_Div_Points_N_Cr
      real(kind=FT) Max_Pres_P_Cr,Factor_Max_Pres
      real(kind=FT) Max_InSitu
      
      A0 = sqrt(THR)
      A1 = -0.156D0
      B1 = 0.0663D0
      E_1 = Material_Para(1,1)/(ONE-Material_Para(1,2)**2)
      u_1 = 12.0D0*Viscosity
      Max_InSitu = max(InSitu_x,InSitu_y)
      
      if (ifra==1) then
          do i_C = 1,num_Crack
              if (Cracks_HF_State(i_C) == 1) then      
                  Num_Div_Points = Cracks_CalP_Num(i_C)
                  call Tool_Length_of_Crack(Inject_Crack_Num,
     &                                             Inject_Crack_L)
                  if(Key_Symm_HF==0)then
                      Inject_Crack_L=Inject_Crack_L*HLF
                  end if
                  call Cal_HF_Theoretic_Pressure(iter,Counter,
     &                     Inject_Crack_L,HF_Theor_Time,HF_Ini_Pressure,
     &                                              HF_Theor_Aper) 
                  Num_Div_Points = Cracks_CalP_Num(i_C)
                  
                  Initial_Cr_CalP_Pres(i_C,1:Num_Div_Points)=
     &                                  HF_Ini_Pressure
     

              endif
          end do
      else
          do i_C=1,num_Crack
              if (Cracks_HF_State(i_C) == 1) then   
                  Num_Div_Points = Cracks_CalP_Num(i_C)
                  Initial_Cr_CalP_Pres(i_C,1:Num_Div_Points)=
     &                                      Last_Inj_Pres
             endif
          enddo
      endif
      
      return 
      end SUBROUTINE Cal_HF_Initial_Pressure             
