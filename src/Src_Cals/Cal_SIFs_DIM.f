 
      subroutine Cal_SIFs_DIM(iter,c_DISP)

      use Global_Crack
      use Global_Crack_Common
      use Global_Elem_Area_Vol
      use Global_Model
      use Global_Common
      use Global_Material
      
      implicit none
      integer,intent(in)::iter
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      integer i_C,i_Tip,Tip_Elem,mat_num
      real(kind=FT) cr_Tip(2,2),omega_Tip(2),delta_L,
     &                 c_E,c_v,c_G,k,offset_delta,Ux,Uy,Dx,Dy,
     &                 Disp_U(2),Disp_D(2),delta_Disp_x,delta_Disp_y,
     &                 delta_Disp_x_CRACK,delta_Disp_y_CRACK
      integer c_Elem
      real(kind=FT) c_Kesi,c_Yita
      real(kind=FT) c_Ave_Length
      real(kind=FT) x,y

      real(kind=FT) delta_Disp_x1,delta_Disp_y1,
     &              delta_Disp_x1_CRACK,delta_Disp_y1_CRACK
      real(kind=FT) delta_Disp_x2,delta_Disp_y2,
     &              delta_Disp_x2_CRACK,delta_Disp_y2_CRACK
      real(kind=FT) delta_Disp_x3,delta_Disp_y3,
     &              delta_Disp_x3_CRACK,delta_Disp_y3_CRACK   

      real(kind=FT) r_1_factor,r_2_factor,r_3_factor
      real(kind=FT) r_1,r_2,r_3
      real(kind=FT) KI_r1,KI_r2,KI_r3,KII_r1,KII_r2,KII_r3
      real(kind=FT) x1,y1,x2,y2,x3,y3
      real(kind=FT) Ux1,Uy1,Dx1,Dy1
      real(kind=FT) Ux2,Uy2,Dx2,Dy2
      real(kind=FT) Ux3,Uy3,Dx3,Dy3
      real(kind=FT) Disp_U1(2),Disp_D1(2)
      real(kind=FT) Disp_U2(2),Disp_D2(2)
      real(kind=FT) Disp_U3(2),Disp_D3(2)
      real(kind=FT) c1,c2,c3,temp
      real(kind=FT) Matrix(3,3),U(3),F(3)
      print *,'    Calculating SIFs of each crack......'   
      
      
      r_1_factor = 0.1D0
      r_2_factor = 0.2D0
      r_3_factor = 0.3D0

      
      
      
      offset_delta = Delta_Factor_Aper*Ave_Elem_L_Enrich
      
      
      KI(1:Max_Num_Cr,1:2) = ZR
      KII(1:Max_Num_Cr,1:2) = ZR
      
      do i_C =1,num_Crack
          cr_Tip(1,1:2) = Cr_First_Tip(i_C,1:2)
          cr_Tip(2,1:2) = Cr_Second_Tip(i_C,1:2)
          omega_Tip(1) = Cr_First_Tip_Ori(i_C)
          omega_Tip(2) = Cr_Second_Tip_Ori(i_C)
          do i_Tip = 1,2
              if ((Crack_Tip_Type(i_C,i_Tip) .ne. 1) .and.
     &            (Crack_Tip_Type(i_C,i_Tip) .ne. -2)) then
                  call Cal_Length_of_ELe_by_Coors_Ave(
     &                 cr_Tip(i_Tip,1),cr_Tip(i_Tip,2),c_Ave_Length)
                  r_1 = r_1_factor*c_Ave_Length
                  r_2 = r_2_factor*c_Ave_Length
                  r_3 = r_3_factor*c_Ave_Length

                  x1 = cr_Tip(i_Tip,1) - r_1*cos(omega_Tip(i_Tip))
                  y1 = cr_Tip(i_Tip,2) - r_1*sin(omega_Tip(i_Tip))
                  x2 = cr_Tip(i_Tip,1) - r_2*cos(omega_Tip(i_Tip))
                  y2 = cr_Tip(i_Tip,2) - r_2*sin(omega_Tip(i_Tip))
                  x3 = cr_Tip(i_Tip,1) - r_3*cos(omega_Tip(i_Tip))
                  y3 = cr_Tip(i_Tip,2) - r_3*sin(omega_Tip(i_Tip))     
                  

                  call Cal_Ele_Num_by_Coors(cr_Tip(i_Tip,1),
     &                                      cr_Tip(i_Tip,2),Tip_Elem)
                  mat_num = Elem_Mat(Tip_Elem)
                  if (Material_Type(mat_num). eq. 1  )then
                      c_E  = Material_Para(mat_num,1)
                    c_v  = Material_Para(mat_num,2)
                  else
                      c_E  = Material_Para(mat_num,1)
                    c_v  = Material_Para(mat_num,2)
                  end if
                  c_G          = c_E/TWO/(ONE+c_v)
                  if (Key_Type_2D == 1) then
                      k = (THR-c_v)/(ONE+c_v)
                  elseif (Key_Type_2D == 2) then
                      k = THR - FOU*c_v
                  end if
                  Ux1 = x1 - offset_delta*sin(omega_Tip(i_Tip))
                  Uy1 = y1 + offset_delta*cos(omega_Tip(i_Tip))
                  Dx1 = x1 + offset_delta*sin(omega_Tip(i_Tip))
                  Dy1 = y1 - offset_delta*cos(omega_Tip(i_Tip))
                  Ux2 = x2 - offset_delta*sin(omega_Tip(i_Tip))
                  Uy2 = y2 + offset_delta*cos(omega_Tip(i_Tip))
                  Dx2 = x2 + offset_delta*sin(omega_Tip(i_Tip))
                  Dy2 = y2 - offset_delta*cos(omega_Tip(i_Tip))     
                  Ux3 = x3 - offset_delta*sin(omega_Tip(i_Tip))
                  Uy3 = y3 + offset_delta*cos(omega_Tip(i_Tip))
                  Dx3 = x3 + offset_delta*sin(omega_Tip(i_Tip))
                  Dy3 = y3 - offset_delta*cos(omega_Tip(i_Tip))  
                  call Cal_Ele_Num_by_Coors(Ux1,Uy1,c_Elem)
                  call Cal_KesiYita_by_Coor([Ux1,Uy1],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_U1)
                  call Cal_Ele_Num_by_Coors(Dx1,Dy1,c_Elem)
                  call Cal_KesiYita_by_Coor([Dx1,Dy1],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_D1)
                 
                  call Cal_Ele_Num_by_Coors(Ux2,Uy2,c_Elem)
                  call Cal_KesiYita_by_Coor([Ux2,Uy2],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_U2)
                  call Cal_Ele_Num_by_Coors(Dx2,Dy2,c_Elem)
                  call Cal_KesiYita_by_Coor([Dx2,Dy2],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_D2)
                  call Cal_Ele_Num_by_Coors(Ux3,Uy3,c_Elem)
                  call Cal_KesiYita_by_Coor([Ux3,Uy3],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_U3)
                  call Cal_Ele_Num_by_Coors(Dx3,Dy3,c_Elem)
                  call Cal_KesiYita_by_Coor([Dx3,Dy3],c_Elem,
     &                                     c_Kesi,c_Yita)
                  call Cal_Any_Point_Disp_KesiYita(c_Elem,
     &                                   c_Kesi,c_Yita,0,c_DISP,Disp_D3)
                  delta_Disp_x1 = Disp_U1(1)-Disp_D1(1)
                  delta_Disp_y1 = Disp_U1(2)-Disp_D1(2)
                  delta_Disp_x1_CRACK =  
     &                         delta_Disp_x1*cos(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y1*sin(omega_Tip(i_Tip))
                  delta_Disp_y1_CRACK = 
     &                        -delta_Disp_x1*sin(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y1*cos(omega_Tip(i_Tip))     
                  KI_r1   = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_y1_CRACK/sqrt(r_1)
                  KII_r1  = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_x1_CRACK/sqrt(r_1)  
                  delta_Disp_x2 = Disp_U2(1)-Disp_D2(1)
                  delta_Disp_y2 = Disp_U2(2)-Disp_D2(2)
                  delta_Disp_x2_CRACK =  
     &                         delta_Disp_x2*cos(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y2*sin(omega_Tip(i_Tip))
                  delta_Disp_y2_CRACK = 
     &                        -delta_Disp_x2*sin(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y2*cos(omega_Tip(i_Tip))     
                  KI_r2   = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_y2_CRACK/sqrt(r_2)
                  KII_r2  = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_x2_CRACK/sqrt(r_2)  
                  delta_Disp_x3 = Disp_U3(1)-Disp_D3(1)
                  delta_Disp_y3 = Disp_U3(2)-Disp_D3(2)
                  delta_Disp_x3_CRACK =  
     &                         delta_Disp_x3*cos(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y3*sin(omega_Tip(i_Tip))
                  delta_Disp_y3_CRACK = 
     &                        -delta_Disp_x3*sin(omega_Tip(i_Tip)) + 
     &                         delta_Disp_y3*cos(omega_Tip(i_Tip))     
                  KI_r3   = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_y3_CRACK/sqrt(r_3)
                  KII_r3  = sqrt(2*pi)*c_G/(1+k)*
     &                   delta_Disp_x3_CRACK/sqrt(r_3)       
     
                  
                  if(i_Tip==1)then
                  KI(i_C,i_Tip)   = KI_r1
                  KII(i_C,i_Tip)  = KII_r1  

                  KI(i_C,i_Tip)  = KI_r1+(r_1/(r_2-r_1))*(KI_r1-KI_r2)
                  KII(i_C,i_Tip) =KII_r1+(r_1/(r_2-r_1))*(KII_r1-KII_r2)

                  c1 = r_2*r_3**2 - r_2**2*r_3
                  c2 = r_3*r_1**2 - r_3**2*r_1
                  c3 = r_1*r_2**2 - r_1**2*r_2
                  temp = c1+c2+c3
                  KI(i_C,i_Tip) = (KI_r1*c1+KI_r2*c2+KI_r3*c3)/temp
                  KII(i_C,i_Tip)= (KII_r1*c1+KII_r2*c2+KII_r3*c3)/temp 
            
                  endif
                  if(Key_SIFs_DIM_Points==1)then
                      KI(i_C,i_Tip)   = KI_r1
                      KII(i_C,i_Tip)  = KII_r1  
                  elseif(Key_SIFs_DIM_Points==2)then
                      KI(i_C,i_Tip)   = KI_r1+
     &                                  (r_1/(r_2-r_1))*(KI_r1-KI_r2)
                      KII(i_C,i_Tip)  = KII_r1+
     &                                  (r_1/(r_2-r_1))*(KII_r1-KII_r2)
                  elseif(Key_SIFs_DIM_Points==3)then
                      c1 = r_2*r_3**2 - r_2**2*r_3
                      c2 = r_3*r_1**2 - r_3**2*r_1
                      c3 = r_1*r_2**2 - r_1**2*r_2
                      temp = c1+c2+c3
                      KI(i_C,i_Tip) = (KI_r1*c1+KI_r2*c2+KI_r3*c3)/temp
                      KII(i_C,i_Tip)= (KII_r1*c1+KII_r2*c2+KII_r3*c3)
     &                                /temp 
                  endif
              end if
          end do
      end do
      
      
      RETURN
      end subroutine Cal_SIFs_DIM