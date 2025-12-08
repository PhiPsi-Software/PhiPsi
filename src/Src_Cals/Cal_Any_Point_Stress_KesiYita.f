 
      subroutine Cal_Any_Point_Stress_KesiYita(c_Elem,in_Kesi,
     &                                       in_Yita,i_G,
     &                                       c_DISP,
     &                                       Sxx_GP, Syy_GP, Sxy_GP)
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_DISP
      use Global_Inclusion
      use Global_Cross
      use Global_Material
      use Global_Stress
      
      implicit none
      
      integer,intent(in)::c_Elem,i_G
      real(kind=FT),intent(in)::in_Kesi,in_Yita,c_DISP(Total_FD)
      real(kind=FT),intent(out)::Sxx_GP, Syy_GP, Sxy_GP
      integer c_NN(4),G_Counter
      real(kind=FT) c_thick,c_D(3,3),U(80)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),
     &                 c_Stress(3),B(3,80),
     &                 tem_B(3,80)
      integer i_C,num_B,num_tem_B
      integer:: Location_ESM(MDOF_2D)
      integer Location_ESM_C_Cr_NoFEM(60)
      integer num_Loc_ESM_C_Cr_NoFEM
      integer num_Loc_ESM
      integer::Location_ESM_C_Crack(80)
      integer num_Loc_ESM_C_Crack
      integer i_H,i_Incl
      real(kind=FT) detJ,c_N(2,8)
      real(kind=FT) c_G_x,c_G_y
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) c_T_Alpha,c_v 
      integer i_Cross
      real(kind=FT) c_TStress(3)
     


      c_thick = thick(Elem_Mat(c_Elem))
      
      c_D     = D(Elem_Mat(c_Elem),:,:)  
      
      if(Flag_Weibull_E)then
          if (Key_Weibull_E(Elem_Mat(c_Elem)) ==1)then
              c_D = Weibull_Elements_D_Matrix(c_Elem,1:3,1:3)
          endif
      endif
          
      c_v     = v(Elem_Mat(c_Elem),1)
      c_T_Alpha  = T_Alpha(Elem_Mat(c_Elem)) 
      c_NN    = G_NN(:,c_Elem)
      c_X_NODES = G_X_NODES(:,c_Elem)
      c_Y_NODES = G_Y_NODES(:,c_Elem)       
      Location_ESM(1:MDOF_2D)  = 0
      num_Loc_ESM = 0
      if(num_Crack/=0)then
        do i_C =1,num_Crack 
          call Location_Element_Stiff_Matrix(c_Elem,i_C,
     &                                    c_POS(:,i_C),
     &                                    Location_ESM_C_Crack,
     &                                    num_Loc_ESM_C_Crack,
     &                                    Location_ESM_C_Cr_NoFEM,
     &                                    num_Loc_ESM_C_Cr_NoFEM)
          Location_ESM(num_Loc_ESM+1:
     &                   num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                   Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack                     
        end do
      endif
      if(num_Hole/=0)then
        do i_H =1,num_Hole
          call Location_Element_Stiff_Matrix_Hl(c_Elem,i_H,
     &                                c_POS_Hl(:,i_H),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
          Location_ESM(num_Loc_ESM+1:
     &              num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &              Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
      endif
      if(num_Cross/=0)then
        do i_Cross =1,num_Cross
          call Location_Element_Stiff_Matrix_Cross(c_Elem,i_Cross,
     &                                c_POS_Cross(:,i_Cross),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
          Location_ESM(num_Loc_ESM+1:
     &              num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &              Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
      endif
      if(num_Inclusion/=0)then
        do i_Incl =1,num_Inclusion
          call Location_Element_Stiff_Matrix_Incl(c_Elem,i_Incl,
     &                                c_POS_Incl(:,i_Incl),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
          Location_ESM(num_Loc_ESM+1:
     &              num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &              Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
      endif
      U(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))
    
      
      c_D     = D(Elem_Mat(c_Elem),:,:) 
      
      if(Flag_Weibull_E)then
          if (Key_Weibull_E(Elem_Mat(c_Elem)) ==1)then
              c_D = Weibull_Elements_D_Matrix(c_Elem,1:3,1:3)
          endif
      endif
      
      c_v     = v(Elem_Mat(c_Elem),1)
      c_T_Alpha = T_Alpha(Elem_Mat(c_Elem))
      G_Counter =  G_Counter +1    
      
      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.
      
      B(1:3,1:80) = ZR
      num_B = 0 
      call Cal_N(in_Kesi,in_Yita,c_N)
      call Cal_detJ(in_Kesi,in_Yita,c_X_NODES,c_Y_NODES,detJ)  
      c_G_x = DOT_PRODUCT(c_N(1,1:7:2),c_X_NODES(1:4))
      c_G_y = DOT_PRODUCT(c_N(1,1:7:2),c_Y_NODES(1:4)) 
      if(num_Crack/=0)then
        do i_C =1,num_Crack 
          call Cal_B_Matrix_Crack(in_Kesi,in_Yita,i_C,c_Elem,i_G,
     &                      c_NN,c_X_NODES,c_Y_NODES,
     &                      tem_B,num_tem_B)                  
          B(1:3,num_B+1:num_B+num_tem_B)=tem_B(1:3,1:num_tem_B)
          num_B = num_B + num_tem_B                      
        end do
      endif
      if(num_Hole/=0)then
          do i_H =1,num_Hole 
              call Cal_B_Matrix_Hl(in_Kesi,in_Yita,
     &                            i_H,c_Elem,i_G,
     &                            c_NN,c_X_NODES,c_Y_NODES,
     &                            tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  
     &                            tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
      endif
      if(num_Cross/=0)then
          do i_Cross =1,num_Cross
              call Cal_B_Matrix_Cross(in_Kesi,in_Yita,
     &                            i_Cross,c_Elem,i_G,
     &                            c_NN,c_X_NODES,c_Y_NODES,
     &                            tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  
     &                            tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
      endif
      if(num_Circ_Incl/=0)then
          do i_Incl =1,num_Circ_Incl 
              call Cal_B_Matrix_Circ_Incl(in_Kesi,in_Yita,
     &                                i_Incl,c_Elem,i_G,
     &                                c_NN,c_X_NODES,c_Y_NODES,
     &                                tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
          call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                        Yes_Gauss_in_Incl,c_Incl_Num)
          if(Yes_Gauss_in_Incl)then
              c_MatNum = Circ_Inclu_Mat_Num(c_Incl_Num)
              c_D     = D(c_MatNum,:,:)   
              c_T_Alpha = T_Alpha(c_MatNum)
              c_v = v(c_MatNum,1)
          endif
      endif
      if(num_Poly_Incl/=0)then
          do i_Incl =1,num_Poly_Incl 
              call Cal_B_Matrix_Poly_Incl(in_Kesi,in_Yita,
     &                                i_Incl,c_Elem,i_G,
     &                                c_NN,c_X_NODES,c_Y_NODES,
     &                                tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
          call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                        Yes_Gauss_in_Incl,c_Incl_Num)
          if(Yes_Gauss_in_Incl)then
              c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
              c_D     = D(c_MatNum,:,:)   
              c_T_Alpha = T_Alpha(c_MatNum)
              c_v = v(c_MatNum,1)
          endif
      endif
      c_Stress=MATMUL(MATMUL(c_D,B(1:3,1:num_Loc_ESM)),
     &                U(1:num_Loc_ESM))
     
      
      if(Key_InSitu_Strategy==2)then
           c_Stress(1)=c_Stress(1)
     &      + sum(InSitu_Strs_Gaus_xx(c_Elem,1:Num_Gauss_P_FEM))/
     &        Num_Gauss_P_FEM
           c_Stress(2)=c_Stress(2)
     &      + sum(InSitu_Strs_Gaus_yy(c_Elem,1:Num_Gauss_P_FEM))/
     &        Num_Gauss_P_FEM
           c_Stress(3)=c_Stress(3)
     &      + sum(InSitu_Strs_Gaus_xy(c_Elem,1:Num_Gauss_P_FEM))/
     &        Num_Gauss_P_FEM 
      endif


      if(Key_Thermal_Stress==1)then
          c_T_Alpha = T_Alpha(Elem_Mat(c_Elem))
          c_TStress = c_T_Alpha*Thermal_Str_Temper(Elem_Mat(c_Elem))*
     &                MATMUL(c_D,[ONE,ONE,ZR]) 
          c_Stress(1) = c_Stress(1) - c_TStress(1)
          c_Stress(2) = c_Stress(2) - c_TStress(2)  
          c_Stress(3) = c_Stress(3) - c_TStress(3)
      endif
      
      Sxx_GP = c_Stress(1)
      Syy_GP = c_Stress(2)
      Sxy_GP = c_Stress(3)

      
      return 
      end SUBROUTINE Cal_Any_Point_Stress_KesiYita        
