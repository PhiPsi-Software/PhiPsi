 
      subroutine Cal_Any_Point_Str_KesiYita_3D(c_Elem,i_N,Key_to_cal,
     &                              Key_CoorSys,
     &                              in_Kesi,
     &                              in_Yita,in_Zeta,i_G,c_DISP,
     &                              Sxx_GP,Syy_GP,Szz_GP,
     &                              Sxy_GP,Syz_GP,Sxz_GP)
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack_Common
      use Global_Crack_3D
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_DISP
      use Global_Cross
      use Global_Material
      use Global_Stress
      
      use Global_Inter_Location_Element_Stiff_Matrix_3D
      
      implicit none
      
      integer,intent(in)::c_Elem,i_N,Key_to_cal,Key_CoorSys,i_G
      real(kind=FT),intent(in)::in_Kesi,in_Yita,in_Zeta,c_DISP(Total_FD)
      real(kind=FT),intent(out)::Sxx_GP,Syy_GP,Szz_GP
      real(kind=FT),intent(out)::Sxy_GP,Syz_GP,Sxz_GP
      integer c_NN(8)
      real(kind=FT) c_D(6,6),U(MDOF_3D)
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),
     &                 c_Str(6),B(6,MDOF_3D),
     &                 tem_B(6,MDOF_3D)
      integer i_C,num_B,num_tem_B
      integer:: Location_ESM(MDOF_3D)
      integer::Location_ESM_C_Crack(MDOF_3D),
     &         Location_ESM_C_Cr_NoFEM(MDOF_3D)
      integer num_Loc_ESM_C_Cr_NoFEM
      integer num_Loc_ESM
      integer num_Loc_ESM_C_Crack
      real(kind=FT) c_T_Alpha,c_TStress(6)     
      integer mat_num
      real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
      real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)      
      real(kind=FT) c_Theta,Tem_Str(6),sTheta,cTheta,c_Stress(6)  
      integer c_POS_3D_c_Ele(8)
      

      
      if(Key_XA/=2) then
          mat_num = Elem_Mat(c_Elem)
          c_D     = D(Elem_Mat(c_Elem),:,:)  
          if (Material_Type(mat_num)==5)then
              Volume_Ratio = Material_Para_Added(mat_num,10)
              c_D_Comp = D_Comp(mat_num,1:6,1:6)
              T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
              TT_Matrix= TRANSPOSE(T_Matrix)
              Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
              Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
              c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
          endif
      else
          c_D     = Elem_D_XA(c_Elem,1:6,1:6)
      endif

      c_NN    = G_NN(1:8,c_Elem)
      c_X_NODES = G_X_NODES(1:8,c_Elem)
      c_Y_NODES = G_Y_NODES(1:8,c_Elem)    
      c_Z_NODES = G_Z_NODES(1:8,c_Elem)        
      
      Location_ESM(1:MDOF_3D)   = 0
      num_Loc_ESM           = 0
      if(num_Crack/=0)then
        do i_C =1,num_Crack
          c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
          call Location_Element_Stiff_Matrix_3D(c_Elem,i_C,
     &                                c_POS_3D_c_Ele(1:8),
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


      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.
      B(1:6,1:MDOF_3D) = ZR
      num_B = 0 
      if(num_Crack/=0)then
        do i_C =1,num_Crack 
          call Cal_B_Matrix_Crack_3D(in_Kesi,in_Yita,in_Zeta,
     &                      i_C,c_Elem,i_G,
     &                      c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,
     &                      tem_B,num_tem_B)     
     
          B(1:6,num_B+1:num_B+num_tem_B)=tem_B(1:6,1:num_tem_B)
          num_B = num_B + num_tem_B     
        end do
      endif
      
      
      
      if(Key_to_cal==1)then
         if (Key_CoorSys==1) then
             c_Str=MATMUL(MATMUL(c_D,B(1:6,1:num_Loc_ESM)),
     &                          U(1:num_Loc_ESM))
         elseif (Key_CoorSys==2) then
             c_Str=MATMUL(MATMUL(c_D,B(1:6,1:num_Loc_ESM)),
     &                          U(1:num_Loc_ESM))
             Tem_Str(1:6)=  c_Str(1:6) 
             c_Theta     =  Theta_Cartesian_to_Cylinder(c_Elem,i_N)
             sTheta      =  sin(c_Theta)
             cTheta      =  cos(c_Theta)
             c_Str(1)    =  Tem_Str(1)*cTheta**2 + Tem_Str(2)*sTheta**2+ 
     &                      TWO*Tem_Str(4)*sTheta*cTheta
             c_Str(2)    =  Tem_Str(1)*sTheta**2 + Tem_Str(2)*cTheta**2- 
     &                      TWO*Tem_Str(4)*sTheta*cTheta     
             c_Str(3)    =  Tem_Str(3)   
             c_Str(4)    =  (Tem_Str(2)-Tem_Str(1))*sTheta*cTheta  + 
     &                      Tem_Str(4)*(cTheta**2 - sTheta**2) 
             c_Str(5)    =  -Tem_Str(6)*sTheta + Tem_Str(5)*cTheta   
             c_Str(6)    =   Tem_Str(6)*cTheta + Tem_Str(5)*sTheta   
         endif         
         
      elseif(Key_to_cal==2)then
         if (Key_CoorSys==1) then
             c_Str=MATMUL(B(1:6,1:num_Loc_ESM),U(1:num_Loc_ESM))
         elseif (Key_CoorSys==2) then            
            Tem_Str=MATMUL(MATMUL(c_D,B(1:6,1:num_Loc_ESM)),
     &                          U(1:num_Loc_ESM))
            c_Theta     =  Theta_Cartesian_to_Cylinder(c_Elem,i_N)
            sTheta      =  sin(c_Theta)
            cTheta      =  cos(c_Theta)
            
            c_Str(1)    =  Tem_Str(1)*cTheta**2 + Tem_Str(2)*sTheta**2 + 
     &                     TWO*Tem_Str(4)*sTheta*cTheta
            c_Str(2)    =  Tem_Str(1)*sTheta**2 + Tem_Str(2)*cTheta**2 - 
     &                     TWO*Tem_Str(4)*sTheta*cTheta     
            c_Str(3)    =  Tem_Str(3)   
            c_Str(4)    =  (Tem_Str(2)-Tem_Str(1))*sTheta*cTheta  + 
     &                     Tem_Str(4)*(cTheta**2 - sTheta**2) 
            c_Str(5)    =  -Tem_Str(6)*sTheta + Tem_Str(5)*cTheta   
            c_Str(6)    =   Tem_Str(6)*cTheta + Tem_Str(5)*sTheta   
            c_Str(1:6) = MATMUL(D_for_cylindrical(mat_num,1:6,1:6),
     &                          c_Str)                
         endif
      endif


      if(Key_Thermal_Stress==1 .and. Key_to_cal==1)then
          if(Key_XA/=2) then
              c_T_Alpha = T_Alpha(Elem_Mat(c_Elem))
          else
              c_T_Alpha =  Elem_TEC_XA(c_Elem)
          endif

          c_TStress = c_T_Alpha*Elem_T_for_Stress(c_Elem)*
     &                MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
          c_Str(1) = c_Str(1) - c_TStress(1)
          c_Str(2) = c_Str(2) - c_TStress(2)  
          c_Str(3) = c_Str(3) - c_TStress(3)
          c_Str(4) = c_Str(4) - c_TStress(4)
          c_Str(5) = c_Str(5) - c_TStress(5)
          c_Str(6) = c_Str(6) - c_TStress(6)
      endif


      if(Key_InSitu_Strategy==4  .and. Key_to_cal==1)then
          c_Stress = MATMUL(c_D,[ InSitu_Strain_Gaus_xx(c_Elem,1),
     &                            InSitu_Strain_Gaus_yy(c_Elem,1),
     &                            InSitu_Strain_Gaus_zz(c_Elem,1),
     &                            InSitu_Strain_Gaus_xy(c_Elem,1),
     &                            InSitu_Strain_Gaus_yz(c_Elem,1),
     &                            InSitu_Strain_Gaus_xz(c_Elem,1)])   
          c_Str(1) = c_Str(1) - c_Stress(1)
          c_Str(2) = c_Str(2) - c_Stress(2)  
          c_Str(3) = c_Str(3) - c_Stress(3)
          c_Str(4) = c_Str(4) - c_Stress(4)
          c_Str(5) = c_Str(5) - c_Stress(5)
          c_Str(6) = c_Str(6) - c_Stress(6)     
      endif
      
      if(Key_InSitu_Strategy==2 .and. Key_to_cal==1)then
           c_Str(1)=c_Str(1)
     &      + sum(InSitu_Strs_Gaus_xx(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
           c_Str(2)=c_Str(2)
     &      + sum(InSitu_Strs_Gaus_yy(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
           c_Str(3)=c_Str(3)
     &      + sum(InSitu_Strs_Gaus_zz(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
           c_Str(4)=c_Str(4)
     &      + sum(InSitu_Strs_Gaus_xy(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
           c_Str(5)=c_Str(5)
     &      + sum(InSitu_Strs_Gaus_yz(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
           c_Str(6)=c_Str(6)
     &      + sum(InSitu_Strs_Gaus_xz(c_Elem,1:Num_Gauss_P_FEM_3D))/
     &        Num_Gauss_P_FEM_3D
      endif

      
      Sxx_GP = c_Str(1)
      Syy_GP = c_Str(2)
      Szz_GP = c_Str(3)
      Sxy_GP = c_Str(4)
      Syz_GP = c_Str(5)
      Sxz_GP = c_Str(6)      
      
      return 
      end SUBROUTINE Cal_Any_Point_Str_KesiYita_3D        
