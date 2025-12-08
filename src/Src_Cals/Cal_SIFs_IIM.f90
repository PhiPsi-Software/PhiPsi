 
subroutine Cal_SIFs_IIM(iter,Display_info,c_DISP)    

use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_Common
use Global_Material
use Global_HF
use Global_Dynamic

implicit none
integer,intent(in)::iter
real(kind=FT),intent(in)::c_DISP(Total_FD)
logical,intent(in)::Display_info
integer i_C,i_Tip,Tip_Elem,mat_num_Tip,i_JE,i_G,j_C
real(kind=FT) cr_Tip(2,2),omega_Tip(2),delta_L,x,y,c_E_Tip,c_v_Tip,c_G_Tip,k_Tip,&
    offset_delta,Ux,Uy,Dx,Dy,Disp_U(2),Disp_D(2),&
    delta_Disp_x,delta_Disp_y,delta_Disp_x_CRACK,delta_Disp_y_CRACK
real(kind=FT) Tip_x,Tip_y,r_JDomain
integer J_Elem(300),num_J_Elem
real(kind=FT) Q_Elem_Nodes(num_Elem,4)
integer c_Ele,c_mat_num
real(kind=FT)  c_E,c_v,c_G,k     
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
real(kind=FT) kesi_Enr(Num_Gauss_Points),yita_Enr(Num_Gauss_Points),weight_Enr(Num_Gauss_Points)
real(kind=FT) kesi_N_Enr(4),yita_N_Enr(4),weight_N_Enr(4)
real(kind=FT) kesi(Num_Gauss_Points),yita(Num_Gauss_Points),weight(Num_Gauss_Points)       
integer c_NN(4),c_Num_Gauss_Point
integer num_B,num_tem_B
integer:: Location_ESM(MDOF_2D)
integer Location_ESM_C_Cr_NoFEM(60)
integer num_Loc_ESM_C_Cr_NoFEM
integer num_Loc_ESM
integer::Location_ESM_C_Crack(80)
integer num_Loc_ESM_C_Crack     
real(kind=FT) c_thick,c_D(3,3),U(80),c_kesi,c_yita,B(3,80),tem_B(3,80),H(2,2)
real(kind=FT) detJ,J(2,2), Inverse_J(2,2)
real(kind=FT) N(2,8),dNdkesi(4,2),dNdx(4,2),Coor_AB(2,2), &
              c_Q(4),Q_Grad(2),QT(2,2),omega,Stress_Gauss(3),&
              Q_Grad_Local(2),H_Local(2,2),CalStress(2,2)
real(kind=FT) tem_loc_Stress(2,2),Local_Gauss_coor(2), &
               Global_Gauss(2),Tip_Coor(2),r,theta,K1,K2
real(kind=FT) SQR,CosT,SinT,CT2,ST2,C3T2,S3T2
real(kind=FT) AuxStress(2,2),AuxGradDisp(2,2),FAC_Disp1, &
               FAC_Disp2,FAC_Stress1,FAC_Stress2, &
               drdx,drdy,dtdx,dtdy,AuxStrain(2,2),du1dr, &
               du1dt,du2dr,du2dt,I1,I2,StrainEnergy,I(2), Kcalc(2),Eeff
integer Num_CalP,i_CalP,Ele_CalP,i_S
real(kind=FT) x_CalP,y_CalP,P_CalP,Distance,Global_CalP(2), &
               Ori_CalP,Local_CalP_coor(2),N_4(4),Q_CalP, &
               c_Pres_x,c_Pres_y,I_Press,L_HF_Elem,x_CalP_temp, &
               y_CalP_temp,crack_p1(2),crack_p2(2)
integer num_Inter,State,i_SS,c_CalP
real(kind=FT)  Inter(2,2),J_Inter(2),L_P,R_P,L_x,L_y,R_x,R_y,SS_x,SS_y,SS_p,SS_Ori,Global_SS(2)     

integer CalP_Num(Max_Num_Cr_CalP),num_Calp_inJ,Ele_SS

real(kind=FT) I_dyn,c_density,c_Accel(2),Accel(80),c_Q_Gauss

integer i_mode,ii,jj

if(Display_info.eqv..True.)then
  print *,'    Calculating SIFs of each crack......'   
end if      

r_JDomain  = 2.0999D0*Ave_Elem_L_Enrich

KI(1:Max_Num_Cr,1:2)  = ZR
KII(1:Max_Num_Cr,1:2) = ZR

if (Key_Integral_Sol  == 2)then
  call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,weight_Enr)
elseif (Key_Integral_Sol  == 3)then
  call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,kesi_Enr,yita_Enr,weight_Enr)
  Num_Gauss_Points = Num_Sub_Quads*4
endif

call Cal_Gauss_Points_QUAD(4,kesi_N_Enr,yita_N_Enr,weight_N_Enr)

LOOP_i_C: do i_C =1,num_Crack
  cr_Tip(1,1:2) = Cr_First_Tip(i_C,1:2)
  cr_Tip(2,1:2) = Cr_Second_Tip(i_C,1:2)
  omega_Tip(1) = Cr_First_Tip_Ori(i_C)
  omega_Tip(2) = Cr_Second_Tip_Ori(i_C)
  LOOP_i_Tip: do i_Tip = 1,2
      if ((Crack_Tip_Type(i_C,i_Tip) .ne. 1) .and.  &
          (Crack_Tip_Type(i_C,i_Tip) .ne. -2)) then
          I1 = ZR; I2 = ZR
          I(1:2) =ZR
          Tip_x = cr_Tip(i_Tip,1)
          Tip_y = cr_Tip(i_Tip,2)
          Tip_Coor = [Tip_x,Tip_y]
          omega = omega_Tip(i_Tip)
          QT(1,1:2) = [ cos(omega),sin(omega)]
          QT(2,1:2) = [-sin(omega),cos(omega)]
          call Cal_Ele_Num_by_Coors(Tip_x,Tip_y,Tip_Elem)
          mat_num_Tip = Elem_Mat(Tip_Elem)
          if (Material_Type(mat_num_Tip) .eq. 1)then
              c_E_Tip  = Material_Para(mat_num_Tip,1)
            c_v_Tip  = Material_Para(mat_num_Tip,2)
          else
              c_E_Tip  = Material_Para(mat_num_Tip,1)
            c_v_Tip  = Material_Para(mat_num_Tip,2)
          end if
          c_G_Tip  = c_E_Tip/TWO/(ONE+c_v_Tip)
          if (Key_Type_2D == 1) then
              k_Tip = (THR-c_v_Tip)/(ONE+c_v_Tip)
          elseif (Key_Type_2D == 2) then
              k_Tip = THR - FOU*c_v_Tip
          end if
          call Cal_JDomain_Element_of_Point(Tip_x,Tip_y,r_JDomain,J_Elem,num_J_Elem,Q_Elem_Nodes)
          LOOP_i_JE : do i_JE = 1,num_J_Elem
            c_Ele = J_Elem(i_JE)
            c_mat_num = Elem_Mat(c_Ele)
            if (Material_Type(c_mat_num) .eq. 1)then
              c_E  = Material_Para(c_mat_num,1)
              c_v  = Material_Para(c_mat_num,2)
              c_density = Material_Para(c_mat_num,3)       
            else
              c_E  = Material_Para(c_mat_num,1)
              c_v  = Material_Para(c_mat_num,2)
              c_density = Material_Para(c_mat_num,3)     
            end if
            c_G  = c_E/TWO/(ONE+c_v)
            if (Key_Type_2D == 1) then
              k = (THR-c_v)/(ONE+c_v)
            elseif (Key_Type_2D == 2) then
              k = THR - FOU*c_v
            end if  
            c_NN    = G_NN(:,c_Ele)
            c_X_NODES = G_X_NODES(:,c_Ele)
            c_Y_NODES = G_Y_NODES(:,c_Ele) 
            c_D     = D(Elem_Mat(c_Ele),:,:) 
            if(Key_Integral_Sol.eq.1)then
            elseif(Key_Integral_Sol.eq.2 .or.Key_Integral_Sol.eq.3)then
              if (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              else 
                  kesi(1:4)    = kesi_N_Enr
                  yita(1:4)    = yita_N_Enr
                  weight(1:4)  = weight_N_Enr
                  c_Num_Gauss_Point = 4
              end if 
            endif
            Location_ESM(1:MDOF_2D)  = 0
            num_Loc_ESM = 0    
            EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.
            do j_C =1,num_Crack 
              call Location_Element_Stiff_Matrix(c_Ele,j_C,c_POS(:,j_C), &
                                        Location_ESM_C_Crack,num_Loc_ESM_C_Crack, &
                                        Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
              Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack                     
            end do
            U(1:num_Loc_ESM)=c_DISP(Location_ESM(1:num_Loc_ESM))                   
            
            LOOP_i_G : do i_G = 1,c_Num_Gauss_Point
              c_kesi = kesi(i_G)                                               
              c_yita = yita(i_G)
              call Cal_Coor_by_KesiYita(c_kesi,c_yita,c_X_NODES,c_Y_NODES,Global_Gauss(1),Global_Gauss(2))
              call Cal_N_dNdkesi_J_detJ(c_kesi,c_yita,c_X_NODES,c_Y_NODES,detJ,J,N,dNdkesi)    
              call Matrix_Inverse_2x2(J,Inverse_J)
              dNdx = MATMUL(dNdkesi,Inverse_J)
              B(1:3,1:80) = ZR
              num_B = 0 
              do j_C =1,num_Crack 
                  call Cal_B_Matrix_Crack(c_kesi,c_yita,j_C,c_Ele,i_G,c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)                  
                  B(1:3,num_B+1:num_B+num_tem_B)=tem_B(1:3,1:num_tem_B)
                  num_B = num_B + num_tem_B                      
              end do
              H(1,1) =DOT_PRODUCT(B(1,1:num_B:2),U(1:num_B:2))
              H(1,2) =DOT_PRODUCT(B(2,2:num_B:2),U(1:num_B:2))
              H(2,1) =DOT_PRODUCT(B(1,1:num_B:2),U(2:num_B:2))
              H(2,2) =DOT_PRODUCT(B(2,2:num_B:2),U(2:num_B:2))
              c_Q = Q_Elem_Nodes(c_Ele,1:4)
              Q_Grad = MATMUL(c_Q,dNdx)
              Stress_Gauss=MATMUL(MATMUL(c_D,B(1:3,1:num_Loc_ESM)),U(1:num_Loc_ESM))
              
              Q_Grad_Local = MATMUL(QT,Q_Grad)
              H_Local = MATMUL(MATMUL(QT,H),transpose(QT))
              tem_loc_Stress(1,1) = Stress_Gauss(1)
              tem_loc_Stress(1,2) = Stress_Gauss(3)
              tem_loc_Stress(2,1) = Stress_Gauss(3)
              tem_loc_Stress(2,2) = Stress_Gauss(2)                     
              CalStress = MATMUL(MATMUL(QT,tem_loc_Stress),transpose(QT))
              Local_Gauss_coor=MATMUL(QT,Global_Gauss-Tip_Coor)
              r = sqrt(Local_Gauss_coor(1)**2 + Local_Gauss_coor(2)**2)
              theta = atan2(Local_Gauss_coor(2),Local_Gauss_coor(1))
              if (Material_Type(c_mat_num) .eq. 1) then
                K1 = ONE; K2 = ONE
                SQR  = sqrt(r)
                CosT = cos(theta)
                SinT = sin(theta)
                CT2  = cos(theta/TWO)
                ST2  = sin(theta/TWO)
                C3T2 = cos(THR*theta/TWO)
                S3T2 = sin(THR*theta/TWO)
                  drdx =  CosT
                  drdy =  SinT
                  dtdx = -SinT/r
                  dtdy =  CosT/r
                  FAC_Stress1 = sqrt(ONE/(2*pi))            
                  FAC_Stress2 = FAC_Stress1 
                  FAC_Disp1 = sqrt(ONE/(TWO*pi))/(TWO*c_G)    
                        FAC_Disp2 = FAC_Disp1
                  AuxStress(1:2,1:2)   = ZR
                AuxGradDisp(1:2,1:2) = ZR
                AuxStrain(1:2,1:2) = ZR   
                LOOP_i_Mode:do i_Mode=1,2
                  if(i_Mode==1)then
                    AuxStress(1,1) = K1*FAC_Stress1/SQR*CT2*(ONE-ST2*S3T2)
                    AuxStress(2,2) = K1*FAC_Stress1/SQR*CT2*(ONE+ST2*S3T2)
                    AuxStress(1,2) = K1*FAC_Stress1/SQR*ST2*CT2*C3T2
                    AuxStress(2,1) = AuxStress(1,2)
                    du1dr = K1*FAC_Disp1*HLF/SQR*CT2*(k-CosT)
                    du1dt = K1*FAC_Disp1*SQR*(-HLF*ST2*(k-CosT)+CT2*SinT)
                    du2dr = K1*FAC_Disp1*HLF/SQR*ST2*(k-CosT)
                    du2dt = K1*FAC_Disp1*SQR*(HLF*CT2*(k-CosT)+ST2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                    AuxStrain(1,1) = AuxGradDisp(1,1)
                    AuxStrain(1,2) = HLF*(AuxGradDisp(1,2)+AuxGradDisp(2,1))
                    AuxStrain(2,1) = AuxStrain(1,2)
                    AuxStrain(2,2) = AuxGradDisp(2,2)                           
                  elseif(i_Mode==2)then
                    AuxStress(1,1) = -K2*FAC_Stress2/SQR*ST2*(TWO+CT2*C3T2)
                    AuxStress(2,2) =  K2*FAC_Stress2/SQR*ST2*CT2*C3T2
                    AuxStress(1,2) =  K2*FAC_Stress2/SQR*CT2*(ONE-ST2*S3T2)
                    AuxStress(2,1) =  AuxStress(1,2)
                    du1dr =K2*FAC_Disp2*HLF/SQR*ST2*(k+TWO+CosT)
                    du1dt =  K2*FAC_Disp2*SQR*(HLF*CT2*(k+2+CosT)-ST2*SinT)
                    du2dr = -K2*FAC_Disp2*HLF/SQR*CT2*(k-TWO+CosT)
                    du2dt = -K2*FAC_Disp2*SQR*(-HLF*ST2*(k-TWO+CosT)-CT2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                    AuxStrain(1,1) = AuxGradDisp(1,1)
                    AuxStrain(1,2) = HLF*(AuxGradDisp(1,2)+AuxGradDisp(2,1))
                    AuxStrain(2,1) = AuxStrain(1,2)
                    AuxStrain(2,2) = AuxGradDisp(2,2)                            
                  end if
                I1 = (CalStress(1,1)*AuxGradDisp(1,1)+ &
                               CalStress(2,1)*AuxGradDisp(2,1))*  &
                                Q_Grad_Local(1)+ (CalStress(1,2)*AuxGradDisp(1,1)+ &
                                 CalStress(2,2)*AuxGradDisp(2,1))*Q_Grad_Local(2)
                I2 = (AuxStress(1,1)*H_Local(1,1)+ &
                                AuxStress(2,1)*H_Local(2,1))*Q_Grad_Local(1)+ &
                               (AuxStress(1,2)*H_Local(1,1)+AuxStress(2,2)*H_Local(2,1))*Q_Grad_Local(2)
                          StrainEnergy = ZR
                do ii = 1,2
                    do jj = 1,2
                      StrainEnergy = StrainEnergy+CalStress(ii,jj)*AuxStrain(ii,jj)
                    end do
                end do
                  I(i_mode) = I(i_mode)+(I1+I2-StrainEnergy*Q_Grad_Local(1))* detJ*weight(i_G)  

                  if (Key_Gravity==1) then
                      call Cal_N(c_Kesi,c_Yita,N)
                      N_4(1:4) =  [N(1,1),N(1,3),N(1,5),N(1,7)]
                      c_Q_Gauss = DOT_PRODUCT(c_Q,N_4)       
                      c_Accel(1)=g_X_Y_Z(1)
                      c_Accel(2)=g_X_Y_Z(2)
                      c_Accel   = MATMUL(QT,c_Accel)
                      I_dyn = c_density*(c_Accel(1)*AuxGradDisp(1,1)+c_Accel(2)*AuxGradDisp(2,1))*c_Q_Gauss
                     I(i_mode)=I(i_mode)-I_dyn*detJ*weight(i_G)  
                  endif
                  

                  if ((Key_Analysis_Type==2 .or. Key_Analysis_Type==6)) then
                      call Cal_N(c_Kesi,c_Yita,N)
                      N_4(1:4) =  [N(1,1),N(1,3),N(1,5),N(1,7)]
                      c_Q_Gauss = DOT_PRODUCT(c_Q,N_4)       
                      if (Key_Analysis_Type==2) then
                        c_Accel(1)=IDy_ACCL(c_NN(1)*2-1)*N_4(1)+ IDy_ACCL(c_NN(2)*2-1)*N_4(2)+  &
                                   IDy_ACCL(c_NN(3)*2-1)*N_4(3)+ IDy_ACCL(c_NN(4)*2-1)*N_4(4) 
                        c_Accel(2)=IDy_ACCL(c_NN(1)*2)*N_4(1) + IDy_ACCL(c_NN(2)*2)*N_4(2) +  &
                                   IDy_ACCL(c_NN(3)*2)*N_4(3) + IDy_ACCL(c_NN(4)*2)*N_4(4)
                      elseif(Key_Analysis_Type==6) then
                        c_Accel(1)=EDy_ACCL(c_NN(1)*2-1)*N_4(1)+ EDy_ACCL(c_NN(2)*2-1)*N_4(2)+ &
                                   EDy_ACCL(c_NN(3)*2-1)*N_4(3)+ EDy_ACCL(c_NN(4)*2-1)*N_4(4) 
                        c_Accel(2)=EDy_ACCL(c_NN(1)*2)*N_4(1) + EDy_ACCL(c_NN(2)*2)*N_4(2) + &
                                   EDy_ACCL(c_NN(3)*2)*N_4(3) + EDy_ACCL(c_NN(4)*2)*N_4(4)      
                      endif
                      c_Accel   = MATMUL(QT,c_Accel)

                      I_dyn = c_density*(c_Accel(1)*AuxGradDisp(1,1)+c_Accel(2)*AuxGradDisp(2,1))*c_Q_Gauss
                      I(i_mode)=I(i_mode)+I_dyn*detJ*weight(i_G)  
                  endif
                  
                end do LOOP_i_Mode
              else
                  print *,'WORK TO BE DONE HERE in Cal_SIFs_IIM'
              endif
            end do LOOP_i_G
          end do LOOP_i_JE

          if (((Key_Analysis_Type==3 .or. Key_Analysis_Type==4 &
                .or.Key_Analysis_Type==5)                &
               .and.(Cracks_HF_State(i_C)==1))     .or.     &
              (Key_Analysis_Type==1 .and.                   &
               Key_Crack_Inner_Pressure==1 .and. &
               Crack_Pressure(i_C) >Tol_10)    ) then  
           do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
              call Tool_Intersection_Line_and_Circle(Tip_x, &
                              Tip_y,r_JDomain,crack_p1,crack_p2,num_Inter,State,Inter)
              if (num_Inter ==1 .and. State==2) then
                  J_Inter(1:2) = Inter(1,1:2)
                  exit
              end if    
           end do
           
           Num_CalP = Cracks_CalP_Num(i_C)
           num_Calp_inJ = 0
           CalP_Num(1:Max_Num_Cr_CalP) = 0
           do i_CalP=1,Num_CalP
               x_CalP = Cracks_CalP_Coors(i_C,i_CalP,1)
               y_CalP = Cracks_CalP_Coors(i_C,i_CalP,2)
               Global_CalP =[x_CalP,y_CalP]
               Distance = sqrt((Tip_x-x_CalP)**2 +(Tip_y-y_CalP)**2)
               if(Distance <= r_JDomain) then
                   num_Calp_inJ = num_Calp_inJ + 1
                   CalP_Num(num_Calp_inJ) = i_CalP
               end if
           end do          
           
           do i_SS=1,num_Calp_inJ
              c_CalP = CalP_Num(i_SS)  
              L_P =  Cracks_CalP_Pres(i_C,c_CalP-1) 
              R_P =  Cracks_CalP_Pres(i_C,c_CalP)   
              if(Key_Analysis_Type ==3 .or. Key_Analysis_Type ==4) then                      
                  L_P=L_P+Cracks_CalP_Remo_Strs(i_C,c_CalP-1)
                  R_P=R_P+Cracks_CalP_Remo_Strs(i_C,c_CalP)
              endif
              L_x =  Cracks_CalP_Coors(i_C,c_CalP-1,1)
              L_y =  Cracks_CalP_Coors(i_C,c_CalP-1,2)
              R_x =  Cracks_CalP_Coors(i_C,c_CalP,1)
              R_y =  Cracks_CalP_Coors(i_C,c_CalP,2)
              L_HF_Elem = sqrt((L_x-R_x)**2+(L_y-R_y)**2)
              SS_x = (L_x+R_x)/TWO
              SS_y = (L_y+R_y)/TWO
              SS_P = (L_P+R_P)/TWO
              SS_Ori = Cracks_CalP_Orient(i_C,c_CalP)                            
              Global_SS =[SS_x,SS_y]
              call Cal_Ele_Num_by_Coors(SS_x,SS_y,Ele_SS)
              c_mat_num = Elem_Mat(Ele_SS)
              if (Material_Type(c_mat_num).eq. 1  )then
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              else
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              end if
              c_G  = c_E/TWO/(ONE+c_v)
              if (Key_Type_2D == 1) then
                k = (THR-c_v)/(ONE+c_v)
              elseif (Key_Type_2D == 2) then
                k = THR - FOU*c_v
              end if                     
              Local_CalP_coor=MATMUL(QT,Global_SS-Tip_Coor)
                r = sqrt(Local_CalP_coor(1)**2 + Local_CalP_coor(2)**2)
              theta =atan2(Local_CalP_coor(2),Local_CalP_coor(1))
              c_NN    = G_NN(:,Ele_SS)
              c_X_NODES = G_X_NODES(:,Ele_SS)
              c_Y_NODES = G_Y_NODES(:,Ele_SS)    
              call Cal_KesiYita_by_Coor(Global_SS,Ele_SS, c_Kesi,c_Yita)
              call Cal_N(c_Kesi,c_Yita,N)
              N_4(1:4) =  [N(1,1),N(1,3),N(1,5),N(1,7)]
              c_Q = Q_Elem_Nodes(Ele_SS,1:4)
              Q_CalP = DOT_PRODUCT(c_Q,N_4)
              c_Pres_x = -sin(Ori_CalP)*SS_P
              c_Pres_y =  cos(Ori_CalP)*SS_P
              K1 = ONE
              K2 = ONE
              SQR  = sqrt(r)
              CosT = cos(theta)
              SinT = sin(theta)
              CT2  = cos(theta/TWO)
              ST2  = sin(theta/TWO)
              C3T2 = cos(THR*theta/TWO)
              S3T2 = sin(THR*theta/TWO)
              drdx =  CosT
              drdy =  SinT
              dtdx = -SinT/r
              dtdy =  CosT/r
              FAC_Disp1 = sqrt(ONE/(TWO*pi))/(TWO*c_G)    
              FAC_Disp2 = FAC_Disp1
              AuxGradDisp(1:2,1:2) = ZR  
              do i_Mode=1,2
                 if(i_Mode==1)then
                    du1dr = K1*FAC_Disp1*HLF/SQR*CT2*(k-CosT)
                    du1dt = K1*FAC_Disp1*SQR* (-HLF*ST2*(k-CosT)+CT2*SinT)
                    du2dr = K1*FAC_Disp1*HLF/SQR*ST2*(k-CosT)
                    du2dt = K1*FAC_Disp1*SQR*(HLF*CT2*(k-CosT)+ST2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                  elseif(i_Mode==2) then
                    du1dr =K2*FAC_Disp2*HLF/SQR*ST2*(k+TWO+CosT)
                    du1dt =  K2*FAC_Disp2*SQR*(HLF*CT2*(k+2+CosT)-ST2*SinT)
                    du2dr = -K2*FAC_Disp2*HLF/SQR*CT2*(k-TWO+CosT)
                    du2dt = -K2*FAC_Disp2*SQR*(-HLF*ST2*(k-TWO+CosT)-CT2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                 endif
               I_Press = TWO*(abs(c_Pres_x)*abs(AuxGradDisp(1,1)) +  &
                      abs(c_Pres_y)*abs(AuxGradDisp(2,1)))*Q_CalP*L_HF_Elem
               I(i_mode) = I(i_mode) + I_Press                      
              end do                  
           end do
          end if    
          
          if (Key_Analysis_Type == 1 .and. Key_TipEnrich     == 4)  then
           do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
              call Tool_Intersection_Line_and_Circle(Tip_x,Tip_y,r_JDomain,crack_p1,crack_p2, &
                              num_Inter,State,Inter)
              if (num_Inter ==1 .and. State==2) then
                  J_Inter(1:2) = Inter(1,1:2)
                  exit
              end if    
           end do
           Num_CalP = Cracks_CalP_Num(i_C)
           num_Calp_inJ = 0
           CalP_Num(1:Max_Num_Cr_CalP) = 0
           do i_CalP=1,Num_CalP
               x_CalP = Cracks_CalP_Coors(i_C,i_CalP,1)
               y_CalP = Cracks_CalP_Coors(i_C,i_CalP,2)
               Global_CalP =[x_CalP,y_CalP]
               Distance = sqrt((Tip_x-x_CalP)**2 +(Tip_y-y_CalP)**2)
               if(Distance <= r_JDomain) then
                   num_Calp_inJ = num_Calp_inJ + 1
                   CalP_Num(num_Calp_inJ) = i_CalP
               end if
           end do          
           do i_SS=1,num_Calp_inJ
              c_CalP = CalP_Num(i_SS)  
              L_P =  Cracks_CalP_Pres(i_C,c_CalP-1) 
              R_P =  Cracks_CalP_Pres(i_C,c_CalP)   
              L_x =  Cracks_CalP_Coors(i_C,c_CalP-1,1)
              L_y =  Cracks_CalP_Coors(i_C,c_CalP-1,2)
              R_x =  Cracks_CalP_Coors(i_C,c_CalP,1)
              R_y =  Cracks_CalP_Coors(i_C,c_CalP,2)
              L_HF_Elem = sqrt((L_x-R_x)**2+(L_y-R_y)**2)
              SS_x = (L_x+R_x)/TWO
              SS_y = (L_y+R_y)/TWO
              SS_P = (L_P+R_P)/TWO
              SS_Ori = Cracks_CalP_Orient(i_C,c_CalP)                            
              Global_SS =[SS_x,SS_y]
              call Cal_Ele_Num_by_Coors(SS_x,SS_y,Ele_SS)
              c_mat_num = Elem_Mat(Ele_SS)
              if (Material_Type(c_mat_num) .eq. 1  )then
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              else
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              end if
              c_G  = c_E/TWO/(ONE+c_v)
              if (Key_Type_2D == 1) then
                k = (THR-c_v)/(ONE+c_v)
              elseif (Key_Type_2D == 2) then
                k = THR - FOU*c_v
              end if                     
              Local_CalP_coor=MATMUL(QT,Global_SS-Tip_Coor)
              r = sqrt(Local_CalP_coor(1)**2 + Local_CalP_coor(2)**2)
              theta =atan2(Local_CalP_coor(2),Local_CalP_coor(1))
              c_NN    = G_NN(:,Ele_SS)
              c_X_NODES = G_X_NODES(:,Ele_SS)
              c_Y_NODES = G_Y_NODES(:,Ele_SS)    
              call Cal_KesiYita_by_Coor(Global_SS,Ele_SS,c_Kesi,c_Yita)
              call Cal_N(c_Kesi,c_Yita,N)
              N_4(1:4) =  [N(1,1),N(1,3),N(1,5),N(1,7)]
              c_Q    = Q_Elem_Nodes(Ele_SS,1:4)
              Q_CalP = DOT_PRODUCT(c_Q,N_4)
              c_Pres_x = Cracks_CalP_Tractions(i_C,c_CalP,1)
              c_Pres_y = Cracks_CalP_Tractions(i_C,c_CalP,2)
              
              
              K1 = ONE
              K2 = ONE
              SQR  = sqrt(r)
              CosT = cos(theta)
              SinT = sin(theta)
              CT2  = cos(theta/TWO)
              ST2  = sin(theta/TWO)
              C3T2 = cos(THR*theta/TWO)
              S3T2 = sin(THR*theta/TWO)
              drdx =  CosT
              drdy =  SinT
              dtdx = -SinT/r
              dtdy =  CosT/r
              FAC_Disp1 = sqrt(ONE/(TWO*pi))/(TWO*c_G)    
              FAC_Disp2 = FAC_Disp1
              AuxGradDisp(1:2,1:2) = ZR  
              do i_Mode=1,2
                 if(i_Mode==1)then
                    du1dr = K1*FAC_Disp1*HLF/SQR*CT2*(k-CosT)
                    du1dt = K1*FAC_Disp1*SQR*(-HLF*ST2*(k-CosT)+CT2*SinT)
                    du2dr = K1*FAC_Disp1*HLF/SQR*ST2*(k-CosT)
                    du2dt = K1*FAC_Disp1*SQR*(HLF*CT2*(k-CosT)+ST2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                 elseif(i_Mode==2) then
                    du1dr =K2*FAC_Disp2*HLF/SQR*ST2*(k+TWO+CosT)
                    du1dt =  K2*FAC_Disp2*SQR*(HLF*CT2*(k+2+CosT)-ST2*SinT)
                    du2dr = -K2*FAC_Disp2*HLF/SQR*CT2*(k-TWO+CosT)
                    du2dt = -K2*FAC_Disp2*SQR*(-HLF*ST2*(k-TWO+CosT)-CT2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                 endif
               I_Press = TWO* (abs(c_Pres_x)*abs(AuxGradDisp(1,1)) + &
                      abs(c_Pres_y)*abs(AuxGradDisp(2,1)))*Q_CalP*L_HF_Elem  
                 I(i_mode) = I(i_mode) + I_Press            
              end do                  
           end do
          end if    
          
          
          if (Key_Analysis_Type == 1 .and.Key_Contact    == 1)  then
           do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
              call Tool_Intersection_Line_and_Circle(Tip_x, &
                              Tip_y,r_JDomain,crack_p1,crack_p2, &
                              num_Inter,State,Inter)
              if (num_Inter ==1 .and. State==2) then
                  J_Inter(1:2) = Inter(1,1:2)
                  exit
              end if    
           end do
           Num_CalP = Cracks_CalP_Num(i_C)
           num_Calp_inJ = 0
           CalP_Num(1:Max_Num_Cr_CalP) = 0
           do i_CalP=1,Num_CalP
               x_CalP = Cracks_CalP_Coors(i_C,i_CalP,1)
               y_CalP = Cracks_CalP_Coors(i_C,i_CalP,2)
               Global_CalP =[x_CalP,y_CalP]
               Distance = sqrt((Tip_x-x_CalP)**2 + (Tip_y-y_CalP)**2)
               if(Distance <= r_JDomain) then
                   num_Calp_inJ = num_Calp_inJ + 1
                   CalP_Num(num_Calp_inJ) = i_CalP
               end if
           end do          
           do i_SS=1,num_Calp_inJ
              c_CalP = CalP_Num(i_SS)  
              L_x =  Cracks_CalP_Coors(i_C,c_CalP-1,1)
              L_y =  Cracks_CalP_Coors(i_C,c_CalP-1,2)
              R_x =  Cracks_CalP_Coors(i_C,c_CalP,1)
              R_y =  Cracks_CalP_Coors(i_C,c_CalP,2)
              L_HF_Elem = sqrt((L_x-R_x)**2+(L_y-R_y)**2)
              SS_x = (L_x+R_x)/TWO
              SS_y = (L_y+R_y)/TWO
              SS_Ori = Cracks_CalP_Orient(i_C,c_CalP)                            
              Global_SS =[SS_x,SS_y]
              call Cal_Ele_Num_by_Coors(SS_x,SS_y,Ele_SS)
              c_mat_num = Elem_Mat(Ele_SS)
              if (Material_Type(c_mat_num) .eq. 1  )then
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              else
                c_E  = Material_Para(c_mat_num,1)
                c_v  = Material_Para(c_mat_num,2)
              end if
              c_G  = c_E/TWO/(ONE+c_v)
              if (Key_Type_2D == 1) then
                k = (THR-c_v)/(ONE+c_v)
              elseif (Key_Type_2D == 2) then
                k = THR - FOU*c_v
              end if                     
              Local_CalP_coor=MATMUL(QT,Global_SS-Tip_Coor)
              r = sqrt(Local_CalP_coor(1)**2 + Local_CalP_coor(2)**2)
              theta =atan2(Local_CalP_coor(2),Local_CalP_coor(1))
              c_NN    = G_NN(:,Ele_SS)
              c_X_NODES = G_X_NODES(:,Ele_SS)
              c_Y_NODES = G_Y_NODES(:,Ele_SS)    
              call Cal_KesiYita_by_Coor(Global_SS,Ele_SS,c_Kesi,c_Yita)
              call Cal_N(c_Kesi,c_Yita,N)
              N_4(1:4) =  [N(1,1),N(1,3),N(1,5),N(1,7)]
              c_Q    = Q_Elem_Nodes(Ele_SS,1:4)
              Q_CalP = DOT_PRODUCT(c_Q,N_4)
              c_Pres_x= sum(Cracks_CalP_Contact_Force_x(i_C,(c_CalP-1):c_CalP))/TWO    
              c_Pres_y= sum(Cracks_CalP_Contact_Force_y(i_C,(c_CalP-1):c_CalP))/TWO                         
              
              
              
              K1 = ONE
              K2 = ONE
              SQR  = sqrt(r)
              CosT = cos(theta)
              SinT = sin(theta)
              CT2  = cos(theta/TWO)
              ST2  = sin(theta/TWO)
              C3T2 = cos(THR*theta/TWO)
              S3T2 = sin(THR*theta/TWO)
              drdx =  CosT
              drdy =  SinT
              dtdx = -SinT/r
              dtdy =  CosT/r
              FAC_Disp1 = sqrt(ONE/(TWO*pi))/(TWO*c_G)    
              FAC_Disp2 = FAC_Disp1
              AuxGradDisp(1:2,1:2) = ZR  
              do i_Mode=1,2
                 if(i_Mode==1)then
                    du1dr = K1*FAC_Disp1*HLF/SQR*CT2*(k-CosT)
                    du1dt = K1*FAC_Disp1*SQR*(-HLF*ST2*(k-CosT)+CT2*SinT)
                    du2dr = K1*FAC_Disp1*HLF/SQR*ST2*(k-CosT)
                    du2dt = K1*FAC_Disp1*SQR*(HLF*CT2*(k-CosT)+ST2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                 elseif(i_Mode==2) then
                    du1dr =K2*FAC_Disp2*HLF/SQR*ST2*(k+TWO+CosT)
                    du1dt =  K2*FAC_Disp2*SQR*(HLF*CT2*(k+2+CosT)-ST2*SinT)
                    du2dr = -K2*FAC_Disp2*HLF/SQR*CT2*(k-TWO+CosT)
                    du2dt = -K2*FAC_Disp2*SQR*(-HLF*ST2*(k-TWO+CosT)-CT2*SinT)
                    AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx
                    AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy
                    AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx
                    AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy
                 endif
               I_Press = TWO*(abs(c_Pres_x)*abs(AuxGradDisp(1,1)) + &
                      abs(c_Pres_y)*abs(AuxGradDisp(2,1)))*Q_CalP*L_HF_Elem   
                 I(i_mode) = I(i_mode) + I_Press            
              end do                  
           end do
          end if   
          if (Material_Type(c_mat_num) .eq. 1) then
                if (Key_Type_2D == 1) then
                  Eeff = c_E_Tip
                elseif (Key_Type_2D == 2)then
                  Eeff = c_E_Tip/(ONE-c_v_Tip**2)                                              
                end if
               Kcalc   = I*Eeff/TWO
               KI(i_C,i_Tip)   = Kcalc(1)
               KII(i_C,i_Tip)  = Kcalc(2)
          else
              print *,'WORK TO BE DONE HERE in Cal_SIFs_IIM'
          endif
      end if
  end do LOOP_i_Tip
end do LOOP_i_C
RETURN
end subroutine Cal_SIFs_IIM