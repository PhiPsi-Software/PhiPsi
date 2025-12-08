 
subroutine Cal_Contact_Resid_3D(iter,ifra,Counter_Iter,i_NR_P,          &
         c_Total_Freedom,c_num_freeDOF,F_U,                             &
         c_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI)
      

use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Material
use Global_XFEM_Elements
use OMP_LIB

use Global_Inter_Location_Element_Stiff_Matrix_3D
      
implicit none
integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
integer, intent(in)::c_Total_Freedom,c_num_freeDOF
real(kind=FT),intent(in)::F_U(c_Total_Freedom),c_DISP(c_Total_Freedom)
integer, intent(in)::freeDOF(c_Total_Freedom)
real(kind=FT),intent(in)::&
            PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D),   &
            PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D),   &
            PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)     
real(kind=FT),intent(out)::R_PSI(c_Total_Freedom)
integer i_E,i_N,i_C,i_G
real(kind=FT) U(MDOF_3D)
integer c_Num_Gauss_Point
real(kind=FT) c_kesi,c_yita,c_eta,c_Stress(6),B(6,MDOF_3D),tem_B(6,MDOF_3D)
integer num_B,num_tem_B
integer:: Location_ESM(MDOF_3D)
integer Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
integer num_Loc_ESM
integer::Location_ESM_C_Crack(MDOF_3D)
integer num_Loc_ESM_C_Crack

real(kind=FT) detJ,c_Inter_Force(MDOF_3D)
integer i_CT_Elem
integer Local_W(60),num_Local_W 
real(kind=FT) N_W(3,60),N(3,24)
real(kind=FT) c_PC_Gauss(3)
integer ref_elem

real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)  
real(kind=FT) kesi_Enr(Num_Gau_Points_3D),      &
              yita_Enr(Num_Gau_Points_3D),&
               zeta_Enr(Num_Gau_Points_3D),&
               weight_Enr(Num_Gau_Points_3D)    
real(kind=FT) kesi_No_Enr(Num_Gauss_P_FEM_3D),     &
              yita_No_Enr(Num_Gauss_P_FEM_3D),&
               zeta_No_Enr(Num_Gauss_P_FEM_3D),&
               weight_No_Enr(Num_Gauss_P_FEM_3D)            
real(kind=FT) kesi(Num_Gau_Points_3D),yita(Num_Gau_Points_3D),&
              zeta(Num_Gau_Points_3D),weight(Num_Gau_Points_3D)
real(kind=FT) kesi_Contact(27),yita_Contact(27),zeta_Contact(27),weight_Contact(27)
real(kind=FT) c_Center(3),ori_n(3),Flu_Ele_Area
integer Solid_El      
real(kind=FT) tem_n(8),N_HF(3)
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3)
real(kind=FT) Tip_T_Matrx(3,3),r_Gauss,theta_Gauss     
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)  
integer mat_num      
real(kind=FT),ALLOCATABLE::R_PSI_thread(:,:)
integer i_Thread,c_Thread,max_threads
integer c_Cr_Location
integer c_POS_3D_c_Ele(8)

      
R_PSI(1:c_Total_Freedom)= ZR
call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_No_Enr,yita_No_Enr,zeta_No_Enr,weight_No_Enr)

call Cal_Gauss_Points_3D_8nodes(27,kesi_Contact,yita_Contact,zeta_Contact,weight_Contact) 

kesi(1:Num_Gau_Points_3D)=ZR
yita(1:Num_Gau_Points_3D)=ZR
zeta(1:Num_Gau_Points_3D)=ZR
weight(1:Num_Gau_Points_3D)=ZR
max_threads = omp_get_max_threads()
ALLOCATE(R_PSI_thread(c_Total_Freedom,max_threads))



EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)=.false.
R_PSI_thread(1:c_Total_Freedom,1:max_threads)= ZR 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_Thread,                        &
!$OMP             i_E,mat_num,c_D,Volume_Ratio,                         &
!$OMP             c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,Rot_c_D_Comp,   &
!$OMP             c_Y_NODES,c_Z_NODES,c_NN,kesi,yita,zeta,weight,       &
!$OMP             Location_ESM_C_Crack,num_Loc_ESM_C_Crack,             &
!$OMP             Location_ESM_C_Cr_NoFEM,                              &
!$OMP             c_Num_Gauss_Point,i_C,Location_ESM,num_Loc_ESM,B,     & 
!$OMP             num_Loc_ESM_C_Cr_NoFEM,i_G,c_Stress,c_Inter_Force,    &
!$OMP             num_B,tem_B,num_tem_B,detJ,U,c_POS_3D_c_Ele)    
c_Thread = omp_get_thread_num()+1
!$OMP do SCHEDULE(STATIC) 
do i_E = 1,Num_Elem
  
  if(Key_XA/=2) then
    mat_num = Elem_Mat(i_E)
    c_D(1:6,1:6) = D(Elem_Mat(i_E),1:6,1:6)     
    if (Material_Type(mat_num)==5)then
          Volume_Ratio = Material_Para_Added(mat_num,10)
          c_D_Comp = D_Comp(mat_num,1:6,1:6)
          T_Matrix = Ele_ComMat_RotMatrix(i_E,1:6,1:6)
          TT_Matrix= TRANSPOSE(T_Matrix)
          Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
          Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
          c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
    endif  
  else
      c_D       = Elem_D_XA(i_E,1:6,1:6)
  endif
  
  c_NN(1:8)      = G_NN(1:8,i_E)
  c_X_NODES(1:8) = G_X_NODES(1:8,i_E)
  c_Y_NODES(1:8) = G_Y_NODES(1:8,i_E)    
  c_Z_NODES(1:8) = G_Z_NODES(1:8,i_E)
  
  Location_ESM(1:MDOF_3D)   = 0
  num_Loc_ESM               = 0
  if(num_Crack/=0)then
    do i_C =1,num_Crack
      c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
      call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, & 
                                num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
      Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) =  &
               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
    end do
  endif
  U(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))    
  
  if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
      

      kesi(1:27)    = kesi_Contact
      yita(1:27)    = yita_Contact
      zeta(1:27)    = zeta_Contact
      weight(1:27)  = weight_Contact
      c_Num_Gauss_Point = 27

  else 
      kesi(1:Num_Gauss_P_FEM_3D)    = kesi_No_Enr
      yita(1:Num_Gauss_P_FEM_3D)    = yita_No_Enr
      zeta(1:Num_Gauss_P_FEM_3D)    = zeta_No_Enr
      weight(1:Num_Gauss_P_FEM_3D)  = weight_No_Enr
      c_Num_Gauss_Point = Num_Gauss_P_FEM_3D
  end if
  
  do i_G = 1,c_Num_Gauss_Point  
      call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)    
      B(1:6,1:MDOF_3D) = ZR
      num_B = 0 
      if(num_Crack/=0)then
        do i_C =1,num_Crack 
        
          c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
          if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
          
          call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),i_C,i_E,i_G,  &
                       c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,tem_B,num_tem_B)
          B(1:6,num_B+1:num_B+num_tem_B) =  tem_B(1:6,1:num_tem_B)
          num_B = num_B + num_tem_B
        end do
      endif
          
      c_Stress = MATMUL(MATMUL(c_D,B(1:6,1:num_Loc_ESM)),U(1:num_Loc_ESM))
      c_Inter_Force(1:num_Loc_ESM) = MATMUL(transpose(B(1:6,1:num_Loc_ESM)),c_Stress)
      R_PSI_thread(Location_ESM(1:num_Loc_ESM),c_Thread) = R_PSI_thread(Location_ESM(1:num_Loc_ESM),c_Thread)  +  &
                                                           c_Inter_Force(1:num_Loc_ESM)*weight(i_G)*detJ     
  end do
enddo
!$omp end do        
!$omp end parallel

DO i_Thread = 1,omp_get_max_threads()
    R_PSI(:)  = R_PSI(:)  + R_PSI_thread(:,i_Thread)
ENDDO

R_PSI =  R_PSI - F_U


R_PSI_thread(1:c_Total_Freedom,1:max_threads)= ZR 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_Thread,                        &
!$OMP             i_C,i_CT_Elem,N_W,c_Center,                           &
!$OMP             ori_n,Solid_El,Flu_Ele_Area,c_PC_Gauss,               &
!$OMP             c_kesi,c_yita,c_eta,c_NN,N,tem_n,N_HF,num_Local_W,    &
!$OMP             i_N,Local_W,ref_elem,BaseLine_A,BaseLine_B,           &
!$OMP             BaseLine_Mid,BaseLine_x_Vec,BaseLine_y_Vec,           &
!$OMP             BaseLine_z_Vec,Tip_T_Matrx,Gauss_Coor_Local,          &
!$OMP             r_Gauss,theta_Gauss,c_Cr_Location)   
c_Thread = omp_get_thread_num()+1
!$OMP do SCHEDULE(STATIC)   
do i_C=1,num_Crack
      do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
          N_W(1:3,1:60)  = ZR   
          c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
          ori_n    = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
          Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)     
          Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_CT_Elem)  
          c_PC_Gauss(1) =  PC_Gauss_x(i_C,i_CT_Elem)
          c_PC_Gauss(2) =  PC_Gauss_y(i_C,i_CT_Elem)
          c_PC_Gauss(3) =  PC_Gauss_z(i_C,i_CT_Elem)          
          
          
          call Cal_KesiYita_by_Coor_3D(c_Center,Solid_El,c_kesi,c_yita,c_eta)
          c_NN(1:8)    = G_NN(1:8,Solid_El)
          call Cal_N_3D(c_kesi,c_yita,c_eta,N)
          tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
          tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
          tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
          tem_n(7) = N(1,19);   tem_n(8) = N(1,22)    
          N_HF(1:3) = ONE/THR      
          num_Local_W = 0
          N_W(1:2,1:60) =ZR
          do i_N = 1,8
              if (Enriched_Node_Type_3D(c_NN(i_N),i_C) ==2) then
                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C)
                  N_W(1,3*num_Local_W-2)  = TWO*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)  = TWO*tem_n(i_N)
                  N_W(3,3*num_Local_W)    = TWO*tem_n(i_N)   
              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==3)then

              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==1)then
                  if ((Elem_Type_3D(Solid_El,i_C).eq.1 )) then
                      ref_elem=Solid_El
                  else 
                      ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(i_C)
                  endif

                  call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
                  
                  BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
                  BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
                  BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
                  BaseLine_x_Vec=Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_y_Vec=Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_z_Vec=Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  Tip_T_Matrx(1:3,1:3)=Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)
                  
                  Gauss_Coor_Local= MATMUL(Tip_T_Matrx,c_Center-BaseLine_Mid)
                  r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
                  theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    
                  
                  
                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C) 


                  N_W(1,3*num_Local_W-2)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(3,3*num_Local_W)  = TWO*sqrt(r_Gauss)*tem_n(i_N)
                  
              endif
          enddo
          
          R_PSI_thread(Local_W(1:3*num_Local_W),c_Thread) = R_PSI_thread(Local_W(1:3*num_Local_W),c_Thread)  +  &
                MATMUL(transpose(N_W(1:3,1:3*num_Local_W)),c_PC_Gauss)*Flu_Ele_Area
                
      enddo
enddo
!$omp end do        
!$omp end parallel


100 continue

DO i_Thread = 1,omp_get_max_threads()
    R_PSI(:)  = R_PSI(:)  + R_PSI_thread(:,i_Thread)
ENDDO      

if (allocated(R_PSI_thread)) deallocate(R_PSI_thread)      
      
      
if(i_NR_P==1)then
  call Vector_Norm2(c_num_freeDOF,R_PSI(freeDOF(1:c_num_freeDOF)),Norm2_Contact_R_PSI_0)   
endif



return 
end SUBROUTINE Cal_Contact_Resid_3D         
