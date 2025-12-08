 
subroutine Cal_InSitu_Stress

use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_DISP
use Global_HF
use Global_Stress
use Global_POST
use Global_Contact
use Global_Material

implicit none
real(kind=FT) F_InSitu(2*Num_Node)
integer Total_Num_G_P_InSitu,num_freeDOF_InSitu
real(kind=FT),ALLOCATABLE::globalK_InSitu(:,:)
integer,ALLOCATABLE::freeDOF_InSitu(:)
real(kind=FT),ALLOCATABLE::tem_DISP_InSitu(:)
real(kind=FT) Max_Disp_x_InSitu,Min_Disp_x_InSitu,Max_Disp_y_InSitu,Min_Disp_y_InSitu
integer Total_Freedom_InSitu
real(kind=FT),ALLOCATABLE:: K_IS_CSR_aa(:)
integer,ALLOCATABLE:: K_IS_CSR_ia(:)
integer,ALLOCATABLE:: K_IS_CSR_ja(:)
integer K_IS_COO_NNZ_Max
integer K_IS_CSR_NNZ
integer IS_Total_FD
integer i_E
real(kind=FT) c_D(3,3),U(8)
real(kind=FT) c_X_NODES(4),c_Y_NODES(4),c_kesi,c_yita,c_Stress(3)
integer c_NN(4)
real(kind=FT) kesi_N(4),yita_N(4),weight_N(4)
integer G_Counter,i_G
real(kind=FT) kesi(4),yita(4)
logical Yes_Add_Insitu
integer c_Key_SLOE

1021 FORMAT(5X,'Range of InSitu displacement x:   ',F13.5,' m to ', F13.5,' m')
1022 FORMAT(5X,'Range of InSitu displacement y:   ',F13.5,' m to ', F13.5,' m')
1121 FORMAT(5X,'Range of InSitu displacement x:   ',F13.5,' mm to ',F13.5,' mm')
1122 FORMAT(5X,'Range of InSitu displacement y:   ',F13.5,' mm to ',F13.5,' mm')
1221 FORMAT(5X,'Value of InSitu stress x:   ',F13.5, ' MPa')
1222 FORMAT(5X,'Value of InSitu stress y:   ',F13.5, ' MPa')

print *,'    ################################'
print *,'         GEO_STATIC CALCULATION     '
print *,'    ################################'
print *,'    Calculating displacement under InSitu stress...'

ALLOCATE(DISP_InSitu(2*Num_Node))
IS_Total_FD = 2*Num_Node


call Force_Vector_no_Killed_Ele(2*Num_Node,1,.False.,ONE,F_InSitu)
ALLOCATE(freeDOF_InSitu(2*Num_Node))
call Boundary_Cond(2*Num_Node,1,freeDOF_InSitu,num_freeDOF_InSitu)

if(Key_K_Sparse==0)then
  ALLOCATE(globalK_InSitu(2*Num_Node,2*Num_Node))
  call Assemble_Stiffness_Matrix_FEM(0,globalK_InSitu,2*Num_Node,&
    Total_Num_G_P_InSitu)
elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
  K_IS_COO_NNZ_Max = CEILING(IS_Total_FD**2*Sparse_Ratio)
  ALLOCATE(K_IS_CSR_aa(K_IS_COO_NNZ_Max))
  ALLOCATE(K_IS_CSR_ja(K_IS_COO_NNZ_Max))
  ALLOCATE(K_IS_CSR_ia(num_freeDOF_InSitu+1))
  call Assemble_Stiffness_Matrix_SPARS_FEM(0,    &
         freeDOF_InSitu(1:num_freeDOF_InSitu),num_freeDOF_InSitu,&
         K_IS_CSR_aa,K_IS_CSR_ja,K_IS_CSR_ia,&
         K_IS_COO_NNZ_Max,K_IS_CSR_NNZ,IS_Total_FD,&
         Total_Num_G_P_InSitu)
#endif
#endif
#endif
endif

DISP_InSitu(1:2*Num_Node) = ZR
ALLOCATE(tem_DISP_InSitu(num_freeDOF_InSitu))
if(Key_K_Sparse==0)then
  if (Key_SLOE==11) then
      
#ifndef Silverfrost
      c_Key_SLOE =  7
#endif

#ifdef Silverfrost
      c_Key_SLOE =  4
#endif

  else
      c_Key_SLOE = Key_SLOE
  endif
  
#ifdef ifort  
    c_Key_SLOE =  3
#endif

#ifdef ifx  
    c_Key_SLOE =  3
#endif

  call Matrix_Solve_LSOE(0,1,c_Key_SLOE,        &
           globalK_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu),&
                          freeDOF_InSitu(1:num_freeDOF_InSitu)),&
           F_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu)),&
           tem_DISP_InSitu,&
           num_freeDOF_InSitu)
elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
  call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
              K_IS_CSR_NNZ,&
              K_IS_CSR_aa(1:K_IS_CSR_NNZ),&
              K_IS_CSR_ja(1:K_IS_CSR_NNZ),    &
              K_IS_CSR_ia(1:num_freeDOF_InSitu+1),  &
              F_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu)),&
              tem_DISP_InSitu,num_freeDOF_InSitu)
#endif
#endif
#endif
endif

DISP_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu))=tem_DISP_InSitu
print *,'    Displacements under InSitu stress obtained.'

#ifndef Silverfrost  
Max_Disp_x_InSitu = maxval(DISP_InSitu(1:2*num_Node:2))
Min_Disp_x_InSitu = minval(DISP_InSitu(1:2*num_Node:2))
Max_Disp_y_InSitu = maxval(DISP_InSitu(2:2*num_Node:2))
Min_Disp_y_InSitu = minval(DISP_InSitu(2:2*num_Node:2))
if(Key_Unit_System==1)then
  WRITE(*,1021) Min_Disp_x_InSitu,Max_Disp_x_InSitu
  WRITE(*,1022) Min_Disp_y_InSitu,Max_Disp_y_InSitu
elseif(Key_Unit_System==2)then
  WRITE(*,1121) Min_Disp_x_InSitu,Max_Disp_x_InSitu
  WRITE(*,1122) Min_Disp_y_InSitu,Max_Disp_y_InSitu
endif
#endif

print *,'    Calculating InSitu stress of nodes...'
ALLOCATE(Str_xx_InSitu(Num_Node))
ALLOCATE(Str_yy_InSitu(Num_Node))
ALLOCATE(Str_xy_InSitu(Num_Node))
ALLOCATE(Str_vm_InSitu(Num_Node))
Yes_Add_Insitu = .False.


call Get_Node_Stress_FEM_IN_OUT(Yes_Add_Insitu,1,DISP_InSitu,&
     Str_xx_InSitu,Str_yy_InSitu,Str_xy_InSitu,Str_vm_InSitu)



InSitu_x = sum(Str_xx_InSitu)/Num_Node
InSitu_y = sum(Str_yy_InSitu)/Num_Node
InSitu_xy = sum(Str_xy_InSitu)/Num_Node

if(Key_Unit_System==1)then
  print *,"    ************************************************"
  WRITE(*,1221) -InSitu_x/1.0D6
  WRITE(*,1222) -InSitu_y/1.0D6
  print *,"    ************************************************"
elseif(Key_Unit_System==2)then
  WRITE(*,1221) -InSitu_x
  WRITE(*,1222) -InSitu_y
endif

print *,'    Calculating InSitu stress of Gauss points...'
ALLOCATE(InSitu_Strs_Gaus_xx(Num_Elem,Num_Gauss_P_FEM))
ALLOCATE(InSitu_Strs_Gaus_yy(Num_Elem,Num_Gauss_P_FEM))
ALLOCATE(InSitu_Strs_Gaus_xy(Num_Elem,Num_Gauss_P_FEM))
call Cal_Gauss_Points_QUAD(4,kesi_N,yita_N,weight_N)
G_Counter = 0
do i_E = 1,Num_Elem
  c_D     = D(Elem_Mat(i_E),:,:)
  
  if(Flag_Weibull_E)then
      if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
          c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
      endif
  endif
  
  c_NN    = G_NN(:,i_E)
  c_X_NODES = G_X_NODES(:,i_E)
  c_Y_NODES = G_Y_NODES(:,i_E)
  U = [DISP_InSitu(c_NN(1)*2-1),DISP_InSitu(c_NN(1)*2),&
         DISP_InSitu(c_NN(2)*2-1),DISP_InSitu(c_NN(2)*2),&
         DISP_InSitu(c_NN(3)*2-1),DISP_InSitu(c_NN(3)*2),&
         DISP_InSitu(c_NN(4)*2-1),DISP_InSitu(c_NN(4)*2)]

  kesi(1:4)   = kesi_N
  yita(1:4)   = yita_N
  do i_G = 1,4
      G_Counter =  G_Counter +1
    c_kesi = kesi(i_G)
    c_yita = yita(i_G)
    call Cal_Ele_Stress_N4(i_E,i_G,c_X_NODES,c_Y_NODES,c_D,c_kesi,c_yita,U,c_Stress)
      InSitu_Strs_Gaus_xx(i_E,i_G) = c_Stress(1)
    InSitu_Strs_Gaus_yy(i_E,i_G) = c_Stress(2)
    InSitu_Strs_Gaus_xy(i_E,i_G) = c_Stress(3)
  end do
end do

if(Key_K_Sparse==0)then
  DEALLOCATE(globalK_InSitu)
elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
  DEALLOCATE(K_IS_CSR_aa); DEALLOCATE(K_IS_CSR_ia)
  DEALLOCATE(K_IS_CSR_ja)
#endif
#endif
endif
DEALLOCATE(freeDOF_InSitu,tem_DISP_InSitu)

State_InSitu = 0
if(InSitu_x < 0.01D6 .or. InSitu_y < 0.01D6)then
  State_InSitu = 1
  InSitu_x = -InSitu_x
  InSitu_y = -InSitu_y
endif

return
end SUBROUTINE Cal_InSitu_Stress
