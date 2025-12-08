 
SUBROUTINE Ele_by_Ele_FEM_PCG_3D(isub,Lambda,c_cg_tol,max_num_PCG,num_FreeD,freeDOF,F,disp)

use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use omp_lib
use Global_Cal_Ele_Stiffness_Matrix_3D_8nodes

implicit none
integer,intent(in)::isub,max_num_PCG,num_FreeD
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
real(kind=FT),intent(out)::disp(Total_FD)
integer i,j,m,i_E,i_G,c_NN(8),mat_num,k
real(kind=FT) c_D(6,6)
real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D),  &
            zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)   
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer ndof
parameter(ndof = 3*8)   
real(kind=FT) localK(ndof,ndof)
integer a_local(ndof),local(ndof),c_All_Local_3D(ndof,Num_Elem)
real(kind=FT) storK(ndof,ndof,Num_Elem)
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:),  &
                         precon_vector(:),pcg_d(:)
integer c_Node,cg_iters,i_PCG,i_Node
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
integer c_loca,c_DOF
real(kind=FT),ALLOCATABLE::u_thread(:,:)
real(kind=FT),ALLOCATABLE::pcg_d_thread(:,:) 
integer i_Thread,c_Thread,max_threads
real(kind=FT) FEM_det
integer c_INFO
real(kind=FT),ALLOCATABLE::Elements_Le(:,:,:),Elements_Le_Inv(:,:,:) 
real(kind=FT),ALLOCATABLE::Elements_De(:,:,:),Elements_De_Inv(:,:,:)
real(kind=FT),ALLOCATABLE::I_Matrx(:,:)
real(kind=FT),ALLOCATABLE::Matrx_For_Crout(:,:)
real(kind=FT),ALLOCATABLE::Crout_L(:,:)
real(kind=FT),ALLOCATABLE::Crout_U(:,:)
real(kind=FT),ALLOCATABLE::Crout_D(:,:)
real(kind=FT),ALLOCATABLE::Ele_Diag_Matrix(:,:,:)
real(kind=FT),ALLOCATABLE::SQRT_Ele_Diag_Matrix(:,:,:)
real(kind=FT),ALLOCATABLE::precon_vector_1_2(:)
integer i_ndof
integer ii,jj
real(kind=FT),ALLOCATABLE::precon_vector_thread(:,:)
real(kind=FT) tem_value
integer date_time(8)
integer(LIT) c_S_Time,F_time
character*10  current_data
5003 FORMAT('    Elapsed CPU time of EBE-PCG solution - ',I8,' s, about ',F10.4,' mins')

print *, "    >>>> Start of element by element PCG solver <<<<"
print *, "    Step 1: prepare data..."
ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD),  &
       x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD),  &
       precon_vector(0:num_FreeD),pcg_d(0:num_FreeD))   
precon_vector= ZR


max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:num_FreeD,max_threads))

print *, "    Step 2: get and store element stiffness matrix..."
call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)    
if(Key_EBE_Precondition > 0)then
  ALLOCATE(precon_vector_thread(0:num_FreeD,max_threads))
  precon_vector_thread(0:num_FreeD,1:max_threads)= ZR 
endif

c_All_Local_3D(1:ndof,1:Num_Elem) = 0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i_E,mat_num,c_D,Volume_Ratio,    &
!$OMP             c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,Rot_c_D_Comp,   &
!$OMP             c_Y_NODES,c_Z_NODES,c_NN, c_thread,                   &
!$OMP             a_local,c_loca,k,j,localK,local) 
c_thread = omp_get_thread_num()+1
!$OMP DO SCHEDULE(static)         
do i_E=1,Num_Elem
  mat_num = Elem_Mat(i_E)
  c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)        
  if (Material_Type(mat_num)==5)then
      Volume_Ratio = Material_Para_Added(mat_num,10)
      c_D_Comp = D_Comp(mat_num,1:6,1:6)
      T_Matrix = Ele_ComMat_RotMatrix(i_E,1:6,1:6)
      TT_Matrix= TRANSPOSE(T_Matrix)
      Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
      Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
      c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
  endif          
  c_NN    = G_NN(1:8,i_E)
  c_X_NODES = G_X_NODES(1:8,i_E)
  c_Y_NODES = G_Y_NODES(1:8,i_E)    
  c_Z_NODES = G_Z_NODES(1:8,i_E)           
  a_local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, &
           c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, &
           c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
           c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3] 
  DO k=1,ndof 
      
      if (Flag_FreeDOF(a_local(k)) == 1 )then
          
          c_loca = Location_FreeDOF(a_local(k))
      else
          c_loca =0
      endif    
      
      local(k) = c_loca
  enddo
  c_All_Local_3D(1:ndof,i_E) = local(1:ndof)
  localK(1:ndof,1:ndof) = ZR
  call Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E,                  &
       Num_Gauss_P_FEM_3D,c_X_NODES,c_Y_NODES,c_Z_NODES,        &
                          c_D,kesi,yita,zeta,weight,localK)  
  storK(1:ndof,1:ndof,i_E)=localK(1:ndof,1:ndof)       
  
  if(Key_EBE_Precondition > 0) then
      DO j=1,ndof 
          c_loca=local(j)
          precon_vector_thread(c_loca,c_thread) = precon_vector_thread(c_loca,c_thread) + localK(j,j) 
      END DO
  endif
enddo
!$omp end do
!$omp end parallel    

DO i_Thread = 1,omp_get_max_threads()
precon_vector(0:num_FreeD)  =  precon_vector(0:num_FreeD)  + precon_vector_thread(0:num_FreeD,i_Thread)
ENDDO


if(Key_EBE_Precondition ==2) then 
#ifndef github 
#ifndef Silverfrost
  if(.not. allocated(Elements_Le_Inv)) allocate(Elements_Le_Inv(ndof,ndof,Num_Elem))
  if(.not. allocated(Elements_De_Inv)) allocate(Elements_De_Inv(ndof,ndof,Num_Elem))
  allocate(I_Matrx(ndof,ndof))
  I_Matrx =  ZR
  forall (i_ndof = 1:ndof) I_Matrx(i_ndof,i_ndof) = ONE

  print *, "            Getting Le of elements for HW preconditioner..." 

  allocate(Matrx_For_Crout(ndof,ndof))
  allocate(Crout_L(ndof,ndof))
  allocate(Crout_U(ndof,ndof))
  allocate(Crout_D(ndof,ndof))
  Crout_L = ZR
  Crout_U = ZR
  Crout_D = ZR
  allocate(Ele_Diag_Matrix(ndof,ndof,Num_Elem))  
  allocate(SQRT_Ele_Diag_Matrix(ndof,ndof,Num_Elem)) 
  Ele_Diag_Matrix = ZR
  SQRT_Ele_Diag_Matrix =  ZR

  do i_E = 1,Num_Elem
      localK(1:ndof,1:ndof) = storK(1:ndof,1:ndof,i_E)
      
      local(1:ndof) = c_All_Local_3D(1:ndof,i_E) 
      do i_ndof = 1,ndof
          c_loca=local(i_ndof)
          Ele_Diag_Matrix(i_ndof,i_ndof,i_E)      = localK(i_ndof,i_ndof)
          SQRT_Ele_Diag_Matrix(i_ndof,i_ndof,i_E) = SQRT(ONE/precon_vector(c_loca))
      enddo          
    
      Matrx_For_Crout(1:ndof,1:ndof) = I_Matrx +  &
           MATMUL(MATMUL(SQRT_Ele_Diag_Matrix(1:ndof,1:ndof,i_E),(localK(1:ndof,1:ndof)-&
                        Ele_Diag_Matrix(1:ndof,1:ndof,i_E))),SQRT_Ele_Diag_Matrix(1:ndof,1:ndof,i_E))
                        
      call Square_Matrix_Crout_Decomposition(ndof,Matrx_For_Crout(1:ndof,1:ndof),Crout_L(1:ndof,1:ndof),&
                                             Crout_U(1:ndof,1:ndof),Crout_D(1:ndof,1:ndof))
      

      
      call Matrix_Inverse_Lapack(Crout_L(1:ndof,1:ndof),ndof)
      Elements_Le_Inv(1:ndof,1:ndof,i_E) = Crout_L(1:ndof,1:ndof)
      call Matrix_Inverse_Lapack(Crout_D(1:ndof,1:ndof),ndof)  
      Elements_De_Inv(1:ndof,1:ndof,i_E) = Crout_D(1:ndof,1:ndof) 
enddo         
deallocate(I_Matrx)
deallocate(Matrx_For_Crout)
deallocate(Crout_L)
deallocate(Crout_U)
#endif 
#endif  
endif


print *, "    Step 3: Get starting loads..."

loads=ZR 
loads(1:num_FreeD) = F(1:num_FreeD)

if(Key_EBE_Precondition ==0) then  
  precon_vector(1:) =  ONE
  precon_vector(0)=ZR 
  pcg_d= precon_vector*loads 
  p    = pcg_d
elseif(Key_EBE_Precondition ==1) then  
  precon_vector(1:)=ONE/precon_vector(1:) 
  precon_vector(0)=ZR 
  pcg_d= precon_vector*loads 
  p    = pcg_d
elseif(Key_EBE_Precondition ==2) then 
  if (allocated(pcg_d_thread)) deallocate(pcg_d_thread)
  ALLOCATE(pcg_d_thread(0:num_FreeD,max_threads))
  
  
  ALLOCATE(precon_vector_1_2(0:num_FreeD))
  precon_vector_1_2(1:num_FreeD) = ONE/SQRT(precon_vector(1:num_FreeD))
  precon_vector_1_2(0)=ZR 
  
  
  pcg_d = precon_vector_1_2*loads
  
  do i_E = 1,Num_Elem
      local     = c_All_Local_3D(1:ndof,i_E) 
      pcg_d(local) = MATMUL(Elements_Le_Inv(1:ndof,1:ndof,i_E),pcg_d(local)) 
  enddo 
  
  do i_E = 1,Num_Elem
      local     = c_All_Local_3D(1:ndof,i_E) 
      pcg_d(local) = MATMUL(Elements_De_Inv(1:ndof,1:ndof,i_E),pcg_d(local))
  enddo 

  do i_E = Num_Elem,1,-1
      local     = c_All_Local_3D(1:ndof,i_E) 
      pcg_d(local) = MATMUL(TRANSPOSE(Elements_Le_Inv(1:ndof,1:ndof,i_E)),pcg_d(local)) 
  enddo 
  
  pcg_d = precon_vector_1_2*pcg_d
  
  p     = pcg_d
  pcg_d(0) = ZR
  p(0)     = ZR
  

endif


x=ZR
cg_iters=0

if(Key_Print_EBEPCG_Solution_Time==1)then
    call Tool_Get_Current_Time(current_data,date_time,c_S_Time)
endif

print *, "    Step 4: pcg equation solution..."
do i_PCG =1,max_num_PCG
  cg_iters=cg_iters+1 
  u=ZR
   
  u_thread(0:num_FreeD,1:max_threads)= ZR  
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,local) 
  c_thread = omp_get_thread_num()+1
  !$OMP DO         
  do i_E=1,Num_Elem
      local=c_All_Local_3D(1:ndof,i_E) 
      u_thread(local,c_thread)=u_thread(local,c_thread)+MATMUL(storK(1:ndof,1:ndof,i_E) ,p(local))
  enddo  
  !$omp end do
  !$omp end parallel       
  
  DO i_Thread = 1,omp_get_max_threads()
      u  =  u  + u_thread(:,i_Thread)
  ENDDO      

  up=DOT_PRODUCT(loads,pcg_d)
  alpha=up/DOT_PRODUCT(p,u) 
  xnew=x+p*alpha 
  loads=loads-u*alpha
  
  if(Key_EBE_Precondition ==0) then  
      pcg_d=precon_vector*loads
      
  elseif(Key_EBE_Precondition ==1) then  
      pcg_d=precon_vector*loads
      
  elseif(Key_EBE_Precondition ==2) then 
      pcg_d = precon_vector_1_2*loads

      do i_E = 1,Num_Elem
          local     = c_All_Local_3D(1:ndof,i_E) 
          pcg_d(local) = MATMUL(Elements_Le_Inv(1:ndof,1:ndof,i_E),pcg_d(local)) 
      enddo           

      do i_E = 1,Num_Elem
          local     = c_All_Local_3D(1:ndof,i_E) 
          pcg_d(local) = MATMUL(Elements_De_Inv(1:ndof,1:ndof,i_E),pcg_d(local)) 
      enddo 
      
      do i_E = Num_Elem,1,-1
          local     = c_All_Local_3D(1:ndof,i_E) 
          pcg_d(local) = MATMUL(TRANSPOSE(Elements_Le_Inv(1:ndof,1:ndof,i_E)),pcg_d(local)) 
      enddo 

      pcg_d = precon_vector_1_2*pcg_d
      pcg_d(0) = ZR
  endif
  
  beta=DOT_PRODUCT(loads,pcg_d)/up 

  p= pcg_d + p*beta
  
  
  tem_value = MAXVAL(ABS(x(0:num_FreeD)))
  if(tem_value==ZR) tem_value = Tol_20
  Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/tem_value
      
  x=xnew
  if(Tol<c_cg_tol.OR.cg_iters==max_num_PCG)then
      exit
  endif
  
    if(cg_iters<=300)then
        if(mod(cg_iters,50) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
    else
        if(mod(cg_iters,100) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
    endif   
    
    if(Tol<c_cg_tol)then  
          exit
    endif
    
    if(cg_iters==max_num_PCG)then
          print *, '    -------------------------------------'
          print *, '        Warning :: PCG-EBE failed!       '  
          print *, '    -------------------------------------'
          exit
    endif     
enddo
101 FORMAT(13X,'Number of iterations to convergence is ',I5)  
102 FORMAT(13X,'Convergence factor is ',E12.5,' / ',E12.5)  
103 FORMAT(13X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)  
write(*,101)  cg_iters      
write(*,102)  Tol,c_cg_tol   

if(Key_Print_EBEPCG_Solution_Time==1)then
    call Tool_Get_Current_Time(current_data,date_time,F_time)
    WRITE(*,5003) F_time-c_S_Time,(dble(F_time)-dble(c_S_Time))/Con_60
endif


print *, "    Step 5: retrive nodal displacement..."         
loads=xnew
disp = ZR
disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)

DEALLOCATE(p,loads,x,xnew,u,precon_vector,pcg_d)  
if(allocated(Elements_Le)) deallocate(Elements_Le)
if(allocated(Elements_Le_Inv)) deallocate(Elements_Le_Inv)

print *, "    >>>>  End of element by element PCG solver  <<<<"
RETURN
END SUBROUTINE Ele_by_Ele_FEM_PCG_3D
