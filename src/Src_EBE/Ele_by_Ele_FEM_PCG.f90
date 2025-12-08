 
SUBROUTINE Ele_by_Ele_FEM_PCG(c_isub,c_Lambda,c_cg_tol, &
                              c_max_num_PCG,c_num_FreeD,c_freeDOF,c_F,c_disp, &
                              c_Total_Num_G_P)

use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use omp_lib

implicit none
integer,intent(in)::c_isub,c_max_num_PCG,c_num_FreeD
integer,intent(in)::c_freeDOF(c_num_FreeD)
real(kind=FT),intent(in)::c_Lambda,c_cg_tol,c_F(c_num_FreeD)
real(kind=FT),intent(out)::c_disp(Total_FD)
integer,intent(out)::c_Total_Num_G_P
integer j,i_E,c_NN(4),mat_num,k
real(kind=FT) c_D(3,3)
real(kind=FT) kesi(Num_Gauss_P_FEM),yita(Num_Gauss_P_FEM),weight(Num_Gauss_P_FEM)   
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer ndof
parameter(ndof = 2*4)   
real(kind=FT) localK(ndof,ndof)
integer a_local(ndof),local(ndof),all_local(ndof,Num_Elem)
real(kind=FT) storK(ndof,ndof,Num_Elem)
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:),diag_precon(:),pcg_d(:)
integer cg_iters,i_PCG
real(kind=FT) alpha,beta,up,Tol,c_thick
integer c_loca
integer i_Thread,c_Thread,max_threads
real(kind=FT),ALLOCATABLE::u_thread(:,:) 
real(kind=FT) tem_value

print *, "    >>>> Start of element by element PCG solver <<<<"
print *, "    Step 1: prepare data..."
ALLOCATE(p(0:c_num_FreeD),loads(0:c_num_FreeD), &
         x(0:c_num_FreeD),xnew(0:c_num_FreeD),u(0:c_num_FreeD), &
         diag_precon(0:c_num_FreeD),pcg_d(0:c_num_FreeD))   
diag_precon= ZR
c_Total_Num_G_P = 0
do i_E = 1,Num_Elem    
  Ele_GP_Start_Num(i_E) = c_Total_Num_G_P + 1
  c_Total_Num_G_P = c_Total_Num_G_P +Num_Gauss_P_FEM
enddo

max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:c_num_FreeD,max_threads))
  
print *, "    Step 2: get and store element stiffness matrix..."  
call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi,yita,weight)


all_local(1:ndof,1:Num_Elem) = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,mat_num,c_D, &
!$OMP            c_X_NODES,c_Y_NODES,c_NN, &
!$OMP            a_local,c_loca,k,j,localK,local,c_thick)    
do i_E=1,Num_Elem
  mat_num = Elem_Mat(i_E)
  c_thick = thick(Elem_Mat(i_E))
  
  c_D(1:3,1:3)     = D(Elem_Mat(i_E),1:3,1:3)         
  
  if(Flag_Weibull_E)then
      if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
          c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
      endif
  endif
          
  c_NN    = G_NN(1:4,i_E)
  c_X_NODES = G_X_NODES(1:4,i_E)
  c_Y_NODES = G_Y_NODES(1:4,i_E) 
  a_local =[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]
  local(1:ndof) = 0
  DO k=1,ndof 
      if (any(c_freeDOF(1:c_num_FreeD)== a_local(k)))then
          c_loca =minloc(c_freeDOF(1:c_num_FreeD),1,MASK=(c_freeDOF(1:c_num_FreeD).eq.a_local(k)))
      else
          c_loca =0
      endif     
      local(k) = c_loca
  enddo
  all_local(1:ndof,i_E) = local(1:ndof)
  localK(1:ndof,1:ndof) = ZR
    call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES,c_Y_NODES,c_thick,c_D,kesi,yita,weight,localK)     
  storK(1:ndof,1:ndof,i_E)=localK(1:ndof,1:ndof)       
  DO j=1,ndof 
      diag_precon(local(j))=diag_precon(local(j))+localK(j,j) 
  END DO
enddo
!$omp end parallel do   


print *, "    Step 3: Invert the preconditioner and get starting loads..."
loads=ZR 
loads(1:c_num_FreeD) = c_F(1:c_num_FreeD)
diag_precon(1:)=ONE/diag_precon(1:) 
diag_precon(0)=ZR 

pcg_d   =diag_precon*loads 
p       =pcg_d
x       =ZR
cg_iters=0

print *, "    Step 4: pcg equation solution..."
do i_PCG =1,c_max_num_PCG
  cg_iters=cg_iters+1 
  u=ZR

  u_thread(0:c_num_FreeD,1:max_threads)= ZR  
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,local,localK) 
  c_thread = omp_get_thread_num()+1
  !$OMP DO         
  do i_E=1,Num_Elem
      local=all_local(1:ndof,i_E) 
      localK=storK(1:ndof,1:ndof,i_E) 
      u_thread(local,c_thread)=u_thread(local,c_thread)+MATMUL(localK,p(local)) 
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
  pcg_d=diag_precon*loads 
  beta=DOT_PRODUCT(loads,pcg_d)/up 
  p=pcg_d+p*beta
  
  
  tem_value = MAXVAL(ABS(x(0:c_num_FreeD)))
  if(tem_value==ZR) tem_value = Tol_20
  Tol = MAXVAL(ABS(x(0:c_num_FreeD)-xnew(0:c_num_FreeD)))/tem_value
  
  x=xnew
  if(Tol<c_cg_tol.OR.cg_iters==c_max_num_PCG)then
      exit
  endif
  
enddo
print *, "    Number of CG iterations to convergence was",cg_iters

loads=xnew
c_disp = ZR

print *, "    Step 5: retrive nodal displacement..." 
     
c_disp(c_freeDOF(1:c_num_FreeD)) =loads(1:c_num_FreeD)

DEALLOCATE(p,loads,x,xnew,u,diag_precon,pcg_d)  

print *, "    >>>>  End of element by element PCG solver  <<<<"
RETURN
END SUBROUTINE Ele_by_Ele_FEM_PCG
