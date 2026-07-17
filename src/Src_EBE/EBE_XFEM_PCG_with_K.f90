!-----------------------------------------------------------
! Brief: Solve the 2D XFEM system via diagonally preconditioned PCG.
!
! Parameters:
!   Input:  isub                  - load step index
!           Lambda                - load factor
!           c_cg_tol              - PCG convergence tolerance
!           max_num_PCG           - maximum PCG iterations
!           num_FreeD             - number of free DOF
!           freeDOF               - list of free DOF indices
!           F                     - global load vector (free DOF)
!           storK                 - pre-stored element XFEM stiffnesses
!           size_local            - per-element local DOF size
!           all_local             - per-element local DOF map
!           diag_precon_no_invert - diagonal preconditioner (uninverted)
!   Output: disp                  - solution displacement vector
!
! Notes:   Uses EBE matrix-vector multiply; the preconditioner is inverted
!          on entry to the iteration loop.
!-----------------------------------------------------------

SUBROUTINE EBE_XFEM_PCG_with_K(isub,Lambda,c_cg_tol,max_num_PCG, num_FreeD,freeDOF,F, storK,size_local,all_local, &
diag_precon_no_invert,disp)
!Element by Element, diagonally preconditioned conjugate gradient solver.
!Ref: Programming the finite element method_2014_Smith_5th, Program 5.6.
!First written on 2020-03-25.

use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_HF
use Global_Contact
use Global_Inclusion
use Global_Cross
use omp_lib

implicit none
integer,intent(in)::isub,max_num_PCG,num_FreeD
real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(out)::disp(Total_FD)
real(kind=FT),intent(in)::storK(MDOF_2D,MDOF_2D,Num_Elem)
real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
integer,intent(in):: size_local(Num_Elem)
integer,intent(in):: all_local(MDOF_2D,Num_Elem)
real(kind=FT)::diag_precon(0:num_FreeD)
integer i,j,m,i_E,i_G,c_NN(4),mat_num,k
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer local(MDOF_2D)
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:), pcg_d(:)
integer c_Node,cg_iters,i_PCG,i_Node
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) kesi_Enr(Num_Gauss_Points), &
yita_Enr(Num_Gauss_Points), weight_Enr(Num_Gauss_Points)
real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM), &
yita_N_Enr(Num_Gauss_P_FEM), weight_N_Enr(Num_Gauss_P_FEM)
integer c_Num_Gauss_Point,i_C
real(kind=FT) kesi(900),yita(900),weight(900)
integer:: Location_ESM(MDOF_2D)
integer::Location_ESM_C_Crack(MDOF_2D),Location_ESM_C_Cr_NoFEM(MDOF_2D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(3,80),tem_B(3,80),detJ,tem_Val
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer num_Bou,c_loca,c_DOF,j_C,i_Incl,i_H,i_Cross
real(kind=FT) max_KK,min_kk
real(kind=FT) c_N(2,8),c_G_x,c_G_y,c_thick
integer c_Incl_Num

integer i_Thread,c_Thread,max_threads
real(kind=FT),ALLOCATABLE::u_thread(:,:) 
real(kind=FT) tem_value
real(kind=FT) current_residual_norm,initial_residual_norm


ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD), x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD), pcg_d(0:num_FreeD))

max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:num_FreeD,max_threads))

diag_precon= ZR

print *, "    EBE-PCG: invert the preconditioner and get starting loads..."
loads=ZR 
loads(1:num_FreeD) = F(1:num_FreeD)
if (Key_EBE_Precondition==0) then
    diag_precon=ONE
    diag_precon(0) =ZR 
elseif(Key_EBE_Precondition==1) then
    diag_precon(1:)=ONE/diag_precon_no_invert(1:)
    diag_precon(0) =ZR 
endif

pcg_d=diag_precon*loads 
p=pcg_d
x=ZR
cg_iters=0

print *, "    EBE-PCG: pcg equation solution..."

initial_residual_norm = SQRT(DOT_PRODUCT(loads(1:num_FreeD), loads(1:num_FreeD)))
if (initial_residual_norm<Tol_20) initial_residual_norm = Tol_20


do i_PCG =1,max_num_PCG
    cg_iters=cg_iters+1 
    u=ZR


    u_thread(0:num_FreeD,1:max_threads)= ZR  
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,num_Loc_ESM,local) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO         
    do i_E=1,Num_Elem
        num_Loc_ESM = size_local(i_E)
        local(1:num_Loc_ESM)=all_local(1:num_Loc_ESM,i_E) 
        u_thread(local(1:num_Loc_ESM),c_thread) = u_thread(local(1:num_Loc_ESM),c_thread)+ &
        MATMUL(storK(1:num_Loc_ESM,1:num_Loc_ESM,i_E),p(local(1:num_Loc_ESM)))
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

    current_residual_norm = SQRT(DOT_PRODUCT(loads(1:num_FreeD), loads(1:num_FreeD)))
    Tol = current_residual_norm / initial_residual_norm

    pcg_d=diag_precon*loads 
    beta=DOT_PRODUCT(loads,pcg_d)/up 
    p=pcg_d+p*beta




    x=xnew

    if(cg_iters<=300)then
        if(mod(cg_iters,50) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
    else
        if(mod(cg_iters,100) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
    endif  

    if(Tol<c_cg_tol.OR.cg_iters==max_num_PCG)then
        exit
    endif
enddo
103 FORMAT(14X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)   
101 FORMAT(5X,'PCG-EBE: number of iterations to convergence is ',I5)  
102 FORMAT(5X,'PCG-EBE: convergence factor is ',E12.5,' / ',E12.5) 
write(*,101)  cg_iters      
write(*,102)  Tol,c_cg_tol      
loads=xnew
disp = ZR
print *, "    EBE-PCG: retrive nodal displacement..."
disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)
DEALLOCATE(p,loads,x,xnew,u,pcg_d)  
RETURN
END SUBROUTINE EBE_XFEM_PCG_with_K
