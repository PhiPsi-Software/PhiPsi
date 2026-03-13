!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
SUBROUTINE Ele_by_Ele_FEM_PCG(c_isub,c_Lambda,c_cg_tol, &
                              c_max_num_PCG,c_num_FreeD,c_freeDOF,c_F,c_disp, &
                              c_Total_Num_G_P)
!Element by Element, diagonally preconditioned conjugate gradient solver.
!Ref: Programming the finite element method_2014_Smith_5th, Program 5.6.
!First written on 2020-03-25.
!Fixed on 2026-02-15: Use true residual norm for convergence criterion.

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
real(kind=FT),ALLOCATABLE::diag_precon_thread(:,:)
real(kind=FT) tem_value

real(kind=FT) initial_residual_norm, current_residual_norm
real(kind=FT),ALLOCATABLE::temp_p(:), temp_u(:)

101 FORMAT(13X,'Number of iterations to convergence is ',I5)  
102 FORMAT(13X,'Convergence factor is ',E12.5,' / ',E12.5)  
103 FORMAT(13X,'CG_Iters:',I5, ' -> Residual:', E12.5,' / ',E12.5)

print *, "    >>>> Start of element by element PCG solver <<<<"

print *, "    Step 1: prepare data..."
ALLOCATE(p(0:c_num_FreeD),loads(0:c_num_FreeD), &
         x(0:c_num_FreeD),xnew(0:c_num_FreeD),u(0:c_num_FreeD), &
         diag_precon(0:c_num_FreeD),pcg_d(0:c_num_FreeD))   
diag_precon = ZR

c_Total_Num_G_P = 0
do i_E = 1, Num_Elem    
  Ele_GP_Start_Num(i_E) = c_Total_Num_G_P + 1
  c_Total_Num_G_P = c_Total_Num_G_P + Num_Gauss_P_FEM
enddo

max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:c_num_FreeD, max_threads))
ALLOCATE(diag_precon_thread(0:c_num_FreeD, max_threads))
diag_precon_thread = ZR

print *, "    Step 2: get and store element stiffness matrix..."  
call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM, kesi, yita, weight)

all_local(1:ndof, 1:Num_Elem) = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i_E,mat_num,c_D, &
!$OMP            c_X_NODES,c_Y_NODES,c_NN,c_thread, &
!$OMP            a_local,c_loca,k,j,localK,local,c_thick)
c_thread = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(static)
do i_E = 1, Num_Elem
  mat_num = Elem_Mat(i_E)
  c_thick = thick(Elem_Mat(i_E))
  
  c_D(1:3, 1:3) = D(Elem_Mat(i_E), 1:3, 1:3)
  
  if(Flag_Weibull_E) then
      if (Key_Weibull_E(Elem_Mat(i_E)) == 1) then
          c_D = Weibull_Elements_D_Matrix(i_E, 1:3, 1:3)
      endif
  endif
          
  c_NN = G_NN(1:4, i_E)
  c_X_NODES = G_X_NODES(1:4, i_E)
  c_Y_NODES = G_Y_NODES(1:4, i_E)
  
  a_local = [c_NN(1)*2-1, c_NN(1)*2, c_NN(2)*2-1, c_NN(2)*2, &
             c_NN(3)*2-1, c_NN(3)*2, c_NN(4)*2-1, c_NN(4)*2]
  
  local(1:ndof) = 0
  DO k = 1, ndof
      if (any(c_freeDOF(1:c_num_FreeD) == a_local(k))) then
          c_loca = minloc(c_freeDOF(1:c_num_FreeD), 1, &
                         MASK=(c_freeDOF(1:c_num_FreeD).eq.a_local(k)))
      else
          c_loca = 0
      endif     
      local(k) = c_loca
  enddo
  
  all_local(1:ndof, i_E) = local(1:ndof)
  localK(1:ndof, 1:ndof) = ZR
  
  call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES, c_Y_NODES, c_thick, &
                                    c_D, kesi, yita, weight, localK)
  storK(1:ndof, 1:ndof, i_E) = localK(1:ndof, 1:ndof)
  
  DO j = 1, ndof
      c_loca = local(j)
      if(c_loca > 0) then
          diag_precon_thread(c_loca, c_thread) = &
              diag_precon_thread(c_loca, c_thread) + localK(j, j)
      endif
  END DO
enddo
!$omp end do
!$omp end parallel

DO i_Thread = 1, max_threads
    diag_precon(0:c_num_FreeD) = diag_precon(0:c_num_FreeD) + &
                                  diag_precon_thread(0:c_num_FreeD, i_Thread)
ENDDO

print *, "    Step 3: Invert the preconditioner and get starting loads..."
loads = ZR 
loads(1:c_num_FreeD) = c_F(1:c_num_FreeD)

diag_precon(1:c_num_FreeD) = ONE / diag_precon(1:c_num_FreeD)
diag_precon(0) = ZR 

tem_value = MINVAL(diag_precon(1:c_num_FreeD))
if(tem_value <= 0.0_FT) then
    write(*,*) '    ERROR: Preconditioner has non-positive entries!'
    write(*,*) '    Min value:', tem_value
    stop
endif

pcg_d = diag_precon * loads 
p = pcg_d
x = ZR
cg_iters = 0

print *, "    Step 4: pcg equation solution..."

initial_residual_norm = SQRT(DOT_PRODUCT(loads(1:c_num_FreeD), loads(1:c_num_FreeD)))
if (initial_residual_norm<Tol_20) initial_residual_norm = Tol_20

ALLOCATE(temp_p(ndof), temp_u(ndof))

do i_PCG = 1, c_max_num_PCG
  cg_iters = cg_iters + 1
  u = ZR

  u_thread(0:c_num_FreeD, 1:max_threads) = ZR
  
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,local,localK,temp_p,temp_u,j)
  c_thread = omp_get_thread_num() + 1
  !$OMP DO         
  do i_E = 1, Num_Elem
      local = all_local(1:ndof, i_E)
      localK = storK(1:ndof, 1:ndof, i_E)
      
      do j = 1, ndof
          if(local(j) > 0) then
              temp_p(j) = p(local(j))
          else
              temp_p(j) = 0.0_FT
          endif
      enddo
      
      temp_u = MATMUL(localK, temp_p)
      
      do j = 1, ndof
          if(local(j) > 0) then
              u_thread(local(j), c_thread) = u_thread(local(j), c_thread) + temp_u(j)
          endif
      enddo
  enddo
  !$omp end do
  !$omp end parallel
  
  DO i_Thread = 1, max_threads
      u = u + u_thread(:, i_Thread)
  ENDDO
  
  up = DOT_PRODUCT(loads(1:c_num_FreeD), pcg_d(1:c_num_FreeD))
  alpha = up / DOT_PRODUCT(p(1:c_num_FreeD), u(1:c_num_FreeD))
  xnew = x + p * alpha
  loads = loads - u * alpha
  
  current_residual_norm = SQRT(DOT_PRODUCT(loads(1:c_num_FreeD), loads(1:c_num_FreeD)))
  Tol = current_residual_norm / initial_residual_norm
  
  pcg_d(1:c_num_FreeD) = diag_precon(1:c_num_FreeD) * loads(1:c_num_FreeD)
  pcg_d(0) = ZR
  
  beta = DOT_PRODUCT(loads(1:c_num_FreeD), pcg_d(1:c_num_FreeD)) / up
  p = pcg_d + p * beta
  x = xnew
  
  if(cg_iters <= 300) then
      if(mod(cg_iters, 50) == 0 .OR. cg_iters == 1) then
          write(*,103) cg_iters, Tol, c_cg_tol
      endif
  else
      if(mod(cg_iters, 100) == 0) then
          write(*,103) cg_iters, Tol, c_cg_tol
      endif
  endif
  
  if(Tol < c_cg_tol) then
      exit
  endif
  
  if(cg_iters == c_max_num_PCG) then
      print *, '    -------------------------------------'
      print *, '        Warning :: PCG-EBE failed!       '
      print *, '    -------------------------------------'
      exit
  endif
enddo

DEALLOCATE(temp_p, temp_u)

write(*,101) cg_iters
write(*,102) Tol, c_cg_tol

print *, "    Step 5: retrive nodal displacement..." 
loads = xnew
c_disp = ZR
c_disp(c_freeDOF(1:c_num_FreeD)) = loads(1:c_num_FreeD)

DEALLOCATE(p, loads, x, xnew, u, diag_precon, pcg_d)
if(allocated(u_thread)) deallocate(u_thread)
if(allocated(diag_precon_thread)) deallocate(diag_precon_thread)

print *, "    >>>>  End of element by element PCG solver  <<<<"
RETURN
END SUBROUTINE Ele_by_Ele_FEM_PCG
