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
 
SUBROUTINE Matrix_Solve_LSOE_Lis(Key_Indent,input_NZZ,input_n,K_Ptr,K_Index,K_Value,input_F,output_D)
! Lis solver.
      

use Global_Float_Type
use Global_Common
use Global_Model

implicit none

#include "lisf.h"

integer,intent(in)::input_NZZ,input_n,Key_Indent
integer,intent(in)::K_Ptr(input_n+1)
integer,intent(in)::K_Index(input_NZZ)
real(kind=FT),intent(in)::K_Value(input_NZZ)

real(kind=FT),intent(in)::input_F(input_n)
real(kind=FT),intent(out)::output_D(input_n)
real(kind=FT) cputime

LIS_INTEGER      my_rank,nprocs
LIS_INTEGER      matrix_type,comm_world

LIS_INTEGER      omp_get_num_procs,omp_get_max_threads

LIS_INTEGER      n,gn,ln,iter,lis_error
LIS_MATRIX       A
LIS_VECTOR       x,u
LIS_SOLVER       solver


LIS_INTEGER  NZ_NUM
LIS_INTEGER A_Lis_Ptr(input_n+1)
LIS_INTEGER A_Lis_Index(input_NZZ)
LIS_SCALAR  A_Lis_Value(input_NZZ)
LIS_REAL residual

integer c_i

801 FORMAT(5X,'Max iteration number of Lis: ',I8,' default(2000)')         
901 FORMAT(5X,'Convergence tolerance of Lis:   ',E12.5)           
1001 FORMAT(5X,'Number of iteration of Lis:   ',I8)        
1002 FORMAT(5X,'CPU time of Lis:          ',F12.5,' mins')      
1003 FORMAT(5X,'Relative residual norm of Lis:',E12.5)     

811 FORMAT(12X,'Max iteration number of Lis: ',I8,' default(2000)')         
911 FORMAT(12X,'Convergence tolerance of Lis:   ',E12.5)           
1011 FORMAT(12X,'Number of iteration of Lis:   ',I8)        
1012 FORMAT(12X,'CPU time of Lis:          ',F12.5,' mins')      
1013 FORMAT(12X,'Relative residual norm of Lis:',E12.5)    

821 FORMAT(9X,'Max iteration number of Lis: ',I8,' default(2000)')         
921 FORMAT(9X,'Convergence tolerance of Lis:   ',E12.5)           
1021 FORMAT(9X,'Number of iteration of Lis:   ',I8)        
1022 FORMAT(9X,'CPU time of Lis:          ',F12.5,' mins')      
1023 FORMAT(9X,'Relative residual norm of Lis:',E12.5)     

831 FORMAT(7X,'Max iteration number of Lis: ',I8,' default(2000)')         
931 FORMAT(7X,'Convergence tolerance of Lis:   ',E12.5)           
1031 FORMAT(7X,'Number of iteration of Lis:   ',I8)        
1032 FORMAT(7X,'CPU time of Lis:          ',F12.5,' mins')      
1033 FORMAT(7X,'Relative residual norm of Lis:',E12.5)  

call lis_initialize(lis_error)


nprocs  = 1
my_rank = 0

n = input_n

ln = 0
matrix_type = LIS_MATRIX_CSR
comm_world = LIS_COMM_WORLD


call lis_matrix_create(comm_world,A,lis_error)
call lis_matrix_set_size(A,ln,n,lis_error)
call lis_matrix_get_size(A,n,gn,lis_error)



A_Lis_Ptr = K_Ptr
A_Lis_Index = K_Index
A_Lis_Value = K_Value


call lis_matrix_set_csr(NZ_NUM, A_Lis_Ptr,  &
            A_Lis_Index, A_Lis_Value, A,lis_error)



call lis_matrix_assemble(A,lis_error)
call lis_vector_duplicate(A,u,lis_error)
call lis_vector_duplicate(A,x,lis_error)
do c_i=1,input_n
  call lis_vector_set_value(LIS_INS_VALUE, c_i, Input_F(c_i),u,lis_error)
enddo
call lis_solver_create(solver,lis_error)
call lis_solver_set_option("-i bicg",solver,lis_error);
call lis_solver_set_option("-p none",solver,lis_error);

if(Lis_Number_Iteration==500)then
  call lis_solver_set_option("-maxiter 500",solver,lis_error) 
elseif(Lis_Number_Iteration==1000)then
  call lis_solver_set_option("-maxiter 1000",solver,lis_error)
elseif (Lis_Number_Iteration==2000)then
  call lis_solver_set_option("-maxiter 2000",solver,lis_error)
elseif (Lis_Number_Iteration==3000)then
  call lis_solver_set_option("-maxiter 3000",solver,lis_error)
elseif (Lis_Number_Iteration==4000)then
  call lis_solver_set_option("-maxiter 4000",solver,lis_error)
elseif (Lis_Number_Iteration==5000)then
  call lis_solver_set_option("-maxiter 5000",solver,lis_error)
elseif (Lis_Number_Iteration==6000)then
  call lis_solver_set_option("-maxiter 6000",solver,lis_error)
elseif (Lis_Number_Iteration==7000)then
  call lis_solver_set_option("-maxiter 7000",solver,lis_error)
elseif (Lis_Number_Iteration==8000)then
  call lis_solver_set_option("-maxiter 8000",solver,lis_error)
elseif (Lis_Number_Iteration==9000)then
  call lis_solver_set_option("-maxiter 9000",solver,lis_error)      
else
  call lis_solver_set_option("-maxiter 2000",solver,lis_error)
endif


if(Key_Indent==0)then
  write(*,801) Lis_Number_Iteration
elseif(Key_Indent==7)then
  write(*,811) Lis_Number_Iteration
elseif(Key_Indent==4)then
  write(*,821) Lis_Number_Iteration  
elseif(Key_Indent==2)then
  write(*,831) Lis_Number_Iteration            
endif      


call lis_solver_set_option("-tol 1.0e-6",solver,lis_error)
if(Key_Indent==0)then
  write(*,901) 1.0e-6
elseif(Key_Indent==7)then
  write(*,911) 1.0e-6
elseif(Key_Indent==4)then
  write(*,921) 1.0e-6      
elseif(Key_Indent==2)then
  write(*,931) 1.0e-6           
endif



if(Key_Num_Process==1) then
call lis_solver_set_option("-omp_num_threads 1",solver,lis_error);           
elseif(Key_Num_Process==2) then
call lis_solver_set_option("-omp_num_threads 2",solver,lis_error);    
elseif(Key_Num_Process==3) then
call lis_solver_set_option("-omp_num_threads 3",solver,lis_error);  
elseif(Key_Num_Process==4) then
call lis_solver_set_option("-omp_num_threads 4",solver,lis_error);  
elseif(Key_Num_Process==5) then
call lis_solver_set_option("-omp_num_threads 5",solver,lis_error);  
elseif(Key_Num_Process==6) then
call lis_solver_set_option("-omp_num_threads 6",solver,lis_error);  
elseif(Key_Num_Process==7) then
call lis_solver_set_option("-omp_num_threads 7",solver,lis_error);  
elseif(Key_Num_Process==8) then
call lis_solver_set_option("-omp_num_threads 8",solver,lis_error);  
elseif(Key_Num_Process==9) then
call lis_solver_set_option("-omp_num_threads 9",solver,lis_error);  
elseif(Key_Num_Process==10) then
call lis_solver_set_option("-omp_num_threads 10",solver,lis_error);  
elseif(Key_Num_Process==11) then
call lis_solver_set_option("-omp_num_threads 11",solver,lis_error);  
elseif(Key_Num_Process==12) then
call lis_solver_set_option("-omp_num_threads 12",solver,lis_error);  
elseif(Key_Num_Process==13) then
call lis_solver_set_option("-omp_num_threads 13",solver,lis_error);  
elseif(Key_Num_Process==14) then
call lis_solver_set_option("-omp_num_threads 14",solver,lis_error);  
elseif(Key_Num_Process==15) then
call lis_solver_set_option("-omp_num_threads 15",solver,lis_error);  
elseif(Key_Num_Process==16) then
call lis_solver_set_option("-omp_num_threads 16",solver,lis_error);  
elseif(Key_Num_Process==17) then
call lis_solver_set_option("-omp_num_threads 17",solver,lis_error);  
elseif(Key_Num_Process==18) then
call lis_solver_set_option("-omp_num_threads 18",solver,lis_error);  
elseif(Key_Num_Process==19) then
call lis_solver_set_option("-omp_num_threads 19",solver,lis_error);  
elseif(Key_Num_Process==20) then
call lis_solver_set_option("-omp_num_threads 20",solver,lis_error);  
elseif(Key_Num_Process==21) then
call lis_solver_set_option("-omp_num_threads 21",solver,lis_error);  
elseif(Key_Num_Process==22) then
call lis_solver_set_option("-omp_num_threads 22",solver,lis_error);  
elseif(Key_Num_Process==23) then
call lis_solver_set_option("-omp_num_threads 23",solver,lis_error);  
elseif(Key_Num_Process==24) then
call lis_solver_set_option("-omp_num_threads 24",solver,lis_error);       
else
call lis_solver_set_option("-omp_num_threads [t]",solver,lis_error);
endif



call lis_solve(A,u,x,solver,lis_error);

call lis_solver_get_iter(solver,iter,lis_error);

if(Key_Indent==0)then
  write(*,1001) iter-1
elseif(Key_Indent==7)then
  write(*,1011) iter-1
elseif(Key_Indent==4)then
  write(*,1021) iter-1    
elseif(Key_Indent==2)then
  write(*,1031) iter-1             
endif      

call lis_solver_get_time(solver,cputime,lis_error)
if(Key_Indent==0)then
  write(*,1002) cputime/60.0 
elseif(Key_Indent==7)then
  write(*,1012) cputime/60.0 
elseif(Key_Indent==4)then
  write(*,1022) cputime/60.0      
elseif(Key_Indent==2)then
  write(*,1032) cputime/60.0           
endif           
  

call lis_solver_get_residualnorm(solver,residual,lis_error)
if(Key_Indent==0)then
  write(*,1003) residual
elseif(Key_Indent==7)then
  write(*,1013) residual
elseif(Key_Indent==4)then
  write(*,1023) residual     
elseif(Key_Indent==2)then
  write(*,1033) residual            
endif           
 

do c_i=1,input_n
  call lis_vector_get_value(x,c_i,output_D(c_i), lis_error)
enddo


END SUBROUTINE Matrix_Solve_LSOE_Lis   
