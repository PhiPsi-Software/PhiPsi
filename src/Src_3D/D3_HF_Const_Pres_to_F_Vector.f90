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
 
SUBROUTINE D3_HF_Const_Pres_to_F_Vector(isub,c_Pres,num_FreeD,  &
                           in_Total_FD,in_num_Tol_CalP_Water,freeDOF,in_F_U,F)

! Convert constant water pressure to F vector.
! It turns out this part of the code is in the main program.
!2022-05-31.    

!#############################
! Read public variable module
!#############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Function_MVMUL_LP 

!###########################
! Variable Type Declaration
!###########################
implicit none
integer,intent(in)::isub,num_FreeD,in_num_Tol_CalP_Water,in_Total_FD
real(kind=FT),intent(in)::c_Pres,in_F_U(in_Total_FD)
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(out)::F(in_Total_FD)
real(kind=FT) temp(num_FreeD)
integer i_Dof
!--------------------------------------------------
real(kind=FT) CalP_Pres(in_num_Tol_CalP_Water)
integer  num_Count,i_C, i_CalP      
!**************************************************************************
! Convert water pressure into a pressure vector for all calculation points
!**************************************************************************
num_Count = 0
do i_C = 1,num_Crack    
    ! If it is an HF fracture
    if(Crack_Type_Status_3D(i_C,1) == 1) then
    do i_CalP=1,Cracks_Real_CalP_Num_3D(i_C)   
        num_Count = num_Count + 1
        CalP_Pres(num_Count) = c_Pres -Cracks_FluidEle_CalP_Glo_Insitu(num_Count)
    end do                     
    endif
end do   

!*****************************************************
! Assemble the load vector considering fluid pressure
!*****************************************************
! KNISUE2022081602. MATMUL may occupy more memory (occurs with the gfortran compiler).
!//////////////////////////////////////////////////////////
! OPTION 1: Minimal memory usage (Intel Fortran compiler).
!//////////////////////////////////////////////////////////
! Using the gfortran compiler, MATMUL consumes a lot of memory. 2022-08-17.
! F(freeDOF(1:num_FreeD)) =  in_F_U(freeDOF(1:num_FreeD))+   &
!            MATMUL(Coupled_Q_3D(freeDOF(1:num_FreeD),1:in_num_Tol_CalP_Water), &
!        CalP_Pres(1:in_num_Tol_CalP_Water))

!///////////
!  OPTION 2
!///////////
!C     F(freeDOF(1:num_FreeD)) =  in_F_U(freeDOF(1:num_FreeD))+  
!C & MVMUL_LP(Coupled_Q_3D(freeDOF(1:num_FreeD), !Using LAPACK, but Intel compiler could not run it
!(fixed on 2022-08-17).
!C    &                     1:in_num_Tol_CalP_Water),
!C    &          CalP_Pres(1:in_num_Tol_CalP_Water)) 

!//////////////////////////////
!  OPTION 3
!2022-08-17. IMPROV2022081702.
!//////////////////////////////
!A(m,k),B(k,n),Output(m,n)      
!C     CALL DGEMM('N','N',num_FreeD,1,in_num_Tol_CalP_Water,ONE,
!C    &       Coupled_Q_3D(freeDOF(1:num_FreeD),1:in_num_Tol_CalP_Water),
!C    &       num_FreeD,CalP_Pres(1:in_num_Tol_CalP_Water),
!C    &       in_num_Tol_CalP_Water,ZR,temp(1:num_FreeD),num_FreeD)
!C     F(freeDOF(1:num_FreeD)) =  in_F_U(freeDOF(1:num_FreeD))+
!C    &                           temp(1:num_FreeD)

!///////////////////////////////////////////////////
! OPTION 4: Use staggered arrays. IMPROV2022091601.
!///////////////////////////////////////////////////
!F(freeDOF(1:num_FreeD)) =  in_F_U(freeDOF(1:num_FreeD))+   &
!            MATMUL(Coupled_Q_3D(freeDOF(1:num_FreeD),1:in_num_Tol_CalP_Water), &
!        CalP_Pres(1:in_num_Tol_CalP_Water))
! Modified with reference to the above Option 1.

do i_Dof =1,num_FreeD
    if(allocated(Coupled_Q_3D(freeDOF(i_Dof))%row))then 
        !BUGFIX2022111101.                        
        F(freeDOF(i_Dof)) =  in_F_U(freeDOF(i_Dof))+   &
                    dot_product(Coupled_Q_3D(freeDOF(i_Dof))%row(1:3), &
                                CalP_Pres(Coupled_Q_3D_Index(freeDOF(i_Dof))%row(1:3)))
    endif
enddo


RETURN
END SUBROUTINE D3_HF_Const_Pres_to_F_Vector
