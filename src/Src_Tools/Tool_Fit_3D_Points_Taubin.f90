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
 
subroutine Tool_Fit_3D_Points_Taubin(i_C,num_Point,In_Points,Points_Flag_Inside,Out_Points)
! Smooth treatment of crack front edge. IMPROV2022121001.
! Taubin algorithm. Key_Smooth_Front = 6.
! Ref: \theory_documents\042.2 Taubin_1995_Curve and surface smoothing without shrinkage.pdf
! Ref: \theory_documents\042 Mesh Smoothing-2022-12-10.pdf, P15-P25.
! Corresponding Matlab test code: \theory_documents\042.1 crack_front_smooth_test.m
!
! Process only the points that are inside the model. Added this feature on 2026-02-01.
! IMPROV-2026020101.
!
!

!......................
! Variable Declaration
!......................
use Global_Float_Type     
use Global_Model

implicit none
integer,intent(in)::i_C,num_Point
real(kind=FT),intent(in)::In_Points(num_Point,3)
logical,intent(in)::Points_Flag_Inside(num_Point)
real(kind=FT),intent(out)::Out_Points(num_Point,3)
real(kind=FT) point(num_Point+1,3)
integer n_Loops,i_Loop,i_point
real(kind=FT) lambda,mu,c_point(3),pre_point(3),nex_point(3)
real(kind=FT) pre_L,nex_L,weight_pre,weight_nex
real(kind=FT) L1_vector(num_Point,3),L2_vector(num_Point,3)


n_Loops =    20


lambda  =  0.330D0
mu      = -0.331D0

point(1:num_Point,1:3) = In_Points(1:num_Point,1:3)
point(num_point+1,1:3) = In_Points(1,1:3)
L1_vector(1:num_Point,1:3) = ZR
L2_vector(1:num_Point,1:3) = ZR

      
do i_Loop=1,n_Loops
    do i_point=1,num_point
        if (.not. Points_Flag_Inside(i_point)) cycle
        
        c_point = point(i_point,1:3)
        nex_point = point(i_point+1,1:3)
        if (i_point>=2) then
            pre_point = point(i_point-1,1:3)
        else
            pre_point = point(num_point,1:3)
        endif
		

        pre_L = sqrt((pre_point(1)-c_point(1))**2+(pre_point(2)-c_point(2))**2+(pre_point(3)-c_point(3))**2)
        nex_L = sqrt((nex_point(1)-c_point(1))**2+(nex_point(2)-c_point(2))**2+(nex_point(3)-c_point(3))**2)
        weight_pre = ONE/pre_L
        weight_nex = ONE/nex_L
        L1_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point
    enddo
    
    do i_point=1,num_point
        if (.not. Points_Flag_Inside(i_point)) cycle
        
        point(i_point,1:3) = point(i_point,1:3) + lambda*L1_vector(i_point,1:3)
    enddo
    point(num_point+1,1:3) = point(1,1:3)
    
    do i_point=1,num_point
        if (.not. Points_Flag_Inside(i_point)) cycle
        
        c_point = point(i_point,1:3)
        nex_point = point(i_point+1,1:3)
        if (i_point>=2) then
            pre_point = point(i_point-1,1:3)
        else
            pre_point = point(num_point,1:3)
        endif
		

        pre_L = sqrt((pre_point(1)-c_point(1))**2+(pre_point(2)-c_point(2))**2+(pre_point(3)-c_point(3))**2)
        nex_L = sqrt((nex_point(1)-c_point(1))**2+(nex_point(2)-c_point(2))**2+(nex_point(3)-c_point(3))**2)
        weight_pre = ONE/pre_L
        weight_nex = ONE/nex_L
        L2_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point
		
        
    enddo
    
    do i_point=1,num_point
        if (.not. Points_Flag_Inside(i_point)) cycle
        
        point(i_point,1:3) = point(i_point,1:3)  + mu*L2_vector(i_point,1:3)
    enddo
    point(num_point+1,1:3) = point(1,1:3)
    
enddo

Out_Points(1:num_Point,1:3)  = point(1:num_Point,1:3) 


return 
end SUBROUTINE Tool_Fit_3D_Points_Taubin             
