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
 
SUBROUTINE D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Overlap_Status)
! Check whether there is an overlap in the coordinate ranges between each crack. NEWFTU2022050402.
! There is an overlap only if x, y, and z all intersect.
!2022-05-04.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use omp_lib      
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
logical,intent(out):: Crack_Coor_Overlap_Status(num_Crack,num_Crack)
integer i_C,j_C
logical x_Logical_Yes,y_Logical_Yes,z_Logical_Yes


!------------
! Initialize
!------------
Crack_Coor_Overlap_Status(1:num_Crack,1:num_Crack) = .False.

!-----------------------------------------------
! If there is only one crack, exit immediately.
!-----------------------------------------------
if(num_Crack==1) then
    return
endif

!-------------
! Crack cycle
!-------------
! IMPROV2022072203. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,j_C,x_Logical_Yes,y_Logical_Yes,z_Logical_Yes) SCHEDULE(static) 
do i_C =1,num_Crack
    do j_C =1,num_Crack
        if(i_C/=j_C)then
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,1,1:2),Crack_Coor_Range(j_C,1,1:2),x_Logical_Yes)  
            if (x_Logical_Yes .eqv. .False.) then
                cycle
            endif
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,2,1:2),Crack_Coor_Range(j_C,2,1:2),y_Logical_Yes)  
            if (y_Logical_Yes .eqv. .False.) then
                cycle
            endif
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,3,1:2),Crack_Coor_Range(j_C,3,1:2),z_Logical_Yes)  
            if (z_Logical_Yes .eqv. .False.) then
                cycle
            endif  
            Crack_Coor_Overlap_Status(i_C,j_C) = .True.
        endif
    enddo    
enddo
!$omp end parallel do         

RETURN
END SUBROUTINE D3_Get_Cracks_Coor_Overlap_Status
