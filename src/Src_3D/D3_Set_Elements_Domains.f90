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
 
subroutine D3_Set_Elements_Domains
! Used for dividing element areas.
! 2022-11-24. NEWFTU2022112402.


!***************
! Public Module
!***************
use Global_Float_Type
use Global_Elem_Area_Vol
use Global_Model
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
real(kind=FT) x_Center,y_Center,z_Center
real(kind=FT) c_E_x,c_E_y,c_E_z,c_Tol
integer,ALLOCATABLE::Tem_Elements(:,:)
integer i_Domain
integer i_E

!********************************************
! Set element partitions: 2^3 = 8 partitions
!********************************************
Domain_Elements_Num(1:8)  = 0
c_Tol    = Ave_Elem_L*TWO
x_Center = Model_Center_x
y_Center = Model_Center_y
z_Center = Model_Center_z
allocate(Tem_Elements(8,Num_Elem))  
Tem_Elements(1:8,1:Num_Elem) = 0
allocate(Ele_Domain_ID(Num_Elem))  
Ele_Domain_ID(1:Num_Elem)    = 0
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E_x,c_E_y,c_E_z) schedule(static) 
! Parallel conflict.
do i_E =1,Num_Elem  
    c_E_x = Elem_Centroid(i_E,1)
    c_E_y = Elem_Centroid(i_E,2)
    c_E_z = Elem_Centroid(i_E,3)
    if( ((c_E_x - x_Center) > -c_Tol) .and. ((c_E_y - y_Center) > -c_Tol) .and. ((c_E_z - z_Center) > -c_Tol)) then
        Ele_Domain_ID(i_E)  = 1
        Domain_Elements_Num(1) = Domain_Elements_Num(1) +1
        Tem_Elements(1,Domain_Elements_Num(1))  = i_E
    endif
    if( ((c_E_x - x_Center) <  c_Tol) .and. ((c_E_y - y_Center) > -c_Tol) .and. ((c_E_z - z_Center) > -c_Tol)) then
        Ele_Domain_ID(i_E)  = 2
        Domain_Elements_Num(2) = Domain_Elements_Num(2) +1
        Tem_Elements(2,Domain_Elements_Num(2))  = i_E
    endif
    if( ((c_E_x - x_Center) <  c_Tol) .and. ((c_E_y - y_Center) < c_Tol) .and. ((c_E_z - z_Center) > -c_Tol)) then
        Ele_Domain_ID(i_E)  = 3
        Domain_Elements_Num(3) = Domain_Elements_Num(3) +1
        Tem_Elements(3,Domain_Elements_Num(3))  = i_E
    endif
    if( ((c_E_x - x_Center) > -c_Tol) .and. ((c_E_y - y_Center) < c_Tol) .and. ((c_E_z - z_Center) > -c_Tol)) then
        Ele_Domain_ID(i_E)  = 4
        Domain_Elements_Num(4) = Domain_Elements_Num(4) +1
        Tem_Elements(4,Domain_Elements_Num(4))  = i_E
    endif
    if( ((c_E_x - x_Center) > -c_Tol) .and. ((c_E_y - y_Center) > -c_Tol) .and. ((c_E_z - z_Center) <  c_Tol)) then
        Ele_Domain_ID(i_E)  = 5
        Domain_Elements_Num(5) = Domain_Elements_Num(5) +1
        Tem_Elements(5,Domain_Elements_Num(5))  = i_E
    endif
    if( ((c_E_x - x_Center) <  c_Tol) .and. ((c_E_y - y_Center) > -c_Tol) .and. ((c_E_z - z_Center) <  c_Tol)) then
        Ele_Domain_ID(i_E)  = 6
        Domain_Elements_Num(6) = Domain_Elements_Num(6) +1
        Tem_Elements(6,Domain_Elements_Num(6))  = i_E
    endif
    if( ((c_E_x - x_Center) <  c_Tol) .and. ((c_E_y - y_Center) < c_Tol) .and. ((c_E_z - z_Center) <  c_Tol)) then
        Ele_Domain_ID(i_E)  = 7
        Domain_Elements_Num(7) = Domain_Elements_Num(7) +1
        Tem_Elements(7,Domain_Elements_Num(7))  = i_E
    endif
    if( ((c_E_x - x_Center) > -c_Tol) .and. ((c_E_y - y_Center) < c_Tol) .and. ((c_E_z - z_Center) <  c_Tol)) then
        Ele_Domain_ID(i_E)  = 8
        Domain_Elements_Num(8) = Domain_Elements_Num(8) +1
        Tem_Elements(8,Domain_Elements_Num(8))  = i_E
    endif

enddo
!!$OMP END PARALLEL DO   

allocate(Domain_Elements(8,maxval(Domain_Elements_Num(1:8))))
do i_Domain = 1,8
    Domain_Elements(i_Domain,1:Domain_Elements_Num(i_Domain)) = Tem_Elements(i_Domain,1:Domain_Elements_Num(i_Domain)) 
enddo

deallocate(Tem_Elements)

    
RETURN
END subroutine D3_Set_Elements_Domains