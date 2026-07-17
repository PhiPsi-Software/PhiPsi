!-----------------------------------------------------------
! Brief: Partition all elements into 8 octant sub-domains.
!
! Parameters:
!   (none)
!
! Notes:   Uses the model center and average element length as a
!          tolerance to assign each element to one of the 2^3=8
!          octants, populating Ele_Domain_ID and Domain_Elements_Num.
!-----------------------------------------------------------

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