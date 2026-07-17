!-----------------------------------------------------------
! Brief: Apply 2D displacement boundary conditions and
!        return both free and fixed DOF lists for output.
!
! Parameters:
!   Input:  c_Total_FD - total DOF count (2D)
!           isub       - current load step index
!   Output: freeDOF, c_num_FD     - free DOFs and count
!           fixedDOF, num_FixedD  - fixed DOFs and count
!
! Notes:   Also constrains nodes inside hole boundaries
!          and supports circular/elliptical hole geometry
!          on the first fracturing step.
!-----------------------------------------------------------

subroutine Boundary_Cond_output_Fixed_Freedom(c_Total_FD,isub, freeDOF,c_num_FD, fixedDOF,num_FixedD)
!     Boundary conditions.
! Total_FD: Total Degrees of Freedom
! c_num_FD: Degrees of freedom after removing constraints

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_Crack_Common


!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none

integer,intent(in)::isub
integer,intent(in)::c_Total_FD
integer,intent(out)::freeDOF(c_Total_FD),c_num_FD
integer,intent(out)::fixedDOF(c_Total_FD),num_FixedD

integer i,i_E,i_N,cur_Node,i_fd
integer :: c_index
integer N1,N2,N3,N4,NN(4),NN_L(5)
integer :: i_H
real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r
logical :: Nodes_4_in(4)
real(kind=FT) :: c_Dis
real(kind=FT) :: Tool_Function_2Point_Dis
integer :: num_in_Nodes
real(kind=FT) x0_Hole,y0_Hole,R0_Hole
real(kind=FT) a_Hole,b_Hole,theta_Hole,Tol  
integer :: Return_Statu

print *,'    Dealing with boundary condition...'   

!.................................
!Initialize vector of fixed DOFs.
!.................................
freeDOF(1:c_Total_FD) = 0
fixedDOF(1:c_Total_FD) = 0
c_index = 0

!.................................
!Fix displacement in x-direction.
!.................................
do i = 1,Num_Bou_x
    cur_Node = Bou_x(i)
    c_index = c_index+1     
    fixedDOF(c_index) = 2*cur_Node-1         
end do   

!.................................
!Fix displacement in y-direction.
!.................................
do i = 1,Num_Bou_y
    cur_Node = Bou_y(i)
    c_index = c_index+1  
    fixedDOF(c_index) = 2*cur_Node        
end do  

!............................................................................................
! Nodes of the hole, elements inside the hole that are not intersected by the hole boundary
! The node is set as a constrained degree of freedom, meaning it does not participate in the
! calculation.
!2017-05-29
! Execute only on the first rupture step
!............................................................................................
! Circular hole
if (num_Hole /= 0 .and. num_Circ_Hole>=1) then        
    ! Cycle between elements
    do i_E = 1,Num_Elem
        N1  = Elem_Node(i_E,1)                                         
        N2  = Elem_Node(i_E,2)                                             
        N3  = Elem_Node(i_E,3)                                             
        N4  = Elem_Node(i_E,4)  
        NN  = [N1,N2,N3,N4]
        NN_L= [N1,N2,N3,N4,N1]
        !X_NODES_L = Coor(NN_L,1)
        !Y_NODES_L = Coor(NN_L,2)
        !**************************
        ! Circulation of each hole
        !**************************
        do i_H = 1,num_Hole
            c_Hole_x = Hole_Coor(i_H,1)
            c_Hole_y = Hole_Coor(i_H,2)
            c_Hole_r = Hole_Coor(i_H,3)
            ! Loop through the 4 nodes of the current element
            Nodes_4_in(1:4) = .False.
            do i_N = 1,4
                c_Dis=Tool_Function_2Point_Dis(Coor(NN(i_N),1:2), &
                Hole_Coor(i_H,1:2))
                if (c_Dis<=c_Hole_r)then
                    Nodes_4_in(i_N) = .True.

                endif
            enddo
            num_in_Nodes = count(Nodes_4_in(1:4))
            if(num_in_Nodes == 4)then
                do i_N=1,4
                    c_index = c_index+1  
                    fixedDOF(c_index) = 2*NN(i_N) -1
                    c_index = c_index+1  
                    fixedDOF(c_index) = 2*NN(i_N)
                enddo
            endif
            !Enriched_Node_Type_HL(1:Num_Node,1:Max_Num_Hl)
        enddo
    enddo
endif
! Elliptical hole, 2020-08-09
if (num_Hole /= 0 .and. num_Ellip_Hole>=1) then   
    ! Cycle between each element
    do i_E = 1,Num_Elem
        N1  = Elem_Node(i_E,1)                                         
        N2  = Elem_Node(i_E,2)                                             
        N3  = Elem_Node(i_E,3)                                             
        N4  = Elem_Node(i_E,4)  
        NN  = [N1,N2,N3,N4]
        NN_L= [N1,N2,N3,N4,N1]
        !X_NODES_L = Coor(NN_L,1)
        !Y_NODES_L = Coor(NN_L,2)
        !**************************
        ! Circulation of each hole
        !**************************
        do i_H = 1,num_Hole
            x0_Hole =  Ellip_Hole_Coor(i_H,1)
            y0_Hole =  Ellip_Hole_Coor(i_H,2)
            a_Hole  =  Ellip_Hole_Coor(i_H,3) 
            b_Hole  =  Ellip_Hole_Coor(i_H,4) 
            theta_Hole =  Ellip_Hole_Coor(i_H,5)
            ! Loop through the 4 nodes of the current element
            Nodes_4_in(1:4) = .False.
            do i_N = 1,4
                Tol = 1.0D-10
                call Tool_Yes_Point_in_Oblique_Ellipse( Coor(NN(i_N),1:2), x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole, &
                Return_Statu,Tol)
                if (Return_Statu==1)then
                    Nodes_4_in(i_N) = .True.

                endif
            enddo
            num_in_Nodes = count(Nodes_4_in(1:4))
            if(num_in_Nodes == 4)then
                do i_N=1,4
                    c_index = c_index+1  
                    fixedDOF(c_index) = 2*NN(i_N) -1
                    c_index = c_index+1  
                    fixedDOF(c_index) = 2*NN(i_N)
                enddo
            endif
            !Enriched_Node_Type_HL(1:Num_Node,1:Max_Num_Hl)
        enddo
    enddo
endif

!NEWFTU2024110701.
if(allocated(Flag_FreeDOF)) deallocate(Flag_FreeDOF)
allocate(Flag_FreeDOF(c_Total_FD))
Flag_FreeDOF(1:c_Total_FD) = 0

!IMPROV2024110703.
if(allocated(Real_FreeDOF_Index)) deallocate(Real_FreeDOF_Index)
allocate(Real_FreeDOF_Index(c_Total_FD)) 
Real_FreeDOF_Index(1:c_Total_FD) = 0

c_num_FD = 0        
num_FixedD = c_index
do i_fd =1,c_Total_FD
    if ( .not.(ANY( fixedDOF(1:c_index) .eq. i_fd)) ) then
        c_num_FD = c_num_FD +1
        freeDOF(c_num_FD) = i_fd
        Flag_FreeDOF(i_fd) = 1
        Real_FreeDOF_Index(i_fd) = c_num_FD
    end if
end do


return
END subroutine Boundary_Cond_output_Fixed_Freedom
