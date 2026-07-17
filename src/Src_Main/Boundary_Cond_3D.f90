!-----------------------------------------------------------
! Brief: Apply 3D displacement boundary conditions and
!        build the list of free DOFs and location map.
!
! Parameters:
!   Input:  c_Total_FD - total DOF count (3D)
!           isub       - current load step index
!   Output: freeDOF    - indices of free DOFs
!           num_FreeD  - number of free DOFs
!           fixedDOF   - indices of constrained DOFs
!           num_FixedD - number of constrained DOFs
!
! Notes:   Populates the global Flag_FreeDOF and
!          Location_FreeDOF arrays used by sparse solvers.
!-----------------------------------------------------------

subroutine Boundary_Cond_3D(c_Total_FD,isub, freeDOF,num_FreeD, fixedDOF,num_FixedD)
!     Boundary Conditions (Three-Dimensional Problem)
! Total_FD: Total Degrees of Freedom
! num_FreeD: Degrees of freedom after removing constraints

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none

integer,intent(in)::isub
integer,intent(in)::c_Total_FD
integer,intent(out)::freeDOF(c_Total_FD),num_FreeD
integer,intent(out)::fixedDOF(c_Total_FD),num_FixedD
integer i,cur_Node,i_fd
integer :: c_index
!fixedDOF(c_Total_FD)
integer :: temp_flag(c_Total_FD)
integer :: Location_Count

print *,'    Dealing with boundary condition...'   

!Initialize vectors.
freeDOF(1:c_Total_FD)  = 0
fixedDOF(1:c_Total_FD) = 0

c_index = 0

temp_flag(1:c_Total_FD) = 0

!Fix displacement in x-direction.
do i = 1,Num_Bou_x
    cur_Node = Bou_x(i)
    c_index = c_index+1     
    fixedDOF(c_index) = 3*cur_Node-2     
    temp_flag(3*cur_Node-2) = 1
end do   

!Fix displacement in y-direction.
do i = 1,Num_Bou_y
    cur_Node = Bou_y(i)
    c_index = c_index+1  
    fixedDOF(c_index) = 3*cur_Node-1     
    temp_flag(3*cur_Node-1) = 1
end do  

!Fix displacement in z-direction.
do i = 1,Num_Bou_z
    cur_Node = Bou_z(i)
    c_index = c_index+1  
    fixedDOF(c_index) = 3*cur_Node  
    temp_flag(3*cur_Node) = 1
end do        

num_FreeD  = 0    
num_FixedD = c_index

! Mark whether the degree of freedom is free. 2022-06-12. IMPROV2022061203.
if (allocated(Flag_FreeDOF)) deallocate(Flag_FreeDOF)
allocate(Flag_FreeDOF(c_Total_FD))
Flag_FreeDOF(1:c_Total_FD) = 0

! Added a global variable. 2022-10-25. Used to mark the position of the degree of freedom
! corresponding to the free degree of freedom. IMPROV2022102501.
if (allocated(Location_FreeDOF)) deallocate(Location_FreeDOF)
allocate(Location_FreeDOF(c_Total_FD))
Location_FreeDOF(1:c_Total_FD) = 0

! Get the list of degrees of freedom.
!=====================================================================
! OPTION 1. The efficiency is low due to the use of the any function.
!=====================================================================
!     do i_fd =1,c_Total_FD
!         if (.not.(ANY(fixedDOF(1:num_FixedD) .eq. i_fd))) then
!             num_FreeD = num_FreeD +1
!             freeDOF(num_FreeD) = i_fd  
! !Marked as free degree of freedom. 2022-06-12. IMPROV2022061203.
!             Flag_FreeDOF(i_fd) = 1 
!         end if
!     end do

!=============================================================================
! OPTION 2. Significant Efficiency Improvement. 2022062501. IMPROV2022062501.
!=============================================================================
!IMPROV2024110703.
!if(allocated(Real_FreeDOF_Index)) deallocate(Real_FreeDOF_Index)
!allocate(Real_FreeDOF_Index(c_Total_FD)) 
!Real_FreeDOF_Index(1:c_Total_FD) = 0

Location_Count = 0
do i_fd =1,c_Total_FD          
    if(temp_flag(i_fd) == 0)then
        num_FreeD = num_FreeD +1
        freeDOF(num_FreeD) = i_fd  
        Flag_FreeDOF(i_fd) = 1 

        !IMPROV2022102501.
        Location_Count = Location_Count +1
        Location_FreeDOF(i_fd) = Location_Count
        !Real_FreeDOF_Index(i_fd) = Location_Count !IMPROV2024110703.

    endif
end do      




return
END subroutine Boundary_Cond_3D
