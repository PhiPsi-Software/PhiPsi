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
 
subroutine Tool_Sort_by_End_to_End(m,m_OP,Input_Outline,Output_Outline,cou)
! Make the outline connect end to end
! m: Number of rows in the Outline matrix
! m_OP: Number of rows to operate on

use Global_Float_Type
implicit none
integer,intent(in)::m,m_OP 
integer,intent(in)::Input_Outline(m,2)
integer,intent(out)::cou 
integer,intent(out)::Output_Outline(m,2)
      
integer i,j,c_Outline,c_node_num_2
logical tem_Yes
integer Location



Output_Outline = 0

Output_Outline(1,:)  = Input_Outline(1,:)


cou = 1
c_Outline    = 1

do i = 1,m_OP
    if (i.eq.1)then
        c_Outline = 1
        c_node_num_2 = Input_Outline(1,2)
    end if
    
    do j = 1,m
        if(j .ne. c_Outline) then
            if (c_node_num_2 .eq. Input_Outline(j,1)) then
                call Vector_belongs_Matrix_Is_Int(cou,2,Output_Outline,[Input_Outline(j,1),Input_Outline(j,2)],Location,tem_Yes)
                if(tem_Yes.eqv..False.) then
                    cou = cou + 1
                    Output_Outline(cou,:) = [Input_Outline(j,1),Input_Outline(j,2)]
                    c_Outline = j
                    c_node_num_2 = Input_Outline(j,2)
                    exit
               end if
            else if(c_node_num_2 .eq. Input_Outline(j,2)) then
                call Vector_belongs_Matrix_Is_Int(cou,2,Output_Outline,&
                    [Input_Outline(j,2),Input_Outline(j,1)],Location,tem_Yes)
                if(tem_Yes.eqv..False.)then
                   cou = cou + 1
                   Output_Outline(cou,:)=[Input_Outline(j,2),Input_Outline(j,1)]
                   c_Outline = j
                   c_node_num_2 = Input_Outline(j,1)
                   exit
                  end if
             end if
        end if
     end do
end do


RETURN
end subroutine Tool_Sort_by_End_to_End
