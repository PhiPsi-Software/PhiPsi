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
 
SUBROUTINE D3_Get_Element_List_of_Natural_Crack(c_Crack,Out_num_Element,Out_Element_List)
! Used to obtain the list of solid elements containing natural fractures. NEWFTU2023011203.
! Only quadrilateral and polygonal planar natural fractures are supported.
! 2023-01-12.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::c_Crack
integer,intent(out)::Out_num_Element
integer,intent(out)::Out_Element_List(10000)
!
real(kind=FT) c_NF_Coor_Ranges(3,2),eight_points(8,3)
integer Domain_Numbers(8),Uniqued_Domain_Numbers(8)
integer,ALLOCATABLE::Domain_List(:)
integer Num_Domain
integer i_Point,i_Domain,c_Domain,i_El,c_Ele
integer c_count
real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min,c_z_max,c_z_min
logical c_Logical_Yes
integer,ALLOCATABLE::Potent_Elems(:),Unique_Potent_Elems(:)
integer Num_Potent_Elems
integer c_NN(8)
integer i_Edge
real(kind=FT) A(3),B(3)
real(kind=FT) c_InterSection_P(3)
logical c_Yes_Inter
integer tem_Out_num_Element
integer tem_Out_Element_List(10000)
integer c_Uniqued_count

!**************************
! Variable initialization.
!**************************
Out_num_Element = 0
Out_Element_List(1:10000) = 0
tem_Out_num_Element = 0
tem_Out_Element_List(1:10000) = 0

!***************************************************
! Obtain the coordinate range of natural fractures.
!***************************************************
call D3_Get_NF_Coor_Ranges(c_Crack,c_NF_Coor_Ranges)

!**********************************************************************************
! The fracture zone number where the natural fracture coordinate range is located.
!**********************************************************************************
! Eight points determined by the range of crack coordinates
eight_points(1,1:3) = [c_NF_Coor_Ranges(1,1),c_NF_Coor_Ranges(2,1),c_NF_Coor_Ranges(3,1)]
eight_points(2,1:3) = [c_NF_Coor_Ranges(1,2),c_NF_Coor_Ranges(2,1),c_NF_Coor_Ranges(3,1)]
eight_points(3,1:3) = [c_NF_Coor_Ranges(1,2),c_NF_Coor_Ranges(2,2),c_NF_Coor_Ranges(3,1)]
eight_points(4,1:3) = [c_NF_Coor_Ranges(1,1),c_NF_Coor_Ranges(2,2),c_NF_Coor_Ranges(3,1)]
eight_points(5,1:3) = [c_NF_Coor_Ranges(1,1),c_NF_Coor_Ranges(2,1),c_NF_Coor_Ranges(3,2)]
eight_points(6,1:3) = [c_NF_Coor_Ranges(1,2),c_NF_Coor_Ranges(2,1),c_NF_Coor_Ranges(3,2)]
eight_points(7,1:3) = [c_NF_Coor_Ranges(1,2),c_NF_Coor_Ranges(2,2),c_NF_Coor_Ranges(3,2)]
eight_points(8,1:3) = [c_NF_Coor_Ranges(1,1),c_NF_Coor_Ranges(2,2),c_NF_Coor_Ranges(3,2)]

do i_Point = 1,8
    call D3_Get_Point_Domain_Number(eight_points(i_Point,1:3),Domain_Numbers(i_Point))
enddo

call Vector_Unique_Int(8,8,Domain_Numbers(1:8),Uniqued_Domain_Numbers(1:8),Num_Domain)
allocate(Domain_List(Num_Domain))
Domain_List(1:Num_Domain) = Uniqued_Domain_Numbers(1:Num_Domain)


if(Num_Domain<=0)then
    print *, "    Error-2023011201 :: Num_Domain<=0 in D3_Get_Element_List_of_Natural_Crack.f90!"       
    call Warning_Message('S',Keywords_Blank) 
endif

!***********************************  
! Get the potential element number.
!***********************************
! Variable memory allocation.
Num_Potent_Elems = 0
do i_Domain=1,Num_Domain
    c_Domain = Domain_List(i_Domain)
    Num_Potent_Elems = Num_Potent_Elems + Domain_Elements_Num(c_Domain)
enddo
allocate(Potent_Elems(Num_Potent_Elems))
Potent_Elems = 0
allocate(Unique_Potent_Elems(Num_Potent_Elems))
Unique_Potent_Elems = 0
! Cycle
c_count = 0
do i_Domain=1,Num_Domain
    c_Domain = Domain_List(i_Domain)
    do i_El =1,Domain_Elements_Num(c_Domain)
      c_Ele = Domain_Elements(c_Domain,i_El)
      c_x_max = x_max_Elements(c_Ele)
      c_x_min = x_min_Elements(c_Ele)
      c_y_max = y_max_Elements(c_Ele)
      c_y_min = y_min_Elements(c_Ele)
      c_z_max = z_max_Elements(c_Ele)
      c_z_min = z_min_Elements(c_Ele)
      ! Check whether the x-coordinate range of the detection element overlaps with the x-coordinate range
      ! of the crack
      call Tool_Yes_Two_Ranges_Overlapped_Double([c_x_min,c_x_max],&
                       [c_NF_Coor_Ranges(1,1),c_NF_Coor_Ranges(1,2)],c_Logical_Yes)   
      if(.not. c_Logical_Yes) then
          cycle
      endif
      ! Check whether the y-coordinate range of the detection element overlaps with the y-coordinate range
      ! of the crack
      call Tool_Yes_Two_Ranges_Overlapped_Double([c_y_min,c_y_max],&
                       [c_NF_Coor_Ranges(2,1),c_NF_Coor_Ranges(2,2)],c_Logical_Yes)   
      if(.not. c_Logical_Yes) then
          cycle
      endif
      ! Check whether the z-coordinate range of the detection element overlaps with the z-coordinate range
      ! of the crack
      call Tool_Yes_Two_Ranges_Overlapped_Double([c_z_min,c_z_max],&
                       [c_NF_Coor_Ranges(3,1),c_NF_Coor_Ranges(3,2)],c_Logical_Yes)   
      if(.not. c_Logical_Yes) then
          cycle
      endif      
      c_count = c_count +1
      Potent_Elems(c_count) = c_Ele
  enddo
enddo

! Delete duplicates
call Vector_Unique_Int(Num_Potent_Elems,c_count,Potent_Elems(1:Num_Potent_Elems),&
                       Unique_Potent_Elems(1:Num_Potent_Elems),c_Uniqued_count)


!**************************************
! Loop search from potential elements.
!**************************************
!///////////////////////////////////////////////////////
! Then accurately determine whether the point (x, y, z) 
! is inside the element or on its boundary.
!///////////////////////////////////////////////////////
do i_El =1,c_Uniqued_count
    c_Ele = Unique_Potent_Elems(i_El)
    c_NN  = G_NN(1:8,c_Ele)
    ! The 12 edges of the element cycle
    do i_Edge=1,12
        A(1:3) = Coor(Element_Edges(i_Edge,1,c_Ele),1:3)
        B(1:3) = Coor(Element_Edges(i_Edge,2,c_Ele),1:3)
        
        ! Using Tool_Intersection_of_AB_and_3D_Plane_Polygon
        call Tool_Intersection_of_AB_and_3D_Plane_Polygon(A,B, &
                  Each_NaCr3D_Poi_Num(c_Crack),Na_Crack3D_Coor(c_Crack,1:Each_NaCr3D_Poi_Num(c_Crack),1:3), &
                  c_Yes_Inter,c_InterSection_P) 
        ! If the edges of a element intersect with a natural fracture, the element is marked as containing a
        ! natural fracture.
        if (c_Yes_Inter) then
            tem_Out_num_Element = tem_Out_num_Element + 1
            tem_Out_Element_List(tem_Out_num_Element) = c_Ele
            goto 100
        endif
    enddo
    100 continue
end do

!********************
! Delete duplicates.
!********************
call Vector_Unique_Int(10000,tem_Out_num_Element,tem_Out_Element_List(1:10000),&
                       Out_Element_List(1:10000),Out_num_Element)

                                        
!*****************
! Clear variable.
!*****************
deallocate(Domain_List)
deallocate(Potent_Elems)
deallocate(Unique_Potent_Elems)

RETURN
END SUBROUTINE D3_Get_Element_List_of_Natural_Crack
