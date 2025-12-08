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
 
SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Point,i_C,&
                               Check_Ball_R,Signed_Dis,                 &
                               c_Yes_Node_PER_in_FS,c_PER_Node_to_FS,   &
                               Yes_Found_Min_Signed_Dis,n_Vector)

! Calculate the distance from a point to a discrete fracture surface 
! (discrete triangles within the detection sphere).
! n_Vector is the outward normal vector of the crack surface at the foot point, 2022-05-12.
! Firstly written by Fang Shi on 2019-04-02.
! Since it is In-Plane Growth, there is no need to look for the minimum value. 2022-09-08.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D

!---------------------------     
! Variable Type Declaration
!---------------------------     
implicit none
integer,intent(in)::i_C
real(kind=FT),intent(in)::Point(3)
real(kind=FT),intent(in)::Check_Ball_R
real(kind=FT),intent(out)::Signed_Dis,c_PER_Node_to_FS(3)
real(kind=FT),intent(out)::n_Vector(3)
logical,intent(out)::c_Yes_Node_PER_in_FS,Yes_Found_Min_Signed_Dis
real(kind=FT) Tri123_Distance,Tri123_PER(3)
logical Tri123_Yes_PER_in,Tri123_Yes_PER_on
integer i_Crack_Ele
real(kind=FT) c_Min_Signed_Dis
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)  
real(kind=FT) c_dis
real(kind=FT) Tool_Function_Dis_Point_to_3D_Tri_v2
real(kind=FT) c_Crack_x_min,c_Crack_x_max
real(kind=FT) c_Crack_y_min,c_Crack_y_max
real(kind=FT) c_Crack_z_min,c_Crack_z_max
real(kind=FT) c_Check_x_min,c_Check_x_max
real(kind=FT) c_Check_y_min,c_Check_y_max
real(kind=FT) c_Check_z_min,c_Check_z_max
integer Num_CrMesh_Outlines
integer P1,P2,P3,i_Tri,Num_Tri
real(kind=FT) In_Points_Three(3,3),Points_Three_n_Vector(3),Points_Three_n_Vector_old(3)
real(kind=FT) minus_result_n_vector(3),sum_abs_minus_result_n_vector



!/////////////////////////
! Variable Initialization
!/////////////////////////
c_Yes_Node_PER_in_FS     = .False.
Yes_Found_Min_Signed_Dis = .False.
c_Min_Signed_Dis         = Con_Big_20
n_Vector(1:3)            = ZR
! Cycles of each discrete fracture element
Signed_Dis               = ZR


!////////////////////////////////////////////////////////
! Check the relationship between the point and the crack 
! coordinate range, considering Check_Ball_R
! IMPROV2022111401.
!////////////////////////////////////////////////////////
c_Crack_x_min = Crack_Coor_Range(i_C,1,1)
c_Check_x_min = c_Crack_x_min - check_Ball_R * TWO
if(Point(1) < c_Check_x_min) then
    return
endif

c_Crack_x_max = Crack_Coor_Range(i_C,1,2)
c_Check_x_max = c_Crack_x_max + Check_Ball_R * TWO
if(Point(1) > c_Check_x_max) then
    return 
endif

c_Crack_y_min = Crack_Coor_Range(i_C,2,1)
c_Check_y_min = c_Crack_y_min - check_Ball_R * TWO
if(Point(2) < c_Check_y_min) then
    return 
endif

c_Crack_y_max = Crack_Coor_Range(i_C,2,2)
c_Check_y_max = c_Crack_y_max + Check_Ball_R * TWO
if(Point(2) > c_Check_y_max) then
    return 
endif

c_Crack_z_min = Crack_Coor_Range(i_C,3,1)
c_Check_z_min = c_Crack_z_min - check_Ball_R * TWO
if(Point(3) < c_Check_z_min) then
    return 
endif

c_Crack_z_max = Crack_Coor_Range(i_C,3,2)
c_Check_z_max = c_Crack_z_max + Check_Ball_R * TWO
if(Point(3) > c_Check_z_max) then
    return 
endif

!///////////////////////////////////////////////////////////////////////
! Standard algorithm: Traverse all discrete fracture surface triangles.
!///////////////////////////////////////////////////////////////////////
if (Key_Scheme_Signed_Dis_InPlane==1) then  
    do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
        Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
        Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
        Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
        
        Point1(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
        Point2(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
        Point3(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
        
        ! Calculate the distance from the point to the discrete triangle. IMPROV2022050501.
        c_dis = Tool_Function_Dis_Point_to_3D_Tri_v2(Point,Point1,Point2,Point3)
        
        !BUGFIX2022122302.
        ! Determine the outward normal based on the element number of the discrete fracture surface element. 
        ! 2022-12-23. Adjusted here to prevent some node values from being zero.
        n_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)
              
        ! If the distance is greater than Check_Ball_R (detection sphere radius), then proceed to the next
        ! loop. IMPROV2022050501.
        if(c_dis > Check_Ball_R) then
            cycle
        endif
        
        ! Calculate the signed distance from the current node to the fracture element and determine the
        ! foot-of-perpendicular location
        call Tool_Dis_Point_to_3D_Tri(Point,Point1,Point2,Point3,Tri123_Distance,Tri123_PER,Tri123_Yes_PER_in,Tri123_Yes_PER_on)
        
        
        ! If the perpendicular foot lies within the element surface
        if (Tri123_Yes_PER_in .or. Tri123_Yes_PER_on)then
              Signed_Dis     = Tri123_Distance
              c_Yes_Node_PER_in_FS  = .True.
              c_PER_Node_to_FS(1:3)  = Tri123_PER
                
              
              !n_Vector = [0.0D0,0.0D0,1.0D0]
              
              call Vector_Normalize(3,n_Vector)
              Yes_Found_Min_Signed_Dis = .True.
              
              !Possi_Crack_Ele_Num = i_Crack_Ele  !2022-09-08.
              
              ! Exit as soon as found. 2022-09-08.
              return
        endif
        

        ! Find the minimum symbolic distance (the foot of the perpendicular does not need to be on the
        ! discrete crack surface). 2022-08-03.
        !if(abs(Tri123_Distance) < c_Min_Signed_Dis)then
              
        !endif
    enddo 
    
!/////////////////////////////////////////////////////////////////////////////
! Fast algorithm: recombine discrete triangles.
!          2024-05-02. IMPROV202450201.
! Note: This algorithm may not be suitable for 3D cracks in concave polygons.
!/////////////////////////////////////////////////////////////////////////////
elseif (Key_Scheme_Signed_Dis_InPlane==2) then  
    Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
    Num_Tri = Num_CrMesh_Outlines
    Point1 = Crack3D_Centroid(i_C,1:3)
    
    do i_Tri =1,Num_Tri
        P2 = Crack3D_Meshed_Outline(i_C)%row(i_Tri,1)
        if (i_Tri < Num_Tri)then
            P3 = Crack3D_Meshed_Outline(i_C)%row(i_Tri+1,1)
        elseif (i_Tri == Num_Tri)then
            P3 = Crack3D_Meshed_Outline(i_C)%row(1,1)
        endif
        
        Point2 = Crack3D_Meshed_Node(i_C)%row(P2,1:3) 
        Point3 = Crack3D_Meshed_Node(i_C)%row(P3,1:3) 
        
        
        ! Check the distance from the sphere center (Point) to the discrete triangle.
        c_dis = Tool_Function_Dis_Point_to_3D_Tri_v2(Point,Point1,Point2,Point3)
        
        
        ! Determine the outward normal based on the element number of the discrete fracture surface element.
        ! Adjust it here to prevent some node values from being zero.
        n_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(1,1:3)
        call Vector_Normalize(3,n_Vector)
        
        
        ! If the distance is greater than Check_Ball_R (detection sphere radius), then proceed to the next
        ! iteration.
        if(c_dis > Check_Ball_R) then
            cycle
        endif
        
        ! Check whether the triangle is abnormal (such as in the case of a concave polygon).
        In_Points_Three(1,1:3) = Point1
        In_Points_Three(2,1:3) = Point2
        In_Points_Three(3,1:3) = Point3
        call Tool_Get_Normal_Vector_of_Points_3D(In_Points_Three,3,Points_Three_n_Vector)        
        call Vector_Normalize(3,Points_Three_n_Vector)
        if (i_Tri>=2) then
            minus_result_n_vector = Points_Three_n_Vector - Points_Three_n_Vector_old
            sum_abs_minus_result_n_vector = sum(abs(minus_result_n_vector))
            ! If the normal vector of the triangle is flipped compared to the previous triangle, it indicates an
            ! error.
            if(sum_abs_minus_result_n_vector>=Tol_3) then
                print *,'    ERROR-2024050201 :: Illegal triangle detected!'
                print *,'                        In D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth.f90!'
                print *,'                        Points_Three_n_Vector_old:',Points_Three_n_Vector_old
                print *,'                        Points_Three_n_Vector:',Points_Three_n_Vector
                print *,'                        Set *Key_Scheme_Signed_Dis_InPlane to 1 and try again!'
                call Warning_Message('S',Keywords_Blank)  
            endif
        endif
        
        ! Calculate the signed distance from the current node to the fracture element and determine the
        ! foot-of-perpendicular location
        call Tool_Dis_Point_to_3D_Tri(Point,Point1,Point2,Point3,Tri123_Distance,Tri123_PER,Tri123_Yes_PER_in,Tri123_Yes_PER_on)
        
        ! If the perpendicular foot lies within the element surface
        if (Tri123_Yes_PER_in .or. Tri123_Yes_PER_on)then
            Signed_Dis     = Tri123_Distance
            c_Yes_Node_PER_in_FS  = .True.
            c_PER_Node_to_FS(1:3)  = Tri123_PER

            
            Yes_Found_Min_Signed_Dis = .True.
              
            ! Exit once found. 2022-09-08.
            return
        endif
        Points_Three_n_Vector_old = Points_Three_n_Vector
    enddo
!/////////////////////////////////////////////////////////////////////////////
! Fast algorithm: Recombine discrete triangles. At the same time, remove 
!                 duplicate points on the same line segment.
!                 2024-05-03. IMPROV202450203.
! Note: This algorithm may not be suitable for 3D cracks in concave polygons.
!/////////////////////////////////////////////////////////////////////////////
elseif (Key_Scheme_Signed_Dis_InPlane==3) then              
    Num_CrMesh_Outlines = simplified_crack_outline_num(i_C)
    Num_Tri = Num_CrMesh_Outlines
      
    Point1 = Crack3D_Centroid(i_C,1:3)
    
    do i_Tri =1,Num_Tri
        !P2 = Crack3D_Meshed_Outline(i_C)%row(i_Tri,1)
        Point2 = simplified_crack_outline_points(i_C,i_Tri,1:3)
        if (i_Tri < Num_Tri)then
            !P3 = Crack3D_Meshed_Outline(i_C)%row(i_Tri+1,1)
            Point3 = simplified_crack_outline_points(i_C,i_Tri+1,1:3)
        elseif (i_Tri == Num_Tri)then
            !P3 = Crack3D_Meshed_Outline(i_C)%row(1,1)
            Point3 = simplified_crack_outline_points(i_C,1,1:3)
        endif
        
        !Point2 = Crack3D_Meshed_Node(i_C)%row(P2,1:3) 
        !Point3 = Crack3D_Meshed_Node(i_C)%row(P3,1:3) 
        
        
        ! Check the distance from the sphere center (Point) to the discrete triangle.
        c_dis = Tool_Function_Dis_Point_to_3D_Tri_v2(Point,Point1,Point2,Point3)
        
        
        ! Determine the outward normal based on the element number of the discrete 
        ! fracture surface element. Adjust it here to prevent some node values from being zero.
        n_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(1,1:3)
        call Vector_Normalize(3,n_Vector)
        
        
        ! If the distance is greater than Check_Ball_R (detection sphere radius), then proceed to the next
        ! iteration.
        if(c_dis > Check_Ball_R) then
            cycle
        endif
        
        ! Check whether the triangle is abnormal (such as in the case of concave polygons).
        In_Points_Three(1,1:3) = Point1
        In_Points_Three(2,1:3) = Point2
        In_Points_Three(3,1:3) = Point3
        call Tool_Get_Normal_Vector_of_Points_3D(In_Points_Three,3,Points_Three_n_Vector)        
        call Vector_Normalize(3,Points_Three_n_Vector)
        if (i_Tri>=2) then
            minus_result_n_vector = Points_Three_n_Vector - Points_Three_n_Vector_old
            sum_abs_minus_result_n_vector = sum(abs(minus_result_n_vector))
            ! If the normal vector of the triangle is flipped compared to the previous triangle, it indicates an
            ! error.
            if(sum_abs_minus_result_n_vector>=Tol_3) then
                print *,'    ERROR-2024050201 :: Illegal triangle detected!'
                print *,'                        In D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth.f90!'
                print *,'                        Points_Three_n_Vector_old:',Points_Three_n_Vector_old
                print *,'                        Points_Three_n_Vector:',Points_Three_n_Vector
                print *,'                        Set *Key_Scheme_Signed_Dis_InPlane to 1 and try again!'
                call Warning_Message('S',Keywords_Blank)  
            endif
        endif
        
        ! Calculate the signed distance from the current node to the fracture element
        ! and determine the foot-of-perpendicular location
        call Tool_Dis_Point_to_3D_Tri(Point,Point1,Point2,Point3,Tri123_Distance,Tri123_PER,Tri123_Yes_PER_in,Tri123_Yes_PER_on)
        
        ! If the perpendicular foot lies within the element surface
        if (Tri123_Yes_PER_in .or. Tri123_Yes_PER_on)then
            Signed_Dis     = Tri123_Distance
            c_Yes_Node_PER_in_FS  = .True.
            c_PER_Node_to_FS(1:3)  = Tri123_PER

            
            Yes_Found_Min_Signed_Dis = .True.
              
            ! Exit once found. 2022-09-08.
            return
        endif
        Points_Three_n_Vector_old = Points_Three_n_Vector
    enddo    
endif

RETURN
END SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth
