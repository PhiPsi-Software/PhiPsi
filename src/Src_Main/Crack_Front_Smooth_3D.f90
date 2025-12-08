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
 
subroutine Crack_Front_Smooth_3D(ifra,iter)
! Smooth the front edge of the crack. For 3D cracks only.
! Called when Key_Smooth_Front=1.
! Called by PhiPsi3D_Static.f, etc.
! 2021-09-08.

!***************
! Public Module
!***************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_HF
use Global_Common
use Global_Material
use Global_Filename

use Global_Cal_Ele_Num_by_Coors_3D  
  
  
!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::ifra,iter
integer i_C
integer Num_CrMesh_Outlines,i_Out_Node,c_Mesh_Node
integer in_Ele_Num
real(kind=FT) c_x,c_y,c_z
integer,ALLOCATABLE:: Vertexes_in_Model(:)
real(kind=FT),ALLOCATABLE:: Points_to_Smooth(:,:)
real(kind=FT),ALLOCATABLE:: Points_Smoothed(:,:)
integer Num_Point_Smooth,c_count,in_Elem_num
integer i_Crack_Ele,Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) Tool_Function_2Point_Dis_3D
real(kind=FT) Vector_V1(3),Vector_V2(3),Vector_Normal(3)
integer i_Crack_Node,c_count_Ele,num_Cr_Nodes
real(kind=FT) c_All_Nor_Vector(3),centroid_x,centroid_y,centroid_z
integer Ele_Num_Cache
integer Smoothed_Cracks(num_Crack),n_Smoothed_Cracks,c_C
real(kind=FT) c_Kesi,c_Yita,c_Zeta 
!******************************************************************************************
! If there is no crack tip (crack tip enhancement node), there is no need to perform crack
! propagation calculations.
! Added on 2021-01-01.
!******************************************************************************************
if(Key_Smooth_Front==6) return
if(Key_Smooth_Front<=0) return

!*********************
! Interactive display
!*********************
print *,'    Smoothing crack front vertex...'   
 
!**********************
! Cycle between cracks
!**********************
Ele_Num_Cache = 1

! IMPROV2022121702. OpenMP Parallelization.
n_Smoothed_Cracks = 0
do i_C =1,num_Crack    
    if(Cracks_Allow_Propa(i_C) /= 0)then  
        n_Smoothed_Cracks = n_Smoothed_Cracks +1
        Smoothed_Cracks(n_Smoothed_Cracks) = i_C
    endif
enddo
  
do i_C =1,n_Smoothed_Cracks   
    c_C = Smoothed_Cracks(i_C) 
    Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(c_C)
    
    ! Clear temporary variables
    if (allocated(Vertexes_in_Model)) DEALLOCATE(Vertexes_in_Model)
    if (allocated(Points_to_Smooth)) deallocate(Points_to_Smooth) 
    if (allocated(Points_Smoothed)) deallocate(Points_Smoothed) 
    
    !---------   STEP   1 ----------
    ! Loop over the front edge points of the crack to obtain the points within the model,
    ! Thus, dealing with the situation of edge cracks.
    !-------------------------------
    allocate(Vertexes_in_Model(Num_CrMesh_Outlines))
    Vertexes_in_Model(1:Num_CrMesh_Outlines) = 0        
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(c_C)%row(i_Out_Node,1)
      c_x  = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,1) 
      c_y  = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,2) 
      c_z  = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,3)
      ! element number
      call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,Ele_Num_Cache,in_Ele_Num)
      if(in_Ele_Num>=1) then
          Vertexes_in_Model(i_Out_Node) =1
      endif
    enddo
    ! Total number of points that need smoothing
    Num_Point_Smooth = sum(Vertexes_in_Model(1:Num_CrMesh_Outlines))
    
    ! If the leading edge of the crack is located outside the model, it is not suitable for the
    ! smoothing algorithm 6, 2022-12-17.
    if(Num_Point_Smooth < Num_CrMesh_Outlines) then
      if(Key_Smooth_Front==6)then
        print *,'    ERROR-2022121702 :: Vertex outside the model!'
        print *,'          Key_Smooth_Front=6 is not available!'
        print *,'          In Crack_Front_Smooth_3D.f!'    
        goto 100 
      endif
    endif
     
    !---------   STEP   2 ------------------------------
    ! Vector of points in the mesh that need smoothing.
    !---------------------------------------------------
    c_count = 0
    allocate(Points_to_Smooth(Num_Point_Smooth,3)) 
    Points_to_Smooth(1:Num_Point_Smooth,1:3) = ZR
    do i_Out_Node = 1,Num_CrMesh_Outlines
      if(Vertexes_in_Model(i_Out_Node)==1)then
          c_Mesh_Node =Crack3D_Meshed_Outline(c_C)%row(i_Out_Node,1)
          c_count = c_count +1
          Points_to_Smooth(c_count,1:3) = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,1:3) 
      endif
    enddo   
    
    !---------   STEP   3 ----------
    ! Smooth finish.
    !-------------------------------
    allocate(Points_Smoothed(Num_Point_Smooth,3)) 
    !//////////////
    ! Fit to line.
    !//////////////
    if(Key_Smooth_Front==1)then
        call Tool_Fit_3D_Points_to_Straight_Line(Num_Point_Smooth, &
                    Points_to_Smooth(1:Num_Point_Smooth,1:3), &
                    Points_Smoothed(1:Num_Point_Smooth,1:3))

    !///////////////////////////////
    !     Fit to circle, 2021-10-30
    !///////////////////////////////
    elseif(Key_Smooth_Front==2)then
        call Tool_Fit_3D_Points_to_Circle(Num_Point_Smooth, &
                    Points_to_Smooth(1:Num_Point_Smooth,1:3), &
                    Points_Smoothed(1:Num_Point_Smooth,1:3))
    
    !//////////////////////////////////////////////////
    !Fit using \Python_Tools\Smooth_3D_Points_CSAPS.py
    ! Cubic spline fitting, calling a Python script
    !date: 2021-11-01
    !//////////////////////////////////////////////////
    elseif(Key_Smooth_Front==3)then
        call Tool_Fit_3D_Points_using_Python( &
                    Key_Smooth_Front_Twice, &
                    Num_Point_Smooth, &
                    Points_to_Smooth(1:Num_Point_Smooth,1:3), &
                    Points_Smoothed(1:Num_Point_Smooth,1:3))   
 
    !//////////////////////////////////////////////////
    !Fit using \Python_Tools\Smooth_3D_Points_Scipy.py
    ! B-spline curve fitting, calling a Python script
    !date: 2021-11-02
    !//////////////////////////////////////////////////
    elseif(Key_Smooth_Front==4)then
        call Tool_Fit_3D_Points_using_Python( &
                    Key_Smooth_Front_Twice, &
                    Num_Point_Smooth, &
                    Points_to_Smooth(1:Num_Point_Smooth,1:3), &
                    Points_Smoothed(1:Num_Point_Smooth,1:3))  
 

    !//////////////////////////////////////
    !    Fit to Ellipse, NEWFTU2022051701.
    !//////////////////////////////////////
    elseif(Key_Smooth_Front==5)then
        call Tool_Fit_3D_Points_to_Ellipse(Num_Point_Smooth, &
                    Points_to_Smooth(1:Num_Point_Smooth,1:3), &
                    Points_Smoothed(1:Num_Point_Smooth,1:3))

    !//////////////////////////////////////////////////////////////////////////////////////////
    ! Taubin algorithm. Key_Smooth_Front=6. IMPROV2022121001.
    !  Ref: \theory_documents\042 Mesh Smoothing-2022-12-10.pdf，P15-P25.
    ! Only applicable to the leading edge of closed cracks.
    ! It should be ensured that the apex of the front edge of all crack surfaces is within the
    ! model.
    !//////////////////////////////////////////////////////////////////////////////////////////
    elseif(Key_Smooth_Front==6)then
        ! Executed in Check_Crack_Grows_3D.f90.
    endif
    
    
    
    
    
    !---------   STEP   4 ----------
    ! Update the front edge coordinates of the crack.
    !-------------------------------
    c_count = 0
    Ele_Num_Cache = 1
    do i_Out_Node = 1,Num_CrMesh_Outlines
      if(Vertexes_in_Model(i_Out_Node)==1)then
          c_Mesh_Node =Crack3D_Meshed_Outline(c_C)%row(i_Out_Node,1)
          c_count = c_count +1
          Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,1:3)  = Points_Smoothed(c_count,1:3)
          ! Update element number
          call Cal_Ele_Num_by_Coors_3D(Points_Smoothed(c_count,1),Points_Smoothed(c_count,2),Points_Smoothed(c_count,3), &
                         Ele_Num_Cache,in_Elem_num)  
          if(in_Elem_num<=0)then
              print *, '    Error-2022121701:: in_Elem_num<=0!'
              print *, '          In Crack_Front_Smooth_3D.f!'
              print *, '          Points_Smoothed:',Points_Smoothed(c_count,1:3)
              call Warning_Message('S',Keywords_Blank)
          else
            Cr3D_Meshed_Node_in_Ele_Num(c_C)%row(c_Mesh_Node)=in_Elem_num 
            !2022-12-18.
            call Cal_KesiYita_by_Coor_3D(Points_Smoothed(c_count,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
            Cr3D_Meshed_Node_in_Ele_Local(c_C)%row(c_Mesh_Node,1:3) = [c_Kesi,c_Yita,c_Zeta] 
          endif
      endif
    enddo   
    
    !---------   STEP   5 ----------
    ! Other variable updates.
    !-------------------------------
    do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(c_C)
      Crack_Node1 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,1)
      Crack_Node2 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,2)
      Crack_Node3 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,3)
      Point1(1:3) = Crack3D_Meshed_Node(c_C)%row(Crack_Node1,1:3)
      Point2(1:3) = Crack3D_Meshed_Node(c_C)%row(Crack_Node2,1:3)
      Point3(1:3) = Crack3D_Meshed_Node(c_C)%row(Crack_Node3,1:3)
      ! Perimeter of a discrete element (spatial triangle)
      Crack3D_Meshed_Ele_Attri(c_C)%row(i_Crack_Ele,1)  =  &
           Tool_Function_2Point_Dis_3D(Point1,Point2) +  &
           Tool_Function_2Point_Dis_3D(Point2,Point3)+ &
           Tool_Function_2Point_Dis_3D(Point1,Point3)
      ! Area of discrete elements (spatial triangles)
      call Tool_Area_Tri_3D(Point1,Point2,Point3,Crack3D_Meshed_Ele_Attri(c_C)%row(i_Crack_Ele,2))
      ! Outer normal vector Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
      Vector_V1 = Point2- Point1
      Vector_V2 = Point3- Point1
      call Vector_Normalize(3,Vector_V1)  
      call Vector_Normalize(3,Vector_V2)  
      call Vector_Cross_Product_3(Vector_V1,Vector_V2,Vector_Normal) 
      call Vector_Normalize(3,Vector_Normal)  
      Crack3D_Meshed_Ele_Nor_Vector(c_C)%row(i_Crack_Ele,1:3)=Vector_Normal
    enddo
  
    ! Outer normal vectors of 3D crack nodes after discretization
    ! Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D,Max_N_Node_3D,3)  
    ! Crack3D_Meshed_Node_Nor_Vector(1:Max_Num_Cr_3D,1:Max_N_Node_3D,1:3)   = ZR
    Crack3D_Meshed_Node_Nor_Vector(c_C)%row(1:Max_N_Node_3D(c_C),1:3) =ZR
    do i_Crack_Node =1,Crack3D_Meshed_Node_num(c_C) 
      c_count_Ele = 0
      c_All_Nor_Vector = [ZR,ZR,ZR]
      do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(c_C)
          Crack_Node1 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,1)
          Crack_Node2 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,2)
          Crack_Node3 = Crack3D_Meshed_Ele(c_C)%row(i_Crack_Ele,3)
          if(i_Crack_Node == Crack_Node1 .or.i_Crack_Node == Crack_Node2 .or. i_Crack_Node == Crack_Node3 )then
             c_All_Nor_Vector  = c_All_Nor_Vector  +Crack3D_Meshed_Ele_Nor_Vector(c_C)%row(i_Crack_Ele,1:3)
              c_count_Ele = c_count_Ele +1
          endif
      enddo
      Crack3D_Meshed_Node_Nor_Vector(c_C)%row(i_Crack_Node,1:3)  = c_All_Nor_Vector/c_count_Ele
    enddo
  
    !----------------------------------------------------------
    ! Centroid of the 3D fracture surface after discretization
    !----------------------------------------------------------
    num_Cr_Nodes = Crack3D_Meshed_Node_num(c_C)
    centroid_x=sum(Crack3D_Meshed_Node(c_C)%row(1:num_Cr_Nodes,1))/dble(num_Cr_Nodes)
    centroid_y=sum(Crack3D_Meshed_Node(c_C)%row(1:num_Cr_Nodes,2))/dble(num_Cr_Nodes)
    centroid_z=sum(Crack3D_Meshed_Node(c_C)%row(1:num_Cr_Nodes,3))/dble(num_Cr_Nodes)
    Crack3D_Centroid(c_C,1) = centroid_x
    Crack3D_Centroid(c_C,2) = centroid_y
    Crack3D_Centroid(c_C,3) = centroid_z
100   continue  
end do
  
RETURN
END SUBROUTINE Crack_Front_Smooth_3D