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
 
      SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh(Point,i_C,
     &                               Check_Ball_R,
     &                               Signed_Dis,Signed_Dis_v2,
     &                               c_Yes_Node_PER_in_FS,
     &                               c_PER_Node_to_FS,
     &                               Yes_Found_Min_Signed_Dis,
     &                               n_Vector)
     
c     Calculate the distance from a point to a discrete fracture surface (discrete triangles within the detection sphere).
c     n_Vector is the outward normal vector of the crack surface at the foot point, 2022-05-12.
c     c_Yes_Node_PER_in_FS: Whether the hanging wall foot is within the fracture surface.
c     c_PER_Node_to_FS: perpendicular foot.
c     Firstly written by Fang Shi on 2019-04-02.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Crack_Common
      use Global_Crack_3D
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::i_C
      real(kind=FT),intent(in)::Point(3)
      real(kind=FT),intent(in)::Check_Ball_R
      real(kind=FT),intent(out)::Signed_Dis,c_PER_Node_to_FS(3)
      real(kind=FT),intent(out)::Signed_Dis_v2
      real(kind=FT),intent(out)::n_Vector(3)
      logical,intent(out)::Yes_Found_Min_Signed_Dis
      logical,intent(out)::c_Yes_Node_PER_in_FS
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
      
c     ------------------------
C     Variable Initialization
c     ------------------------
      c_Yes_Node_PER_in_FS     = .False.
      c_Min_Signed_Dis         = Con_Big_20
      Yes_Found_Min_Signed_Dis = .False.
      n_Vector(1:3)            = ZR
      ! Cycles of each discrete fracture element
      Signed_Dis               = ZR
      Signed_Dis_v2            = ZR
      
c     --------------------------------------------------------
C      Check the relationship between the point and the crack 
c      coordinate range, considering Check_Ball_R
c      IMPROV2022111401.
c     --------------------------------------------------------
      c_Crack_x_min = Crack_Coor_Range(i_C,1,1)
      c_Check_x_min = c_Crack_x_min - check_Ball_R * TWO
      if(Point(1) < c_Check_x_min) return
      
      c_Crack_x_max = Crack_Coor_Range(i_C,1,2)
      c_Check_x_max = c_Crack_x_max + Check_Ball_R * TWO
      if(Point(1) > c_Check_x_max) return
      
      c_Crack_y_min = Crack_Coor_Range(i_C,2,1)
      c_Check_y_min = c_Crack_y_min - check_Ball_R * TWO
      if(Point(2) < c_Check_y_min) return
      
      c_Crack_y_max = Crack_Coor_Range(i_C,2,2)
      c_Check_y_max = c_Crack_y_max + Check_Ball_R * TWO
      if(Point(2) > c_Check_y_max) return
      
      c_Crack_z_min = Crack_Coor_Range(i_C,3,1)
      c_Check_z_min = c_Crack_z_min - check_Ball_R * TWO
      if(Point(3) < c_Check_z_min) return
      
      c_Crack_z_max = Crack_Coor_Range(i_C,3,2)
      c_Check_z_max = c_Crack_z_max + Check_Ball_R * TWO
      if(Point(3) > c_Check_z_max) return
      
c     ------------
C     Loop search
c     ------------
      do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
          Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
          Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
          Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
          
          Point1(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
          Point2(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
          Point3(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
          
          ! Calculate the distance from the point to the discrete triangle. IMPROV2022050501.
          c_dis = Tool_Function_Dis_Point_to_3D_Tri_v2(Point,
     &                            Point1,Point2,Point3)
          
          !BUGFIX2022122303.
          ! Determine the outward normal based on the element number of the discrete fracture surface element.
          ! 2022-12-23. Adjusted here to prevent some node values from being 0.
          n_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector
     &                    (i_C)%row(i_Crack_Ele,1:3)

          ! If the distance is greater than Check_Ball_R (detection sphere radius), then proceed to the next
          ! loop. IMPROV2022050501.
          if(c_dis > Check_Ball_R) then
              cycle
          endif
          
          ! Calculate the signed distance from the current node to the fracture element and determine the
          ! foot-of-perpendicular location
          call Tool_Dis_Point_to_3D_Tri
     &           (Point,Point1,Point2,Point3,
     &            Tri123_Distance,Tri123_PER,
     &            Tri123_Yes_PER_in,Tri123_Yes_PER_on)
          
          ! If the perpendicular foot lies within the element surface
          if (Tri123_Yes_PER_in .or. Tri123_Yes_PER_on)then
            ! Find the minimum symbol distance
            if(abs(Tri123_Distance) < c_Min_Signed_Dis)then
              Signed_Dis     = Tri123_Distance
              c_Min_Signed_Dis = abs(Tri123_Distance) 
              c_Yes_Node_PER_in_FS  = .True.
              c_PER_Node_to_FS(1:3)  = Tri123_PER
              Yes_Found_Min_Signed_Dis = .True.
              
              call Vector_Normalize(3,n_Vector)
              
            endif
          endif
          
          ! Find the minimum symbolic distance (the foot of the perpendicular does not need to be on the
          ! discrete crack surface). 2022-08-03.
          if(abs(Tri123_Distance) < c_Min_Signed_Dis)then
              Signed_Dis_v2     = Tri123_Distance
          endif
          
          
      enddo 
      
      
      RETURN
      END SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh
