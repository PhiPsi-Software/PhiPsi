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
 
      subroutine Tool_Intersections_of_Two_Triangles_3D(Tri_1,Tri_2,
     &                            Logical_Inter,num_Inters,
     &                            Inter_Points,Logical_Parallel)   
     
c     Used to determine whether two spatial triangles intersect and to calculate the intersection point (if they intersect). NEWFTU2022042901.
c     Added on 2022-04-29. 
c     Ref: \theory_documents\031 Do Two Triangles Intersect and Intersection Point Calculation-2022-04-29.pdf or https://stackoverflow.com/questions/7113344/find-whether-two-triangles-intersect-or-not
c      There are four situations in total. In the first situation, two sides of one triangle intersect with another triangle. In the second situation, the two triangles intersect each other.
c     Case 3: If they are coplanar, there may be multiple intersection points. See changelog.src/2022-04-29-04.png, Ref: \theory_documents\031.1 Whether two triangles intersect and calculation of intersection points (including coplanar or other situations)-2022-04-29.pdf
c     Case 4: One edge of a triangle lies within the plane of another triangle. See changelog.src/2022-04-29-05.png. Ref: \theory_documents\031.1 Determining whether two triangles intersect and calculating intersection points (including coplanar or other cases) - 2022-04-29.pdf

c     Test data: https://math.stackexchange.com/questions/1220102/how-do-i-find-the-intersection-of-two-3d-triangles
C     Tri_1(1,1:3) = [6.0,  8.0,  3.0]
C     Tri_1(2,1:3) = [6.0,  8.0, -2.0]
C     Tri_1(3,1:3) = [6.0, -4.0, -2.0]
C     Tri_2(1,1:3) = [0.0,  5.0,  0.0]
C     Tri_2(2,1:3) = [0.0,  0.0,  0.0]
C     Tri_2(3,1:3) = [8.0,  0.0,  0.0]
c     Two points of intersection
C     Inter_Points(1,1:3):   6.0       0.80        0.00
C     Inter_Points(2,1:3):   6.0       1.25        0.00  
c     The corresponding Matlab plotting verification code: changelog.src/2022-04-29-02.m
      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      use Global_Common
      implicit none
      real(kind=FT),intent(in)::Tri_1(3,3),Tri_2(3,3)
      integer,intent(out)::num_Inters
      logical,intent(out)::Logical_Inter
      real(kind=FT),intent(out)::Inter_Points(10,3)
      integer Tri_Edges(3,3)
      integer i_Edge
      real(kind=FT) P1(3),P2(3),c_InterSection_P(3)
      logical c_Yes_Inter,Logical_Parallel
      real(kind=FT) Dis_two_P,Tool_Function_2Point_Dis_3D
      real(kind=FT) Np_Tri_1(3),Np_Tri_2(3)
      logical Logical_Yes_x,Logical_Yes_y,Logical_Yes_z
      real(kind=FT) max_x_Tri_1,max_y_Tri_1,max_z_Tri_1
      real(kind=FT) min_x_Tri_1,min_y_Tri_1,min_z_Tri_1
      real(kind=FT) max_x_Tri_2,max_y_Tri_2,max_z_Tri_2
      real(kind=FT) min_x_Tri_2,min_y_Tri_2,min_z_Tri_2
      integer c_Location
      logical c_Yes_Exist      
 
      num_Inters    = 0
      Logical_Inter =  .False.
      Logical_Parallel = .False.
      Inter_Points(1:10,1:3) = -TEN_15
      
      Tri_Edges(1,1:2) = [1,2] 
      Tri_Edges(2,1:2) = [2,3] 
      Tri_Edges(3,1:2) = [3,1] 
      
      max_x_Tri_1 = maxval(Tri_1(1:3,1))
      max_y_Tri_1 = maxval(Tri_1(1:3,2))
      max_z_Tri_1 = maxval(Tri_1(1:3,3))
      min_x_Tri_1 = minval(Tri_1(1:3,1))
      min_y_Tri_1 = minval(Tri_1(1:3,2))
      min_z_Tri_1 = minval(Tri_1(1:3,3))      
      max_x_Tri_2 = maxval(Tri_2(1:3,1))
      max_y_Tri_2 = maxval(Tri_2(1:3,2))
      max_z_Tri_2 = maxval(Tri_2(1:3,3))
      min_x_Tri_2 = minval(Tri_2(1:3,1))
      min_y_Tri_2 = minval(Tri_2(1:3,2))
      min_z_Tri_2 = minval(Tri_2(1:3,3))         
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_x_Tri_1,max_x_Tri_1],
     &                            [min_x_Tri_2,max_x_Tri_2],
     &                            Logical_Yes_x) 
      if(Logical_Yes_x .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_y_Tri_1,max_y_Tri_1],
     &                            [min_y_Tri_2,max_y_Tri_2],
     &                            Logical_Yes_y) 
      if(Logical_Yes_y .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_z_Tri_1,max_z_Tri_1],
     &                            [min_z_Tri_2,max_z_Tri_2],
     &                            Logical_Yes_z) 
      if(Logical_Yes_z .eqv. .False.)then
          return
      endif      
      
      call Tool_Normal_Unit_vector_of_3D_Tri(Tri_1(1,1:3),
     &                         Tri_1(2,1:3),Tri_1(3,1:3),Np_Tri_1(1:3))
      call Tool_Normal_Unit_vector_of_3D_Tri(Tri_2(1,1:3),
     &                         Tri_2(2,1:3),Tri_2(3,1:3),Np_Tri_2(1:3))     
      call Vectors_Equal_Is_Dou_with_Tol(Np_Tri_1,Np_Tri_2,3,
     &                                   Logical_Parallel)   
      if(Logical_Parallel) then
          return
      endif

      do i_Edge = 1,3
          P1(1:3) = Tri_1(Tri_Edges(i_Edge,1),1:3)
          P2(1:3) = Tri_1(Tri_Edges(i_Edge,2),1:3)
          call Tool_Intersection_of_AB_and_Triangle_3D(P1(1:3),P2(1:3),
     &                 Tri_2(1,1:3),Tri_2(2,1:3),Tri_2(3,1:3),
     &                 c_Yes_Inter,c_InterSection_P) 
         if (c_Yes_Inter) then
             call Vector_belongs_Matrix_Is_Dou(num_Inters,3,
     &                Inter_Points(1:num_Inters,1:3),c_InterSection_P,
     &                c_Location,c_Yes_Exist)  
             if(c_Yes_Exist .eqv. .False.) then
                 num_Inters = num_Inters +1
                 Inter_Points(num_Inters,1:3) = c_InterSection_P
                 Logical_Inter = .True.
             endif
         endif
      enddo
      
      do i_Edge = 1,3
          P1(1:3) = Tri_2(Tri_Edges(i_Edge,1),1:3)
          P2(1:3) = Tri_2(Tri_Edges(i_Edge,2),1:3)
          call Tool_Intersection_of_AB_and_Triangle_3D(P1(1:3),P2(1:3),
     &                 Tri_1(1,1:3),Tri_1(2,1:3),Tri_1(3,1:3),
     &                 c_Yes_Inter,c_InterSection_P) 
         if (c_Yes_Inter) then
             call Vector_belongs_Matrix_Is_Dou(num_Inters,3,
     &                Inter_Points(1:num_Inters,1:3),c_InterSection_P,
     &                c_Location,c_Yes_Exist)  
             if(c_Yes_Exist .eqv. .False.) then     
                 num_Inters = num_Inters +1
                 Inter_Points(num_Inters,1:3) = c_InterSection_P
                 Logical_Inter = .True.
             endif
         endif
      enddo
          
      if(num_Inters==2) then
          Dis_two_P = Tool_Function_2Point_Dis_3D(Inter_Points(1,1:3),
     &                                            Inter_Points(2,1:3))
          if(Dis_two_P<=Tol_11)then
              num_Inters = 1
          endif
      endif
      
      if(num_Inters>=3) then
          print *, '    Warning :: More than 2 intersections found!'
          print *, '         in Tool_Intersections_of_Two_Triangles_3D!'
      endif      
      return 
      end SUBROUTINE Tool_Intersections_of_Two_Triangles_3D          
