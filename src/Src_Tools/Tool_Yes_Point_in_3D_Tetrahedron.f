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
 
      subroutine Tool_Yes_Point_in_3D_Tetrahedron(Point,
     &                            A,B,C,D,
     &                            Yes_in,Yes_on)
C     Determine whether a point in space is inside a spatial tetrahedron. Reference: http://blog.csdn.net/bugrunner/article/details/7423727
c     Or My Notes V3-P176

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::Point(3),A(3),B(3),C(3),D(3)
      logical,intent(out):: Yes_in,Yes_on
      real(kind=FT) alpha,beta,gama,delta
      real(kind=FT) c_Distance1,c_Distance2
      real(kind=FT) c_max,c_min
      real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min
      real(kind=FT) coor_z_max,coor_z_min
      real(kind=FT) Tool_Function_2Point_Dis_3D
      
      Yes_in = .False.
      Yes_on = .False.
      
      
      
      coor_x_max = max(A(1),B(1),C(1),D(1))
      if(Point(1) > (coor_x_max+Tol_8))then
          return
      endif
      coor_x_min = min(A(1),B(1),C(1),D(1))
      if(Point(1) < (coor_x_min-Tol_8))then
          return
      endif
      coor_y_max = max(A(2),B(2),C(2),D(2))
      if(Point(2) > (coor_y_max+Tol_8))then
          return
      endif
      coor_y_min = min(A(2),B(2),C(2),D(2))
      if(Point(2) < (coor_y_min-Tol_8))then
          return
      endif
      coor_z_max = max(A(3),B(3),C(3),D(3))
      if(Point(3) > (coor_z_max+Tol_8))then
          return
      endif
      coor_z_min = min(A(3),B(3),C(3),D(3))
      if(Point(3) < (coor_z_min-Tol_8))then
          return
      endif

      if(Tool_Function_2Point_Dis_3D(A,B)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(A,C)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(A,D)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(B,C)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(B,D)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(C,D)<Tol_20) return
      

      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,B,C,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(A,B,C,D,c_Distance2)
       
      
      alpha = c_Distance1/c_Distance2
      
      if(alpha > (ONE+Tol_10)) then
          return
      endif
      if(alpha < (ZR-Tol_10)) then
          return
      endif
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,C,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(B,A,C,D,c_Distance2)   
      
      beta = c_Distance1/c_Distance2
      
      
      if(beta > (ONE+Tol_10)) then
          return
      endif
      if(beta < (ZR-Tol_10)) then
          return
      endif
      
      
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,B,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(C,A,B,D,c_Distance2)  
      
      gama = c_Distance1/c_Distance2     
      if(gama > (ONE+Tol_10)) then
          return
      endif
      if(gama < (ZR-Tol_10)) then
          return
      endif
      
      
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,B,C,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(D,A,B,C,c_Distance2)   
      
      delta = c_Distance1/c_Distance2
      if(delta > (ONE+Tol_10)) then
          return
      endif
      if(delta < (ZR-Tol_10)) then
          return
      endif
      
      
      c_max = max(alpha,beta,gama,delta)
      c_min = min(alpha,beta,gama,delta)
      
      
      if (c_max <= (ONE+Tol_10) .and. c_min>=(ZR-Tol_10)) then
          Yes_in = .True.
      endif

      if ((abs(c_max-ONE) <Tol_10) .and. (abs(c_min-ONE) <Tol_10))then
          Yes_on = .True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_in_3D_Tetrahedron                   
