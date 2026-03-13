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
 
      subroutine Tool_Yes_Point_in_Inclusions(x,y,
     &                        Yes_Gauss_in_Incl,c_Incl_Num)
C     Determine whether a point is inside a certain inclusion.

      use Global_Float_Type
      use Global_Inclusion
      
      implicit none
      integer c_Incl_Num
      real(kind=FT) x,y
      logical Yes_Gauss_in_Incl
      
      integer i_Incl,n_Vertex
      real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r,c_Dis
      real(kind=FT) Tool_Function_2Point_Dis
      
      Yes_Gauss_in_Incl = .False.
      c_Incl_Num = 0
      
      if(num_Circ_Incl/=0)then
          do i_Incl = 1,num_Circ_Incl
              c_Incl_x = Circ_Inclu_Coor(i_Incl,1)
              c_Incl_y = Circ_Inclu_Coor(i_Incl,2)
              c_Incl_r = Circ_Inclu_Coor(i_Incl,3)
              c_Dis    = Tool_Function_2Point_Dis([x,y],
     &                              [c_Incl_x,c_Incl_y])
              if(c_Dis < c_Incl_r)then
                  Yes_Gauss_in_Incl = .True.
                  c_Incl_Num = i_Incl
                  exit
              endif
          enddo
      endif
      
      if(num_Poly_Incl/=0)then
          do i_Incl = 1,num_Poly_Incl
              n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
              call Tool_Yes_In_Poly(x,y,
     &                        Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),
     &                        Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1),
     &                        n_Vertex+1,
     &                        Yes_Gauss_in_Incl)
     
              if(Yes_Gauss_in_Incl)then
                  c_Incl_Num = i_Incl
                  exit
              endif
          enddo
      endif
  
      return 
      end SUBROUTINE Tool_Yes_Point_in_Inclusions                          
