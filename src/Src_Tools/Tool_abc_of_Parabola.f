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
 
      subroutine Tool_abc_of_Parabola(Point_1,Point_2,Point_3,a,b,c)
C     Calculate the coefficients a, b, c of the parabola
C     The input variables are the x and y values of three points arranged sequentially on a parabolic curve.
C     Calculate the equation of a parabola: ax**2 + bx + c = 0
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: Point_1(2),Point_2(2),Point_3(2)
      real(kind=FT),intent(out):: a,b,c
      
      real(kind=FT) D(3),K(3,3),inv_K(3,3),F(3),x1,x2,x3,y1,y2,y3
      
      x1 = Point_1(1); y1 = Point_1(2)
      x2 = Point_2(1); y2 = Point_2(2)
      x3 = Point_3(1); y3 = Point_3(2)
      
      K(1,1:3) = [x1**2, x1, ONE]
      K(2,1:3) = [x2**2, x2, ONE]
      K(3,1:3) = [x3**2, x3, ONE]

      
      F(1:3)   = [y1,y2,y3]
      
      call Matrix_Inverse(K,inv_K,3) 
      
      D = MATMUL(inv_K,F) 
      
      a = D(1)
      b = D(2)
      c = D(3)

      
      return 
      end SUBROUTINE Tool_abc_of_Parabola                          
