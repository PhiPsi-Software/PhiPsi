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
 
      subroutine Tool_Area_Tri_3D(A,B,C,area)
       
      ! Given the coordinates of the vertices, calculate the area of a triangle in three-dimensional
      ! space.
      ! See my notes V3-P139 for details.
      
      use Global_Float_Type
      IMPLICIT NONE
      real(kind=FT),intent(in)::A(3),B(3),C(3)
      real(kind=FT),intent(out)::area
      real(kind=FT) Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
      real(kind=FT) x1,x2,x3,y1,y2,y3
      
      Ax = A(1)
      Ay = A(2)
      Az = A(3)
      Bx = B(1)
      By = B(2)
      Bz = B(3)
      Cx = C(1)
      Cy = C(2)
      Cz = C(3)
      
      x1 = Bx - Ax;x2 = By - Ay;x3 = Bz - Az
      
      y1 = Cx - Ax;y2 = Cy - Ay;y3 = Cz - Az
      
      area = HLF*sqrt((x2*y3-x3*y2)**2 + 
     &                (x3*y1-x1*y3)**2 + 
     &                (x1*y2-x2*y1)**2)
      
      RETURN
      END subroutine Tool_Area_Tri_3D
