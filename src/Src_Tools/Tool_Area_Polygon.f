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
 
      subroutine Tool_Area_Polygon(x,y,nb,area)
       
      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2000-07-04  Time: 12:24:06
      use Global_Float_Type
      IMPLICIT NONE
      integer,intent(in)::nb
      real(kind=FT) ,intent(in)::x(nb),y(nb)
      real(kind=FT) ,intent(out)::area
      INTEGER i, n, nm1
      real(kind=FT) a


      n = nb
      
      
      IF((abs(x(1)-x(n)) <=Tol_15) .AND. (abs(y(1)-y(n))<=Tol_15))then 
          n = n - 1
      endif
      
      SELECT CASE (n)
      CASE (:2)
          area = ZR
      CASE (3)
          area=HLF*((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)))
      CASE DEFAULT
          nm1 = n - 1
          a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
          DO  i = 2, nm1
              a = a + x(i)*(y(i+1) - y(i-1))
          END DO
          area = HLF*a
      END SELECT

      RETURN
      END subroutine Tool_Area_Polygon
