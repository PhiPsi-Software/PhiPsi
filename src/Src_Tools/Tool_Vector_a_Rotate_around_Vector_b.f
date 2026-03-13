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
 
      SUBROUTINE Tool_Vector_a_Rotate_around_Vector_b(a,b,theta,output)
C     Vector a rotate around vector b for angle theta.
c     2020-03-13 
c     https://www.cnblogs.com/wubugui/p/3734627.html

      use Global_Float_Type 
      
      implicit none
      real(kind=FT),intent(in)::a(3),b(3),theta
      real(kind=FT),intent(out)::output(3)
      real(kind=FT) tem(3),ax,ay,az,bx,by,bz
      
      tem(1:3) = b(1:3)
      call Vector_Normalize(3,tem(1:3))   
      
      ax = a(1)
      ay = a(2)
      az = a(3)
      bx = tem(1)
      by = tem(2)
      bz = tem(3)
      
      output(1)= ax*cos(theta)+(by*az-bz*ay)*sin(theta) + 
     &           bx*(bx*ax + ay*ay + az*az)*(ONE - cos(theta))
      output(2)= ay*cos(theta)+(bx*az-bz*ax)*sin(theta) + 
     &           by*(bx*ax + ay*ay + az*az)*(ONE  - cos(theta))
      output(3)= az*cos(theta)+(bx*ay-by*ax)*sin(theta) + 
     &           bz*(bx*ax + ay*ay + az*az)*(ONE  - cos(theta))
      
      return
      END SUBROUTINE Tool_Vector_a_Rotate_around_Vector_b