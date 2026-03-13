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
 
      subroutine Tool_Get_3D_Rect_by_n_Vector_and_Center(n_Vector,
     &                        c_Crack_Center,Crack_L,
     &                        Out_Coors)
c     Obtain a 3D rectangle based on the normal vector, center coordinates, and rectangle side lengths.
c     NEWFTU2022053103.
c     2022-05-31.
c     Algorithm: Implemented using a 3D circle.

      use Global_Float_Type 
      implicit none
      real(kind=FT),intent(in) :: n_Vector(3),c_Crack_Center(3),Crack_L
      real(kind=FT),intent(out) :: Out_Coors(4,3)
      integer num_Circle_Divison
      integer i_theta
      real(kind=FT) theta,a(3),b(3),norm_a,norm_b
      real(kind=FT) c_xyz(4,3),c_Radius
      
      c_Radius = sqrt(TWO)/TWO*Crack_L
      
       call Vector_Cross_Product_3(n_Vector,[ONE,ZR,ZR],a)   
       if(sum(abs(a))<=Tol_20) then
          call Vector_Cross_Product_3(n_Vector,[ZR,ONE,ZR],a)   
       endif
       call Vector_Cross_Product_3(n_Vector,a,b)   
       call Vector_Norm2(3,a,norm_a)   
       call Vector_Norm2(3,b,norm_b)  
       a=a/norm_a
       b=b/norm_b
       
           
      num_Circle_Divison = 4
      do i_theta = 1,num_Circle_Divison
           theta = (i_theta-1)*TWO*pi/num_Circle_Divison - pi/FOU
           c_xyz(i_theta,1)=c_Crack_Center(1)+
     &                 c_Radius*a(1)*cos(theta)
     &                +c_Radius*b(1)*sin(theta)
           c_xyz(i_theta,2)=c_Crack_Center(2)+
     &                 c_Radius*a(2)*cos(theta)
     &                +c_Radius*b(2)*sin(theta)
           c_xyz(i_theta,3)=c_Crack_Center(3)+
     &                 c_Radius*a(3)*cos(theta)
     &                +c_Radius*b(3)*sin(theta)
      enddo
      

      Out_Coors(1,1:3) =  c_xyz(1,1:3)
      Out_Coors(2,1:3) =  c_xyz(2,1:3)
      Out_Coors(3,1:3) =  c_xyz(3,1:3)
      Out_Coors(4,1:3) =  c_xyz(4,1:3)
      
      return 
      end SUBROUTINE Tool_Get_3D_Rect_by_n_Vector_and_Center                      
