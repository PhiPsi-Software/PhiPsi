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
 
      subroutine Tool_Dis_Point_to_Line_AB_3D(P,A,B,Dis,foot_point,Key)
C     Calculate the normal distance from point P to the spatial line segment AB (output the foot of the perpendicular coordinates)
c     Key =1 Treat AB as a straight line; Key =2 Treat AB as a line segment
c     2020-01-04
c     Reference: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html or
c 	  theoretical_documents\006 Distance from a Point to a Straight Line in Space-2020-01-04.pdf

      use Global_Float_Type 
      use Global_Common
      implicit none
      real(kind=FT), intent(in) :: P(3),A(3),B(3)
      integer, intent(in) :: Key
      real(kind=FT), intent(out) :: Dis,foot_point(3)
      real(kind=FT) tem1(3),Norm_1,Norm_2,c_Dis(3)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      real(kind=FT) c_dot_product,para_t

      
      call Vector_Cross_Product_3(P-A,P-B,tem1)  
      call Vector_Norm2(3,tem1,Norm_1)   
      call Vector_Norm2(3,B-A,Norm_2)  
      c_dot_product = DOT_PRODUCT(A-P,B-A)
      para_t = -c_dot_product/(Norm_2**2)
      
      if(Norm_2<=Tol_11) then
          print *, '    Error :: in Tool_Dis_Point_to_Line_AB_3D.f!'
          print *, '             AB is too short!'
          call Warning_Message('S',Keywords_Blank) 
          return
      endif
      
      if (Key==1) then
          Dis = Norm_1/Norm_2
          foot_point(1) = A(1)+(B(1)-A(1))*para_t
          foot_point(2) = A(2)+(B(2)-A(2))*para_t
          foot_point(3) = A(3)+(B(3)-A(3))*para_t
      elseif(Key==2)then
          c_Dis(1) = Norm_1/Norm_2
          c_Dis(2) = Tool_Function_2Point_Dis_3D(P,A)
          c_Dis(3) = Tool_Function_2Point_Dis_3D(P,B)
          Dis = minval(c_Dis)
          foot_point(1) = A(1)+(B(1)-A(1))*para_t
          foot_point(2) = A(2)+(B(2)-A(2))*para_t
          foot_point(3) = A(3)+(B(3)-A(3))*para_t
      endif
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_Line_AB_3D                     
