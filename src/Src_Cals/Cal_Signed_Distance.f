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
 
      subroutine Cal_Signed_Distance(Line_AB,Point_C,S_Distance)
C     Calculate Symbol Distance
c     This function calculates the signed distance from the Point_C to the Line_AB.
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in)::Line_AB(2,2),Point_C(2)
      real(kind=FT),intent(out)::S_Distance
      
      real(kind=FT) tem_1(2,2), tem_2(2)
      real(kind=FT) tem_Det,tem_Norm
      
      
      tem_1(1,:) = Line_AB(2,:)-  Line_AB(1,:)
      tem_1(2,:) = Point_C     -  Line_AB(1,:)
      
      
      tem_2(:)   = Line_AB(2,:)-Line_AB(1,:)
      
      
      tem_Det = tem_1(1,1)*tem_1(2,2) - tem_1(2,1)*tem_1(1,2)
      
      call Vector_Norm2(2,tem_2,tem_Norm) 
      
      S_Distance = tem_Det / tem_Norm
      
      return 
      end SUBROUTINE Cal_Signed_Distance                          
