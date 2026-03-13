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
 
      subroutine Tool_Get_Farthest_Two_Points_3D(num_Point,Points,
     &                                       Point1_Num,Point2_Num) 
c     Find the two points that are farthest apart in a 3D point set.
c     2022-04-30.
c     NEWFTU2022043004.
c     The simplest algorithm, algorithm complexity: O(n^2)/2
 
      use Global_Float_Type
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::Points(num_Point,3)       
      integer,intent(out)::Point1_Num,Point2_Num 
      real(kind=FT) Dis_Matrix(num_Point,num_Point)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      integer i_P,j_P
      
      if(num_Point<=1)then
        print *, '    Error :: num_Point<=1!'
        print *, '             in Tool_Get_Farthest_Points_3D.f!'
        call Warning_Message('S',Keywords_Blank)      
      endif
      Dis_Matrix = -TEN_15
      
      
      do i_P=1,num_Point
          do j_P=i_P+1,num_Point
              Dis_Matrix(i_P,j_P) = Tool_Function_2Point_Dis_3D(
     &                              Points(i_P,1:3),
     &                              Points(j_P,1:3))
              print *,i_P,j_P,Dis_Matrix(i_P,j_P) 
          enddo
      enddo
      
      call Matrix_Max_Location_Dou(num_Point,num_Point,Dis_Matrix,
     &                             Point1_Num,Point2_Num)
      if(Point1_Num== Point2_Num)then
        print *, '    Error :: Point1_Num = Point2_Num!'
        print *, '             in Tool_Get_Farthest_Points_3D.f!'
        call Warning_Message('S',Keywords_Blank)  
      endif
      
      return 
      end SUBROUTINE Tool_Get_Farthest_Two_Points_3D                    
