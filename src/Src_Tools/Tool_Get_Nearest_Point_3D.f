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
 
      subroutine Tool_Get_Nearest_Point_3D(num_Point,Points,A,
     &                                     Point_Num) 
c     Find the point closest to A from a set of 3D points.
c     2022-04-30.
c     NEWFTU2022043007.
 
      use Global_Float_Type
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::Points(num_Point,3),A(3)       
      integer,intent(out)::Point_Num
      real(kind=FT) Dis_Vector(num_Point)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      integer i_P
      
      if(num_Point<=1)then
        print *, '    Error :: num_Point<=1!'
        print *, '             in Tool_Get_Nearest_Point_3D.f!'
        call Warning_Message('S',Keywords_Blank)      
      endif
      
      
      do i_P=1,num_Point
          Dis_Vector(i_P) = Tool_Function_2Point_Dis_3D(
     &                              A,Points(i_P,1:3))
      enddo
      
      Point_Num = MINLOC(Dis_Vector(1:num_Point),1)
      
      return 
      end SUBROUTINE Tool_Get_Nearest_Point_3D                     
