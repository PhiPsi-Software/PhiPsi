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
 
      subroutine Tool_Length_of_Crack(i_Crack,Length_Crack)
C     Calculate the length of the specified cracks, which is the sum of the lengths of all crack segments.
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      
      implicit none
      integer,intent(in)::i_Crack
      real(kind=FT),intent(out)::Length_Crack
      integer c_C,i_S
      real(kind=FT) crack_p1(2),crack_p2(2),delta_L
      
      if(i_Crack>num_Crack)then
          print *,'     Error:: wrong crack number for '
     &            // 'Tool_Length_of_Crack.f!'
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      Length_Crack = ZR
      c_C = i_Crack   
      
      do i_S = 1,Each_Cr_Poi_Num(c_C)-1
          crack_p1 = [Crack_Coor(c_C,i_S,1),Crack_Coor(c_C,i_S,2)]
          crack_p2 = [Crack_Coor(c_C,i_S+1,1),Crack_Coor(c_C,i_S+1,2)]
          delta_L  =  sqrt((crack_p2(2)-crack_p1(2))**2 +
     &                     (crack_p2(1)-crack_p1(1))**2)
          Length_Crack = Length_Crack +  delta_L
      end do
      
      return 
      end SUBROUTINE Tool_Length_of_Crack                          
