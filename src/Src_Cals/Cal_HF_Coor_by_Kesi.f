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
 
      subroutine Cal_HF_Coor_by_Kesi(kesi,x1,y1,x2,y2,
     &                                Out_x,Out_y)
C     Calculate the global coordinates of the hydraulic fracture Gauss point based on the local coordinates kesi
c     x1 and y1 are the coordinates of the first node of the hydraulic fracture element. The node is also called the calculation point.
c     x2 and y2 are the coordinates of the first node of the hydraulic fracture element.
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::kesi,x1,y1,x2,y2
      real(kind=FT),intent(out)::Out_x,Out_y
      real(kind=FT) x0,y0
      real(kind=FT) angle,length
      
      x0     = HLF*(x1+x2)
      y0     = HLF*(y1+y2)
      angle  = atan2(y2-y1,x2-x1)
      length = sqrt((x2-x1)**2 + (y2-y1)**2)
      
      Out_x  = x0  + HLF*length*kesi*cos(angle)
      Out_y  = y0  + HLF*length*kesi*sin(angle)
     
      return 
      end SUBROUTINE Cal_HF_Coor_by_Kesi                         
