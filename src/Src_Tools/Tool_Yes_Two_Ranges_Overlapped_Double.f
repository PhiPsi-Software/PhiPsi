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
 
      subroutine Tool_Yes_Two_Ranges_Overlapped_Double(Range_1,Range_2,
     &                            Logical_Yes)   
     
c     Used to check whether the ranges of two double-precision numbers intersect (overlap). Can be used to check if coordinate ranges overlap. NEWFTU2022043001.
c     Added on 2022-04-30. 
c     Ref: https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-if-two-ranges-overlap
c     or
c     \theory_documents\032 Check Whether Two Data Ranges Intersect (Overlap) - 2022-04-30.pdf

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in)::Range_1(2),Range_2(2)
      logical,intent(out)::Logical_Yes
      real(kind=FT) max_x1_y1,min_x2_y2
      
      Logical_Yes =.False.
      
      max_x1_y1 = max(Range_1(1),Range_2(1))
      min_x2_y2 = min(Range_1(2),Range_2(2))
      
      if(max_x1_y1 <= min_x2_y2)then
          Logical_Yes =.True.
      endif
      
      if(abs(max_x1_y1-min_x2_y2) <=Tol_11) then
          Logical_Yes =.True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Two_Ranges_Overlapped_Double          
