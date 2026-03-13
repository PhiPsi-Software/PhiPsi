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
 
      subroutine Vector_Range_Number_of_Value_Dou(value_to_be_find,n,
     &               values,Yes_In,
     &               Range_Num)   
C     Get the position of the value in the vector (which segment it is in). The vector should be ensured to be in ascending order.
c     For example, value_to_be_find: 7.0 is at position 2 in values: 1.0, 5.0, 9.0, 11.1. Range_Num = 2.
C     2022-07-06.
    
      use Global_Float_Type
      implicit none

      integer,intent(in)::n
      real(kind=FT),intent(in)::value_to_be_find,values(n)
      integer,intent(out)::Range_Num
      logical,intent(out)::Yes_In
      integer i               
      
      Yes_In = .False.

      do i=1,n-1
          if(value_to_be_find >= values(i)   .and. 
     &       value_to_be_find <= values(i+1)  )then
              Range_Num = i
              Yes_In = .True.
              return
          endif
      end do
      
      return 
      end subroutine Vector_Range_Number_of_Value_Dou
    


