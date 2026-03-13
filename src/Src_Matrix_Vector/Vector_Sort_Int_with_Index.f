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
 
      SUBROUTINE Vector_Sort_Int_with_Index(n,Vector,V_Index)   
C     One-dimensional array sorting, integer.
c     The V_Index array stores the positions of the sorted elements in the original array.
c     2022-07-27.

      use Global_Float_Type     
      implicit none
      integer,intent(in)::n
      integer,intent(inout):: Vector(n)
      integer,intent(out):: V_Index(n)
      integer Vector_Old(n)
      integer i,j,a
      logical c_Yes_In
      
      Vector_Old = Vector
      
      do j=2, n
          a=Vector(j)
          do i=j-1,1,-1
              if (Vector(i).le.a) goto 10
              Vector(i+1)=Vector(i)
              
              
          end do
          i=0
          
   10     Vector(i+1)=a
   
      end do
      
      
      do i=1,n
          call Vector_Location_Int(n,Vector_Old,Vector(i),
     &                             V_Index(i),c_Yes_In)   
      enddo
      
      
      return
      END SUBROUTINE Vector_Sort_Int_with_Index
    


