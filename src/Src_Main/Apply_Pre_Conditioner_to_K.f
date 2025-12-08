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
 
      SUBROUTINE Apply_Pre_Conditioner_to_K(n,
     &                   K,Vector_Dc)    
 
c     Apply the Pre-Conditioner to K.
c     Written first by Fang Shi on 2020-02-12.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      
      integer,intent(in)::n
      real(kind=FT),intent(inout)::K(n,n)
      real(kind=FT),intent(out)::Vector_Dc(n)
      real(kind=FT) diagona(n)
      integer i,j
      real(kind=FT) max_diagona,min_diagona,tem_diagona
      
      !Get diagonal elements of K.
      call Matrix_Get_Diagonal_Elements(n,K,diagona)   
      max_diagona = maxval(diagona)
      min_diagona = minval(diagona)
      print *,'    Max value of  diagonal element of K:',max_diagona
      print *,'    Min value of  diagonal element of K:',min_diagona      
      
      do i =1,n
          if(diagona(i)==ZR)then
              print *,i
          endif
      enddo
      
      !The min diagonal element must be positive.
      if (min_diagona<=ZR)then
          Print *,'    ERROR :: in Apply_Pre_Conditioner_to_K.f!'
          Print *,'    Error code: 777'
          call Warning_Message('S',Keywords_Blank)  
      endif

      tem_diagona = sqrt((max_diagona+min_diagona)/TWO)
      !Loop to get Vector_Dc.
      do i =1,n
          Vector_Dc(i)=tem_diagona/sqrt(diagona(i))
      enddo

      !Apply Vector_Dc to K.
      do i =1,n
        do j =1,n
          K(i,j) = K(i,j)*Vector_Dc(i)*Vector_Dc(j) 
        enddo                      
      enddo
      
      RETURN
      END SUBROUTINE Apply_Pre_Conditioner_to_K