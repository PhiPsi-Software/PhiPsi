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
 
      subroutine Cal_HF_Jacobian_PenaltyFunc(Counter,NR_Deri,Total_FD,
     &                              num_FreeD,
     &                              num_free_CalP,
     &                              freeDOF,freeDOF_HF)   
      ! Penalty function method for handling junction points with the same water pressure
      ! For details on the penalty function method, see 'The Finite Element Method in Engineering'
      ! translated by Zeng Pan, 3rd edition, Section 3.8
      
c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common

      implicit none
      integer, intent(in)::Counter,Total_FD,num_FreeD,
     &                     num_free_CalP
      integer, intent(in)::freeDOF(Total_FD),
     &                     freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(inout)::
     &  NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
      real(kind=FT) k,max_K(3)
      integer i_Pair,CalP_1,CalP_2
      print *,'    Modifing the Jacobian matrix by the penalty function' 
     &             //' method......'
      
      max_K(1) = maxval(abs(NR_Deri(freeDOF(1:num_FreeD),
     &                              freeDOF_HF(1:num_free_CalP))))
      max_K(2) = maxval(abs(NR_Deri(freeDOF_HF(1:num_free_CalP),
     &                              freeDOF(1:num_FreeD))))
      max_K(3) = maxval(abs(NR_Deri(freeDOF_HF(1:num_free_CalP),
     &                              freeDOF_HF(1:num_free_CalP))))
      k = 1.0D4*maxval(max_K)
      do i_Pair = 1,Num_JunPair
          CalP_1 = Cracks_JunPair(i_Pair,1)
          CalP_2 = Cracks_JunPair(i_Pair,2)
          
          
          NR_Deri(Total_FD+CalP_1,Total_FD+CalP_1) = 
     &                  NR_Deri(Total_FD+CalP_1,Total_FD+CalP_1) + k
          NR_Deri(Total_FD+CalP_1,Total_FD+CalP_2) = 
     &                  NR_Deri(Total_FD+CalP_1,Total_FD+CalP_2) - k
          NR_Deri(Total_FD+CalP_2,Total_FD+CalP_1) = 
     &                  NR_Deri(Total_FD+CalP_2,Total_FD+CalP_1) - k
          NR_Deri(Total_FD+CalP_2,Total_FD+CalP_2) = 
     &                  NR_Deri(Total_FD+CalP_2,Total_FD+CalP_2) + k
      end do
      return 
      end SUBROUTINE Cal_HF_Jacobian_PenaltyFunc     
