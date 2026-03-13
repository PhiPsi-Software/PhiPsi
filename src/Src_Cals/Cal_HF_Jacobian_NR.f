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
 
      subroutine Cal_HF_Jacobian_NR(Counter,NR_Deri,Total_FD,
     &                              num_FreeD,
     &                              num_free_CalP,
     &                              globalK,Coupled_Q,
     &                              freeDOF,freeDOF_HF,
     &                              Local_freeDOF_HF,delta_Time,H)   
      ! Hydraulic Fracturing Newton-Raphson Iterative Algorithm Assembled Derivative Matrix (Jacobian
      ! Matrix)
      
c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common

      implicit none
      integer, intent(in)::Counter,Total_FD,num_FreeD,
     &                     num_free_CalP
      real(kind=FT),intent(in)::
     &              Coupled_Q(Total_FD,num_Tol_CalP_Water),
     &              globalK(Total_FD,Total_FD),delta_Time,
     &              H(num_Tol_CalP_Water,num_Tol_CalP_Water)
      integer, intent(in)::freeDOF(Total_FD),
     &                     freeDOF_HF(num_Tol_CalP_Water),
     &                     Local_freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::
     &  NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
      print *,'    Calculating the Jacobian matrix of ' 
     &             //'NR interation......'
      
      NR_Deri(1:Total_FD+num_Tol_CalP_Water,
     &        1:Total_FD+num_Tol_CalP_Water)=ZR
      
      NR_Deri(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)) = 
     &          globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD))
      NR_Deri(freeDOF(1:num_FreeD),freeDOF_HF(1:num_free_CalP)) = 
     &           -Coupled_Q(freeDOF(1:num_FreeD),
     &                      Local_freeDOF_HF(1:num_free_CalP))
      NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF(1:num_FreeD)) = 
     &         -transpose(Coupled_Q(freeDOF(1:num_FreeD),
     &                             Local_freeDOF_HF(1:num_free_CalP)))
      NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF_HF(1:num_free_CalP))=
     &          -H(Local_freeDOF_HF(1:num_free_CalP),
     &             Local_freeDOF_HF(1:num_free_CalP))*delta_Time
      
      return 
      end SUBROUTINE Cal_HF_Jacobian_NR          
