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
 
      subroutine Cal_HF_Resid(ifra,iter,Counter_Iter,R,Total_FD,
     &                              num_FreeD,
     &                              num_free_CalP,
     &                              globalK,Coupled_Q,F_U,
     &                              Last_DISP,Last_Last_DISP,
     &                              freeDOF,freeDOF_HF,Local_freeDOF_HF,
     &                              Last_CalP_Pres,delta_Time,
     &                              H,c_S)
      ! Hydraulic Fracturing Newton-Raphson Iterative Algorithm for Calculating Residual Vectors
      
c     ----------------------------
c     Read Public Variable Module
c     ----------------------------
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Material
      
      implicit none
      integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,
     &                     num_free_CalP
      real(kind=FT),intent(in)::
     &              Coupled_Q(Total_FD,num_Tol_CalP_Water),
     &              globalK(Total_FD,Total_FD),
     &              F_U(Total_FD),Last_DISP(Total_FD),
     &              Last_Last_DISP(Total_FD),
     &              Last_CalP_Pres(num_Tol_CalP_Water),delta_Time,
     &              H(num_Tol_CalP_Water,num_Tol_CalP_Water),
     &              c_S(num_Tol_CalP_Water)
      integer, intent(in)::freeDOF(Total_FD),
     &                     freeDOF_HF(num_Tol_CalP_Water),
     &                     Local_freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
      real(kind=FT) delta_DISP(Total_FD)
      
      R(1:Total_FD+num_Tol_CalP_Water) = ZR

      R(freeDOF(1:num_FreeD))  =
     &    MATMUL(globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)),
     &            Last_DISP(freeDOF(1:num_FreeD)))
     &   -MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                     Local_freeDOF_HF(1:num_free_CalP)),
     &  thick(1)*Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))     
     &   -F_U(freeDOF(1:num_FreeD))

      if(Key_IniPre_PassOn==0) then
          if (iter == 1) then  
              delta_DISP = Last_DISP
          elseif (iter > 1) then
              delta_DISP(freeDOF(1:num_FreeD)) = 
     &                       Last_DISP(freeDOF(1:num_FreeD))- 
     &                       Last_Last_DISP(freeDOF(1:num_FreeD))
          end if
      elseif(Key_IniPre_PassOn==1) then
          if (iter == 1) then  
              if (ifra == 1) then
                  delta_DISP = Last_DISP
              elseif(ifra > 1) then
                  delta_DISP = Last_DISP
              endif          
          elseif (iter > 1) then
              delta_DISP(freeDOF(1:num_FreeD)) = 
     &                       Last_DISP(freeDOF(1:num_FreeD))- 
     &                       Last_Last_DISP(freeDOF(1:num_FreeD))
          end if              
      endif      
      R(freeDOF_HF(1:num_free_CalP)) = -MATMUL(transpose(
     &    Coupled_Q(freeDOF(1:num_FreeD),
     &              Local_freeDOF_HF(1:num_free_CalP))),
     &    delta_DISP(freeDOF(1:num_FreeD)))
     &    - delta_Time*MATMUL(H(Local_freeDOF_HF(1:num_free_CalP),
     &               Local_freeDOF_HF(1:num_free_CalP)),
     &      Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))
     &    + c_S(Local_freeDOF_HF(1:num_free_CalP))*delta_Time
      return 
      end SUBROUTINE Cal_HF_Resid          
