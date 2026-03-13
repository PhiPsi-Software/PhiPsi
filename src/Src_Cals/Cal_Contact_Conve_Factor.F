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
 
      SUBROUTINE Cal_Contact_Conve_Factor(
     &                 iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance,
     &                 c_Total_Freedom,c_freeDOF,c_num_freeDOF,
     &                 c_F,R_PSI,Last_R_PSI,
     &                 delta_U,U,Contact_DISP_0,
     &                 Yes_Conve,Conve_Factor)
c     Check for convergence and calculate the convergence factor

c     The following procedures together make up the penalty function contact detection algorithm:
c     Cal_Contact_Jacobin.f
c     Cal_Contact_Contact_State_Gauss.f
c     Cal_Contact_PN_and_PT.f
c     Cal_Contact_Resid.f
c     Cal_Contact_Conve_Factor.f

      !**********************
      ! Read public variable
      !**********************
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Common
      
      implicit none
      integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P,c_Total_Freedom
     &                    ,c_num_freeDOF
      real(kind=FT),intent(in)::Conve_Tolerance
      integer, intent(in)::c_freeDOF(c_Total_Freedom)
      real(kind=FT),intent(in)::Contact_DISP_0(c_Total_Freedom)
      real(kind=FT),intent(in)::
     &          c_F(c_Total_Freedom),R_PSI(c_Total_Freedom),
     &          Last_R_PSI(c_Total_Freedom)
      real(kind=FT),intent(in)::
     &          delta_U(c_Total_Freedom),U(c_Total_Freedom)
      logical,intent(out)::Yes_Conve
      real(kind=FT),intent(out)::Conve_Factor
      real(kind=FT) Norm_R_PSI,Norm_U,Norm_delta_U
      real(kind=FT) Old_Max_Aperture(num_Crack)
      real(kind=FT) Old_Min_Aperture(num_Crack)
      real(kind=FT) New_Max_Aperture(num_Crack)
      real(kind=FT) New_Min_Aperture(num_Crack)
      real(kind=FT) Delta_Max_Aperture(num_Crack)
      real(kind=FT) Delta_Min_Aperture(num_Crack)
      real(kind=FT) Norm_Delta_Max_Aperture
      real(kind=FT) Norm_Delta_Min_Aperture
      real(kind=FT) Norm_New_Max_Aperture
      real(kind=FT) Norm_New_Min_Aperture
      real(kind=FT) Conve_Factor_1,Conve_Factor_2
      
      real(kind=FT) Old_Ave_Aperture(num_Crack)
      real(kind=FT) New_Ave_Aperture(num_Crack)
      real(kind=FT) Delta_Ave_Aperture(num_Crack)
      real(kind=FT) Norm_Delta_Ave_Aperture
      real(kind=FT) Norm_New_Ave_Aperture
      

      if(Key_Conta_ConCrit==1)then
          call Vector_Norm2(c_num_freeDOF,
     &                  R_PSI(c_freeDOF(1:c_num_freeDOF)),Norm_R_PSI)   
          Conve_Factor = Norm_R_PSI / Norm2_Contact_R_PSI_0
      endif
      
      if(Key_Conta_ConCrit==2)then
          call Vector_Norm2(c_Total_Freedom,delta_U,Norm_delta_U)   
          call Vector_Norm2(c_Total_Freedom,Contact_DISP_0,Norm_U) 
          Conve_Factor = Norm_delta_U/Norm_U
      endif
      
      if(Key_Conta_ConCrit==3)then
#ifndef Silverfrost
          Old_Ave_Aperture(1:num_Crack)=
     &                        Crack_Max_Min_Aperture(1:num_Crack,3)  
          call Cal_Crack_Aperture_3D(iter,U)
          New_Ave_Aperture(1:num_Crack) =
     &               Crack_Max_Min_Aperture(1:num_Crack,3)   
          Delta_Ave_Aperture(1:num_Crack)=New_Ave_Aperture(1:num_Crack)
     &                        - Old_Ave_Aperture(1:num_Crack)
          call Vector_Norm2(num_Crack,Delta_Ave_Aperture(1:num_Crack),
     &                    Norm_Delta_Ave_Aperture)    
          call Vector_Norm2(num_Crack,New_Ave_Aperture(1:num_Crack),
     &                    Norm_New_Ave_Aperture)    
          Conve_Factor = Norm_Delta_Ave_Aperture/Norm_New_Ave_Aperture
#endif
      endif
      
      Yes_Conve = .False.
      if(Conve_Factor <= Conve_Tolerance)then
          Yes_Conve = .True.
      endif
      
      RETURN
      END SUBROUTINE Cal_Contact_Conve_Factor
