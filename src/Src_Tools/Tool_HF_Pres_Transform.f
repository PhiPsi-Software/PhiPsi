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
 
      subroutine Tool_HF_Pres_Transform(iter,CalP_Pres,
     &                                  temp_Pres)
      ! The water pressure from all calculation points obtained from the fluid-solid coupling calculation
      ! is transferred to each crack for post-processing purposes.
      use Global_Float_Type      
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      
      integer,intent(in)::iter
      real(kind=FT),intent(in)::CalP_Pres(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::temp_Pres(Max_Num_Cr,
     &                                        Max_Num_Cr_CalP)
      
      integer i_C,i_CalP,num_Count
      
      temp_Pres = ZR
      
      num_Count = 0
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then         
              do i_CalP=1,Cracks_CalP_Num(i_C)
                  num_Count = num_Count + 1
                  temp_Pres(i_C,i_CalP) = CalP_Pres(num_Count)
              end do
          end if
      end do

      return 
      end SUBROUTINE Tool_HF_Pres_Transform            
