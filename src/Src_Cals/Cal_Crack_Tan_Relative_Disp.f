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
 
      subroutine Cal_Crack_Tan_Relative_Disp(isub,c_DISP,
     &                                       Cracks_CalP_Tan_Aper)
      ! Calculate tangential crack width and tangential relative sliding displacement
      
      !+++++++++++++++++++++++++++++
      ! Read public variable module
      !+++++++++++++++++++++++++++++
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Common
      use Global_HF
      
      implicit none
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      real(kind=FT),intent(out)::Cracks_CalP_Tan_Aper(num_Crack,
     &                                                Max_Num_Cr_CalP)
      integer,intent(in)::isub
      integer i_C,i_P
      integer Num_CalP
      real(kind=FT) CalP_Points(Max_Num_Cr_CalP,2)   
      real(kind=FT) P_Aperture,Cr_Omega
      real(kind=FT) c_Tangent_disp
      
      Cracks_CalP_Tan_Aper(1:num_Crack,1:Max_Num_Cr_CalP) = ZR
      
      do i_C = 1,num_Crack
          Num_CalP = Cracks_CalP_Num(i_C)
          
          CalP_Points(1:Num_CalP,1:2)  =
     &         Cracks_CalP_Coors(i_C,1:Num_CalP,1:2)  
          
          do i_P =1,Num_CalP
              Cr_Omega = Cracks_CalP_Orient(i_C,i_P)
              call Cal_Point_Aperture(i_C,CalP_Points(i_P,1:2),
     &                                c_DISP,Cr_Omega,P_Aperture)
                  call Cal_Point_Aperture_and_Tangent_disp(i_C,
     &                             CalP_Points(i_P,1:2),
     &                              c_DISP,Cr_Omega,P_Aperture,
     &                              c_Tangent_disp)  
              Cracks_CalP_Tan_Aper(i_C,i_P) = c_Tangent_disp
          enddo
          

      end do
      
      return
      end SUBROUTINE Cal_Crack_Tan_Relative_Disp               
