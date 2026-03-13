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
 
      subroutine Cal_HF_Initial_Aperture(iter,ifra,Counter,
     &                                   Last_Cracks_CalP_Num_ifra,
     &                                   Last_Cr_CalP_Aper_ifra,
     &                                   Initial_Cr_CalP_Aper) 
      ! Option1: If data is not inherited between each rupture step, set the opening to 0 at each rupture
      ! step.
      ! Option 2: If data is inherited between each fracture step, then read the data at the end of the
      ! previous iteration step and assign it to the current fracture step.
      ! The first iteration step
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      
      integer,intent(in)::iter,ifra,Counter
      real(kind=FT),intent(in)::Last_Cr_CalP_Aper_ifra(Max_Num_Cr,
     &                                   Max_Num_Cr_CalP)
      integer,intent(in)::Last_Cracks_CalP_Num_ifra(Max_Num_Cr)
      real(kind=FT),intent(out)::Initial_Cr_CalP_Aper(Max_Num_Cr,
     &                                   Max_Num_Cr_CalP) 
      
      integer i_C,Num_Div_Points
      integer Last_Num_Div_Points
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then
              if (Key_IniPre_PassOn==0)then
                  Num_Div_Points = Cracks_CalP_Num(i_C)
                  Initial_Cr_CalP_Aper(i_C,1:Num_Div_Points)=ZR
              elseif (Key_IniPre_PassOn==1)then
                  if (ifra==1) then
                      Num_Div_Points = Cracks_CalP_Num(i_C)
                      Initial_Cr_CalP_Aper(i_C,1:Num_Div_Points)=ZR  
                  else
                      
                      Num_Div_Points = Cracks_CalP_Num(i_C)
                      Last_Num_Div_Points=Last_Cracks_CalP_Num_ifra(i_C)
                      Initial_Cr_CalP_Aper(i_C,1:Last_Num_Div_Points)=
     &                Last_Cr_CalP_Aper_ifra(i_C,1:Last_Num_Div_Points)
                      Initial_Cr_CalP_Aper(i_C,
     &                         Last_Num_Div_Points+1:Num_Div_Points)=
     &                                                            ZR
                  endif
              end if
          end if
      end do
      
      
      return 
      end SUBROUTINE Cal_HF_Initial_Aperture             
