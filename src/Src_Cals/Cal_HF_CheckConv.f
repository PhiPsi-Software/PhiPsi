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
 
      subroutine Cal_HF_CheckConv(iter,
     &                    Last_Cr_CalP_Aper,
     &                    Last_Cr_CalP_Pres,
     &                    Yes_Conv)
      ! Hydraulic Fracturing NR Iterative Convergence Check
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      
      integer,intent(in)::iter
      real(kind=FT),intent(in):: Last_Cr_CalP_Aper(Max_Num_Cr,
     &                                              Max_Num_Cr_CalP)
      real(kind=FT),intent(in):: Last_Cr_CalP_Pres(Max_Num_Cr,
     &                                              Max_Num_Cr_CalP)
      logical,intent(out):: Yes_Conv
      
      integer i_C,num_CalP
      real(kind=FT) Tol,
     &                 Conv_Aper(Max_Num_Cr),
     &                 Max_Conv_Aper,
     &                 Conv_Pres(Max_Num_Cr),
     &                 Max_Conv_Pres,tem1(Max_Num_Cr_CalP),
     &                 tem2(Max_Num_Cr_CalP),norm1,norm2
      
      Tol = MNR_Tol
      
      do i_C = 1,num_Crack
          num_CalP = Cracks_CalP_Num(i_C)
          tem1(1:num_CalP) = Cracks_CalP_Aper(i_C,1:num_CalP)-
     &                       Last_Cr_CalP_Aper(i_C,1:num_CalP)
          tem2(1:num_CalP) = Last_Cr_CalP_Aper(i_C,1:num_CalP)     
          call Vector_Norm2(num_CalP,tem1(1:num_CalP),norm1) 
          call Vector_Norm2(num_CalP,tem2(1:num_CalP),norm2) 
          Conv_Aper(i_C)=  norm1/norm2
          tem1(1:num_CalP) = Cracks_CalP_Pres(i_C,1:num_CalP)-
     &                       Last_Cr_CalP_Pres(i_C,1:num_CalP)
          tem2(1:num_CalP) = Last_Cr_CalP_Pres(i_C,1:num_CalP)     
          call Vector_Norm2(num_CalP,tem1(1:num_CalP),norm1) 
          call Vector_Norm2(num_CalP,tem2(1:num_CalP),norm2) 
          Conv_Pres(i_C)=  norm1/norm2
      end do
      Max_Conv_Aper = maxval(Conv_Aper(1:num_Crack))
      Max_Conv_Pres = maxval(Conv_Pres(1:num_Crack))      
      
      
      
      
      if (Key_HF_Conv_Crite == 1) then
          if (Max_Conv_Aper <= Tol) then                                    
              Yes_Conv =.True.
              write(*,1000) 
              write(*,1002) Max_Conv_Aper 
              write(*,1005) Max_Conv_Pres
              write(*,1003)
              write(*,1004) 
              write(*,1000) 
          else
              Yes_Conv = .False.
              write(*,1000) 
              write(*,1002) Max_Conv_Aper
              write(*,1005) Max_Conv_Pres
              write(*,1000)
          end if
      elseif (Key_HF_Conv_Crite == 2) then
          if ((Max_Conv_Aper <= Tol) .and .
     &        (Max_Conv_Pres <= Tol)) then                                    
              Yes_Conv =.True.
              write(*,1000) 
              write(*,1002) Max_Conv_Aper 
              write(*,1005) Max_Conv_Pres
              write(*,1003)
              write(*,1004) 
              write(*,1000) 
          else
              Yes_Conv = .False.
              write(*,1000) 
              write(*,1002) Max_Conv_Aper
              write(*,1005) Max_Conv_Pres
              write(*,1000)
          end if          
      end if
      
 1000 FORMAT('     -#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#-')   
 1002 FORMAT('         Aperture EOC of iteration is ',E12.5)
 1005 FORMAT('         Pressure EOC of iteration is ',E12.5)
 1003 FORMAT('         ----------------------------------') 
 1004 FORMAT('         Convergence condition is satisfied.')
 
      return 
      end SUBROUTINE Cal_HF_CheckConv             
