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
 
      subroutine Cal_F(r,theta,omega,c_mat_type,F)
      
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Common
      implicit none
      
      real(kind=FT),intent(in)::r,theta,omega
      integer,intent(in)::c_mat_type
      real(kind=FT),intent(out)::F(4)
      
      real(kind=FT) fac,st2,ct2,s3t2,c3t2,st,ct,
     &                 r2,c_theta
     
      c_theta = theta
      if (r.ne.ZR) then
          r2 = sqrt(r)
      else
          r2 = sqrt(Ave_Elem_Area)*0.1D-4
            c_theta = ZR
      end if
      
      
      if (c_mat_type ==1 ) then
            fac = HLF/r2
            st2 = sin(c_theta/TWO)
            ct2 = cos(c_theta/TWO)
            s3t2= sin(1.5D0*c_theta)
            c3t2= cos(1.5D0*c_theta)
            st  = sin(c_theta)
            ct  = cos(c_theta)

        F(1) = r2*st2
        F(2) = r2*ct2
        F(3) = r2*st2*st
        F(4) = r2*ct2*st
      elseif (c_mat_type ==2 .or. c_mat_type ==3) then
      
      end if
      
      if(Key_TipEnrich==4)then
          F(1) = r*st2
          F(2) = ZR
          F(3) = ZR
          F(4) = ZR
      endif
      return 
      end SUBROUTINE Cal_F                
