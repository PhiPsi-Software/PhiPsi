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
 
      subroutine Tool_Cos_Value_of_Vectors_a_and_b_3D(a,b,cos_value)
c     cos=a*b/(|a|*|b|) 

      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::a(3),b(3)
      real(kind=FT),intent(out)::cos_value
      real(kind=FT) x1,y1,z1,x2,y2,z2
      real(kind=FT) a_plus_b,norm_a,norm_b
      
      x1 = a(1)
      y1 = a(2)
      z1 = a(3)
      x2 = b(1)
      y2 = b(2)
      z2 = b(3)
      
      if(abs(x1+x2)< Tol_11  .and.         
     &   abs(y1+y2)< Tol_11  .and.
     &   abs(z1+z2)< Tol_11) then
         cos_value = -ONE
         goto 100
      endif
      
      a_plus_b = x1*x2+y1*y2+z1*z2
      norm_a   = sqrt(x1**2+y1**2+z1**2)
      norm_b   = sqrt(x2**2+y2**2+z2**2)
      cos_value = a_plus_b/norm_a/norm_b  

  100 continue
  
      return 
      end SUBROUTINE Tool_Cos_Value_of_Vectors_a_and_b_3D        
