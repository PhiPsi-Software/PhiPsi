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
 
      subroutine Cal_HF_Theoretic_Pressure(iter,Counter_Iter,
     &            Inject_Crack_L,Theor_Time,Theor_Press,Theor_Aper)
      ! Calculate theoretical water pressure and theoretical practice based on crack length and other
      ! parameters, for viscosity-dominated cases only
      
c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Material
      use Global_HF

      implicit none
      integer, intent(in)::iter,Counter_Iter
      real(kind=FT),intent(in)::Inject_Crack_L
      real(kind=FT),intent(out)::Theor_Time,Theor_Press,Theor_Aper
      real(kind=FT)  A0,A1,B1,E_1,u_1
      real(kind=FT) Q_0 
      real(kind=FT) tem,tem1,tem2,tem3,ttt1,ttt2,ttt3,Kesi
      
      A0 = sqrt(THR)
      A1 = -0.156D0
      B1 = 0.0663D0
      E_1 = Material_Para(1,1)/(ONE-Material_Para(1,2)**2)
      u_1 = 12.0D0*Viscosity
      Q_0 = Inject_Q_Val(1)
      
      Theor_Time=((Inject_Crack_L/0.616D0)**SIX*u_1/E_1/(Q_0**THR))
     &           **(ONE/FOU)
      
      Kesi = ZR
      call Tool_Function_BETA(ONE/TWO,TWO/THR,ttt1)
      tem1 = ttt1/(THR*pi)
      ttt2 = ONE
      ttt3 = ONE
      tem2 =  A0*ttt2+TEN/SEV*A1*ttt3
      tem3 =  B1*(TWO-pi*abs(Kesi))
      tem  = tem1*tem2+tem3
      Theor_Press = ((u_1/E_1/Theor_Time)**(ONE/THR))*E_1*tem
      print *,'    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      print *,'    Theoretic pressure(MPa):  ',Theor_Press/1.0D6
      print *,'    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      tem1 = A0+A1
        tem2 = FOU
      tem = tem1 + B1*(tem2+tem3)
      Theor_Aper = ((u_1/E_1/Theor_Time)**(ONE/THR))*
     &              Inject_Crack_L*tem
      
      return 
      end subroutine Cal_HF_Theoretic_Pressure