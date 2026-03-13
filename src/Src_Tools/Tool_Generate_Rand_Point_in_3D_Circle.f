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
 
      subroutine Tool_Generate_Rand_Point_in_3D_Circle(
     &                      Center,n_Vector,Radius,
     &                      Rand_Point)
C     Randomly generate points inside a 3D circle.
c     2022-06-16.

      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::Center(3),n_Vector(3),Radius
      real(kind=FT),intent(out)::Rand_Point(3)
      
      real(kind=FT) Tem_1,Ran_theta
      real(kind=FT) a(3),b(3),norm_a,norm_b,Tem_2,Ran_R

      
      if(key_Random==1)then
          call Init_Random_Seed()
      endif      
      if(key_Random==1 .or. key_Random==0)then
          call random_number(Tem_1)
      elseif(key_Random==2 .or. key_Random==3)then
          call Tool_Generate_Random_Number(Tem_1)
      endif     
      Ran_theta = Tem_1*pi
      
      if(key_Random==1 .or. key_Random==0)then
          call random_number (Tem_2)
      elseif(key_Random==2 .or. key_Random==3)then
          call Tool_Generate_Random_Number(Tem_2)
      endif    
      Ran_R = Tem_2*Radius
      
      call Vector_Cross_Product_3(n_Vector,[ONE,ZR,ZR],a)   
      if(sum(abs(a))<=Tol_20) then
          call Vector_Cross_Product_3(n_Vector,[ZR,ONE,ZR],a)   
      endif
      call Vector_Cross_Product_3(n_Vector,a,b)   
      call Vector_Norm2(3,a,norm_a)   
      call Vector_Norm2(3,b,norm_b)  
      a=a/norm_a
      b=b/norm_b


      Rand_Point(1)=Center(1)+Ran_R*a(1)*cos(Ran_theta)+
     &              Ran_R*b(1)*sin(Ran_theta)
      Rand_Point(2)=Center(2)+Ran_R*a(2)*cos(Ran_theta)+
     &              Ran_R*b(2)*sin(Ran_theta)
      Rand_Point(3)=Center(3)+Ran_R*a(3)*cos(Ran_theta)+
     &              Ran_R*b(3)*sin(Ran_theta)
     
      
      return 
      end SUBROUTINE Tool_Generate_Rand_Point_in_3D_Circle                       
