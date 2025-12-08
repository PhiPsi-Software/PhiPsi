 
      subroutine Tool_Generate_Rand_Point_in_3D_Circle(
     &                      Center,n_Vector,Radius,
     &                      Rand_Point)

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
