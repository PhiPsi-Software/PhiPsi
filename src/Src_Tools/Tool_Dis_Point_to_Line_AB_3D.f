 
      subroutine Tool_Dis_Point_to_Line_AB_3D(P,A,B,Dis,foot_point,Key)

      use Global_Float_Type 
      use Global_Common
      implicit none
      real(kind=FT), intent(in) :: P(3),A(3),B(3)
      integer, intent(in) :: Key
      real(kind=FT), intent(out) :: Dis,foot_point(3)
      real(kind=FT) tem1(3),Norm_1,Norm_2,c_Dis(3)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      real(kind=FT) c_dot_product,para_t

      
      call Vector_Cross_Product_3(P-A,P-B,tem1)  
      call Vector_Norm2(3,tem1,Norm_1)   
      call Vector_Norm2(3,B-A,Norm_2)  
      c_dot_product = DOT_PRODUCT(A-P,B-A)
      para_t = -c_dot_product/(Norm_2**2)
      
      if(Norm_2<=Tol_11) then
          print *, '    Error :: in Tool_Dis_Point_to_Line_AB_3D.f!'
          print *, '             AB is too short!'
          call Warning_Message('S',Keywords_Blank) 
          return
      endif
      
      if (Key==1) then
          Dis = Norm_1/Norm_2
          foot_point(1) = A(1)+(B(1)-A(1))*para_t
          foot_point(2) = A(2)+(B(2)-A(2))*para_t
          foot_point(3) = A(3)+(B(3)-A(3))*para_t
      elseif(Key==2)then
          c_Dis(1) = Norm_1/Norm_2
          c_Dis(2) = Tool_Function_2Point_Dis_3D(P,A)
          c_Dis(3) = Tool_Function_2Point_Dis_3D(P,B)
          Dis = minval(c_Dis)
          foot_point(1) = A(1)+(B(1)-A(1))*para_t
          foot_point(2) = A(2)+(B(2)-A(2))*para_t
          foot_point(3) = A(3)+(B(3)-A(3))*para_t
      endif
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_Line_AB_3D                     
