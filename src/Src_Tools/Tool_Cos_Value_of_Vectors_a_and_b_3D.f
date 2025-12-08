 
      subroutine Tool_Cos_Value_of_Vectors_a_and_b_3D(a,b,cos_value)

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
