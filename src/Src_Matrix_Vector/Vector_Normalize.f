 
      SUBROUTINE Vector_Normalize(n,Vector)   

      use Global_Float_Type
      use Global_Common  
      
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(inout):: Vector(n)
      real(kind=FT) Norm_2
     
      call Vector_Norm2(n,Vector,Norm_2)   
      if (Norm_2==ZR)then
          Norm_2 = Tol_30
      endif
      
      Vector = Vector/Norm_2
      
      return
      END SUBROUTINE Vector_Normalize
    


