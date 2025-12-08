 
      SUBROUTINE Matrix_Condition_Number_v2(norm_Type,n,Matrix,
     &                                      Condition_Num_Norm)

      use Global_Float_Type
      
      implicit none
      integer,intent(in)::n,norm_Type
      real(kind=FT),intent(in)::Matrix(n,n)
      real(kind=FT),intent(out)::Condition_Num_Norm
      
      real(kind=FT) Norm_2_K,Norm_2_K_inv
      real(kind=FT) Norm_1_K,Norm_1_K_inv
      real(kind=FT) Matrix_inv(n,n)
      
      
      call Matrix_Inverse(Matrix,Matrix_inv,n)   
      
      if (norm_Type==2) then
          call Matrix_Norm2(n,Matrix,Norm_2_K) 

          call Matrix_Norm2(n,Matrix_inv,Norm_2_K_inv) 
          
          Condition_Num_Norm =Norm_2_K* Norm_2_K_inv
      endif

      if (norm_Type==1) then
          call Matrix_Norm1(n,Matrix,Norm_1_K) 

          call Matrix_Norm1(n,Matrix_inv,Norm_1_K_inv) 
          
          Condition_Num_Norm =Norm_1_K* Norm_1_K_inv
      endif      
      
      return
      END SUBROUTINE Matrix_Condition_Number_v2
      
