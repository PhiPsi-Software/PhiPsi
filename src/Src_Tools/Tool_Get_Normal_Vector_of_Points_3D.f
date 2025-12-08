 
      subroutine Tool_Get_Normal_Vector_of_Points_3D(
     &                 In_Points,num_Point,n_Vector)
      
      use Global_Float_Type
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::n_Vector(3)
      real(kind=FT) M(num_Point,3),MTM_Inv(3,3),MTM(3,3),L(num_Point)
      real(kind=FT) norm_n_Vector
      
      M = In_Points
      MTM = matmul(TRANSPOSE(M),M)
      call Matrix_Inverse_3x3(MTM,MTM_Inv)  
      L(1:num_Point) = ONE
      n_Vector =  matmul(matmul(MTM_Inv,TRANSPOSE(M)),L)
      call Vector_Norm2(3,n_Vector,norm_n_Vector) 
      n_Vector =  n_Vector /norm_n_Vector
      
      
      return 
      end SUBROUTINE Tool_Get_Normal_Vector_of_Points_3D                        
