 
      subroutine Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor(
     &  Stress,n,Normal_Stress)
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::Stress(6),n(3)
      real(kind=FT),intent(out)::Normal_Stress
      real(kind=FT) Stress_Matrx(3,3),t_Vector(3)
      
      real(kind=FT) n_vector(3)
      
      
      n_vector = n
      
      call Vector_Normalize(3,n_vector)
      
      Stress_Matrx(1,1:3) = [Stress(1),Stress(4),Stress(6)]
      Stress_Matrx(2,1:3) = [Stress(4),Stress(2),Stress(5)]
      Stress_Matrx(3,1:3) = [Stress(6),Stress(5),Stress(3)]
      
      t_Vector = MATMUL(Stress_Matrx,n_vector)
      
      Normal_Stress = DOT_PRODUCT(t_Vector,n_vector)
      
      return 
      end SUBROUTINE Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor          
