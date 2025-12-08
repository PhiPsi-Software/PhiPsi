 
      subroutine Tool_Normal_Unit_vector_of_3D_Tri(
     &                            Tri_P1,Tri_P2,Tri_P3,Np)

      use Global_Float_Type     
      implicit none
      real(kind=FT),intent(in)::Tri_P1(3),Tri_P2(3),Tri_P3(3)
      real(kind=FT),intent(out):: Np(3)
      real(kind=FT) P1P2(3),P1P3(3),vectorNorm3d
      
      Np(1:3) = ZR
      
      P1P2 = Tri_P2 - Tri_P1
      P1P3 = Tri_P3 - Tri_P1
      call Vector_Cross_Product_3(P1P2,P1P3,Np)  
      call Vector_Norm2(3,Np,vectorNorm3d)   
      Np = Np/vectorNorm3d
      
      return 
      end SUBROUTINE Tool_Normal_Unit_vector_of_3D_Tri                   
