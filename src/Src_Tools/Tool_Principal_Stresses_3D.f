 
      subroutine Tool_Principal_Stresses_3D(S_xx,S_yy,S_zz,
     &                            S_yz,S_xz,S_xy,
     &                            S_1,S_2,S_3,
     &                            Vector_S1,Vector_S2,Vector_S3)
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::S_xx,S_yy,S_zz,S_xy,S_yz,S_xz
      real(kind=FT),intent(out)::S_1,S_2,S_3
      real(kind=FT),intent(out)::Vector_S1(3),Vector_S2(3),Vector_S3(3)
      real(kind=FT) S_m
      real(kind=FT) pian_Stress(3,3)
      real(kind=FT) pian_J1,pian_J2,pian_J3
      real(kind=FT) omega_sigma
      real(kind=FT) tem_1,tem_2
      real(kind=FT) tem_3,tem_4,tem_5,tem_6
      real(kind=FT) A,B,C
      
      S_1=ZR
      S_2=ZR
      S_3=ZR
      Vector_S1(1:3)=ZR
      Vector_S2(1:3)=ZR
      Vector_S3(1:3)=ZR
      S_m = (S_xx+S_yy+S_zz)/THR
      pian_Stress(1,1) = S_xx - S_m
      pian_Stress(1,2) = S_xy
      pian_Stress(1,3) = S_xz
      pian_Stress(2,1) = S_xy
      pian_Stress(2,2) = S_yy - S_m
      pian_Stress(2,3) = S_yz
      pian_Stress(3,1) = S_xz
      pian_Stress(3,2) = S_yz
      pian_Stress(3,3) = S_zz - S_m
      pian_J2= -((S_xx - S_m)*(S_yy - S_m) + 
     &           (S_yy - S_m)*(S_zz - S_m) +
     &           (S_xx - S_m)*(S_zz - S_m)) +
     &            S_xy**2 + S_yz**2 + S_xz**2  
      
      if(pian_J2<Tol_10)then
          pian_J2 = Tol_10
      endif
      
      call Matrix_Det_3x3(pian_Stress,pian_J3) 
      tem_1 = -THR*SQRT(THR)/TWO*(pian_J3/pian_J2/sqrt(pian_J2))
      omega_sigma =  ONE/THR*ACOS(tem_1)
      
      tem_2 = TWO/SQRT(THR)
      S_1 =  tem_2*sqrt(pian_J2)*cos(omega_sigma - pi/THR) + S_m
      S_2 =  tem_2*sqrt(pian_J2)*cos(omega_sigma + pi/THR) + S_m
      S_3 = -tem_2*sqrt(pian_J2)*cos(omega_sigma) + S_m
      
      A = S_xy*S_yz-(S_yy-S_1)*S_xz
      B = S_xy*S_xz-(S_xx-S_1)*S_yz
      C = (S_xx-S_1)*(S_yy-S_1)-S_xy**2
      tem_6 = sqrt(A**2 + B**2 + C**2)
      Vector_S1(1) = A/tem_6
      Vector_S1(2) = B/tem_6
      Vector_S1(3) = C/tem_6
      
      A = S_xy*S_yz-(S_yy-S_2)*S_xz
      B = S_xy*S_xz-(S_xx-S_2)*S_yz
      C = (S_xx-S_2)*(S_yy-S_2)-S_xy**2
      tem_6 = sqrt(A**2+B**2+C**2)
      Vector_S2(1) = A/tem_6
      Vector_S2(2) = B/tem_6
      Vector_S2(3) = C/tem_6
      
      A = S_xy*S_yz-(S_yy-S_3)*S_xz
      B = S_xy*S_xz-(S_xx-S_3)*S_yz
      C = (S_xx-S_3)*(S_yy-S_3)-S_xy**2
      tem_6 = sqrt(A**2+B**2+C**2)
      Vector_S3(1) = A/tem_6
      Vector_S3(2) = B/tem_6
      Vector_S3(3) = C/tem_6  

      RETURN
      END SUBROUTINE Tool_Principal_Stresses_3D