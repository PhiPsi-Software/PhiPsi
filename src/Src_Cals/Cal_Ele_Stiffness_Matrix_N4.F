 
      subroutine Cal_Ele_Stiffness_Matrix_N4(X_NODES,Y_NODES,
     &                           thick,c_D,kesi,yita,weight,
     &                           localK)
      use Global_Float_Type
      use Global_Common
      implicit none
      
      real(kind=FT),intent(in)::thick,c_D(3,3),
     &                     kesi(Num_Gauss_P_FEM),yita(Num_Gauss_P_FEM),
     &                     weight(Num_Gauss_P_FEM)
      real(kind=FT),intent(in):: X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out):: localK(8,8)
      
      integer i
      real(kind=FT) X(4),Y(4),JM(4,4),a,b,c,d,
     &                Nkesi(4),Nyita(4),detJ,
     &                B1(3,2),B2(3,2),B3(3,2),B4(3,2),ToTal_B(3,8)
#ifdef Silverfrost
      real(kind=FT) :: tmp1(8,3), tmp2(8,8)
#endif
      localK(1:8,1:8)=ZR

      X=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Y=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]

      do i=1,Num_Gauss_P_FEM
            JM(1,:)=[ZR,ONE-yita(i),yita(i)-kesi(i),  kesi(i)-ONE]
            JM(2,:)=[yita(i)-ONE,ZR,ONE+kesi(i), -kesi(i)-yita(i)]
            JM(3,:)=[kesi(i)-yita(i), -kesi(i)-ONE,ZR,yita(i)+ONE]
            JM(4,:)=[ONE-kesi(i), yita(i)+kesi(i),-yita(i)-ONE,ZR]
      
            Nkesi=[(yita(i)-ONE)/FOU,(-yita(i)+ONE)/FOU,
     &           (yita(i)+ONE)/FOU,(-yita(i)-ONE)/FOU]
            Nyita=[(kesi(i)-ONE)/FOU,(-kesi(i)-ONE)/FOU,
     &           (kesi(i)+ONE)/FOU,(-kesi(i)+ONE)/FOU]
      
            a=(Y(1)*(kesi(i)-ONE)+Y(2)*(-ONE-kesi(i))+
     *                     Y(3)*(ONE+kesi(i))+Y(4)*( ONE-kesi(i)))/FOU
            b=(Y(1)*(yita(i)-ONE)+Y(2)*(ONE-yita(i))+ 
     *                     Y(3)*(ONE+yita(i))+Y(4)*(-ONE-yita(i)))/FOU
            c=(X(1)*(yita(i)-ONE)+X(2)*(ONE-yita(i))+ 
     *                     X(3)*(ONE+yita(i))+X(4)*(-ONE-yita(i)))/FOU
            d=(X(1)*(kesi(i)-ONE)+X(2)*(-ONE-kesi(i))+
     *                     X(3)*(ONE+kesi(i))+X(4)*( ONE-kesi(i)))/FOU

            detJ  =   dot_product(MATMUL(X,JM),Y)/EIG

            B1(1,:)=[a*Nkesi(1)-b*Nyita(1),ZR ]
            B1(2,:)=[ZR, c*Nyita(1)-d*Nkesi(1)]
            B1(3,:)=[c*Nyita(1)-d*Nkesi(1), a*Nkesi(1)-b*Nyita(1)]
            B2(1,:)=[a*Nkesi(2)-b*Nyita(2),ZR ]
            B2(2,:)=[ZR ,c*Nyita(2)-d*Nkesi(2)]
            B2(3,:)=[c*Nyita(2)-d*Nkesi(2), a*Nkesi(2)-b*Nyita(2)]
            B3(1,:)=[a*Nkesi(3)-b*Nyita(3),ZR]
            B3(2,:)=[ZR , c*Nyita(3)-d*Nkesi(3)]
            B3(3,:)=[c*Nyita(3)-d*Nkesi(3), a*Nkesi(3)-b*Nyita(3)]
            B4(1,:)=[a*Nkesi(4)-b*Nyita(4),ZR]
            B4(2,:)=[ZR ,c*Nyita(4)-d*Nkesi(4)]
            B4(3,:)=[c*Nyita(4)-d*Nkesi(4),a*Nkesi(4)-b*Nyita(4)]
         
            ToTal_B(1:3,1:2)= B1
            ToTal_B(1:3,3:4)= B2
            ToTal_B(1:3,5:6)= B3
            ToTal_B(1:3,7:8)= B4
            ToTal_B=ToTal_B/detJ
            
#ifndef Silverfrost
            localK = localK + thick*weight(i)*detJ*
     &             MATMUL(MATMUL(transpose(ToTal_B),c_D),ToTal_B)
#endif
#ifdef Silverfrost
            tmp1 = MATMUL( transpose(ToTal_B), c_D ) 
            tmp2 = MATMUL( tmp1, ToTal_B ) 
            localK = localK + thick*weight(i)*detJ*tmp2
#endif      
      end do
      
      return 
      end SUBROUTINE Cal_Ele_Stiffness_Matrix_N4               
