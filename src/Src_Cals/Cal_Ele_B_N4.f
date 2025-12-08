 
      subroutine Cal_Ele_B_N4(i_E,X_NODES,Y_NODES,
     &                             kesi,yita,ToTal_B) 
      
      use Global_Float_Type
      implicit none
      
      integer,intent(in)::i_E
      real(kind=FT),intent(in):: kesi,yita
      real(kind=FT),intent(in):: X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out):: ToTal_B(3,8)
      
      real(kind=FT) X(4),Y(4),JM(4,4),a,b,c,d,
     &                Nkesi(4),Nyita(4),detJ,
     &                B1(3,2),B2(3,2),B3(3,2),B4(3,2)
      
      
      X=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Y=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]

      JM(1,:)=[ZR,     ONE-yita,  yita-kesi,  kesi-ONE ]
      JM(2,:)=[yita-ONE,  ZR,     ONE+kesi,  -kesi-yita]
      JM(3,:)=[kesi-yita, -kesi-ONE, ZR,      yita+ONE ]
      JM(4,:)=[ONE-kesi,  yita+kesi, -yita-ONE,  ZR    ]

      Nkesi=[(yita-ONE)/FOU,(-yita+ONE)/FOU,
     &       (yita+ONE)/FOU,(-yita-ONE)/FOU]
      Nyita=[(kesi-ONE)/FOU,(-kesi-ONE)/FOU,
     &       (kesi+ONE)/FOU,(-kesi+ONE)/FOU]

      a=(Y(1)*(kesi-ONE)+Y(2)*(-ONE-kesi)+
     *   Y(3)*(ONE+kesi)+Y(4)*( ONE-kesi))/FOU
      b=(Y(1)*(yita-ONE)+Y(2)*(ONE-yita)+ 
     *   Y(3)*(ONE+yita)+Y(4)*(-ONE-yita))/FOU
      c=(X(1)*(yita-ONE)+X(2)*(ONE-yita)+ 
     *   X(3)*(ONE+yita)+X(4)*(-ONE-yita))/FOU
      d=(X(1)*(kesi-ONE)+X(2)*(-ONE-kesi)+
     *   X(3)*(ONE+kesi)+X(4)*( ONE-kesi))/FOU

      detJ  =   dot_product(MATMUL(X,JM),Y)/EIG

      B1(1,:)=[a*Nkesi(1)-b*Nyita(1),                ZR ]
      B1(2,:)=[ZR,                 c*Nyita(1)-d*Nkesi(1)]
      B1(3,:)=[c*Nyita(1)-d*Nkesi(1), a*Nkesi(1)-b*Nyita(1)]
      B2(1,:)=[a*Nkesi(2)-b*Nyita(2),                ZR ]
      B2(2,:)=[ZR ,                c*Nyita(2)-d*Nkesi(2)]
      B2(3,:)=[c*Nyita(2)-d*Nkesi(2), a*Nkesi(2)-b*Nyita(2)]
      B3(1,:)=[a*Nkesi(3)-b*Nyita(3),                 ZR]
      B3(2,:)=[ZR ,                c*Nyita(3)-d*Nkesi(3)]
      B3(3,:)=[c*Nyita(3)-d*Nkesi(3), a*Nkesi(3)-b*Nyita(3)]
      B4(1,:)=[a*Nkesi(4)-b*Nyita(4),                 ZR]
      B4(2,:)=[ZR ,                c*Nyita(4)-d*Nkesi(4)]
      B4(3,:)=[c*Nyita(4)-d*Nkesi(4), a*Nkesi(4)-b*Nyita(4)]

      ToTal_B(1:3,1:2)= B1
      ToTal_B(1:3,3:4)= B2
      ToTal_B(1:3,5:6)= B3
      ToTal_B(1:3,7:8)= B4
      
      ToTal_B = ToTal_B /detJ    
      
      return 
      end SUBROUTINE Cal_Ele_B_N4           
