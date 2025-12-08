 
      subroutine Cal_N(kesi,yita,N)
     
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita
      real(kind=FT),intent(out):: N(2,8)
      real(kind=FT) N1,N2,N3,N4
      
      
        N1 = (ONE-kesi)*(ONE-yita)/FOU
        N2 = (ONE+kesi)*(ONE-yita)/FOU
        N3 = (ONE+kesi)*(ONE+yita)/FOU
        N4 = (ONE-kesi)*(ONE+yita)/FOU
        N(1,:) = [N1,ZR,N2,ZR,N3,ZR,N4,ZR]
        N(2,:) = [ZR,N1,ZR,N2,ZR,N3,ZR,N4]


      
      return 
      end SUBROUTINE Cal_N               
