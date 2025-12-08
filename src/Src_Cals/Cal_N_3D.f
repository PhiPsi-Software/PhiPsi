 
      subroutine Cal_N_3D(kesi,yita,zeta,N)
     
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,zeta
      real(kind=FT),intent(out):: N(3,24)
      
      real(kind=FT) N1,N2,N3,N4,N5,N6,N7,N8
      real(kind=FT) o
      
      o      = ZR
      
      N1 = (ONE-kesi)*(ONE-yita)*(ONE-zeta)/EIG
      N2 = (ONE+kesi)*(ONE-yita)*(ONE-zeta)/EIG
      N3 = (ONE+kesi)*(ONE+yita)*(ONE-zeta)/EIG
      N4 = (ONE-kesi)*(ONE+yita)*(ONE-zeta)/EIG
      N5 = (ONE-kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N6 = (ONE+kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N7 = (ONE+kesi)*(ONE+yita)*(ONE+zeta)/EIG
      N8 = (ONE-kesi)*(ONE+yita)*(ONE+zeta)/EIG
      
      N(1,1:24) =
     &        [N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o,o]
      N(2,1:24) =
     &        [o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o]
      N(3,1:24) =
     &        [o,o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8]
      
      return 
      end SUBROUTINE Cal_N_3D               
