 
      subroutine Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_QUAD,
     &                                             kesi,yita,weight)

      use Global_Float_Type      
      implicit none
      integer,intent(in)::Num_QUAD
      real(kind=FT),intent(out)::kesi(Num_QUAD*4),
     &                           yita(Num_QUAD*4),
     &                           weight(Num_QUAD*4)
      integer i,j
      real(kind=FT) np(30),nw(30)
      
      
      real(kind=FT) np_QUAD(2),nw_QUAD(2)
      real(kind=FT) kesi_QUAD(4),yita_QUAD(4),weight_QUAD(4)
      integer pt
      integer i_Gauss
      real(kind=FT) NN(4)
      real(kind=FT) coord(4,2)
      integer num_Box_Edge
      real(kind=FT) Length_Box
      real(kind=FT) N1,N2,N3,N4
      real(kind=FT) c_kesi,c_yita
      real(kind=FT) tem(2)
      real(kind=FT) Box_Area
      np(1:20) = ZR
      nw(1:20) = ZR
      
      kesi(1:Num_QUAD*4)   = ZR
      yita(1:Num_QUAD*4)   = ZR
      weight(1:Num_QUAD*4) = ZR
      
      num_Box_Edge = int(sqrt(dble(Num_QUAD)))
      
      Length_Box = ONE/dble(num_Box_Edge)
      
      np_QUAD(1) = -0.577350269189626D0
      np_QUAD(2) =  0.577350269189626D0
      nw_QUAD(1) =  ONE
      nw_QUAD(2) =  ONE
      do i = 1,2
          do j = 1,2
              kesi_QUAD((i-1)*2+j)   = np_QUAD(i)
              yita_QUAD((i-1)*2+j)   = np_QUAD(j)
              weight_QUAD((i-1)*2+j) = nw_QUAD(i)*nw_QUAD(j)
          end do
      end do   

      Box_Area = ONE/Num_QUAD
      
      
      pt = 0
      do i=1,int(sqrt(dble(Num_QUAD)))
          do j=1,int(sqrt(dble(Num_QUAD)))
              coord(1,1) =  ((i-1)*Length_Box)*TWO - ONE
              coord(2,1) =  (i*Length_Box)*TWO     - ONE
              coord(3,1) =  (i*Length_Box)*TWO     - ONE
              coord(4,1) =  ((i-1)*Length_Box)*TWO - ONE
              coord(1,2) =  ((j-1)*Length_Box)*TWO - ONE
              coord(2,2) =  ((j-1)*Length_Box)*TWO - ONE
              coord(3,2) =  (j*Length_Box)*TWO     - ONE
              coord(4,2) =  (j*Length_Box)*TWO     - ONE
              
              do i_Gauss =1,4
                  c_kesi = kesi_QUAD(i_Gauss)
                  c_yita = yita_QUAD(i_Gauss)
                  N1 = (ONE-c_kesi)*(ONE-c_yita)/FOU
                  N2 = (ONE+c_kesi)*(ONE-c_yita)/FOU
                  N3 = (ONE+c_kesi)*(ONE+c_yita)/FOU
                  N4 = (ONE-c_kesi)*(ONE+c_yita)/FOU
                  NN(1:4)   = [N1,N2,N3,N4]

                    tem(1:2) = MATMUL(NN(1:4),coord(1:4,1:2))   
                            
                  pt = pt +1
 
                  kesi(pt)   = tem(1)
                  yita(pt)   = tem(2)
                  
                  weight(pt) = weight_QUAD(i_Gauss)*Box_Area
              enddo
          enddo
      enddo
      
      return 
      end SUBROUTINE Cal_Gauss_Points_QUAD_for_SUBQUAD                        
