 
subroutine Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_Cubes,Num_Gau_Points_Cube,&
                                            kesi_Cubes,yita_Cubes,zeta_Cubes,weight_Cubes)       

use Global_Float_Type      
implicit none
integer,intent(in)::Num_Sub_Cubes,Num_Gau_Points_Cube
real(kind=FT),intent(out)::kesi_Cubes(Num_Sub_Cubes*Num_Gau_Points_Cube),&
                           yita_Cubes(Num_Sub_Cubes*Num_Gau_Points_Cube),&
                           zeta_Cubes(Num_Sub_Cubes*Num_Gau_Points_Cube),&
                           weight_Cubes(Num_Sub_Cubes*Num_Gau_Points_Cube)
integer i,j,k
real(kind=FT) kesi_CUBE(Num_Gau_Points_Cube),yita_CUBE(Num_Gau_Points_Cube),&
              zeta_CUBE(Num_Gau_Points_Cube),weight_CUBE(Num_Gau_Points_Cube)
real(kind=FT) coord(8,3)
real(kind=FT) N1,N2,N3,N4,N5,N6,N7,N8
real(kind=FT) c_kesi,c_yita,c_zeta
real(kind=FT) NN(8)
integer pt
integer i_Gauss
integer Total_Gauss_Num
integer num_Gauss_Cube_1D
real(kind=FT) Cube_Volume
integer num_Box_Edge
real(kind=FT) Length_Box
real(kind=FT) tem(3)


Total_Gauss_Num = Num_Sub_Cubes*Num_Gau_Points_Cube
kesi_Cubes(1:Total_Gauss_Num)  =  ZR
yita_Cubes(1:Total_Gauss_Num)  =  ZR
zeta_Cubes(1:Total_Gauss_Num)  =  ZR
weight_Cubes(1:Total_Gauss_Num)=  ZR

num_Box_Edge = CEILING((dble(Num_Sub_Cubes))**(ONE/THR))
Length_Box   = ONE/dble(num_Box_Edge)
num_Gauss_Cube_1D =  int((dble(Num_Gau_Points_Cube))**(ONE/THR))

Cube_Volume  = ONE/Num_Sub_Cubes


call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_Cube,kesi_CUBE,yita_CUBE,zeta_CUBE,weight_CUBE)      


pt = 0
do i=1, num_Box_Edge
  do j=1, num_Box_Edge
    do k=1, num_Box_Edge
      coord(1,1) =  ((i-1)*Length_Box)*TWO - ONE
      coord(4,1) =  ((i-1)*Length_Box)*TWO - ONE
      coord(5,1) =  ((i-1)*Length_Box)*TWO - ONE  
      coord(8,1) =  ((i-1)*Length_Box)*TWO - ONE      
      
      coord(2,1) =  (i*Length_Box)*TWO     - ONE
      coord(3,1) =  (i*Length_Box)*TWO     - ONE
      coord(6,1) =  (i*Length_Box)*TWO     - ONE
      coord(7,1) =  (i*Length_Box)*TWO     - ONE

      coord(1,2) =  ((j-1)*Length_Box)*TWO - ONE
      coord(2,2) =  ((j-1)*Length_Box)*TWO - ONE
      coord(5,2) =  ((j-1)*Length_Box)*TWO - ONE
      coord(6,2) =  ((j-1)*Length_Box)*TWO - ONE
      
      coord(3,2) =  (j*Length_Box)*TWO     - ONE
      coord(4,2) =  (j*Length_Box)*TWO     - ONE
      coord(7,2) =  (j*Length_Box)*TWO     - ONE
      coord(8,2) =  (j*Length_Box)*TWO     - ONE

      coord(1,3) =  ((k-1)*Length_Box)*TWO - ONE
      coord(2,3) =  ((k-1)*Length_Box)*TWO - ONE
      coord(3,3) =  ((k-1)*Length_Box)*TWO - ONE
      coord(4,3) =  ((k-1)*Length_Box)*TWO - ONE
      
      coord(5,3) =  (k*Length_Box)*TWO     - ONE
      coord(6,3) =  (k*Length_Box)*TWO     - ONE
      coord(7,3) =  (k*Length_Box)*TWO     - ONE
      coord(8,3) =  (k*Length_Box)*TWO     - ONE
      
      do i_Gauss =1,Num_Gau_Points_Cube
          c_kesi = kesi_CUBE(i_Gauss)
          c_yita = yita_CUBE(i_Gauss)
          c_zeta = zeta_CUBE(i_Gauss)
          N1 = (ONE-c_kesi)*(ONE-c_yita)*(ONE-c_zeta)/EIG
          N2 = (ONE+c_kesi)*(ONE-c_yita)*(ONE-c_zeta)/EIG
          N3 = (ONE+c_kesi)*(ONE+c_yita)*(ONE-c_zeta)/EIG
          N4 = (ONE-c_kesi)*(ONE+c_yita)*(ONE-c_zeta)/EIG
          N5 = (ONE-c_kesi)*(ONE-c_yita)*(ONE+c_zeta)/EIG
          N6 = (ONE+c_kesi)*(ONE-c_yita)*(ONE+c_zeta)/EIG
          N7 = (ONE+c_kesi)*(ONE+c_yita)*(ONE+c_zeta)/EIG
          N8 = (ONE-c_kesi)*(ONE+c_yita)*(ONE+c_zeta)/EIG
          
          NN(1:8)   = [N1,N2,N3,N4,N5,N6,N7,N8]

          tem(1:3) = MATMUL(NN(1:8),coord(1:8,1:3))   
                    
          pt = pt +1

          kesi_Cubes(pt)   = tem(1)
          yita_Cubes(pt)   = tem(2)
          zeta_Cubes(pt)   = tem(3)
          
          weight_Cubes(pt) = weight_CUBE(i_Gauss)*Cube_Volume
      enddo
    enddo
  enddo
enddo



return 
end SUBROUTINE Cal_Gauss_Points_3D_for_SUBCUBES                       
