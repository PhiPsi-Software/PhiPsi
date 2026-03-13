!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
subroutine Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_Cubes,Num_Gau_Points_Cube,&
                                            kesi_Cubes,yita_Cubes,zeta_Cubes,weight_Cubes)       
! Integration Scheme 3: Local Coordinates and Weight Coefficients of Gaussian Points for
! Quadrilateral Subdivision Rules
! Theoretical references for 2D problems: Extended Finite Element Method: Theory and
! Applications_2015 XFEM new book_Page 58_Figure 2.14
! 3D problems refer to 2D.

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
