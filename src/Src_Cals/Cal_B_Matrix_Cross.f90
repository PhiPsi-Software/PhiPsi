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
 
subroutine Cal_B_Matrix_Cross(kesi,yita,i_Cross,i_E,i_G,c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)  

! Calculate the B matrix (for cross).
! References: Belytschko_2001_Arbitrary discontinuities in finite elements_Equation-12
!..........................................
! Variables related to cross-shaped cracks
!..........................................
! Cross_Point_Cr_num(num_Cross,1:2) !The primary and secondary crack numbers corresponding to each
! cross intersection point
! Cross_Point_Ele_num(num_Cross,1:2) !The element numbers corresponding to each cross intersection
! point
! Cross_Point_RABCD(num_Cross,1:10,1:2) !Coordinates of intersection and A, B, C, D points
! corresponding to each cross intersection
! Where R is the intersection point
! For details about this variable, see the notes on V5-PhiPsi cross-shaped fracture handling.
! Elem_Type_Cross(Num_Elem,num_Cross);          ! Element type: 1--cross crack reinforced element
! Enriched_Node_Type_Cross(Num_Node, num_Cross); ! Type of enriched nodes: 1 -- Cross-shaped crack
! enriched node

use Global_Float_Type
use Global_Crack
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_Cross

implicit none

integer,intent(in)::i_Cross,i_E,i_G
real(kind=FT),intent(in)::c_X_NODES(4),c_Y_NODES(4)
integer,intent(in)::c_NN(4)
real(kind=FT),intent(in)::kesi,yita
real(kind=FT),intent(out)::tem_B(3,80)
integer,intent(out)::num_tem_B        
real(kind=FT) detJ, J(2,2), Inverse_J(2,2)
real(kind=FT) N(2,8),dNdkesi(4,2),dNdx(4,2)
integer mat_num,c_mat_type
real(kind=FT) B_FEM(3,8),B_XFEM(3,60),BI_enr(3,2)
integer num_B_XFEM
integer i_N
real(kind=FT) c_G_x,c_G_y
real(kind=FT) Cross_gp,Cross_i
integer c_Cross_elem
real(kind=FT) Coor_AB(2,2),Coor_CD(2,2)
real(kind=FT) distance_Node_M,distance_Gauss_M
real(kind=FT) distance_Node_sm,distance_Gauss_sm
real(kind=FT) H_i_M,H_i_sm

real(kind=FT) H_x0

call Cal_N_dNdkesi_J_detJ(kesi,yita,c_X_NODES,c_Y_NODES,detJ,J,N,dNdkesi)    

call Matrix_Inverse_2x2(J,Inverse_J)

dNdx = MATMUL(dNdkesi,Inverse_J)

c_G_x = DOT_PRODUCT(N(1,1:7:2),c_X_NODES(1:4))
c_G_y = DOT_PRODUCT(N(1,1:7:2),c_Y_NODES(1:4))      
 
mat_num    = Elem_Mat(i_E)

c_mat_type = Material_Type(mat_num)     

B_FEM(1:3,1:8) = ZR

if (EleGaus_yes_FEM_asemd(i_E,i_G) .eqv. .False.) then
  B_FEM(1,1:8:2)   =  dNdx(:,1)
  B_FEM(2,2:8:2)   =  dNdx(:,2)
  B_FEM(3,1:8:2)   =  dNdx(:,2)
  B_FEM(3,2:8:2)   =  dNdx(:,1)
  EleGaus_yes_FEM_asemd(i_E,i_G)= .True.
end if

if(maxval(Enriched_Node_Type_Cross(c_NN,i_Cross)).eq.0 .and. &
 minval(Enriched_Node_Type_Cross(c_NN,i_Cross)).eq.0) then
  if (i_Cross.eq.1) then
      tem_B(1:3,1:8) = B_FEM
      num_tem_B = 8
  else
      num_tem_B = 0
  end if
elseif(maxval(Enriched_Node_Type_Cross(c_NN,i_Cross)).gt.0 .or. &
     minval(Enriched_Node_Type_Cross(c_NN,i_Cross)).gt.0) then 
  B_XFEM(1:3,1:60) = ZR
  num_B_XFEM = 0
  do i_N = 1,4
      if (Enriched_Node_Type_Cross(c_NN(i_N),i_Cross).eq.1)then     
          if (Elem_Type_Cross(i_E,i_Cross).eq.1) then
          
              c_Cross_elem = Node_Cross_elem(c_NN(i_N),i_Cross)
              Coor_AB(1,1:2)=Cross_Point_RABCD(i_Cross,2,1:2)
              Coor_AB(2,1:2)=Cross_Point_RABCD(i_Cross,3,1:2)
              Coor_CD(1,1:2)=Cross_Point_RABCD(i_Cross,4,1:2)
              Coor_CD(2,1:2)=Cross_Point_RABCD(i_Cross,5,1:2)

              call Cal_Signed_Distance(Coor_AB,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_M)    
              call Cal_Signed_Distance(Coor_AB,[c_G_x,c_G_y],distance_Gauss_M)                      
              call Cal_Signed_Distance(Coor_CD,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_sm)    
              call Cal_Signed_Distance(Coor_CD,[c_G_x,c_G_y],distance_Gauss_sm) 
              if(Key_Heaviside_Value==-1)then
                  call Cal_Sign(distance_Gauss_M*distance_Gauss_sm,Cross_gp)
                  call Cal_Sign(distance_Node_M*distance_Node_sm,Cross_i)
              elseif(Key_Heaviside_Value==0)then
                  call Cal_Sign_1_and_0(distance_Gauss_M*distance_Gauss_sm,Cross_gp)
                  call Cal_Sign_1_and_0(distance_Node_M*distance_Node_sm,Cross_i)
               endif

              if (H_i_M*H_x0 > ZR) then
                  if (H_i_sm >= ZR)then
                      Cross_i=  ONE
                  else
                      Cross_i= -ONE
                  end if
              else
                  Cross_i=  ZR
              end if      
          
              BI_enr(1,:) = [dNdx(i_N,1)*(Cross_gp-Cross_i), ZR]
              BI_enr(2,:) = [ZR,dNdx(i_N,2)*(Cross_gp-Cross_i)]
              BI_enr(3,:) = [dNdx(i_N,2)*(Cross_gp-Cross_i), dNdx(i_N,1)*(Cross_gp-Cross_i)]
          else
              BI_enr(1:3,1:2)=ZR
          end if
          B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
          num_B_XFEM = num_B_XFEM + 2                  
      end if
  end do
  if (i_Cross.eq.1) then
      tem_B(1:3,1:8)            = B_FEM
      tem_B(1:3,9:8+num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
      num_tem_B = 8 + num_B_XFEM
  else
      tem_B(1:3,1:num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
      num_tem_B = num_B_XFEM
  end if
end if   

RETURN
END SUBROUTINE Cal_B_Matrix_Cross