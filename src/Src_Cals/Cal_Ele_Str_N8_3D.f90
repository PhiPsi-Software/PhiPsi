!-----------------------------------------------------------
! Brief: Compute nodal stress or strain of an 8-node 3D hex element in Cartesian or cylindrical system.
!
! Parameters:
!   Input:  i_E          - element index
!           i_N          - integration/local point index
!           Key_to_cal   - 1=compute stress, 2=compute strain
!           Key_CoorSys  - 1=Cartesian, 2=cylindrical
!           X_NODES      - x-coordinates of the 8 element nodes
!           Y_NODES      - y-coordinates of the 8 element nodes
!           Z_NODES      - z-coordinates of the 8 element nodes
!           c_D          - 6x6 constitutive (material) matrix
!           kesi         - parametric coordinate of the evaluation point
!           yita         - parametric coordinate of the evaluation point
!           zeta         - parametric coordinate of the evaluation point
!           U            - element nodal displacement vector (24)
!   Output: c_Str        - computed stress or strain (6 components)
!
! Notes:   When Key_CoorSys=2 the result is rotated from Cartesian to cylindrical
!          using the element/node angle Theta_Cartesian_to_Cylinder.
!-----------------------------------------------------------

subroutine Cal_Ele_Str_N8_3D(i_E,i_N,Key_to_cal,Key_CoorSys, X_NODES, Y_NODES,Z_NODES, c_D,kesi,yita,zeta,U, c_Str)
! Used to calculate the FEM nodal stress or strain of an 8-node three-dimensional element.
! Calculate stress when Key_to_cal=1, and calculate strain when Key_to_cal=2.
! Key_CoorSys=1, Cartesian coordinate system, Key_CoorSys==2, cylindrical coordinate system
! Modified on 2021-09-10.
use Global_Float_Type
use Global_Material
use Global_Common
use Global_Material
use Global_Model
implicit none

integer,intent(in)::i_E,i_N,Key_to_cal,Key_CoorSys
real(kind=FT),intent(in):: c_D(6,6),kesi,yita,zeta
real(kind=FT),intent(in):: X_NODES(8),Y_NODES(8),Z_NODES(8),U(24)
real(kind=FT),intent(out):: c_Str(6)
real(kind=FT) detJ, J(3,3), Inverse_J(3,3)
real(kind=FT) N(3,24),dNdkesi(8,3),dNdx(8,3)      
real(kind=FT) X(8),Y(8),Z(8)
real(kind=FT) B_FEM(6,24)
real(kind=FT) c_Theta,Tem_Str(6),sTheta,cTheta,c_E,c_v
integer :: mat_num

X=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4), X_NODES(5),X_NODES(6),X_NODES(7),X_NODES(8)]
Y=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4), Y_NODES(5),Y_NODES(6),Y_NODES(7),Y_NODES(8)]
Z=[Z_NODES(1),Z_NODES(2),Z_NODES(3),Z_NODES(4), Z_NODES(5),Z_NODES(6),Z_NODES(7),Z_NODES(8)]


call Cal_N_dNdkesi_J_detJ_3D(kesi,yita,zeta, X,Y,Z,detJ,J,N,dNdkesi)


call Matrix_Inverse_3x3(J,Inverse_J)

dNdx = MATMUL(dNdkesi,Inverse_J)

B_FEM(1:6,1:24) = ZR

B_FEM(1,1:24:3)   =  dNdx(1:8,1)
B_FEM(2,2:24:3)   =  dNdx(1:8,2)
B_FEM(3,3:24:3)   =  dNdx(1:8,3)
B_FEM(4,1:24:3)   =  dNdx(1:8,2)
B_FEM(4,2:24:3)   =  dNdx(1:8,1)
B_FEM(5,2:24:3)   =  dNdx(1:8,3)
B_FEM(5,3:24:3)   =  dNdx(1:8,2)
B_FEM(6,1:24:3)   =  dNdx(1:8,3)
B_FEM(6,3:24:3)   =  dNdx(1:8,1)  

if(Key_to_cal==1)then
    if (Key_CoorSys==1) then
        c_Str = MATMUL(MATMUL(c_D,B_FEM),U)
    elseif (Key_CoorSys==2) then
        Tem_Str(1:6) = MATMUL(MATMUL(c_D,B_FEM),U)
        c_Theta     =  Theta_Cartesian_to_Cylinder(i_E,i_N)
        sTheta      =  sin(c_Theta)
        cTheta      =  cos(c_Theta)

        c_Str(1)    =  Tem_Str(1)*cTheta**2 + Tem_Str(2)*sTheta**2 + TWO*Tem_Str(4)*sTheta*cTheta
        c_Str(2)    =  Tem_Str(1)*sTheta**2 + Tem_Str(2)*cTheta**2 - TWO*Tem_Str(4)*sTheta*cTheta
        c_Str(3)    =  Tem_Str(3)   
        c_Str(4)    =  (Tem_Str(2)-Tem_Str(1))*sTheta*cTheta  + Tem_Str(4)*(cTheta**2 - sTheta**2)
        c_Str(5)    =  -Tem_Str(6)*sTheta + Tem_Str(5)*cTheta   
        c_Str(6)    =   Tem_Str(6)*cTheta + Tem_Str(5)*sTheta   
    endif             
elseif(Key_to_cal==2)then
    if (Key_CoorSys==1) then
        c_Str = MATMUL(B_FEM,U)
    elseif (Key_CoorSys==2) then
        Tem_Str(1:6) = MATMUL(MATMUL(c_D,B_FEM),U)

        mat_num = Elem_Mat(i_E)

        c_Theta     =  Theta_Cartesian_to_Cylinder(i_E,i_N)
        sTheta      =  sin(c_Theta)
        cTheta      =  cos(c_Theta)

        c_Str(1)    =  Tem_Str(1)*cTheta**2 + Tem_Str(2)*sTheta**2 + TWO*Tem_Str(4)*sTheta*cTheta
        c_Str(2)    =  Tem_Str(1)*sTheta**2 + Tem_Str(2)*cTheta**2 - TWO*Tem_Str(4)*sTheta*cTheta
        c_Str(3)    =  Tem_Str(3)   
        c_Str(4)    =  (Tem_Str(2)-Tem_Str(1))*sTheta*cTheta  + Tem_Str(4)*(cTheta**2 - sTheta**2)
        c_Str(5)    =  -Tem_Str(6)*sTheta + Tem_Str(5)*cTheta   
        c_Str(6)    =   Tem_Str(6)*cTheta + Tem_Str(5)*sTheta   

        c_Str(1:6) = MATMUL(D_for_cylindrical(mat_num,1:6,1:6), c_Str)
    endif     

endif


return 
end subroutine Cal_Ele_Str_N8_3D           
