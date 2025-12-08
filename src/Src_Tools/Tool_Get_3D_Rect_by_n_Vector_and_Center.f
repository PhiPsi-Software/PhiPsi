 
      subroutine Tool_Get_3D_Rect_by_n_Vector_and_Center(n_Vector,
     &                        c_Crack_Center,Crack_L,
     &                        Out_Coors)

      use Global_Float_Type 
      implicit none
      real(kind=FT),intent(in) :: n_Vector(3),c_Crack_Center(3),Crack_L
      real(kind=FT),intent(out) :: Out_Coors(4,3)
      integer num_Circle_Divison
      integer i_theta
      real(kind=FT) theta,a(3),b(3),norm_a,norm_b
      real(kind=FT) c_xyz(4,3),c_Radius
      
      c_Radius = sqrt(TWO)/TWO*Crack_L
      
       call Vector_Cross_Product_3(n_Vector,[ONE,ZR,ZR],a)   
       if(sum(abs(a))<=Tol_20) then
          call Vector_Cross_Product_3(n_Vector,[ZR,ONE,ZR],a)   
       endif
       call Vector_Cross_Product_3(n_Vector,a,b)   
       call Vector_Norm2(3,a,norm_a)   
       call Vector_Norm2(3,b,norm_b)  
       a=a/norm_a
       b=b/norm_b
       
           
      num_Circle_Divison = 4
      do i_theta = 1,num_Circle_Divison
           theta = (i_theta-1)*TWO*pi/num_Circle_Divison - pi/FOU
           c_xyz(i_theta,1)=c_Crack_Center(1)+
     &                 c_Radius*a(1)*cos(theta)
     &                +c_Radius*b(1)*sin(theta)
           c_xyz(i_theta,2)=c_Crack_Center(2)+
     &                 c_Radius*a(2)*cos(theta)
     &                +c_Radius*b(2)*sin(theta)
           c_xyz(i_theta,3)=c_Crack_Center(3)+
     &                 c_Radius*a(3)*cos(theta)
     &                +c_Radius*b(3)*sin(theta)
      enddo
      

      Out_Coors(1,1:3) =  c_xyz(1,1:3)
      Out_Coors(2,1:3) =  c_xyz(2,1:3)
      Out_Coors(3,1:3) =  c_xyz(3,1:3)
      Out_Coors(4,1:3) =  c_xyz(4,1:3)
      
      return 
      end SUBROUTINE Tool_Get_3D_Rect_by_n_Vector_and_Center                      
