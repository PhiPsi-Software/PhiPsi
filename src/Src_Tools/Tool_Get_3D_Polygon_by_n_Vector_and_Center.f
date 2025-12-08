 
      subroutine Tool_Get_3D_Polygon_by_n_Vector_and_Center(n_Vector,
     &                        c_Num_Poly_Edges_NaCr,
     &                        c_Crack_Center,Crack_L,
     &                        Out_Coors)

      use Global_Float_Type 
      implicit none
      integer,intent(in) :: c_Num_Poly_Edges_NaCr
      real(kind=FT),intent(in) :: n_Vector(3),c_Crack_Center(3),Crack_L
      real(kind=FT),intent(out) :: Out_Coors(c_Num_Poly_Edges_NaCr,3)
      
      integer num_Circle_Divison
      integer i_theta
      real(kind=FT) theta,a(3),b(3),norm_a,norm_b
      real(kind=FT) c_xyz(1000,3),c_Radius
      
      
      
      
      c_Radius = Crack_L/TWO
      
      num_Circle_Divison = c_Num_Poly_Edges_NaCr
      
      c_xyz(1:1000,1:3)   = ZR
      
      do i_theta = 1,num_Circle_Divison
           theta = real(i_theta-1)*TWO*pi/real(num_Circle_Divison) -
     &        pi/FOU
     
           call Vector_Cross_Product_3(n_Vector,[ONE,ZR,ZR],a)   
           if(sum(abs(a))<=Tol_20) then
              call Vector_Cross_Product_3(n_Vector,[ZR,ONE,ZR],a)   
           endif
           call Vector_Cross_Product_3(n_Vector,a,b)   
           call Vector_Norm2(3,a,norm_a)   
           call Vector_Norm2(3,b,norm_b)  
           a=a/norm_a
           b=b/norm_b
 
           c_xyz(i_theta,1)=c_Crack_Center(1)+
     &                 c_Radius*a(1)*cos(theta)
     &                +c_Radius*b(1)*sin(theta)
           c_xyz(i_theta,2)=c_Crack_Center(2)+
     &                 c_Radius*a(2)*cos(theta)
     &                +c_Radius*b(2)*sin(theta)
           c_xyz(i_theta,3)=c_Crack_Center(3)+
     &                 c_Radius*a(3)*cos(theta)
     &                +c_Radius*b(3)*sin(theta)
     
           Out_Coors(i_theta,1:3) =  c_xyz(i_theta,1:3)
      enddo

      
      return 
      end SUBROUTINE Tool_Get_3D_Polygon_by_n_Vector_and_Center                      
