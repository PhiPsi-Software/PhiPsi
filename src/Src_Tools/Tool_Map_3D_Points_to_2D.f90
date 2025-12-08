 
subroutine Tool_Map_3D_Points_to_2D(num_Points,X_3D,Y_3D,Z_3D,X_2D,Y_2D)
use Global_Float_Type     
implicit none
integer,intent(in):: num_Points
real(kind=FT),intent(in)::X_3D(num_Points),Y_3D(num_Points),Z_3D(num_Points)
real(kind=FT),intent(out)::X_2D(num_Points),Y_2D(num_Points)

real(kind=FT) center(3),normal_vector(3),Vector_U(3),Vector_V(3)
real(kind=FT) divisors(3),divisor
integer i_Point,divisor_type

center = [sum(X_3D)/dble(num_Points),sum(Y_3D)/dble(num_Points),sum(Z_3D)/dble(num_Points)]


call Tool_Normal_vector_of_3D_Tri([X_3D(1),Y_3D(1),Z_3D(1)],&
                                [X_3D(2),Y_3D(2),Z_3D(2)],&
                                [X_3D(3),Y_3D(3),Z_3D(3)],normal_vector)
call Vector_Normalize(3,normal_vector)   


Vector_U = [X_3D(2),Y_3D(2),Z_3D(2)]-[X_3D(1),Y_3D(1),Z_3D(1)]
call Vector_Normalize(3,Vector_U) 

call Vector_Cross_Product_3(normal_vector,Vector_U,Vector_V)   
call Vector_Normalize(3,Vector_V)  

divisors(1) =  Vector_U(1)*Vector_V(2)-Vector_V(1)*Vector_U(2) 
divisors(2) =  Vector_U(1)*Vector_V(3)-Vector_V(1)*Vector_U(3) 
divisors(3) =  Vector_U(2)*Vector_V(3)-Vector_V(2)*Vector_U(3) 
divisor_type = maxloc(abs(divisors),1)
divisor      = divisors(divisor_type)
do i_Point =1,num_Points
  if(divisor_type==1) then
      X_2D(i_Point) = ((X_3D(i_Point) - center(1))*Vector_V(2) - (Y_3D(i_Point) - center(2))*Vector_V(1))/divisor
      Y_2D(i_Point) = ((Y_3D(i_Point) - center(2))*Vector_U(1) - (X_3D(i_Point) - center(1))*Vector_U(2))/divisor
  elseif(divisor_type==2) then
      X_2D(i_Point) = ((X_3D(i_Point) - center(1))*Vector_V(3) - (Z_3D(i_Point) - center(3))*Vector_V(1))/divisor
      Y_2D(i_Point) = ((Z_3D(i_Point) - center(3))*Vector_U(1) - (X_3D(i_Point) - center(1))*Vector_U(3))/divisor
  elseif(divisor_type==3) then
      X_2D(i_Point) = ((Y_3D(i_Point) - center(2))*Vector_V(3) - (Z_3D(i_Point) - center(3))*Vector_V(2))/divisor
      Y_2D(i_Point) = ((Z_3D(i_Point) - center(3))*Vector_U(2) - (Y_3D(i_Point) - center(2))*Vector_U(3))/divisor
  endif
enddo

return 
end SUBROUTINE Tool_Map_3D_Points_to_2D                 
