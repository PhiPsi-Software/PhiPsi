 
subroutine Cal_Ele_Num_by_Coors_3D(x,y,z,Ele_Num_Cache,OUT_Elem)

use Global_Float_Type
use Global_Model
use Global_Common
use Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron
use Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron_with_Tol

implicit none
real(kind=FT),intent(in):: x,y,z
integer,intent(out)::OUT_Elem
integer,intent(inout)::Ele_Num_Cache
real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min,c_z_max,c_z_min
integer c_count,Potent_Elem(100)
integer i
integer c_NN(8)
logical c_Yes_in,c_Yes_on
real(kind=FT) In_Hexahedron_Tol
integer c_Ele_Domain_ID,c_Ele
real(kind=FT) Point_A(3),Point_B(3),Point_C(3),Point_D(3),Point_E(3),Point_F(3),Point_G(3),Point_H(3)  

#ifndef Silverfrost
if (isnan(x) .or. isnan(y) .or. isnan(z)) then
  print *, '    ERROR2023040401 :: x or y or z is NaN in Cal_Ele_Num_by_Coors_3D.f90! '
  print *, '                       x, y, z:', x, y, z
  call Warning_Message('S',Keywords_Blank)
end if
#endif
c_count = 0
OUT_Elem= 0

if (x > Max_X_Coor .or. x < Min_X_Coor)then
  return
endif
if (y > Max_Y_Coor .or. y < Min_Y_Coor)then
  return
endif
if (z > Max_Z_Coor .or. z < Min_Z_Coor)then
  return
endif


if(Ele_Num_Cache>Num_Elem) then
    Ele_Num_Cache = Num_Elem 
endif
if(Ele_Num_Cache<1) then
    Ele_Num_Cache = 1 
endif

c_NN  = G_NN(1:8,Ele_Num_Cache)
Point_A = Coor(c_NN(1),1:3)
Point_B = Coor(c_NN(2),1:3)
Point_C = Coor(c_NN(3),1:3)
Point_D = Coor(c_NN(4),1:3)
Point_E = Coor(c_NN(5),1:3)
Point_F = Coor(c_NN(6),1:3)
Point_G = Coor(c_NN(7),1:3)
Point_H = Coor(c_NN(8),1:3)
call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
         Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
         c_Yes_in,c_Yes_on)

if(c_Yes_in)then
 OUT_Elem = Ele_Num_Cache
 Ele_Num_Cache = OUT_Elem
 return
endif

if((Ele_Num_Cache-1) >=1)then
  c_NN  = G_NN(1:8,Ele_Num_Cache-1)
  Point_A = Coor(c_NN(1),1:3)
  Point_B = Coor(c_NN(2),1:3)
  Point_C = Coor(c_NN(3),1:3)
  Point_D = Coor(c_NN(4),1:3)
  Point_E = Coor(c_NN(5),1:3)
  Point_F = Coor(c_NN(6),1:3)
  Point_G = Coor(c_NN(7),1:3)
  Point_H = Coor(c_NN(8),1:3)
  call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             c_Yes_in,c_Yes_on)
  if(c_Yes_in)then
     OUT_Elem = Ele_Num_Cache - 1
     Ele_Num_Cache = OUT_Elem
     return
  endif
endif

if((Ele_Num_Cache+1) <= num_Elem)then
  c_NN  = G_NN(1:8,Ele_Num_Cache+1)
  Point_A = Coor(c_NN(1),1:3)
  Point_B = Coor(c_NN(2),1:3)
  Point_C = Coor(c_NN(3),1:3)
  Point_D = Coor(c_NN(4),1:3)
  Point_E = Coor(c_NN(5),1:3)
  Point_F = Coor(c_NN(6),1:3)
  Point_G = Coor(c_NN(7),1:3)
  Point_H = Coor(c_NN(8),1:3)
  call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             c_Yes_in,c_Yes_on)
  if(c_Yes_in)then
     OUT_Elem = Ele_Num_Cache + 1
     Ele_Num_Cache = OUT_Elem
     return
  endif
endif      


call D3_Get_Point_Domain_Number([x,y,z],c_Ele_Domain_ID)
Potent_Elem =0
do i=1,Domain_Elements_Num(c_Ele_Domain_ID)
  c_Ele = Domain_Elements(c_Ele_Domain_ID,i)
  c_x_max = x_max_Elements(c_Ele)
  c_x_min = x_min_Elements(c_Ele)
  c_y_max = y_max_Elements(c_Ele)
  c_y_min = y_min_Elements(c_Ele)
  c_z_max = z_max_Elements(c_Ele)
  c_z_min = z_min_Elements(c_Ele)
  if(x.lt.(c_x_min-Tol_8)) then
      cycle
  endif
  if(y.lt.(c_y_min-Tol_8)) then
      cycle
  endif
  if(z.lt.(c_z_min-Tol_8)) then
      cycle
  endif
  if(x.gt.(c_x_max+Tol_8)) then
      cycle
  endif
  if(y.gt.(c_y_max+Tol_8)) then
      cycle
  endif
  if(z.gt.(c_z_max+Tol_8)) then
      cycle
  endif
  c_count = c_count +1

  if(c_count>100)then
      print *,'    ERROR-2023040501 :: c_count>100 in Cal_Ele_Num_by_Coors_3D.f90!'
      call Warning_Message('S',Keywords_Blank)
  endif
  Potent_Elem(c_count) = c_Ele
enddo


do i =1,c_count
  c_NN  = G_NN(1:8,Potent_Elem(i))
  Point_A = Coor(c_NN(1),1:3)
  Point_B = Coor(c_NN(2),1:3)
  Point_C = Coor(c_NN(3),1:3)
  Point_D = Coor(c_NN(4),1:3)
  Point_E = Coor(c_NN(5),1:3)
  Point_F = Coor(c_NN(6),1:3)
  Point_G = Coor(c_NN(7),1:3)
  Point_H = Coor(c_NN(8),1:3)
  call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             c_Yes_in,c_Yes_on)
 if(c_Yes_in)then
     OUT_Elem = Potent_Elem(i)
     Ele_Num_Cache = OUT_Elem 
     return
 endif
end do


if(c_count>1 .and. OUT_Elem==0)then
  In_Hexahedron_Tol =  1.0D-5
  do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         Ele_Num_Cache = OUT_Elem 
         return
     endif
  end do
endif

if(c_count>1 .and. OUT_Elem==0)then
In_Hexahedron_Tol =  5.0D-5
do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do
endif

if(c_count>1 .and. OUT_Elem==0)then
In_Hexahedron_Tol =  1.0D-3
do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do
endif  



c_count = 0
OUT_Elem= 0
Potent_Elem = 0
do i=1,Num_Elem
  c_x_max = x_max_Elements(i)
  c_x_min = x_min_Elements(i)
  c_y_max = y_max_Elements(i)
  c_y_min = y_min_Elements(i)
  c_z_max = z_max_Elements(i)
  c_z_min = z_min_Elements(i)
  if(x.lt.(c_x_min-Tol_8)) then
      cycle
  endif
  if(y.lt.(c_y_min-Tol_8)) then
      cycle
  endif
  if(z.lt.(c_z_min-Tol_8)) then
      cycle
  endif
  if(x.gt.(c_x_max+Tol_8)) then
      cycle
  endif
  if(y.gt.(c_y_max+Tol_8)) then
      cycle
  endif
  if(z.gt.(c_z_max+Tol_8)) then
      cycle
  endif
  c_count = c_count +1

  if(c_count>100)then
      print *,'    ERROR-2023040502 :: c_count>100 in Cal_Ele_Num_by_Coors_3D.f90!'
      call Warning_Message('S',Keywords_Blank)
  endif
  
  Potent_Elem(c_count) = i
end do  


do i =1,c_count
  c_NN  = G_NN(1:8,Potent_Elem(i))
  Point_A = Coor(c_NN(1),1:3)
  Point_B = Coor(c_NN(2),1:3)
  Point_C = Coor(c_NN(3),1:3)
  Point_D = Coor(c_NN(4),1:3)
  Point_E = Coor(c_NN(5),1:3)
  Point_F = Coor(c_NN(6),1:3)
  Point_G = Coor(c_NN(7),1:3)
  Point_H = Coor(c_NN(8),1:3)
  call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             c_Yes_in,c_Yes_on)
 if(c_Yes_in)then
     OUT_Elem = Potent_Elem(i)
     Ele_Num_Cache = OUT_Elem 
     return
 endif
end do


if(c_count>1 .and. OUT_Elem==0)then
  In_Hexahedron_Tol =  1.0D-5
  do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         Ele_Num_Cache = OUT_Elem 
         return
     endif
  end do
endif

if(c_count>1 .and. OUT_Elem==0)then
In_Hexahedron_Tol =  5.0D-5
do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do
endif  

if(c_count>1 .and. OUT_Elem==0)then
In_Hexahedron_Tol =  1.0D-3
do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do
endif  

if(c_count>1 .and. OUT_Elem==0)then
In_Hexahedron_Tol =  1.0D-2
do i =1,c_count
     c_NN  = G_NN(1:8,Potent_Elem(i))
     Point_A = Coor(c_NN(1),1:3)
     Point_B = Coor(c_NN(2),1:3)
     Point_C = Coor(c_NN(3),1:3)
     Point_D = Coor(c_NN(4),1:3)
     Point_E = Coor(c_NN(5),1:3)
     Point_F = Coor(c_NN(6),1:3)
     Point_G = Coor(c_NN(7),1:3)
     Point_H = Coor(c_NN(8),1:3)
     call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
             Point_A,Point_B,Point_C,Point_D,Point_E,Point_F,Point_G,Point_H,& 
             In_Hexahedron_Tol,c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do
endif  

return
end SUBROUTINE Cal_Ele_Num_by_Coors_3D
