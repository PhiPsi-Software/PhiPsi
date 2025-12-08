
subroutine Cal_Ele_Num_by_Coors_3D_v2(x,y,z,OUT_Elem)


use Global_Float_Type
use Global_Model
use Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron

implicit none

real(kind=FT),intent(in):: x,y,z
integer,intent(out)::OUT_Elem
real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min,c_z_max,c_z_min

integer c_count,Potent_Elem(500)
integer i
integer c_NN(8)
logical c_Yes_in,c_Yes_on
real(kind=FT) In_Hexahedron_Tol

c_count = 0
OUT_Elem= 0
Potent_Elem = 0

if (x > Max_X_Coor .or. x < Min_X_Coor)then
      return
endif
if (y > Max_Y_Coor .or. y < Min_Y_Coor)then
      return
endif
if (z > Max_Z_Coor .or. z < Min_Z_Coor)then
      return
endif


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
      
      Potent_Elem(c_count) = i
    
end do



do i =1,c_count
      c_NN  = G_NN(1:8,Potent_Elem(i))
      call Tool_Yes_Point_in_3D_Hexahedron([x,y,z], &
                 Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),&
                 Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),&
                 c_Yes_in,c_Yes_on)
     if(c_Yes_in)then
         OUT_Elem = Potent_Elem(i)
         return
     endif
end do

if(c_count>1 .and. OUT_Elem==0)then
    In_Hexahedron_Tol =  1.0D-5
    print *,'    WARNING2022100702 :: Retry 1 in Cal_Ele_Num_by_Coors_3D_v2.f90!'
    do i =1,c_count
        c_NN  = G_NN(1:8,Potent_Elem(i))
        call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
                     Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),&
                     Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),&
                     In_Hexahedron_Tol,   &
                     c_Yes_in,c_Yes_on)
         if(c_Yes_in)then
             OUT_Elem = Potent_Elem(i)
             return
         endif
    end do
endif

if(c_count>1 .and. OUT_Elem==0)then
    In_Hexahedron_Tol =  5.0D-5
    print *,'    WARNING2022100801 :: Retry 2 in Cal_Ele_Num_by_Coors_3D_v2.f90!'
    do i =1,c_count
        c_NN  = G_NN(1:8,Potent_Elem(i))
        call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
                     Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),&
                     Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),&
                     In_Hexahedron_Tol,   &
                     c_Yes_in,c_Yes_on)
         if(c_Yes_in)then
             OUT_Elem = Potent_Elem(i)
             return
         endif
    end do
endif


if(c_count>1 .and. OUT_Elem==0)then
    In_Hexahedron_Tol =  1.0D-3
    print *,'    WARNING2022100804 :: Retry 3 in Cal_Ele_Num_by_Coors_3D_v2.f90!'
    do i =1,c_count
        c_NN  = G_NN(1:8,Potent_Elem(i))
        call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
                     Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),&
                     Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),&
                     In_Hexahedron_Tol,   &
                     c_Yes_in,c_Yes_on)
         if(c_Yes_in)then
             OUT_Elem = Potent_Elem(i)
             return
         endif
    end do
endif  

if(c_count>1 .and. OUT_Elem==0)then
    In_Hexahedron_Tol =  1.0D-2
    print *,'    WARNING2022100805 :: Retry 4 in Cal_Ele_Num_by_Coors_3D_v2.f90!'
    do i =1,c_count
        c_NN  = G_NN(1:8,Potent_Elem(i))
        call Tool_Yes_Point_in_3D_Hexahedron_with_Tol([x,y,z], &
                     Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),&
                     Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),&
                     In_Hexahedron_Tol,   &
                     c_Yes_in,c_Yes_on)
         if(c_Yes_in)then
             OUT_Elem = Potent_Elem(i)
             return
         endif
    end do
endif  

return
end SUBROUTINE Cal_Ele_Num_by_Coors_3D_v2
