 
subroutine Tool_Nearest_Point_from_Point_to_Segment_3D(point,line,nearestPoint,distance)

use Global_Float_Type

implicit none
real(kind=FT),intent(in) :: point(3)
real(kind=FT),intent(in) :: line(2,3)
real(kind=FT),intent(out) :: nearestPoint(3),distance
real(kind=FT) :: P1(3),P2(3),t

real(kind=FT) :: P2_P1(3)
real(kind=FT) :: dis_to_P1,dis_to_P2
real(kind=FT) :: direction_vector(3),v(3)

logical Yes

P1(1:3) = line(1,1:3)
P2(1:3) = line(2,1:3)

P2_P1 = P2-P1

direction_vector = (P2-P1)/(sqrt(P2_P1(1)**2 + P2_P1(2)**2 + P2_P1(3)**2))
v = Point-P1

t = dot_product(v,direction_vector)

nearestPoint = P1+t*direction_vector
distance = sqrt((point(1)-nearestPoint(1))**2 +&
                (point(2)-nearestPoint(2))**2 +&
                (point(3)-nearestPoint(3))**2)

dis_to_P1 = sqrt((point(1)-P1(1))**2 +(point(2)-P1(2))**2 +(point(3)-P1(3))**2)
dis_to_P2 = sqrt((point(1)-P2(1))**2 +(point(2)-P2(2))**2 +(point(3)-P2(3))**2)

call Tool_Yes_Point_on_Line_Segment_3D(P1,P2,nearestPoint,Yes)

if(Yes .eqv. .false.) then
    if(dis_to_P1<=dis_to_P2)then
        nearestPoint = P1
        distance = dis_to_P1
    else
        nearestPoint = P2
        distance = dis_to_P2
    endif
endif 

return 
end SUBROUTINE Tool_Nearest_Point_from_Point_to_Segment_3D                         
