 
subroutine Tool_Fit_3D_Points_Check_Rollback(i_C,num_Point,Pre_Fronts,Points_to_Smooth,Points_Smoothed)

use Global_Float_Type     
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_Common
      
implicit none
integer,intent(in)::i_C,num_Point
real(kind=FT),intent(in)::Pre_Fronts(num_Point,3),Points_to_Smooth(num_Point,3)
real(kind=FT),intent(inout)::Points_Smoothed(num_Point,3)
integer i_V,i_Try,num_Try
real(kind=FT) c_point_after_Smooth(3)
real(kind=FT) c_point_before_Smooth(3)
real(kind=FT) c_point_after_Smooth_Updated(3)
logical Yes_Smoothed_Point_in_Pre_Fronts
real(kind=FT) c_ratio
real(kind=FT) Pre_Point(3),Vector_V1(3),Vector_V2(3),c_angle


num_Try= 100

if(Key_InPlane_Growth==1) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_V,c_point_after_Smooth,c_point_before_Smooth,i_Try, &
    !$OMP              c_ratio,c_point_after_Smooth_Updated,Yes_Smoothed_Point_in_Pre_Fronts)   
    do i_V = 1,num_Point
        
        c_point_after_Smooth = Points_Smoothed(i_V,1:3)
        
        c_point_before_Smooth = Points_to_Smooth(i_V,1:3)
        
        
        do i_Try= 1,num_Try
            c_ratio = dble(num_Try-i_Try+1)/dble(num_Try)
            c_point_after_Smooth_Updated = c_point_before_Smooth + c_ratio*(c_point_after_Smooth - c_point_before_Smooth)
            call Tool_Yes_Point_in_3D_Plane_Polygon(num_Point,Pre_Fronts,c_point_after_Smooth_Updated,&
                                                    Yes_Smoothed_Point_in_Pre_Fronts)
                                                
            if(Yes_Smoothed_Point_in_Pre_Fronts .eqv. .False.) then  
                Points_Smoothed(i_V,1:3) = c_point_after_Smooth_Updated
                goto 100
            endif
        enddo
        
        Points_Smoothed(i_V,1:3) = c_point_before_Smooth
        
        100 continue
        
    enddo
    !$OMP END PARALLEL DO
    
else

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_V,c_point_after_Smooth,c_point_before_Smooth,i_Try, &
    !$OMP              Pre_Point,Vector_V1,Vector_V2,c_ratio,c_point_after_Smooth_Updated,c_angle)   
    do i_V = 1,num_Point
        
        c_point_after_Smooth = Points_Smoothed(i_V,1:3)
        
        c_point_before_Smooth = Points_to_Smooth(i_V,1:3)
        
        
        Pre_Point = Pre_Fronts(i_V,1:3)
        
        Vector_V1 = c_point_before_Smooth - Pre_Point
        
        do i_Try= 1,num_Try
            c_ratio = dble(num_Try-i_Try+1)/dble(num_Try)
            c_point_after_Smooth_Updated = c_point_before_Smooth + c_ratio*(c_point_after_Smooth - c_point_before_Smooth)
            
            Vector_V2 = c_point_after_Smooth_Updated - Pre_Point
            
            call Tool_Angle_of_Vectors_a_and_b_3D(Vector_V1,Vector_V2,c_angle,2)

                                                    
                                                
            if(c_angle <= (pi/TWO)) then  
                Points_Smoothed(i_V,1:3) = c_point_after_Smooth_Updated
                goto 200
            endif
        enddo
        
        Points_Smoothed(i_V,1:3) = c_point_before_Smooth
        
        200 continue
        
    enddo
    !$OMP END PARALLEL DO
endif
return 
end SUBROUTINE Tool_Fit_3D_Points_Check_Rollback           
