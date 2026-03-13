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
 
subroutine Tool_Fit_3D_Points_Check_Rollback(i_C,num_Point,Pre_Fronts,Points_to_Smooth,Points_Smoothed)
!Prevent the fracture front vertex from falling back into the original fracture front polygon after
!optimization.
!IMPROV2022122001.
!Ref: PhiPsi V1,P91.
!Note that Points Smoothed is an input/output variable.

!.....................
!Variable declaration
!.....................
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
