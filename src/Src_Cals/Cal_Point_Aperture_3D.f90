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
 
subroutine Cal_Point_Aperture_3D(c_Crack,c_Point,c_DISP,ori_n,T,Aperture,print_tab)
! Calculate the opening at a point on the crack (3D).
! 2020-01-18.

!+++++++++++++++++++++++++++++
! Read public variable module
!+++++++++++++++++++++++++++++
use Global_Float_Type
use Global_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol

use Global_Cal_Ele_Num_by_Coors_3D  

implicit none
real(kind=FT),intent(in)::c_Point(3),c_DISP(Total_FD),ori_n(3),T(3,3)
integer,intent(in)::print_tab
integer,intent(in)::c_Crack
real(kind=FT),intent(out)::Aperture
real(kind=FT) delta_L
real(kind=FT) up_offset_P(3),low_offset_P(3)
integer up_Elem_num,low_Elem_num
real(kind=FT) up_Kesi,up_Yita,up_Zeta
real(kind=FT) low_Kesi,low_Yita,low_Zeta
real(kind=FT) up_P_Disp(3),low_P_Disp(3),Relative_Disp(3)
real(kind=FT) norm_Flu_Ele_Vector
real(kind=FT) delta_L_Factor
real(kind=FT) offset_dis
integer Ele_Num_Cache

delta_L_Factor = 0.02D0
Aperture       = ZR

Ele_Num_Cache = 1
delta_L = delta_L_Factor*Ave_Elem_L_Enrich
up_offset_P(1:3) = c_Point + delta_L*ori_n(1:3)
low_offset_P(1:3)= c_Point - delta_L*ori_n(1:3)
call Cal_Ele_Num_by_Coors_3D(up_offset_P(1),up_offset_P(2),up_offset_P(3),Ele_Num_Cache,up_Elem_num)        
call Cal_Ele_Num_by_Coors_3D(low_offset_P(1),low_offset_P(2),low_offset_P(3),Ele_Num_Cache,low_Elem_num)  
call Cal_KesiYita_by_Coor_3D(up_offset_P,up_Elem_num,up_Kesi,up_Yita,up_Zeta)
call Cal_KesiYita_by_Coor_3D(low_offset_P,low_Elem_num,low_Kesi,low_Yita,low_Zeta)  
if(up_Elem_num<=0)then
      offset_dis = delta_L/TWO
      call Cal_Ele_Num_by_Coors_3D_Try(up_offset_P(1), &
                  up_offset_P(2),up_offset_P(3),offset_dis,     &
                  up_offset_P(1),up_offset_P(2),up_offset_P(3),up_Elem_num)
      if(up_Elem_num<=0)then
        if (print_tab==5)then
          print *,"     WARNING :: illegal up_Elem_num in Cal_Point_Aperture_3D.f"
        elseif(print_tab==10)then
          print *,"           WARNING :: illegal up_Elem_num in Cal_Point_Aperture_3D.f"          
        endif
        return            
      endif
endif
if(low_Elem_num<=0)then
      offset_dis = delta_L/TWO
      call Cal_Ele_Num_by_Coors_3D_Try(low_offset_P(1),  &
                                       low_offset_P(2),low_offset_P(3),offset_dis,  &
                                       low_offset_P(1),low_offset_P(2),low_offset_P(3),low_Elem_num)
      if(low_Elem_num<=0)then                      
          if (print_tab==5)then
            print *,"     WARNING :: illegal low_Elem_num in Cal_Point_Aperture_3D.f"
          elseif(print_tab==10)then
            print *,"           WARNING :: illegal low_Elem_num in Cal_Point_Aperture_3D.f"         
          endif     
          return
      endif
endif

call Cal_Any_Point_Disp_KesiYita_3D(up_Elem_num,up_Kesi,up_Yita,up_Zeta,c_DISP,up_P_Disp)
call Cal_Any_Point_Disp_KesiYita_3D(low_Elem_num,low_Kesi,low_Yita,low_Zeta,c_DISP,low_P_Disp)    
    
Relative_Disp(1:3) = up_P_Disp- low_P_Disp      

norm_Flu_Ele_Vector = ONE
Aperture   = (Relative_Disp(1)*ori_n(1)+ Relative_Disp(2)*ori_n(2)+ Relative_Disp(3)*ori_n(3))/norm_Flu_Ele_Vector
  

return 
end SUBROUTINE Cal_Point_Aperture_3D           
