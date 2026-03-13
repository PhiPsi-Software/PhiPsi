subroutine Cal_Shaped_Cracks_2D(isub,c_DISP)
! Calculate and save shaped cracks for post-processing display. 
!
! Save to sccx, sccy and scdx, scdy files.
!
! sccx: store the x coodinates of shaped crack points. 
!       Each line store 1 crack x cooridnate
! sccy: store the y coodinates of shaped crack points. 
!       Each line store 1 crack y cooridnate
! sccx: store the x coodinates of shaped crack points. 
!       Each line store 1 crack x cooridnate
! sccy: store the y coodinates of shaped crack points. 
!       Each line store 1 crack y cooridnate
!
! NEWFTU-2025122602. 2025-12-26.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Common
use Global_HF
use Global_Filename

implicit none
real(kind=FT),intent(in)::c_DISP(Total_FD)
integer,intent(in)::isub
integer i_C,i_S
integer c_Num_Divis_Elment
real(kind=FT) offset_delta_Factor
real(kind=FT) crack_p1(2),crack_p2(2),x(2),y(2)
real(kind=FT) l_seg,l_elem,l_division
integer Num_Division
real(kind=FT) Tool_Function_2Point_Dis
real(kind=FT) offset_delta
integer Flag_Edge_Tip
real(kind=FT) first_Tip(2),second_Tip(2),l_AB(2,2)
real(kind=FT) Offs_Edge_Tip_Up(2), Offs_Edge_Tip_Down(2)
integer Num_S_Div_Points
real(kind=FT) Offsetted_UP(2000,2)
real(kind=FT) Offsetted_DOWN(2000,2)
real(kind=FT) Offsetted_DOWN_Flipped(2000,2)
real(kind=FT) Line_AB(2,2)
real(kind=FT),ALLOCATABLE::Div_Points(:,:) 
real(kind=FT),ALLOCATABLE::Offsetted_D_P_Up(:,:) 
real(kind=FT),ALLOCATABLE::Offsetted_D_P_Down(:,:) 
integer p_count
real(kind=FT),ALLOCATABLE::Shaped_Points(:,:)
integer total_shaped_points, idx, i
character(200) c_File_name_1,c_File_name_2,c_File_name_3,c_File_name_4
character(5) temp  
integer file_unit_1,file_unit_2,file_unit_3,file_unit_4
integer j,c_Elem
real(kind=FT) c_point(2),c_Kesi,c_Yita
real(kind=FT),ALLOCATABLE::shaped_points_disp_x(:) 
real(kind=FT),ALLOCATABLE::shaped_points_disp_y(:) 
real(kind=FT) temp_Disp(2)
real(kind=FT) first_Tip_UP(2),first_Tip_DOWN(2)
real(kind=FT) second_Tip_UP(2),second_Tip_DOWN(2)
real(kind=FT) second_Tip2(2)
real(kind=FT) Offsetted_UP2(2000,2)
real(kind=FT) Offsetted_DOWN2(2000,2)
real(kind=FT) Fir_Poit_Offsetted_UP2(2),Signed_Distance_1,Signed_Distance_2
integer :: n_up, n_down
real(kind=FT) first_Tip3(2)
real(kind=FT), dimension(:,:), allocatable :: Offsetted_DOWN2_Flipped
real(kind=FT) Offsetted_UP3(2000,2)
real(kind=FT) Offsetted_DOWN3(2000,2)
real(kind=FT) Last_Poit_Offsetted_UP3(2)
integer n_down3,n_up3
real(kind=FT), dimension(:,:), allocatable :: Offsetted_DOWN3_Flipped



c_Num_Divis_Elment  = 6
offset_delta_Factor = 0.001


write(temp,'(I5)') isub
c_File_name_1   =  trim(Full_Pathname)//'.sccx'//'_'//ADJUSTL(temp)   
c_File_name_2   =  trim(Full_Pathname)//'.sccy'//'_'//ADJUSTL(temp)   
c_File_name_3   =  trim(Full_Pathname)//'.scdx'//'_'//ADJUSTL(temp)   
c_File_name_4   =  trim(Full_Pathname)//'.scdy'//'_'//ADJUSTL(temp)   
file_unit_1 = 1001
file_unit_2 = 1002
file_unit_3 = 1003
file_unit_4 = 1004
open(unit=file_unit_1, file=c_File_name_1, status='replace', action='write')
close(file_unit_1)
open(unit=file_unit_2, file=c_File_name_2, status='replace', action='write')
close(file_unit_2)
open(unit=file_unit_3, file=c_File_name_3, status='replace', action='write')
close(file_unit_3)
open(unit=file_unit_4, file=c_File_name_4, status='replace', action='write')
close(file_unit_4)

do i_C = 1,num_Crack
    Flag_Edge_Tip      = 0
    Offs_Edge_Tip_Up   = ZR
    Offs_Edge_Tip_Down = ZR
    first_Tip          = ZR
    second_Tip         = ZR
    Offsetted_UP       = ZR
    Offsetted_DOWN     = ZR
    p_count            = 0
    first_Tip_UP       = ZR
    first_Tip_DOWN     = ZR
    second_Tip2        = ZR
    Offsetted_UP2      = ZR
    Offsetted_DOWN2    = ZR
    Offsetted_UP3      = ZR
    Offsetted_DOWN3    = ZR
    if (sum(Crack_Tip_Type(i_C,:)) <= 0) then 
        do i_S = 1,Each_Cr_Poi_Num(i_C)-1
            crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
            crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
            
            Line_AB(1,:) = crack_p1
            Line_AB(2,:) = crack_p2
            
            
            x = [crack_p1(1),crack_p2(1)]
            y = [crack_p1(2),crack_p2(2)]
            l_seg = Tool_Function_2Point_Dis(crack_p1,crack_p2)
            l_elem= Ave_Elem_L_Enrich
            l_division = l_elem/c_Num_Divis_Elment
            Num_Division = nint(l_seg / l_division)
            if (Num_Division<=1) then
                Num_Division = 2;
            endif
            
            
            offset_delta = offset_delta_Factor*l_elem
            if (i_S == 1) then
                if (Crack_Tip_Type(i_C, 1) /= -2) then
                    first_Tip(1) = x(1)
                    first_Tip(2) = y(1)
                else
                    Flag_Edge_Tip = 1
                    l_AB(1, 1) = x(1)
                    l_AB(1, 2) = y(1)
                    l_AB(2, 1) = x(2)
                    l_AB(2, 2) = y(2)
                    call Cal_Offseted_Single_Point([x(1), y(1)], l_AB, offset_delta, Offs_Edge_Tip_Up, Offs_Edge_Tip_Down)
                                                                              
                end if
            end if

            if (i_S == Each_Cr_Poi_Num(i_C)-1) then
                if (Crack_Tip_Type(i_C, 2) /= -2) then
                    second_Tip(1) = x(2)
                    second_Tip(2) = y(2)
                else
                    Flag_Edge_Tip = 2
                    l_AB(1, 1) = x(1)
                    l_AB(1, 2) = y(1)
                    l_AB(2, 1) = x(2)
                    l_AB(2, 2) = y(2)
                    call Cal_Offseted_Single_Point([x(2), y(2)], l_AB, offset_delta, Offs_Edge_Tip_Up, Offs_Edge_Tip_Down)
                end if
            end if
            
            
            
            allocate(Div_Points(Num_Division-1,2)) 
            allocate(Offsetted_D_P_Up(Num_Division-1,2)) 
            allocate(Offsetted_D_P_Down(Num_Division-1,2)) 
            call Cal_Equal_Division_Points_and_Offset(Num_Division, Line_AB,&
                    offset_delta,Div_Points, &
                    Offsetted_D_P_Up, Offsetted_D_P_Down)
                                      
            
            
            Offsetted_UP(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Up(1:Num_Division-1,1)
            Offsetted_UP(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Up(1:Num_Division-1,2)
            
            Offsetted_DOWN(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Down(1:Num_Division-1,1)
            Offsetted_DOWN(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Down(1:Num_Division-1,2)
            
            p_count = p_count+(Num_Division-1)
            
            deallocate(Div_Points)
            deallocate(Offsetted_D_P_Up)
            deallocate(Offsetted_D_P_Down)
            
            
        enddo
        
        
              
        do i = 1, p_count
            Offsetted_DOWN_Flipped(i, :) = Offsetted_DOWN(p_count - i + 1, :)
        end do
        
        
        
        if (Flag_Edge_Tip == 0) then
            total_shaped_points = 1 + p_count + 1 + p_count
            allocate(Shaped_Points(total_shaped_points, 2))
            
            idx = 1
            Shaped_Points(idx, :) = first_Tip
            idx = idx + 1
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_UP(1:p_count, :)
            idx = idx + p_count
            
            Shaped_Points(idx, :) = second_Tip
            idx = idx + 1
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_DOWN_Flipped(1:p_count, :)
            
        elseif (Flag_Edge_Tip == 1) then
            total_shaped_points = 1 + p_count + 1 + p_count + 1
            allocate(Shaped_Points(total_shaped_points, 2))
            
            idx = 1
            Shaped_Points(idx, :) = Offs_Edge_Tip_Up
            idx = idx + 1
            
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_UP(1:p_count, :)
            idx = idx + p_count
            
            Shaped_Points(idx, :) = second_Tip
            idx = idx + 1
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_DOWN_Flipped(1:p_count, :)
            idx = idx + p_count
            
            Shaped_Points(idx, :) = Offs_Edge_Tip_Down
            
        elseif (Flag_Edge_Tip == 2) then
            total_shaped_points = 1 + p_count + 1 + 1 + p_count
            allocate(Shaped_Points(total_shaped_points, 2))
            
            idx = 1
            
            Shaped_Points(idx, :) = first_Tip
            idx = idx + 1
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_UP(1:p_count, :)
            idx = idx + p_count
            
            Shaped_Points(idx, :) = Offs_Edge_Tip_Up
            idx = idx + 1
            
            
            Shaped_Points(idx, :) = Offs_Edge_Tip_Down
            idx = idx + 1
            
            Shaped_Points(idx:idx+p_count-1, :) = Offsetted_DOWN_Flipped(1:p_count, :)
            
        endif
        
        
    else
        if (Crack_Tip_Type(i_C,1)==1) then
            do i_S = 1,Each_Cr_Poi_Num(i_C)-1
                crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
                crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
                
                Line_AB(1,:) = crack_p1
                Line_AB(2,:) = crack_p2
                
                
                x = [crack_p1(1),crack_p2(1)]
                y = [crack_p1(2),crack_p2(2)]
                l_seg = Tool_Function_2Point_Dis(crack_p1,crack_p2)
                l_elem= Ave_Elem_L_Enrich
                l_division = l_elem/c_Num_Divis_Elment
                Num_Division = nint(l_seg / l_division)
                if (Num_Division<=1) then
                    Num_Division = 2;
                endif
                
                
                offset_delta = offset_delta_Factor*l_elem
                if (i_S == 1) then
                    l_AB(1, 1) = x(1)
                    l_AB(1, 2) = y(1)
                    l_AB(2, 1) = x(2)
                    l_AB(2, 2) = y(2)
                    call Cal_Offseted_Single_Point([x(1), y(1)], l_AB, offset_delta, first_Tip_UP,first_Tip_DOWN)
                end if

                if (i_S == Each_Cr_Poi_Num(i_C)-1) then
                    second_Tip2(1) = x(2)
                    second_Tip2(2) = y(2)
                end if
                
                
                
                allocate(Div_Points(Num_Division-1,2)) 
                allocate(Offsetted_D_P_Up(Num_Division-1,2)) 
                allocate(Offsetted_D_P_Down(Num_Division-1,2)) 
                call Cal_Equal_Division_Points_and_Offset(Num_Division, &
                     Line_AB, offset_delta,Div_Points, Offsetted_D_P_Up, &
                     Offsetted_D_P_Down)
                                          
                
                
                Offsetted_UP2(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Up(1:Num_Division-1,1)
                Offsetted_UP2(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Up(1:Num_Division-1,2)
                
                Offsetted_DOWN2(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Down(1:Num_Division-1,1)
                Offsetted_DOWN2(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Down(1:Num_Division-1,2)
                
                p_count = p_count+(Num_Division-1)
                
                deallocate(Div_Points)
                deallocate(Offsetted_D_P_Up)
                deallocate(Offsetted_D_P_Down)
                
            enddo
            
            Fir_Poit_Offsetted_UP2 = Offsetted_UP2(1,:)
            
            call Cal_Signed_Distance(l_AB,Fir_Poit_Offsetted_UP2,Signed_Distance_1)
            call Cal_Signed_Distance(l_AB,first_Tip_UP,Signed_Distance_2)
            


            n_down = p_count
            allocate(Offsetted_DOWN2_Flipped(n_down, 2))
            do i = 1, n_down
                Offsetted_DOWN2_Flipped(i, :) = Offsetted_DOWN2(n_down - i + 1, :)
            end do


            n_up = p_count

            if (Signed_Distance_1 * Signed_Distance_2 > ZR) then
                total_Shaped_Points = 1 + 1 + n_up + 1 + n_down
                allocate(Shaped_Points(total_Shaped_Points, 2))
                
                idx = 1
                Shaped_Points(idx, :) = first_Tip_DOWN
                idx = idx + 1
                
                Shaped_Points(idx, :) = first_Tip_UP
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_up-1, :) = Offsetted_UP2(1:n_up, :)
                idx = idx + n_up
                
                Shaped_Points(idx, :) = second_Tip2
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_down-1, :) = Offsetted_DOWN2_Flipped(1:n_down, :)
                
            else
                total_Shaped_Points = 1 + 1 + n_up + 1 + n_down
                allocate(Shaped_Points(total_Shaped_Points, 2))
                
                idx = 1
                Shaped_Points(idx, :) = first_Tip_UP
                idx = idx + 1
                
                Shaped_Points(idx, :) = first_Tip_DOWN
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_up-1, :) = Offsetted_UP2(1:n_up, :)
                idx = idx + n_up
                
                Shaped_Points(idx, :) = second_Tip2
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_down-1, :) = Offsetted_DOWN2_Flipped(1:n_down, :)
                
            end if
            
            deallocate(Offsetted_DOWN2_Flipped)
        
        endif
        
        if (Crack_Tip_Type(i_C,2)==1) then
            p_count=0
            Offsetted_UP3 = ZR
            Offsetted_DOWN3 = ZR
            do i_S = 1,Each_Cr_Poi_Num(i_C)-1
                crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
                crack_p2 = [Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
                
                Line_AB(1,:) = crack_p1
                Line_AB(2,:) = crack_p2
                
                
                x = [crack_p1(1),crack_p2(1)]
                y = [crack_p1(2),crack_p2(2)]
                l_seg = Tool_Function_2Point_Dis(crack_p1,crack_p2)
                l_elem= Ave_Elem_L_Enrich
                l_division = l_elem/c_Num_Divis_Elment
                Num_Division = nint(l_seg / l_division)
                if (Num_Division<=1) then
                    Num_Division = 2;
                endif
                
                
                offset_delta = offset_delta_Factor*l_elem
                if (i_S == 1) then
                    first_Tip3(1) = x(1)
                    first_Tip3(2) = y(1)

                end if

                if (i_S == Each_Cr_Poi_Num(i_C)-1) then
                    l_AB(1, 1) = x(1)
                    l_AB(1, 2) = y(1)
                    l_AB(2, 1) = x(2)
                    l_AB(2, 2) = y(2)
                    call Cal_Offseted_Single_Point([x(2), y(2)], l_AB, offset_delta, second_Tip_UP,second_Tip_DOWN)
                end if
                
                
                
                allocate(Div_Points(Num_Division-1,2)) 
                allocate(Offsetted_D_P_Up(Num_Division-1,2)) 
                allocate(Offsetted_D_P_Down(Num_Division-1,2)) 
                call Cal_Equal_Division_Points_and_Offset(Num_Division, &
                    Line_AB, offset_delta,Div_Points,&
                    Offsetted_D_P_Up, Offsetted_D_P_Down)
                
                Offsetted_UP3(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Up(1:Num_Division-1,1)
                Offsetted_UP3(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Up(1:Num_Division-1,2)
                
                Offsetted_DOWN3(p_count+1:p_count+(Num_Division-1),1) = Offsetted_D_P_Down(1:Num_Division-1,1)
                Offsetted_DOWN3(p_count+1:p_count+(Num_Division-1),2) = Offsetted_D_P_Down(1:Num_Division-1,2)
                
                p_count = p_count+(Num_Division-1)
                
                deallocate(Div_Points)
                deallocate(Offsetted_D_P_Up)
                deallocate(Offsetted_D_P_Down)
                
            enddo
            
            Last_Poit_Offsetted_UP3 = Offsetted_UP3(p_count,:)
            
            call Cal_Signed_Distance(l_AB,Last_Poit_Offsetted_UP3,Signed_Distance_1)
            call Cal_Signed_Distance(l_AB,second_Tip_UP,Signed_Distance_2)
            


            n_down3 = p_count
            allocate(Offsetted_DOWN3_Flipped(n_down3, 2))
            do i = 1, n_down3
                Offsetted_DOWN3_Flipped(i, :) = Offsetted_DOWN3(n_down3 - i + 1, :)
            end do

            n_up3 = p_count

            if (Signed_Distance_1 * Signed_Distance_2 > ZR) then
                total_shaped_points = 1 + n_up3 + 1 + 1 + n_down3
                if (allocated(Shaped_Points)) deallocate(Shaped_Points)
                allocate(Shaped_Points(total_shaped_points, 2))
                
                idx = 1
                Shaped_Points(idx, :) = first_Tip3
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_up3-1, :) = Offsetted_UP3(1:n_up3, :)
                idx = idx + n_up3
                
                Shaped_Points(idx, :) = second_Tip_UP
                idx = idx + 1
                
                Shaped_Points(idx, :) = second_Tip_DOWN
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_down3-1, :) = Offsetted_DOWN3_Flipped(1:n_down3, :)
                
            else
                total_shaped_points = 1 + n_up3 + 1 + 1 + n_down3
                allocate(Shaped_Points(total_shaped_points, 2))
                
                idx = 1
                Shaped_Points(idx, :) = first_Tip3
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_up3-1, :) = Offsetted_UP3(1:n_up3, :)
                idx = idx + n_up3
                
                Shaped_Points(idx, :) = second_Tip_DOWN
                idx = idx + 1
                
                Shaped_Points(idx, :) = second_Tip_UP
                idx = idx + 1
                
                Shaped_Points(idx:idx+n_down3-1, :) = Offsetted_DOWN3_Flipped(1:n_down3, :)
                
            end if

            deallocate(Offsetted_DOWN3_Flipped)
            
            
        endif
        
                
    endif
    

    
    open(unit=file_unit_1, file=c_File_name_1, status='old', position='append', action='write')
    open(unit=file_unit_2, file=c_File_name_2, status='old', position='append', action='write')
    
    do j = 1, total_shaped_points
        if (j < total_shaped_points) then
            write(file_unit_1, '(ES25.16E3,1X)', advance='no') Shaped_Points(j, 1)
        else
            write(file_unit_1, '(ES25.16E3)') Shaped_Points(j, 1)
        endif
    end do
    
    do j = 1, total_shaped_points
        if (j < total_shaped_points) then
            write(file_unit_2, '(ES25.16E3,1X)', advance='no') Shaped_Points(j, 2)
        else
            write(file_unit_2, '(ES25.16E3)') Shaped_Points(j, 2)
        endif
    end do
    
    close(file_unit_1)
    close(file_unit_2)
    
    allocate(shaped_points_disp_x(total_shaped_points)) 
    allocate(shaped_points_disp_y(total_shaped_points)) 
    shaped_points_disp_x =  ZR
    shaped_points_disp_y =  ZR
    
    do j =1,total_shaped_points
        c_point = Shaped_Points(j,1:2)
        call Cal_Ele_Num_by_Coors(c_point(1),c_point(2),c_Elem)
        if (c_Elem >=1) then
            call Cal_KesiYita_by_Coor(c_point,c_Elem,c_Kesi,c_Yita)
            call Cal_Any_Point_Disp_KesiYita(c_Elem,c_Kesi,c_Yita,0,c_DISP,temp_Disp)
            shaped_points_disp_x(j) = temp_Disp(1)
            shaped_points_disp_y(j) = temp_Disp(2)
        elseif (c_Elem ==0) then
        
            shaped_points_disp_x(j) = ZR
            shaped_points_disp_y(j) = ZR
        endif
    enddo
    
    open(unit=file_unit_3, file=c_File_name_3, status='old', position='append', action='write')
    open(unit=file_unit_4, file=c_File_name_4, status='old', position='append', action='write')
    
    do j = 1, total_shaped_points
        if (j < total_shaped_points) then
            write(file_unit_3, '(ES25.16E3,1X)', advance='no') shaped_points_disp_x(j)
        else
            write(file_unit_3, '(ES25.16E3)') shaped_points_disp_x(j)
        endif
    end do
    
    do j = 1, total_shaped_points
        if (j < total_shaped_points) then
            write(file_unit_4, '(ES25.16E3,1X)', advance='no') shaped_points_disp_y(j)
        else
            write(file_unit_4, '(ES25.16E3)') shaped_points_disp_y(j)
        endif
    end do
    
    close(file_unit_3)
    close(file_unit_4)
    
    deallocate(Shaped_Points)
    deallocate(shaped_points_disp_x) 
    deallocate(shaped_points_disp_y) 
    
    


enddo


return 
end SUBROUTINE Cal_Shaped_Cracks_2D              
