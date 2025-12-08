 
subroutine Tool_Sort_by_End_to_End(m,m_OP,Input_Outline,Output_Outline,cou)

use Global_Float_Type
implicit none
integer,intent(in)::m,m_OP 
integer,intent(in)::Input_Outline(m,2)
integer,intent(out)::cou 
integer,intent(out)::Output_Outline(m,2)
      
integer i,j,c_Outline,c_node_num_2
logical tem_Yes
integer Location



Output_Outline = 0

Output_Outline(1,:)  = Input_Outline(1,:)


cou = 1
c_Outline    = 1

do i = 1,m_OP
    if (i.eq.1)then
        c_Outline = 1
        c_node_num_2 = Input_Outline(1,2)
    end if
    
    do j = 1,m
        if(j .ne. c_Outline) then
            if (c_node_num_2 .eq. Input_Outline(j,1)) then
                call Vector_belongs_Matrix_Is_Int(cou,2,Output_Outline,[Input_Outline(j,1),Input_Outline(j,2)],Location,tem_Yes)
                if(tem_Yes.eqv..False.) then
                    cou = cou + 1
                    Output_Outline(cou,:) = [Input_Outline(j,1),Input_Outline(j,2)]
                    c_Outline = j
                    c_node_num_2 = Input_Outline(j,2)
                    exit
               end if
            else if(c_node_num_2 .eq. Input_Outline(j,2)) then
                call Vector_belongs_Matrix_Is_Int(cou,2,Output_Outline,&
                    [Input_Outline(j,2),Input_Outline(j,1)],Location,tem_Yes)
                if(tem_Yes.eqv..False.)then
                   cou = cou + 1
                   Output_Outline(cou,:)=[Input_Outline(j,2),Input_Outline(j,1)]
                   c_Outline = j
                   c_node_num_2 = Input_Outline(j,1)
                   exit
                  end if
             end if
        end if
     end do
end do


RETURN
end subroutine Tool_Sort_by_End_to_End
