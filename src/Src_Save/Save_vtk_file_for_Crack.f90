 
SUBROUTINE Save_vtk_file_for_Crack(isub)

use Global_Float_Type      
use Global_Common
use Global_Model
use Global_Filename
use Global_Stress
use Global_DISP
use Global_POST
use Global_Crack_Common
use Global_Crack
use Global_Crack_3D
      
implicit none

integer isub
integer I
character(256) DATFIL,VTKFIL, TEMP
integer IUNIT,i_Node,ind,i_E
integer num_Crack_Elements
integer num_Crack_Points
real(kind=FT),ALLOCATABLE::Global_Crack_Points(:,:)
integer,ALLOCATABLE::Global_Crack_Element_Points(:,:)
integer,allocatable::Temp_Global_Index(:,:)
integer,ALLOCATABLE::Global_Crack_Element_Crack_ID(:)
integer,ALLOCATABLE::Global_Crack_Node_Crack_ID(:,:)
integer c_count,i_C
integer i_Crack_Node,i_Crack_Ele
integer c_crack_number,c_local_node
real(kind=FT),allocatable::Crack_nodes_coor_2d(:,:)
integer,allocatable::Crack_nodes_index_2d(:,:)
integer num_2d_crack_nodes,num_2d_crack_elements
integer,allocatable::num_Elem_Eack_2d_Crack(:)
integer i_C_2D,i_Node_2d,counter_2d_node,i_Elem_2d
integer max_num_Elem_Eack_2d_Crack
real(kind=FT) c_node_point_2d(2)

if (Key_Save_Nothing==1) return 

if (Key_Save_vtk/=1) return

if (Key_Dimension == 2) then
    if (num_Crack==0) return
    print *,'    Saving vtk file for crack...'
    
    IUNIT = 101
    DATFIL= trim(Full_Pathname)
    VTKFIL=''
    I=len_trim(DATFIL)+1
    VTKFIL(1:I)=DATFIL(1:I)
    WRITE(TEMP,'(I5.5)')isub
    VTKFIL(I:I+5)='_CRACK'
    VTKFIL(I+6:I+15)='_'//TRIM(TEMP)//'.vtk'
    
    open(IUNIT,file=VTKFIL,status='unknown')
    WRITE(IUNIT,'(A)')'# vtk DataFile Version 4.0'
    WRITE(IUNIT,'(A)')DATFIL(1:I)//' - results from increment '//TRIM(TEMP)
    
    WRITE(IUNIT,'(A)')'ASCII'
    WRITE(IUNIT,'(A)')'DATASET UNSTRUCTURED_GRID'

    allocate(num_Elem_Eack_2d_Crack(num_Crack))
    do i_C_2D=1,num_Crack
        num_Elem_Eack_2d_Crack(i_C_2D) = Each_Cr_Poi_Num(i_C_2D)-1
    end do
    max_num_Elem_Eack_2d_Crack = maxval(num_Elem_Eack_2d_Crack(:))
    allocate(Crack_nodes_index_2d(num_Crack*max_num_Elem_Eack_2d_Crack,2))
    Crack_nodes_index_2d = 0
    num_2d_crack_nodes = sum(Each_Cr_Poi_Num(1:num_Crack))
    allocate(Crack_nodes_coor_2d(num_2d_crack_nodes,2))
    Crack_nodes_coor_2d =  ZR

    counter_2d_node = 0
    num_2d_crack_elements = 0
    do i_C_2D=1,num_Crack
        do i_Node_2d = 1,Each_Cr_Poi_Num(i_C_2D)
            c_node_point_2d(1:2) = Crack_Coor(i_C_2D,i_Node_2d,1:2)
            counter_2d_node = counter_2d_node + 1
            Crack_nodes_coor_2d(counter_2d_node,1:2) = c_node_point_2d(1:2)
            if(i_Node_2d < Each_Cr_Poi_Num(i_C_2D)) then
                num_2d_crack_elements = num_2d_crack_elements + 1  
                Crack_nodes_index_2d(num_2d_crack_elements,1) = counter_2d_node
                Crack_nodes_index_2d(num_2d_crack_elements,2) = counter_2d_node + 1

            endif
        end do
    end do

    if(counter_2d_node/=num_2d_crack_nodes) then
        print *, '    Error :: counter_2d_node /= num_2d_crack_nodes!'
        print *, '             In Save_vtk_file_for_Crack.f90!'
        call Warning_Message('S',Keywords_Blank)
    endif
    
    WRITE(IUNIT,'(/A,I0,A)')'POINTS ',num_2d_crack_nodes,' double'
    if(num_Crack >0) then
          do i_Node=1,num_2d_crack_nodes
              WRITE(IUNIT,'(3F12.6)') Crack_nodes_coor_2d(i_Node,1:2), ZR 
          enddo 
    endif    



    ind = 0
    ind = ind + num_2d_crack_elements*(2+1)
    WRITE(IUNIT,'(/A,I0,A,I0)')'CELLS ',num_2d_crack_elements,' ',ind
    do i_Elem_2d = 1,num_2d_crack_elements
        write(IUNIT,'(3I10)') 2,Crack_nodes_index_2d(i_Elem_2d,1:2) -1
    enddo

    write(IUNIT,'(/A,I10)') "CELL_TYPES", num_2d_crack_elements
    DO i_Elem_2d=1,num_2d_crack_elements
        write(IUNIT,'(I3)') 3
    END DO
    
    write(IUNIT,'(/A,I10)') "CELL_DATA", num_2d_crack_elements


    WRITE(IUNIT,'(/A,/A)')'SCALARS Crack_Element_ID int','LOOKUP_TABLE default'      
    do i_Elem_2d=1,num_2d_crack_elements
        write(IUNIT,'(I8)') i_Elem_2d
    enddo
    
    write(IUNIT,'(/A,I10)') "POINT_DATA", num_2d_crack_nodes

    WRITE(IUNIT,'(/A,/A)')'SCALARS Crack_Node_Number integer','LOOKUP_TABLE default'
    do i_Node=1,num_2d_crack_nodes
        write(IUNIT,'(I8)') i_Node
    enddo

    
    close(IUNIT)      
elseif (Key_Dimension == 3) then
    if (num_Crack==0) return

    print *,'    Saving vtk file for crack...'
    IUNIT = 101
    DATFIL= trim(Full_Pathname)
    VTKFIL=''
    I=len_trim(DATFIL)+1
    VTKFIL(1:I)=DATFIL(1:I)
    WRITE(TEMP,'(I5.5)')isub
    VTKFIL(I:I+5)='_CRACK'
    VTKFIL(I+6:I+15)='_'//TRIM(TEMP)//'.vtk'

    open(IUNIT,file=VTKFIL,status='unknown')
    WRITE(IUNIT,'(A)')'# vtk DataFile Version 4.0'
    WRITE(IUNIT,'(A)')DATFIL(1:I)//' - results from increment '//TRIM(TEMP)

    WRITE(IUNIT,'(A)')'ASCII'
    WRITE(IUNIT,'(A)')'DATASET UNSTRUCTURED_GRID'
          
    num_Crack_Points = 0
    num_Crack_Elements = 0

    if(num_Crack >0 .and. Key_Dimension==3) then
    do i_C = 1,num_Crack
          num_Crack_Points = num_Crack_Points + Crack3D_Meshed_Node_num(i_C)  
          num_Crack_Elements = num_Crack_Elements + Crack3D_Meshed_Ele_num(i_C)
    enddo
    
    if (allocated(Global_Crack_Points)) deallocate(Global_Crack_Points)
    allocate(Global_Crack_Points(num_Crack_Points,3))
    Global_Crack_Points(1:num_Crack_Points,1:3) = ZR
    if (allocated(Global_Crack_Element_Points)) deallocate(Global_Crack_Element_Points)
    allocate(Global_Crack_Element_Points(num_Crack_Elements,3))  
    Global_Crack_Element_Points(1:num_Crack_Elements,1:3) = 0 

    if (allocated(Temp_Global_Index)) deallocate(Temp_Global_Index)
    allocate(Temp_Global_Index(num_Crack,maxval(Crack3D_Meshed_Node_num(1:num_Crack)))) 

    if (allocated(Global_Crack_Element_Crack_ID)) deallocate(Global_Crack_Element_Crack_ID)
    allocate(Global_Crack_Element_Crack_ID(num_Crack_Elements)) 
    Global_Crack_Element_Crack_ID(1:num_Crack_Elements) = 0

    if (allocated(Global_Crack_Node_Crack_ID)) deallocate(Global_Crack_Node_Crack_ID)
    allocate(Global_Crack_Node_Crack_ID(num_Crack_Points,2)) 
    Global_Crack_Node_Crack_ID(num_Crack_Points,2) = 0
           
    Temp_Global_Index(1:num_Crack,1:maxval(Crack3D_Meshed_Node_num(1:num_Crack))) = 0
    c_count = 0
    do i_C = 1,num_Crack
          do i_Crack_Node = 1,Crack3D_Meshed_Node_num(i_C) 
              c_count = c_count + 1
              Global_Crack_Points(c_count,1:3) = Crack3D_Meshed_Node(i_C)%row(i_Crack_Node,1:3) 
              Temp_Global_Index(i_C,i_Crack_Node) = c_count         
              Global_Crack_Node_Crack_ID(c_count,1) = i_C            
              Global_Crack_Node_Crack_ID(c_count,2) = i_Crack_Node  
          enddo
    enddo
    
    c_count = 0
    do i_C = 1,num_Crack
          do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C) 
              c_count = c_count + 1
              Global_Crack_Element_Points(c_count,1) = Temp_Global_Index(i_C,Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1))  
              Global_Crack_Element_Points(c_count,2) = Temp_Global_Index(i_C,Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2))  
              Global_Crack_Element_Points(c_count,3) = Temp_Global_Index(i_C,Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3))  
              Global_Crack_Element_Crack_ID(c_count) = i_C  
          enddo
    enddo
    endif
          

    WRITE(IUNIT,'(/A,I0,A)')'POINTS ',num_Crack_Points,' double'
    if(num_Crack >0) then
          do i_Node=1,num_Crack_Points
              WRITE(IUNIT,'(3F12.6)') Global_Crack_Points(i_Node,1:3) 
          enddo 
    endif    

    ind = 0

    if(num_Crack >0 .and. Key_Dimension==3) then
      ind = ind + num_Crack_Elements*(3+1)
    endif
          
    WRITE(IUNIT,'(/A,I0,A,I0)')'CELLS ',num_Crack_Elements,' ',ind

    IF(Key_Dimension == 2) THEN
    ELSEIF(Key_Dimension == 3) then
        if(num_Crack >0) then
            do i_Crack_Ele = 1,num_Crack_Elements
                write(IUNIT,'(4I10)') 3,Global_Crack_Element_Points(i_Crack_Ele,1:3) -1
            enddo
        endif
    END IF
          
    write(IUNIT,'(/A,I10)') "CELL_TYPES", num_Crack_Elements
    DO i_E=1,num_Crack_Elements
        write(IUNIT,'(I3)') 5
    END DO

    write(IUNIT,'(/A,I10)') "CELL_DATA", num_Crack_Elements

    WRITE(IUNIT,'(/A,/A)')'SCALARS Crack_ID int','LOOKUP_TABLE default'
    do i_E=1,num_Crack_Elements
        write(IUNIT,'(I8)') Global_Crack_Element_Crack_ID(i_E)
    enddo

    WRITE(IUNIT,'(/A,/A)')'SCALARS Crack_Element_ID int','LOOKUP_TABLE default'      
    do i_E=1,num_Crack_Elements
        write(IUNIT,'(I8)') i_E
    enddo
          
    write(IUNIT,'(/A,I10)') "POINT_DATA", num_Crack_Points

    WRITE(IUNIT,'(/A,/A)')'SCALARS Crack_Node_Number integer','LOOKUP_TABLE default'
    do i_Node=1,num_Crack_Points
        write(IUNIT,'(I8)') i_Node
    enddo

    write(IUNIT,'(/A)') "SCALARS Crack_Node_Aperture double"
    write(IUNIT,'(A)')  "LOOKUP_TABLE default"
    if(Key_Dimension==3 .and. allocated(Crack3D_Meshed_Node_Value))then
        do i_Node=1,num_Crack_Points
            c_crack_number = Global_Crack_Node_Crack_ID(i_Node,1) 
            c_local_node   = Global_Crack_Node_Crack_ID(i_Node,2) 
            write(IUNIT,'(E20.12)') Crack3D_Meshed_Node_Value(c_crack_number)%row(c_local_node,1)
        END DO   
    endif

    
    close(IUNIT)      
endif


RETURN
END SUBROUTINE Save_vtk_file_for_Crack
