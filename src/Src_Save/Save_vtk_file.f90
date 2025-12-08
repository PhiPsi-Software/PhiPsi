 
SUBROUTINE Save_vtk_file(isub)

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
integer IUNIT,i_Node,ind,i_E,n1,npElem
integer i_N
      
if (Key_Save_Nothing==1) return

if (Key_Save_vtk/=1) then
  return
endif

print *,'    Saving vtk file...'


IUNIT = 101
DATFIL= trim(Full_Pathname)
VTKFIL=''


I=len_trim(DATFIL)+1

VTKFIL(1:I)=DATFIL(1:I)
WRITE(TEMP,'(I5.5)')isub
VTKFIL(I:I+9)='_'//TRIM(TEMP)//'.vtk'

open(IUNIT,file=VTKFIL,status='unknown')
WRITE(IUNIT,'(A)')'# vtk DataFile Version 4.0'
WRITE(IUNIT,'(A)')DATFIL(1:I)//' - results from increment '//TRIM(TEMP)

WRITE(IUNIT,'(A)')'ASCII'
WRITE(IUNIT,'(A)')'DATASET UNSTRUCTURED_GRID'
      

WRITE(IUNIT,'(/A,I0,A)')'POINTS ',Num_Node,' double'
if(Key_Dimension==2)then
  do i_Node=1,num_Node
      WRITE(IUNIT,'(3F12.6)') Coor(i_Node,1:2),ZR
  enddo
elseif(Key_Dimension==3)then
  do i_Node=1,num_Node
      WRITE(IUNIT,'(3F12.6)') Coor(i_Node,1:3)
  enddo      
endif


if(Key_Dimension==2)then
  npElem = 4
elseif(Key_Dimension==3)then
  npElem = 8     
endif
      
ind = Num_Elem*(npElem+1)
WRITE(IUNIT,'(/A,I0,A,I0)')'CELLS ',Num_Elem,' ',ind

IF(Key_Dimension == 2) THEN
    IF(npElem == 3) THEN
    n1 = 5
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1
    END DO
    ELSE IF(npElem == 6) THEN
    n1 = 22
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1,&
                                                            Elem_Node(i_E,4)-1,Elem_Node(i_E,5)-1,Elem_Node(i_E,6)-1
    END DO
    ELSE IF(npElem == 4) THEN
    n1 = 9
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1,Elem_Node(i_E,4)-1
    END DO
    END IF
    ELSEIF(Key_Dimension == 3) then
    IF(npElem == 4) THEN
    n1 = 10
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1,Elem_Node(i_E,4)-1
    END DO
    ELSE IF(npElem == 6) THEN
    n1 = 13
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1,&
                                                            Elem_Node(i_E,4)-1,Elem_Node(i_E,5)-1,Elem_Node(i_E,6)-1
    END DO
    ELSE IF(npElem == 8) THEN
    n1 = 12
    DO i_E=1,Num_Elem
        write(IUNIT,'(I10,I10,I10,I10,I10,I10,I10,I10,I10)') npElem,Elem_Node(i_E,1)-1,Elem_Node(i_E,2)-1,Elem_Node(i_E,3)-1,&
                Elem_Node(i_E,4)-1,Elem_Node(i_E,5)-1,Elem_Node(i_E,6)-1,&
                Elem_Node(i_E,7)-1,Elem_Node(i_E,8)-1     
    END DO          
    END IF
END IF
      
write(IUNIT,'(/A,I10)') "CELL_TYPES", Num_Elem
DO i_E=1,Num_Elem
    write(IUNIT,'(I3)') n1
END DO

write(IUNIT,'(/A,I10)') "CELL_DATA", Num_Elem

WRITE(IUNIT,'(/A,/A)')'SCALARS Material_ID int','LOOKUP_TABLE default'
do i_E=1,Num_Elem
      write(IUNIT,'(I8)') Elem_Mat(i_E)
enddo

WRITE(IUNIT,'(/A,/A)')'SCALARS Element_ID int','LOOKUP_TABLE default'      
do i_E=1,Num_Elem
      write(IUNIT,'(I8)') i_E
enddo
      
if(allocated(Ele_Permeability_3D))then
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_xx double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,1)
  END DO
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_yy double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,2)
  END DO
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_zz double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,3)
  END DO
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_xy double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,4)
  END DO
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_yz double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,5)
  END DO
  write(IUNIT,'(/A,/A)') 'SCALARS Element_Permeability_xz double','LOOKUP_TABLE default'
  do i_E=1,Num_Elem
      write(IUNIT,'(E20.12)') Ele_Permeability_3D(i_E,6)
  END DO          
endif
      
      
if (Enrich_Freedom >0)then
  WRITE(IUNIT,'(/A)')'SCALARS Enriched_Element_Type int'
  WRITE(IUNIT,'(/A)')'LOOKUP_TABLE default'   
  if (num_crack >=1 .and. Key_Dimension == 2) then
      do i_E=1,Num_Elem
          if (any(Elem_Type(i_E,1:num_crack) == 1))then
              write(IUNIT,'(I1)') 1
          elseif (any(Elem_Type(i_E,1:num_crack)== 2))then
              write(IUNIT,'(I1)') 2
          elseif (any(Elem_Type(i_E,1:num_crack)== 3))then
              write(IUNIT,'(I1)') 3
          elseif (any(Elem_Type(i_E,1:num_crack)== 4))then
              write(IUNIT,'(I1)') 4
          elseif (any(Elem_Type(i_E,1:num_crack)== 5))then
              write(IUNIT,'(I1)') 5
          elseif (any(Elem_Type(i_E,1:num_crack)== 6))then
              write(IUNIT,'(I1)') 6 
          else
              write(IUNIT,'(I1)') 0
          endif
      enddo  
  endif
  if (num_crack >=1 .and. Key_Dimension == 3) then
      do i_E=1,Num_Elem
          if (any(Elem_Type_3D(i_E,1:num_crack) == 1))then
              write(IUNIT,'(I1)') 1
          elseif (any(Elem_Type_3D(i_E,1:num_crack)== 2))then
              write(IUNIT,'(I1)') 2
          else
              write(IUNIT,'(I1)') 0
          endif
      enddo  
  endif          
endif
      
      
write(IUNIT,'(/A,I10)') "POINT_DATA", Num_Node

WRITE(IUNIT,'(/A)')'VECTORS Displacement double'   
if (Key_Dimension == 2) then   
  do i_Node=1,num_Node
      write(IUNIT,'(3E20.12)')DISP(2*i_Node-1),DISP(2*i_Node),ZR
  END DO        
elseif(Key_Dimension == 3)then
  do i_Node=1,num_Node
      write(IUNIT,'(3E20.12)') DISP(3*i_Node-2),DISP(3*i_Node-1),DISP(3*i_Node)
  END DO         
endif      

if (Key_Dimension == 2 .and. allocated(Stress_xx_Node)) then   
  write(IUNIT,'(/A)') "SCALARS stress_xx double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_xx_Node(i_Node)
  END DO
  write(IUNIT,'(/A)') "SCALARS stress_yy double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_yy_Node(i_Node)
  END DO  
  write(IUNIT,'(/A)') "SCALARS stress_xy double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_xy_Node(i_Node)
  END DO    
  write(IUNIT,'(/A)') "SCALARS stress_vm double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_vm_Node(i_Node)
  END DO          
elseif(Key_Dimension == 3 .and. allocated(Stress_xx_Node) )then
  write(IUNIT,'(/A)') "SCALARS stress_xx double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_xx_Node(i_Node)
  END DO
  write(IUNIT,'(/A)') "SCALARS stress_yy double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_yy_Node(i_Node)
  END DO  
  write(IUNIT,'(/A)') "SCALARS stress_zz double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_zz_Node(i_Node)
  END DO            
  write(IUNIT,'(/A)') "SCALARS stress_xy double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_xy_Node(i_Node)
  END DO   
  write(IUNIT,'(/A)') "SCALARS stress_yz double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_yz_Node(i_Node)
  END DO     
  write(IUNIT,'(/A)') "SCALARS stress_xz double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_xz_Node(i_Node)
  END DO           
  write(IUNIT,'(/A)') "SCALARS stress_vm double"
  write(IUNIT,'(A)') "LOOKUP_TABLE default"
  do i_Node=1,num_Node
      write(IUNIT,'(E20.12)') Stress_vm_Node(i_Node)
  END DO               
endif     
      
WRITE(IUNIT,'(/A,/A)')'SCALARS Node_Number integer','LOOKUP_TABLE default'
do i_Node=1,num_Node
  write(IUNIT,'(I8)') i_Node
enddo

if (Enrich_Freedom > 0 )then
  WRITE(IUNIT,'(/A,/A)')'SCALARS Enriched_Node_Type int','LOOKUP_TABLE default'   
  if (num_crack >=1 .and. Key_Dimension == 2) then
    do i_N=1,num_Node
      if (any(Enriched_Node_Type(i_N,1:num_crack)== 1))then
          write(IUNIT,'(I1)') 1
      elseif (any(Enriched_Node_Type(i_N,1:num_crack)==2))then
          write(IUNIT,'(I1)') 2
      elseif (any(Enriched_Node_Type(i_N,1:num_crack)==3))then
          write(IUNIT,'(I1)') 3
      elseif (any(Enriched_Node_Type(i_N,1:num_crack)==6))then
          write(IUNIT,'(I1)') 6 
      else
          write(IUNIT,'(I1)') 0
      endif
    enddo  
  endif
  if (num_crack >=1 .and. Key_Dimension == 3) then
    do i_N=1,num_Node
      if (any(Enriched_Node_Type_3D(i_N,1:num_crack)== 1))then
          write(IUNIT,'(I1)') 1
      elseif (any(Enriched_Node_Type_3D(i_N,1:num_crack)==2))then
          write(IUNIT,'(I1)') 2
      else
          write(IUNIT,'(I1)') 0
      endif
    enddo  
  endif          
endif

close(IUNIT)      



RETURN
END SUBROUTINE Save_vtk_file
