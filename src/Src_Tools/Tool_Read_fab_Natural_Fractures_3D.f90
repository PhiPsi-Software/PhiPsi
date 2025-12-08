 
SUBROUTINE Tool_Read_fab_Natural_Fractures_3D

use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_Filename

implicit none
character(len=127) c_File_name
character(len=127) c_line_string,c_m_line_string
character(len=127),ALLOCATABLE::Input_Info(:)
integer stat_read,ls1,ls2,i_Char,num_Input_Info,i_Info
integer num_Data
logical Yes_Even,alive
integer tem_num_Crack

integer,ALLOCATABLE::tem_Crack_Vertex(:)
real(kind=FT),ALLOCATABLE::tem_Crack3D_Coor(:,:,:)
real(kind=FT),ALLOCATABLE::tem_Crack3D_n_Vector(:,:)
real(kind=FT) tem_Data(100)
integer j_Info,c_Info
integer c_len_trim,i_V,i_C
real(kind=FT) M_x_min,M_x_max,M_y_min,M_y_max,M_z_min,M_z_max 
logical Yes_Points_In_Model
real(kind=FT) c_Point(3)
integer c_num_Crack,max_Vextex
integer num_possibe_cracks

print *, "    Reading FracMan natural fracture *.fab file..."

c_num_Crack = 0

c_File_name = trim(trim(Full_Pathname)//'.fab')
inquire(file=c_File_name, exist=alive)
if(alive.EQV..FALSE.)then
  print *, "    Error(2022072403) :: Can not find *.fab file!"    
  call Warning_Message('S',Keywords_Blank) 
endif

open(1102,file=c_File_name,status='old')
num_Input_Info = 0
do
  read(1102,"(A127)",iostat=stat_read) c_line_string
  if(stat_read/=0) exit
  if (len_trim(c_line_string) /= 0)then

      c_line_string = adjustl(c_line_string)
      c_line_string = trim(c_line_string)
      ls1 = len_trim(c_line_string)
      ls2 = 0
      do i_Char = 1,ls1
        if(c_line_string(i_Char:i_Char).ne.' ') then
          ls2 = ls2 + 1
          c_m_line_string(ls2:ls2)=c_line_string(i_Char:i_Char)
        endif
      enddo
      num_Input_Info = num_Input_Info +1
  endif
end do
print *,'    Number of non-empty lines of *.fab file:      ',num_Input_Info
close(1102)

allocate(Input_Info(num_Input_Info))
open(1502,file=c_File_name,status='old')
num_Input_Info = 0
do
  read(1502,"(A127)",iostat=stat_read) c_line_string
  if(stat_read/=0) exit
  if (len_trim(c_line_string) /= 0)then

      c_line_string = adjustl(c_line_string)
      c_line_string = trim(c_line_string)
      ls1 = len_trim(c_line_string)
      ls2 = 0
      do i_Char = 1,ls1
        if(c_line_string(i_Char:i_Char).ne.' ') then
          ls2 = ls2 + 1
          c_m_line_string(ls2:ls2)=c_line_string(i_Char:i_Char)
        endif
      enddo
      num_Input_Info = num_Input_Info +1
      Input_Info(num_Input_Info)=c_line_string
  endif
end do
close(1502)


num_possibe_cracks = int(num_Input_Info/6)
allocate(tem_Crack_Vertex(num_possibe_cracks))
allocate(tem_Crack3D_Coor(num_possibe_cracks,10,3))
allocate(tem_Crack3D_n_Vector(num_possibe_cracks,3))


tem_num_Crack = 0
tem_Crack_Vertex(1:num_possibe_cracks)  = 0
tem_Crack3D_Coor(1:num_possibe_cracks,1:10,1:3) = ZR 
tem_Crack3D_n_Vector(1:num_possibe_cracks,1:3)  = ZR
if (Key_Fracture_Zone==0)then
  M_x_min = Min_X_Coor
  M_x_max = Max_X_Coor
  M_y_min = Min_Y_Coor
  M_y_max = Max_Y_Coor
  M_z_min = Min_Z_Coor
  M_z_max = Max_Z_Coor
elseif(Key_Fracture_Zone==1)then
  M_x_min = Frac_Zone_Minx
  M_x_max = Frac_Zone_Maxx
  M_y_min = Frac_Zone_Miny
  M_y_max = Frac_Zone_Maxy
  M_z_min = Frac_Zone_Minz
  M_z_max = Frac_Zone_Maxz     
endif

do i_Info=1,num_Input_Info-1
  c_info = i_Info
  if (Input_Info(i_Info)(1:14)=='BEGIN FRACTURE')then
      
      
      100 continue
      c_info = c_info +1
      
      if(c_info>=num_Input_Info-1) then
          EXIT
      endif
      
      
      if (Input_Info(c_info)(1:12)=='END FRACTURE')then
          exit
      endif   
      

      c_len_trim = len_trim(Input_Info(c_info))
      
      call Tool_Get_num_Data_from_a_String_v2(Input_Info(c_info),c_len_trim,num_Data,Yes_Even)     
      
      Read(Input_Info(c_info) ,*) tem_Data(1:num_Data)
      
      tem_num_Crack = tem_num_Crack + 1
      tem_Crack_Vertex(tem_num_Crack) =  int(tem_Data(2))
      if(int(tem_Data(2))>10)then
          print *, "    ERROR(2022072401) :: vertex number>10!"    
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      do j_Info = 1,int(tem_Data(2))
        c_info = c_info + 1
        Read(Input_Info(c_info) ,*) tem_Data(1:4)
        tem_Crack3D_Coor(tem_num_Crack,j_Info,1:3)=tem_Data(2:4)
      enddo
      
      c_info = c_info + 1
      Read(Input_Info(c_info) ,*) tem_Data(1:4)
      tem_Crack3D_n_Vector(tem_num_Crack,1:3)=tem_Data(2:4)
      
      goto 100
      
  endif       
enddo

print *,'    Number of cracks in FracMan *.fab file:       ',tem_num_Crack


num_Rand_Na_Crack = 0
max_Vextex = -1
do i_C=1,tem_num_Crack
  Yes_Points_In_Model = .True.
  do i_V = 1,tem_Crack_Vertex(i_C) 
      c_Point =  tem_Crack3D_Coor(i_C,i_V,1:3)  
      if(c_Point(1) >= M_x_min .and. c_Point(1) <= M_x_max .and. &
         c_Point(2) >= M_y_min .and. c_Point(2) <= M_y_max .and. &
         c_Point(3) >= M_z_min .and. c_Point(3) <= M_z_max) then
         
      else
          Yes_Points_In_Model = .False.
          exit
      endif
  enddo
  if(Yes_Points_In_Model) then
      num_Rand_Na_Crack= num_Rand_Na_Crack+ 1
      if(tem_Crack_Vertex(i_C) > max_Vextex) then
          max_Vextex = tem_Crack_Vertex(i_C) 
      endif
  endif
enddo

allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,max_Vextex,3))   
allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  
num_Rand_Na_Crack = 0
do i_C=1,tem_num_Crack
  Yes_Points_In_Model = .True.
  do i_V = 1,tem_Crack_Vertex(i_C) 
      c_Point =  tem_Crack3D_Coor(i_C,i_V,1:3)  
      if(c_Point(1) >= M_x_min .and. c_Point(1) <= M_x_max .and. &
         c_Point(2) >= M_y_min .and. c_Point(2) <= M_y_max .and. &
         c_Point(3) >= M_z_min .and. c_Point(3) <= M_z_max) then
      else
          Yes_Points_In_Model = .False.
          exit
      endif
  enddo
  if(Yes_Points_In_Model) then
      num_Rand_Na_Crack= num_Rand_Na_Crack+ 1
      Each_NaCr3D_Poi_Num(num_Rand_Na_Crack)   = tem_Crack_Vertex(i_C)
      do i_V = 1,tem_Crack_Vertex(i_C) 
          Na_Crack3D_Coor(num_Rand_Na_Crack,i_V,1:3) = tem_Crack3D_Coor(i_C,i_V,1:3)  
      enddo
  endif
enddo

call Save_Files_NF_3D(1,1)


if(Key_NaCr_Active_Scheme_3D == 1) then
  do i_C = 1,num_Rand_Na_Crack
    Crack3D_Coor(num_Crack+i_C,1:Each_NaCr3D_Poi_Num(i_C) ,1:3) = Na_Crack3D_Coor(i_C,1:Each_NaCr3D_Poi_Num(i_C),1:3)
    Crack_Type_Status_3D(num_Crack+i_C,1) = 2
    Key_CS_Crack(num_Crack+i_C) = 1
    if(Key_NaCr_Growth==0)then
        Crack_Type_Status_3D(num_Crack+i_C,3) = 0
    elseif(Key_NaCr_Growth==1) then
        Crack_Type_Status_3D(num_Crack+i_C,3) = 1
    endif
    Each_Cr_Poi_Num(num_Crack+i_C)   = Each_NaCr3D_Poi_Num(i_C)   
  enddo

  num_Crack = num_Crack + num_Rand_Na_Crack
  NaCr3D_Status(1:num_Rand_Na_Crack,1) = 1
  
elseif(Key_NaCr_Active_Scheme_3D == 2) then

elseif(Key_NaCr_Active_Scheme_3D == 3) then  
  
  allocate(Na_Crack3D_Ele_List(num_Rand_Na_Crack))
endif  

print *,'    Number of imported natural cracks from *.fab: ',  num_Rand_Na_Crack

if(allocated(Input_Info)) deallocate(Input_Info)
if(allocated(tem_Crack_Vertex)) deallocate(tem_Crack_Vertex)
if(allocated(tem_Crack3D_Coor)) deallocate(tem_Crack3D_Coor)
if(allocated(tem_Crack3D_n_Vector)) deallocate(tem_Crack3D_n_Vector)

RETURN
END SUBROUTINE Tool_Read_fab_Natural_Fractures_3D