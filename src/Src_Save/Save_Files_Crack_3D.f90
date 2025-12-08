 
SUBROUTINE Save_Files_Crack_3D(isub)
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Filename
use Global_POST

implicit none
integer isub
integer i,j
character(200) Filename_1,Filename_2,Filename_3
character(200) Filename_4,Filename_5,Filename_6
character(200) Filename_7,Filename_8,Filename_9
character(5) temp    
real(kind=FT) tem_1
integer i_E,i_C,i_Node
logical c_write
integer c_Cr_Location
real(kind=FT) tem_crack_radius
integer :: temp_Solid_El_Crs(Solid_El_Max_num_Crs)

if (Key_Save_Nothing==1) return

tem_1 = 1000.0D0

print *,'    Saving coordinates of cracks...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.crax'//'_'//ADJUSTL(temp)  
Filename_2=trim(Full_Pathname)//'.cray'//'_'//ADJUSTL(temp)      
Filename_3=trim(Full_Pathname)//'.craz'//'_'//ADJUSTL(temp)   
Filename_4=trim(Full_Pathname)//'.cnox'//'_'//ADJUSTL(temp)      
Filename_5=trim(Full_Pathname)//'.cnoy'//'_'//ADJUSTL(temp)   
Filename_6=trim(Full_Pathname)//'.cnoz'//'_'//ADJUSTL(temp)        
open(101,file=Filename_1,status='unknown')    
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown') 
open(104,file=Filename_4,status='unknown') 
open(105,file=Filename_5,status='unknown') 
open(106,file=Filename_6,status='unknown')
      
do i=1,num_Crack
      write(101, '(2000E20.12)') (Crack3D_Coor(i,j,1),j=1,Each_Cr_Poi_Num(i))
      write(102, '(2000E20.12)') (Crack3D_Coor(i,j,2),j=1,Each_Cr_Poi_Num(i))
      write(103, '(2000E20.12)') (Crack3D_Coor(i,j,3),j=1,Each_Cr_Poi_Num(i))
end do
close(101)
close(102)   
close(103)

do i=1,num_Crack
  write(104, '(50000E20.12)') (Crack3D_Meshed_Node(i)%row(j,1),j=1,Crack3D_Meshed_Node_num(i))       
  write(105, '(50000E20.12)') (Crack3D_Meshed_Node(i)%row(j,2),j=1,Crack3D_Meshed_Node_num(i))   
  write(106, '(50000E20.12)') (Crack3D_Meshed_Node(i)%row(j,3),j=1,Crack3D_Meshed_Node_num(i))        
end do
close(104)          
close(105)    
close(106)    

print *,'    Saving apertures of crack nodes...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.cmap'//'_'//ADJUSTL(temp)    
open(101,file=Filename_1,status='unknown') 
do i=1,num_Crack
write(101,'(50000E20.12)')(Crack3D_Meshed_Node_Value(i)%row(j,1),j=1,Crack3D_Meshed_Node_num(i))
end do
close(101) 

print *,'    Saving mesh of cracks...'
write(temp,'(I5)') isub
Filename_1   =  trim(Full_Pathname)//'.cms1_'//ADJUSTL(temp)  
Filename_2   =  trim(Full_Pathname)//'.cms2_'//ADJUSTL(temp)      
Filename_3   =  trim(Full_Pathname)//'.cms3_'//ADJUSTL(temp)         
open(101,file=Filename_1,status='unknown')    
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown')     
do i=1,num_Crack
  write(101, '(50000I10)') (Crack3D_Meshed_Ele(i)%row(j,1),j=1,Crack3D_Meshed_Ele_num(i))
  write(102, '(50000I10)') (Crack3D_Meshed_Ele(i)%row(j,2),j=1,Crack3D_Meshed_Ele_num(i))
  write(103, '(50000I10)') (Crack3D_Meshed_Ele(i)%row(j,3),j=1,Crack3D_Meshed_Ele_num(i))
end do
close(101)
close(102)   
close(103)  


print *,'    Saving ennd file...'
Filename_1  = trim(Full_Pathname)//'.ennd'//'_'//ADJUSTL(temp)      
if (Key_Data_Format==1) then
   open(103,file=Filename_1,status='unknown') 
elseif(Key_Data_Format==2)then
   open(103,file=Filename_1,status='unknown',form='unformatted',access='stream') 
endif
do i=1,Num_Node
  if (Key_Data_Format==1) then
      write(103, '(200I10)') (Enriched_Node_Type_3D(i,j),j=1,num_Crack)
  elseif(Key_Data_Format==2)then
      write(103) (Enriched_Node_Type_3D(i,j),j=1,num_Crack)
  endif     
end do
close(103)     



print *,'    Saving elty file...'
Filename_1  = trim(Full_Pathname)//'.elty'//'_'//ADJUSTL(temp)      
if (Key_Data_Format==1) then
   open(104,file=Filename_1,status='unknown')    
elseif(Key_Data_Format==2)then
   open(104,file=Filename_1,status='unknown', form='unformatted',access='stream') 
endif 
do i=1,Num_Elem
  if (Key_Data_Format==1) then
      write(104, '(200I10)') (Elem_Type_3D(i,j),j=1,num_Crack)
  elseif(Key_Data_Format==2)then
      write(104) (Elem_Type_3D(i,j),j=1,num_Crack)
  endif                  
end do
close(104)   



print *,'    Saving posi file...'
Filename_1  = trim(Full_Pathname)//'.posi'//'_'//ADJUSTL(temp)    
if (Key_Data_Format==1) then
   open(110,file=Filename_1,status='unknown')  
elseif(Key_Data_Format==2)then
   open(110,file=Filename_1,status='unknown',form='unformatted',access='stream') 
endif   
do i=1,Num_Node
  if (Key_Data_Format==1) then
      write(110, '(200I10)') (c_POS_3D(i,j),j=1,num_Crack)
  elseif(Key_Data_Format==2)then
      write(110) (c_POS_3D(i,j),j=1,num_Crack)
  endif            
end do
close(110) 

print *,'    Saving number of fluid elements of cracks...'
write(temp,'(I5)') isub
Filename_1 =trim(Full_Pathname)//'.cpfn'//'_'//ADJUSTL(temp)         
open(101,file=Filename_1,status='unknown')     
do i=1,num_Crack
  write(101, '(I10)') Cracks_FluidEle_num_3D(i)  
end do
close(101)

print *,'    Saving fluid node number of cracks...'
write(temp,'(I5)') isub
Filename_1 = trim(Full_Pathname)//'.cpno'//'_'//ADJUSTL(temp)        
if (Key_Data_Format==1) then
open(101,file=Filename_1,status='unknown')     
do i=1,num_Crack
  do j=1,Cracks_FluidEle_num_3D(i)
      write(101, '(3I10)') Cracks_FluidEle_CalP_3D(i)%row(j,1:3)
  enddo
end do
close(101)
elseif(Key_Data_Format==2)then
open(101,file=Filename_1,status='unknown',form='unformatted',access='stream')     
do i=1,num_Crack
  do j=1,Cracks_FluidEle_num_3D(i)
      write(101) Cracks_FluidEle_CalP_3D(i)%row(j,1:3)
  enddo
end do
close(101)          
endif



print *,'    Saving fluid node coordinates of cracks...'
Filename_4   = trim(Full_Pathname)//'.ccpx_'//ADJUSTL(temp)     
Filename_5   = trim(Full_Pathname)//'.ccpy_'//ADJUSTL(temp)    
Filename_6   = trim(Full_Pathname)//'.ccpz_'//ADJUSTL(temp)  
open(701,file=Filename_4,status='unknown') 
open(702,file=Filename_5,status='unknown') 
open(703,file=Filename_6,status='unknown')            
do i=1,num_Crack
  write(701, '(50000E20.12)') (Cracks_CalP_Coors_3D(i)%row(j,1),j=1,Cracks_CalP_Num_3D(i))       
  write(702, '(50000E20.12)') (Cracks_CalP_Coors_3D(i)%row(j,2),j=1,Cracks_CalP_Num_3D(i))   
  write(703, '(50000E20.12)') (Cracks_CalP_Coors_3D(i)%row(j,3),j=1,Cracks_CalP_Num_3D(i))        
end do
close(701)          
close(702)    
close(703)     


print *,'    Saving vertex of cracks...'
write(temp,'(I5)') isub
Filename_1  =trim(Full_Pathname)//'.cmso'//'_'//ADJUSTL(temp)         
open(101,file=Filename_1,status='unknown')     
do i=1,num_Crack
  write(101, '(50000I10)') (Crack3D_Meshed_Outline(i)%row(j,1),j=1,Crack3D_Meshed_Outline_num(i))
end do
close(101)     

if(Key_Simple_Post==1) then
  return
endif

print *,'    Saving element number containing crack node...'
write(temp,'(I5)') isub
Filename_1 =trim(Full_Pathname)//'.cmse'//'_'//ADJUSTL(temp)         
open(101,file=Filename_1,status='unknown')    
do i=1,num_Crack
  write(101, '(50000I10)')  (Cr3D_Meshed_Node_in_Ele_Num(i)%row(j),j=1,Crack3D_Meshed_Node_num(i))
end do
close(101)          

print *,'    Saving local coordinates of crack node...'
Filename_1  = trim(Full_Pathname)//'.cnlx'//'_'//ADJUSTL(temp)      
Filename_2  = trim(Full_Pathname)//'.cnly'//'_'//ADJUSTL(temp)   
Filename_3  = trim(Full_Pathname)//'.cnlz'//'_'//ADJUSTL(temp)        
open(101,file=Filename_1,status='unknown')  
do i=1,num_Crack
  write(101, '(50000E20.12)')(Cr3D_Meshed_Node_in_Ele_Local(i)%row(j,1),j=1,Crack3D_Meshed_Node_num(i))              
end do
close(101)
open(102,file=Filename_2,status='unknown') 
do i=1,num_Crack    
  write(102, '(50000E20.12)')(Cr3D_Meshed_Node_in_Ele_Local(i)%row(j,2),j=1,Crack3D_Meshed_Node_num(i))         
end do
close(102)  
open(103,file=Filename_3,status='unknown')  
do i=1,num_Crack  
  write(103, '(50000E20.12)')(Cr3D_Meshed_Node_in_Ele_Local(i)%row(j,3),j=1,Crack3D_Meshed_Node_num(i))        
end do
close(103) 

print *,'    Saving *.cvxx and other files...'
Filename_1=trim(Full_Pathname)//'.cvxx_'//ADJUSTL(temp)     
Filename_2=trim(Full_Pathname)//'.cvxy_'//ADJUSTL(temp)    
Filename_3=trim(Full_Pathname)//'.cvxz_'//ADJUSTL(temp)  
Filename_4=trim(Full_Pathname)//'.cvyx_'//ADJUSTL(temp)     
Filename_5=trim(Full_Pathname)//'.cvyy_'//ADJUSTL(temp)    
Filename_6=trim(Full_Pathname)//'.cvyz_'//ADJUSTL(temp)  
Filename_7=trim(Full_Pathname)//'.cvzx_'//ADJUSTL(temp)     
Filename_8=trim(Full_Pathname)//'.cvzy_'//ADJUSTL(temp)    
Filename_9=trim(Full_Pathname)//'.cvzz_'//ADJUSTL(temp)            
open(701,file=Filename_1,status='unknown') 
do i=1,num_Crack
write(701, '(50000E20.12)') (Crack3D_Meshed_Vertex_x_Vector(i)%row(j,1),j=1,Crack3D_Meshed_Outline_num(i))        
end do
close(701)             
open(702,file=Filename_2,status='unknown') 
do i=1,num_Crack     
write(702, '(50000E20.12)') (Crack3D_Meshed_Vertex_x_Vector(i)%row(j,2),j=1,Crack3D_Meshed_Outline_num(i))     
end do          
close(702)  
open(703,file=Filename_3,status='unknown')  
do i=1,num_Crack
write(703, '(50000E20.12)') (Crack3D_Meshed_Vertex_x_Vector(i)%row(j,3),j=1,Crack3D_Meshed_Outline_num(i))   
end do
close(703)          
open(704,file=Filename_4,status='unknown') 
do i=1,num_Crack
write(704, '(50000E20.12)') (Crack3D_Meshed_Vertex_y_Vector(i)%row(j,1),j=1,Crack3D_Meshed_Outline_num(i))            
end do
close(704)  
open(705,file=Filename_5,status='unknown') 
do i=1,num_Crack     
write(705, '(50000E20.12)') (Crack3D_Meshed_Vertex_y_Vector(i)%row(j,2),j=1,Crack3D_Meshed_Outline_num(i))      
end do  
close(705)   
open(706,file=Filename_6,status='unknown')   
do i=1,num_Crack     
write(706, '(50000E20.12)') (Crack3D_Meshed_Vertex_y_Vector(i)%row(j,3),j=1,Crack3D_Meshed_Outline_num(i))   
end do          
close(706) 
open(707,file=Filename_7,status='unknown') 
do i=1,num_Crack     
write(707, '(50000E20.12)') (Crack3D_Meshed_Vertex_z_Vector(i)%row(j,1),j=1,Crack3D_Meshed_Outline_num(i))           
end do          
close(707) 
open(708,file=Filename_8,status='unknown') 
do i=1,num_Crack       
write(708, '(50000E20.12)') (Crack3D_Meshed_Vertex_z_Vector(i)%row(j,2),j=1,Crack3D_Meshed_Outline_num(i))       
end do          
close(708) 
open(709,file=Filename_9,status='unknown') 
do i=1,num_Crack     
write(709, '(50000E20.12)') (Crack3D_Meshed_Vertex_z_Vector(i)%row(j,3),j=1,Crack3D_Meshed_Outline_num(i))      
end do
close(709)     

print *,'    Saving *.blab file...'
Filename_1=trim(Full_Pathname)//'.blab_'//ADJUSTL(temp)
if (Key_Data_Format==1) then
  open(701,file=Filename_1,status='unknown')   
elseif(Key_Data_Format==2)then
  open(701,file=Filename_1,status='unknown',form='unformatted',access='stream')
endif
do i_E=1,Num_Elem
  c_write = .False.
  do i_C=1,num_Crack
      temp_Solid_El_Crs = Solid_El_Crs(i_E, 1:Solid_El_Max_num_Crs)
      call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_Solid_El_Crs, i_C, c_Cr_Location)
      if(c_Cr_Location>=1)then
          if (Key_Data_Format==1) then
             if(allocated(Solid_El_Tip_BaseLine(i_E)%row)) then
                 write(701, '(6E20.12)') Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3), &
                                         Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3)
             else
                 write(701, '(6E20.12)') ZR,ZR,ZR,ZR,ZR,ZR
             endif
          elseif(Key_Data_Format==2)then
             if(allocated(Solid_El_Tip_BaseLine(i_E)%row)) then
               write(701) Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3),&
                          Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3)
             else
                 write(701) ZR,ZR,ZR,ZR,ZR,ZR
             endif                       
          endif

          c_write = .True.
          exit         
      endif
  end do
  if(c_write .eqv. .False.) then
      if (Key_Data_Format==1) then
          write(701, '(6E20.12)') ZR,ZR,ZR,ZR,ZR,ZR   
      elseif(Key_Data_Format==2)then
          write(701) ZR,ZR,ZR,ZR,ZR,ZR   
      endif                  
  endif
end do
close(701)     

print *,'    Saving *.tere file...'
Filename_1=trim(Full_Pathname)//'.tere_'//ADJUSTL(temp)
if (Key_Data_Format==1) then
  open(701,file=Filename_1,status='unknown')   
  do i_Node=1,Num_Node
      c_write = .False.
      do i_C=1,num_Crack
          if(allocated(Ele_Num_Tip_Enriched_Node_3D(i_Node)%row)) then
              if (Ele_Num_Tip_Enriched_Node_3D(i_Node)%row(i_C)>=1)then
                  write(701, '(2I10)') Ele_Num_Tip_Enriched_Node_3D(i_Node)%row(i_C),i_C  
                  c_write = .True.
                  exit
              endif 
          endif
      end do
      if(c_write .eqv. .False.) then
          write(701, '(2I10)') 0,0  
      endif
  end do
  close(701)   
elseif(Key_Data_Format==2)then
  open(701,file=Filename_1,status='unknown',form='unformatted',access='stream') 
  do i_Node=1,Num_Node
      c_write = .False.
      do i_C=1,num_Crack
          if(allocated(Ele_Num_Tip_Enriched_Node_3D(i_Node)%row)) then
              if (Ele_Num_Tip_Enriched_Node_3D(i_Node)%row(i_C)>=1)then
                  write(701) Ele_Num_Tip_Enriched_Node_3D(i_Node)%row(i_C),i_C  
                  c_write = .True.
                  exit
              endif  
          endif
      end do
      if(c_write .eqv. .False.) then
          write(701) 0,0  
      endif
  end do
  close(701)  
endif

print *,'    Saving *.blvx and other files...'
Filename_1=trim(Full_Pathname)//'.blvx_'//ADJUSTL(temp)
Filename_2=trim(Full_Pathname)//'.blvy_'//ADJUSTL(temp)
Filename_3=trim(Full_Pathname)//'.blvz_'//ADJUSTL(temp)

if (Key_Data_Format==1) then
  open(701,file=Filename_1,status='unknown') 
  open(702,file=Filename_2,status='unknown')
  open(703,file=Filename_3,status='unknown') 
elseif(Key_Data_Format==2)then
  open(701,file=Filename_1,status='unknown',form='unformatted',access='stream') 
  open(702,file=Filename_2,status='unknown',form='unformatted',access='stream')
  open(703,file=Filename_3,status='unknown',form='unformatted',access='stream')
endif

do i_E=1,Num_Elem
c_write = .False.
do i_C=1,num_Crack     
  temp_Solid_El_Crs = Solid_El_Crs(i_E, 1:Solid_El_Max_num_Crs)
  call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_Solid_El_Crs, i_C, c_Cr_Location)
  
  if(c_Cr_Location>=1)then
        if (Key_Data_Format==1) then
          if(allocated(Solid_El_Tip_BaseLine_x_Vec(i_E)%row)) then                
              write(701, '(3E20.12)') Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(701, '(3E20.12)') ZR,ZR,ZR
          endif
        elseif(Key_Data_Format==2)then
          if(allocated(Solid_El_Tip_BaseLine_x_Vec(i_E)%row)) then    
              write(701)Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(701) ZR,ZR,ZR
          endif
        endif
        c_write = .True.
        exit   
  endif  
end do  
if(c_write .eqv. .False.) then
    if (Key_Data_Format==1) then
        write(701, '(3E20.12)') ZR,ZR,ZR 
    elseif(Key_Data_Format==2)then
        write(701) ZR,ZR,ZR 
    endif
endif
end do

do i_E=1,Num_Elem
c_write = .False.
do i_C=1,num_Crack      
  temp_Solid_El_Crs = Solid_El_Crs(i_E, 1:Solid_El_Max_num_Crs)
  call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_Solid_El_Crs, i_C, c_Cr_Location)
  
  if(c_Cr_Location>=1)then
        if (Key_Data_Format==1) then
          if(allocated(Solid_El_Tip_BaseLine_y_Vec(i_E)%row)) then                
              write(702, '(3E20.12)') Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(702, '(3E20.12)') ZR,ZR,ZR
          endif
        elseif(Key_Data_Format==2)then
          if(allocated(Solid_El_Tip_BaseLine_y_Vec(i_E)%row)) then    
              write(702)Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(702) ZR,ZR,ZR
          endif
        endif
        
        c_write = .True.
        exit     
  endif  
end do  
if(c_write .eqv. .False.) then
    if (Key_Data_Format==1) then
        write(702, '(3E20.12)') ZR,ZR,ZR 
    elseif(Key_Data_Format==2)then
        write(702) ZR,ZR,ZR 
    endif
endif
end do

do i_E=1,Num_Elem
c_write = .False.
do i_C=1,num_Crack            
  temp_Solid_El_Crs = Solid_El_Crs(i_E, 1:Solid_El_Max_num_Crs)
  call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_Solid_El_Crs, i_C, c_Cr_Location)

  if(c_Cr_Location>=1)then
        if (Key_Data_Format==1) then
          if(allocated(Solid_El_Tip_BaseLine_z_Vec(i_E)%row)) then                
              write(703, '(3E20.12)') Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(703, '(3E20.12)') ZR,ZR,ZR
          endif
        elseif(Key_Data_Format==2)then
          if(allocated(Solid_El_Tip_BaseLine_z_Vec(i_E)%row)) then    
              write(703)Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3)
          else
              write(703) ZR,ZR,ZR
          endif
        endif
        c_write = .True.
        exit   
          c_write = .True.
          exit   
  endif  
end do  
if(c_write .eqv. .False.) then
    if (Key_Data_Format==1) then
        write(703, '(3E20.12)') ZR,ZR,ZR 
    elseif(Key_Data_Format==2)then
        write(703) ZR,ZR,ZR 
    endif
endif
end do 

close(701)   
close(702)           
close(703) 

print *,'    Saving fluid element normal vector of cracks...'
Filename_1   = trim(Full_Pathname)//'.fenx_'//ADJUSTL(temp)     
Filename_2   = trim(Full_Pathname)//'.feny_'//ADJUSTL(temp)    
Filename_3   = trim(Full_Pathname)//'.fenz_'//ADJUSTL(temp)  
open(801,file=Filename_1,status='unknown') 
open(802,file=Filename_2,status='unknown') 
open(803,file=Filename_3,status='unknown')            
do i=1,num_Crack
  write(801, '(50000E20.12)') (Cracks_FluidEle_Vector_3D(i)%row(j,1),j=1,Cracks_FluidEle_num_3D(i))       
  write(802, '(50000E20.12)') (Cracks_FluidEle_Vector_3D(i)%row(j,2),j=1,Cracks_FluidEle_num_3D(i))   
  write(803, '(50000E20.12)') (Cracks_FluidEle_Vector_3D(i)%row(j,3),j=1,Cracks_FluidEle_num_3D(i))        
end do
close(801)          
close(802)    
close(803)    


print *,'    Saving up dis vectors of crack nodes...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.fnux_'//ADJUSTL(temp)    
Filename_2=trim(Full_Pathname)//'.fnuy_'//ADJUSTL(temp)   
Filename_3=trim(Full_Pathname)//'.fnuz_'//ADJUSTL(temp)        
open(101,file=Filename_1,status='unknown') 
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown') 
do i=1,num_Crack
write(101,'(50000E20.12)') (Cracks_CalP_UpDis_3D(i)%row(j,1),j=1,Cracks_CalP_Num_3D(i))
write(102,'(50000E20.12)') (Cracks_CalP_UpDis_3D(i)%row(j,2),j=1,Cracks_CalP_Num_3D(i))
write(103,'(50000E20.12)') (Cracks_CalP_UpDis_3D(i)%row(j,3),j=1,Cracks_CalP_Num_3D(i))     
end do
close(101) 
close(102)
close(103)


print *,'    Saving low dis vectors of crack nodes...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.fnlx_'//ADJUSTL(temp)    
Filename_2=trim(Full_Pathname)//'.fnly_'//ADJUSTL(temp)   
Filename_3=trim(Full_Pathname)//'.fnlz_'//ADJUSTL(temp)        
open(101,file=Filename_1,status='unknown') 
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown') 
do i=1,num_Crack
write(101,'(50000E20.12)') (Cracks_CalP_LowDis_3D(i)%row(j,1),j=1,Cracks_CalP_Num_3D(i))
write(102,'(50000E20.12)') (Cracks_CalP_LowDis_3D(i)%row(j,2),j=1,Cracks_CalP_Num_3D(i))
write(103,'(50000E20.12)') (Cracks_CalP_LowDis_3D(i)%row(j,3),j=1,Cracks_CalP_Num_3D(i))     
end do
close(101) 
close(102)
close(103)


print *,'    Saving local CS of fluid elements Gauss points of cracks...'
Filename_1=trim(Full_Pathname)//'.fexx_'//ADJUSTL(temp)     
Filename_2=trim(Full_Pathname)//'.fexy_'//ADJUSTL(temp)    
Filename_3=trim(Full_Pathname)//'.fexz_'//ADJUSTL(temp)  
open(801,file=Filename_1,status='unknown') 
open(802,file=Filename_2,status='unknown') 
open(803,file=Filename_3,status='unknown') 
do i=1,num_Crack
  write(801, '(50000E20.12)') (Cracks_FluidEle_LCS_x_3D(i)%row(j,1),j=1,Cracks_FluidEle_num_3D(i))       
  write(802, '(50000E20.12)') (Cracks_FluidEle_LCS_x_3D(i)%row(j,2),j=1,Cracks_FluidEle_num_3D(i))   
  write(803, '(50000E20.12)') (Cracks_FluidEle_LCS_x_3D(i)%row(j,3),j=1,Cracks_FluidEle_num_3D(i))        
end do
close(801)          
close(802)    
close(803)            
Filename_1=trim(Full_Pathname)//'.feyx_'//ADJUSTL(temp)     
Filename_2=trim(Full_Pathname)//'.feyy_'//ADJUSTL(temp)    
Filename_3=trim(Full_Pathname)//'.feyz_'//ADJUSTL(temp)  
open(801,file=Filename_1,status='unknown') 
open(802,file=Filename_2,status='unknown') 
open(803,file=Filename_3,status='unknown') 
do i=1,num_Crack
  write(801, '(50000E20.12)') (Cracks_FluidEle_LCS_y_3D(i)%row(j,1),j=1,Cracks_FluidEle_num_3D(i))       
  write(802, '(50000E20.12)') (Cracks_FluidEle_LCS_y_3D(i)%row(j,2),j=1,Cracks_FluidEle_num_3D(i))   
  write(803, '(50000E20.12)') (Cracks_FluidEle_LCS_y_3D(i)%row(j,3),j=1,Cracks_FluidEle_num_3D(i))        
end do
close(801)          
close(802)    
close(803)   
Filename_1=trim(Full_Pathname)//'.fezx_'//ADJUSTL(temp)     
Filename_2=trim(Full_Pathname)//'.fezy_'//ADJUSTL(temp)    
Filename_3=trim(Full_Pathname)//'.fezz_'//ADJUSTL(temp)  
open(801,file=Filename_1,status='unknown') 
open(802,file=Filename_2,status='unknown') 
open(803,file=Filename_3,status='unknown') 
do i=1,num_Crack
  write(801, '(50000E20.12)') (Cracks_FluidEle_LCS_z_3D(i)%row(j,1),j=1,Cracks_FluidEle_num_3D(i))       
  write(802, '(50000E20.12)') (Cracks_FluidEle_LCS_z_3D(i)%row(j,2),j=1,Cracks_FluidEle_num_3D(i))   
  write(803, '(50000E20.12)') (Cracks_FluidEle_LCS_z_3D(i)%row(j,3),j=1,Cracks_FluidEle_num_3D(i))        
end do
close(801)          
close(802)    
close(803) 

print *,'    Saving fluid node normal vector of cracks...'
Filename_1   = trim(Full_Pathname)//'.fnnx_'//ADJUSTL(temp)     
Filename_2   = trim(Full_Pathname)//'.fnny_'//ADJUSTL(temp)    
Filename_3   = trim(Full_Pathname)//'.fnnz_'//ADJUSTL(temp)   
open(901,file=Filename_1,status='unknown') 
open(902,file=Filename_2,status='unknown') 
open(903,file=Filename_3,status='unknown')            
do i=1,num_Crack
  write(901, '(50000E20.12)') (Cracks_CalP_Orient_3D(i)%row(j,1),j=1,Cracks_CalP_Num_3D(i))       
  write(902, '(50000E20.12)') (Cracks_CalP_Orient_3D(i)%row(j,2),j=1,Cracks_CalP_Num_3D(i))   
  write(903, '(50000E20.12)') (Cracks_CalP_Orient_3D(i)%row(j,3),j=1,Cracks_CalP_Num_3D(i))        
end do
close(901)          
close(902)    
close(903)  

if(Key_Save_Crack_Radius==1) then
    Filename_1=trim(Full_Pathname)//'.crrd_'//ADJUSTL(temp)     
    open(801,file=Filename_1,status='unknown') 
    do i_C=1,num_Crack
        call Cal_Crack_Avarage_Radius(i_C,tem_crack_radius)
        write(801, '(E20.12)') tem_crack_radius           
    end do
    close(801)          
endif
     
RETURN
END SUBROUTINE Save_Files_Crack_3D
