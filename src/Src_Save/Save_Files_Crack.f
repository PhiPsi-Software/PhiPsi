 
      SUBROUTINE Save_Files_Crack(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_Cross
      use Global_Cross
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j,k
      integer i_Crack
      character(200) c_File_name_1,c_File_name_2
      character(5) temp  
      real(kind=FT) tem_1
      integer i_Line
      
      if (Key_Save_Nothing==1) return
      
      tem_1 = 1000.0D0

      print *,'    Saving coordinates of cracks...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crax'//'_'//ADJUSTL(temp)  
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &             //'cray'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown') 
          
      if(Key_Unit_System==1)then    
          do i=1,num_Crack
              write(101, '(200E20.12)') (Edge_Disposed_Crack(i,
     &                                j,1),j=1,Each_Cr_Poi_Num(i))
              write(102, '(200E20.12)') (Edge_Disposed_Crack(i,
     &                                j,2),j=1,Each_Cr_Poi_Num(i))
          end do
      elseif(Key_Unit_System==2)then 
          do i=1,num_Crack
              write(101, '(200E20.12)') (Edge_Disposed_Crack(i,
     &                      j,1)/tem_1,j=1,Each_Cr_Poi_Num(i))
              write(102, '(200E20.12)') (Edge_Disposed_Crack(i,
     &                      j,2)/tem_1,j=1,Each_Cr_Poi_Num(i))
          end do
      endif
      close(101)
      close(102)    

      print *,'    Saving coordinates of cracks for checking...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crcn'//'_'//ADJUSTL(temp)    
      open(101,file=c_File_name_1,status='unknown')    
      if(Key_Unit_System==1)then     
          do i=1,num_Crack
              do j=1,Each_Cr_Poi_Num(i)
                write(101, '(200E20.12)') Edge_Disposed_Crack(i,j,1),
     &                                    Edge_Disposed_Crack(i,j,2)
              enddo
          end do
      endif
      close(101)  
      

      if(Yes_Arc_Crack)then
          print *,'    Saving coordinates of arc cracks...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'arcc'//'_'//ADJUSTL(temp)  
          open(101,file=c_File_name_1,status='unknown')    
              
          if(Key_Unit_System==1)then     
              do i=1,num_Crack
                  do i_Line=1,11
                      write(101, '(200E20.12)') 
     &                 (Arc_Crack_Coor(i,
     &                         j,i_Line),j=1,Each_Cr_Poi_Num(i)-1)
                  enddo
              end do
          elseif(Key_Unit_System==2)then 
          endif
          close(101)
      endif
      

      print *,'    Saving coordinates of cracks (original)...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crxo'//'_'//ADJUSTL(temp)  
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &             //'cryo'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown') 
          
      if(Key_Unit_System==1)then     
          do i=1,num_Crack
              write(101, '(200E20.12)') (Crack_Coor(i,
     &                                j,1),j=1,Each_Cr_Poi_Num(i))
              write(102, '(200E20.12)') (Crack_Coor(i,
     &                                j,2),j=1,Each_Cr_Poi_Num(i))
          end do
      elseif(Key_Unit_System==2)then 
          do i=1,num_Crack
              write(101, '(200E20.12)') (Crack_Coor(i,
     &                      j,1)/tem_1,j=1,Each_Cr_Poi_Num(i))
              write(102, '(200E20.12)') (Crack_Coor(i,
     &                      j,2)/tem_1,j=1,Each_Cr_Poi_Num(i))
          end do
      endif
      close(101)
      close(102)   
      
      print *,'    Saving ennd file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'ennd'//'_'//ADJUSTL(temp)      
      if (Key_Data_Format ==1)then
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Enriched_Node_Type(i,
     &                                    j),j=1,num_Crack)
          end do
          close(103)   
      elseif (Key_Data_Format ==2)then
          open(103,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')           
          do i=1,Num_Node
              write(103) (Enriched_Node_Type(i,j),j=1,num_Crack)
          end do
          close(103)  
      endif      
      
      print *,'    Saving elty file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'elty'//'_'//ADJUSTL(temp)      
      if (Key_Data_Format ==1)then
          open(104,file=c_File_name_1,status='unknown')         
          do i=1,Num_Elem
              write(104, '(200I10)') (Elem_Type(i,j),j=1,num_Crack)
          end do
          close(104)   
      elseif (Key_Data_Format ==2)then
          open(104,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Elem
              write(104) (Elem_Type(i,j),j=1,num_Crack)
          end do
          close(104)   

      endif
     
      print *,'    Saving celc file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celc'//'_'//ADJUSTL(temp)      
      if (Key_Data_Format ==1)then
          open(105,file=c_File_name_1,status='unknown')  
          do i_Crack = 1,num_Crack
              do i=1,Num_Elem
                write(105,'(200E20.12)')
     &                     (Coors_Element_Crack(i,i_Crack,j),j=1,4)
              end do
          end do
          close(105)   
      elseif (Key_Data_Format ==2)then
          open(105,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i_Crack = 1,num_Crack
              do i=1,Num_Elem
                write(105) (Coors_Element_Crack(i,i_Crack,j),j=1,4)
              end do
          end do
          close(105)   
      endif
      
      print *,'    Saving crab file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crab'//'_'//ADJUSTL(temp)      
      open(105,file=c_File_name_1,status='unknown')  
      do i=1,num_Cross
        write(105,'(200E20.12)')
     &            (Cross_Point_RABCD(i,1:10,1))
        write(105,'(200E20.12)')
     &            (Cross_Point_RABCD(i,1:10,2))
      end do
      close(105) 
      
      print *,'    Saving celt file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celt'//'_'//ADJUSTL(temp)    
      if (Key_Data_Format ==1)then
          open(108,file=c_File_name_1,status='unknown')    
          do i=1,Num_Elem
              write(108, '(200E20.12)') (Coors_Tip(i,j),j=1,2)
          end do
          close(108)  
      elseif (Key_Data_Format ==2)then
          open(108,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Elem
              write(108) (Coors_Tip(i,j),j=1,2)
          end do
          close(108)  

      endif
      
      print *,'    Saving celv file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celv'//'_'//ADJUSTL(temp)      
      if (Key_Data_Format ==1)then     
          open(106,file=c_File_name_1,status='unknown')
          do i=1,Num_Elem
              write(106, '(200E20.12)') (Coors_Vertex(i,j),j=1,2)
          end do
          close(106)   
      elseif (Key_Data_Format ==2)then
          open(106,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Elem
              write(106) (Coors_Vertex(i,j),j=1,2)
          end do
          close(106)   
      endif
      
      print *,'    Saving celj file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celj'//'_'//ADJUSTL(temp)      
      if (Key_Data_Format ==1)then
          open(107,file=c_File_name_1,status='unknown')  
          do i_Crack=1,num_Crack
            do i=1,Num_Elem
              write(107, '(200E20.12)') 
     &                      (Coors_Junction(i,i_Crack,j),j=1,4)
            end do
          end do
          close(107)   
      elseif (Key_Data_Format ==2)then
          open(107,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i_Crack=1,num_Crack
            do i=1,Num_Elem
              write(107) (Coors_Junction(i,i_Crack,j),j=1,4)
            end do
          end do
          close(107)  
      endif
      
      print *,'    Saving ctty file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'ctty'//'_'//ADJUSTL(temp)      
      open(109,file=c_File_name_1,status='unknown')         
      do i=1,num_Crack
          write(109, '(200I10)') (Crack_Tip_Type(i,j),j=1,2)
      end do
      close(109)      
      
      print *,'    Saving posi file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'posi'//'_'//ADJUSTL(temp)   
      if (Key_Data_Format ==1)then
          open(110,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(110, '(200I10)') (c_POS(i,j),j=1,num_Crack)
          end do
          close(110) 
      elseif (Key_Data_Format ==2)then
          open(110,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Node
              write(110) (c_POS(i,j),j=1,num_Crack)
          end do
          close(110)
      endif
      
      print *,'    Saving njel file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'njel'//'_'//ADJUSTL(temp)     
      if (Key_Data_Format ==1)then
          open(111,file=c_File_name_1,status='unknown') 
          do i=1,Num_Node
              write(111, '(200I10)') (Node_Jun_elem(i,
     &                                    j),j=1,num_Crack)
          end do
          close(111) 
      elseif (Key_Data_Format ==2)then
          open(111,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Node
              write(111) (Node_Jun_elem(i,j),j=1,num_Crack)
          end do
          close(111) 
      endif
      
      print *,'    Saving njhl file...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'njhl'//'_'//ADJUSTL(temp)    
      if (Key_Data_Format ==1)then
          open(111,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(111, '(200I10)') (Node_Jun_Hole(i,
     &                                    j),j=1,num_Crack)
          end do
          close(111)    
      elseif (Key_Data_Format ==2)then
          open(111,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')  
          do i=1,Num_Node
              write(111) (Node_Jun_Hole(i,j),j=1,num_Crack)
          end do
          close(111)  
          
      endif
      RETURN
      END SUBROUTINE Save_Files_Crack
