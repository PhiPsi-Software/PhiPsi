 
      SUBROUTINE Save_Files_Holes(isub)

      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp  
      
      if (Key_Save_Nothing==1) return
      
      select case(Key_Data_Format)
      

      case(1:2)
          print *,'    Saving coordinates of holes...'
          write(temp,'(I5)') isub
          if (num_Circ_Hole>=1)then
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'hlcr'   
              open(101,file=c_File_name_1,status='unknown')    
              do i=1,num_Circ_Hole
                  write(101, '(3E20.12)') Hole_Coor(i,1:3)
              end do
              close(101) 
          endif
          if (num_Ellip_Hole>=1)then 
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'ehcr'             
              open(102,file=c_File_name_1,status='unknown')    
              do i=1,num_Ellip_Hole
                  write(102, '(5E20.12)') Ellip_Hole_Coor(i,1:5)
              end do
              close(102)  
          endif
          
          
          print *,'    Saving ennh file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'ennh'//'_'//ADJUSTL(temp)       
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Enriched_Node_Type_Hl(i,
     &                                    j),j=1,num_Hole)
          end do
          close(103)           
          print *,'    Saving elth file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elth'//'_'//ADJUSTL(temp)       
          open(104,file=c_File_name_1,status='unknown')         
          do i=1,Num_Elem
              write(104, '(200I10)') (Elem_Type_Hl(i,j),j=1,num_Hole)
          end do
          close(104)   
          print *,'    Saving posh file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'posh'//'_'//ADJUSTL(temp)         
          open(110,file=c_File_name_1,status='unknown')      
          do i=1,Num_Node
              write(110, '(200I10)') (c_POS_Hl(i,j),j=1,num_Hole)
          end do
          close(110) 
      end select
      RETURN
      END SUBROUTINE Save_Files_Holes
