 
      SUBROUTINE Save_Files_Cross(isub)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Cross
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      write(temp,'(I5)') isub
      
      select case(Key_Data_Format)

      case(1:2)
          print *,'    Saving coordinates of crosses...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'cscr'   
          open(101,file=c_File_name_1,status='unknown')    
          
          do i=1,num_Cross
              write(101, '(2E20.12)') Cross_Point_RABCD(i,1,1:2)
          end do
          close(101)   
          
          print *,'    Saving enns file of cross...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'enns'//'_'//ADJUSTL(temp)       
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Enriched_Node_Type_Cross(i,
     &                                    j),j=1,num_Cross)
          end do
          close(103)     
          
          print *,'    Saving elts file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elts'//'_'//ADJUSTL(temp)       
          open(104,file=c_File_name_1,status='unknown')         
          do i=1,Num_Elem
              write(104, '(200I10)') 
     &                    (Elem_Type_Cross(i,j),j=1,num_Cross)
          end do
          close(104)  
          
          print *,'    Saving poss file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'poss'//'_'//ADJUSTL(temp)         
          open(110,file=c_File_name_1,status='unknown')      
          do i=1,Num_Node
              write(110, '(200I10)') (c_POS_Cross(i,j),j=1,num_Cross)
          end do
          close(110) 
          
          print *,'    Saving nods file of cross...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'nods'//'_'//ADJUSTL(temp)       
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Node_Cross_elem(i,
     &                                    j),j=1,num_Cross)
          end do
          close(103)   
      end select
      RETURN
      END SUBROUTINE Save_Files_Cross
