 
      SUBROUTINE Save_Disp(isub,Key_CoorSys)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer isub,Key_CoorSys
      integer i
      character(200) c_File_name_1,c_File_name_2
      character(5) temp   
      real(kind=FT) tem_1
      integer c_Total_FD 
      
      if (Key_Save_Nothing==1) return
      
      c_Total_FD = Total_FD 
      
      
      if (Key_CoorSys == 1)then
          print *,'    Saving disp...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'disp'//'_'//ADJUSTL(temp)  
          c_File_name_2   =  trim(Full_Pathname)//'.'
     &                     //'disn'//'_'//ADJUSTL(temp)    
          tem_1 = 1000.0D0
          select case(Key_Data_Format)
          case(1)  
              open(201,file=c_File_name_1,status='unknown')     
              do i=1,c_Total_FD
                  write(201, '(1E20.12)') DISP(i)
              end do
              close(201)    
              open(202,file=c_File_name_2,status='unknown')
              if (Key_Dimension == 2) then   
                  if(Key_Unit_System==1)then     
                      do i=1,int(c_Total_FD/2)
                          write(202, '(I8,2(A,E20.12))') 
     &                           i,',',DISP(2*i-1),',',DISP(2*i)
                      end do
                  elseif(Key_Unit_System==2)then 
                      do i=1,int(c_Total_FD/2)
                          write(202, '(I8,2(A,E20.12))') 
     &                       i,',',DISP(2*i-1)/tem_1,',',DISP(2*i)/tem_1
                      end do
                  endif      
              else if(Key_Dimension == 3)then
                  do i=1,int(c_Total_FD/3)
                      write(202, '(I8,3E20.12)') i,
     &                                 DISP(3*i-2),DISP(3*i-1),DISP(3*i)
                  end do          
              end if          
              close(202)  
          case(2)  
              open(201,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              write(201) (DISP(i),i=1,c_Total_FD)
              close(201)    

              if (Key_Dimension == 2) then 
                  open(203,file=c_File_name_2,status='unknown',
     &                     form='unformatted',access='stream')     
                write(203)(DISP(2*i-1),DISP(2*i),i=1,int(c_Total_FD/2))      
                  close(203)   
              else if(Key_Dimension == 3)then
                  open(203,file=c_File_name_2,status='unknown',
     &                     form='unformatted',access='stream') 
                  do i=1,int(c_Total_FD/3)
                      write(203) DISP(3*i-2),DISP(3*i-1),DISP(3*i)
                  end do    
                  close(203)  
              endif
          end select
      
      elseif (Key_CoorSys == 2)then  
          print *,'    Saving disp in cylindrical coodinates...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'dipc'//'_'//ADJUSTL(temp)  
          c_File_name_2   =  trim(Full_Pathname)//'.'
     &                     //'dinc'//'_'//ADJUSTL(temp)    
          tem_1 = 1000.0D0
          select case(Key_Data_Format)   
          case(1)  
              open(201,file=c_File_name_1,status='unknown')     
              do i=1,c_Total_FD
                  write(201, '(1E20.12)') DISP_Cylinder(i)
              end do
              close(201)    
              open(202,file=c_File_name_2,status='unknown')
              if (Key_Dimension == 2) then   
                  if(Key_Unit_System==1)then     
                      do i=1,int(c_Total_FD/2)
                          write(202, '(I8,2(A,E20.12))') 
     &                        i,',',DISP_Cylinder(2*i-1),',',
     &                       DISP_Cylinder(2*i)
                      end do
                  elseif(Key_Unit_System==2)then 
                      do i=1,int(c_Total_FD/2)
                          write(202, '(I8,2(A,E20.12))') 
     &                       i,',',DISP_Cylinder(2*i-1)/tem_1,',',
     &                             DISP_Cylinder(2*i)/tem_1
                      end do
                  endif      
              else if(Key_Dimension == 3)then
                  do i=1,int(c_Total_FD/3)
                      write(202, '(I8,3E20.12)') i,
     &                          DISP_Cylinder(3*i-2),
     &                          DISP_Cylinder(3*i-1),DISP_Cylinder(3*i)
                  end do          
              end if          
              close(202)       
          case(2)  
              open(201,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              write(201) (DISP_Cylinder(i),i=1,c_Total_FD)
              close(201)    

              if (Key_Dimension == 2) then 
                  open(203,file=c_File_name_2,status='unknown',
     &                     form='unformatted',access='stream')     
                  write(203)(DISP_Cylinder(2*i-1),
     &                     DISP_Cylinder(2*i),i=1,int(c_Total_FD/2))      
                  close(203)            
              else if(Key_Dimension == 3)then
                  open(203,file=c_File_name_2,status='unknown',
     &                     form='unformatted',access='stream') 
                  do i=1,int(c_Total_FD/3)
                      write(203) DISP_Cylinder(3*i-2),
     &                 DISP_Cylinder(3*i-1),DISP_Cylinder(3*i)
                  end do    
                  close(203)  
              endif
          end select
      endif
      
      RETURN
      END SUBROUTINE Save_Disp
