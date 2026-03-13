!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      SUBROUTINE Save_Strain_Node(isub,Key_CoorSys)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Strain
      use Global_POST
      
      implicit none
      
      integer isub,Key_CoorSys
      integer i
      character(200) c_File_name_1
      character(5) temp   
      real(kind=FT) tem_2  
      
      if (Key_Save_Nothing==1) return 
      
      tem_2 = 1.0D6

      if (Key_CoorSys==1) then
          print *,'    Saving strain of all nodes...'      
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'sran'//'_'//ADJUSTL(temp)  
           
          select case(Key_Data_Format)
          case(1)  
              open(202,file=c_File_name_1,status='unknown')
              if (Key_Dimension == 2) then 
                  if(Key_Unit_System==1)then    
                      do i=1,num_Node
                          write(202, '(I8,4E20.12)') i,
     &                                Strain_xx_Node(i),
     &                                Strain_yy_Node(i),
     &                                Strain_xy_Node(i),    
     &                                Strain_vm_Node(i)
                      end do
                  elseif(Key_Unit_System==2)then 
                      do i=1,num_Node
                          write(202, '(I8,4E20.12)') i,
     &                                Strain_xx_Node(i)*tem_2,
     &                                Strain_yy_Node(i)*tem_2,
     &                                Strain_xy_Node(i)*tem_2,    
     &                                Strain_vm_Node(i)*tem_2
                      end do
                  endif        
              else if(Key_Dimension == 3)then
                  do i=1,num_Node
                      write(202, '(I8,7E20.12)') i,
     &                            Strain_xx_Node(i),
     &                            Strain_yy_Node(i),
     &                            Strain_zz_Node(i),
     &                            Strain_xy_Node(i),
     &                            Strain_yz_Node(i),
     &                            Strain_xz_Node(i),       
     &                            Strain_vm_Node(i)   
                  end do          
              end if          
              close(202)  
          case(2) 
              open(203,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              if (Key_Dimension == 2) then        
                  write(203) (Strain_xx_Node(i),
     &                        Strain_yy_Node(i),
     &                        Strain_xy_Node(i),    
     &                        Strain_vm_Node(i),i=1,num_Node)
        
              else if(Key_Dimension == 3)then
                  write(203)     (Strain_xx_Node(i),
     &                            Strain_yy_Node(i),
     &                            Strain_zz_Node(i),
     &                            Strain_xy_Node(i),
     &                            Strain_yz_Node(i),
     &                            Strain_xz_Node(i),       
     &                            Strain_vm_Node(i),i=1,num_Node)

              end if          
              close(203)      
          end select
          
      elseif(Key_CoorSys==2) then
          print *,'    Saving strain of all nodes in cylindrical CS...'      
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'srac'//'_'//ADJUSTL(temp)  
     
          select case(Key_Data_Format)
          case(1) 
              open(202,file=c_File_name_1,status='unknown')
              if (Key_Dimension == 2) then 

              else if(Key_Dimension == 3)then
                  do i=1,num_Node
                      write(202, '(I8,7E20.12)') i,
     &                            Strain_Crr_Node(i),
     &                            Strain_Ctt_Node(i),
     &                            Strain_Czz_Node(i),
     &                            Strain_Crt_Node(i),
     &                            Strain_Ctz_Node(i),
     &                            Strain_Crz_Node(i),       
     &                            Strain_Cvm_Node(i)   
                  end do          
              end if          
              close(202)  
          case(2)  
              open(203,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              if (Key_Dimension == 2) then           
              else if(Key_Dimension == 3)then
                  write(203)     (Strain_Crr_Node(i),
     &                            Strain_Ctt_Node(i),
     &                            Strain_Czz_Node(i),
     &                            Strain_Crt_Node(i),
     &                            Strain_Ctz_Node(i),
     &                            Strain_Crz_Node(i),       
     &                            Strain_Cvm_Node(i),i=1,num_Node)

              end if          
              close(203)      
          end select      
      endif

      RETURN
      END SUBROUTINE Save_Strain_Node
