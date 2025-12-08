 
      SUBROUTINE Save_HF_wbfp_File

      use Global_Float_Type      
      use Global_Common
      use Global_Filename
      use Global_POST
      use Global_HF
      
      implicit none
      
      character(200) Filename_1
      real(kind=FT) c_Point(3)
      integer c_i_WB,c_i_Stage
      integer num_Stages
      real(kind=FT) Start_Point(3),End_Point(3)
      real(kind=FT) L_Stage
      real(kind=FT) L_WB_Fracturing
      real(kind=FT) Tool_Function_2Point_Dis_3D
      
      if (Key_Save_Nothing==1) return 
      
      Filename_1=trim(Full_Pathname)//'.wbfp'
      open(101,file=Filename_1,status='unknown')  
      do c_i_WB = 1,num_Wellbore
          num_Stages = num_Stages_Wellbores(c_i_WB)          
          Start_Point = Wellbores_Start_Point(c_i_WB,1:3)  
          if(abs(Start_Point(1))>=1.0D10) then
              close(101)
              return
          endif
          End_Point   = Wellbores_End_Point(c_i_WB,1:3)   
          L_WB_Fracturing=Tool_Function_2Point_Dis_3D(Start_Point,
     &                                                End_Point)              
          L_Stage = L_WB_Fracturing/num_Stages
          write(101, '(3E20.12)') Start_Point(1:3)
          do c_i_Stage = 1,num_Stages_Wellbores(c_i_WB)-1
              call Tool_Offset_Point_A_to_Point_B_3D(Start_Point,
     &                            End_Point,
     &                            L_Stage*c_i_Stage,c_Point(1:3))
              write(101, '(3E20.12)') c_Point(1:3)
          enddo  
          write(101, '(3E20.12)') End_Point(1:3)
      enddo
      close(101)
 
      RETURN
      END SUBROUTINE Save_HF_wbfp_File
