 
      subroutine Cal_InSitu_Stress_for_Hole_and_Inclusion

      use Global_Float_Type
      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_DISP
      use Global_HF
      use Global_Stress
      use Global_POST
      use Global_Contact
      use Global_Material
      use Global_Inclusion
      
      implicit none
      real(kind=FT) ,ALLOCATABLE::F_InSitu(:)
      integer Total_Num_G_P_InSitu
      real(kind=FT),ALLOCATABLE::globalK_InSitu(:,:)
      integer,ALLOCATABLE::freeDOF_InSitu(:)
      integer num_FreeD_InSitu,num_FixedD_InSitu
      integer,ALLOCATABLE::fixedDOF_InSitu(:)
      real(kind=FT),ALLOCATABLE::tem_DISP_InSitu(:)
      real(kind=FT) Max_Disp_x_InSitu,Min_Disp_x_InSitu,
     &                 Max_Disp_y_InSitu,Min_Disp_y_InSitu
      integer Total_Freedom_InSitu
      integer c_Total_Num_G_P
      integer i_E
      real(kind=FT) c_D(3,3)
      integer c_NN(4)
      real(kind=FT) kesi_N(4),
     &              yita_N(4),
     &              weight_N(4)    
      integer G_Counter
      logical Yes_Add_Insitu
      
      integer i
 1021 FORMAT(5X,'Range of InSitu displacement x:   ',F13.5,
     &                 ' m to ', F13.5,' m')  
 1022 FORMAT(5X,'Range of InSitu displacement y:   ',F13.5,
     &                 ' m to ', F13.5,' m')  
 1221 FORMAT(5X,'Value of InSitu stress x:   ',F13.5, ' MPa')
 1222 FORMAT(5X,'Value of InSitu stress y:   ',F13.5, ' MPa')     
 
      print *,'    ################################'
      print *,'         GEO_STATIC CALCULATION     '
      print *,'    ################################'
      print *,'    Calculating displacement under InSitu stress...'
      
      call Determine_Enriched_Nodes(1,1)

      call Number_Enriched_Nodes(1)
      print *,'    Total_FD:',Total_FD
      print *,'    Enrich_Freedom:',Enrich_Freedom

      ALLOCATE(freeDOF_InSitu(Total_FD))
      ALLOCATE(fixedDOF_InSitu(Total_FD))
      call Boundary_Cond_output_Fixed_Freedom(Total_FD,1,
     &                              freeDOF_InSitu,num_FreeD_InSitu,
     &                              fixedDOF_InSitu,num_FixedD_InSitu) 


      ALLOCATE(F_InSitu(Total_FD))
      call Force_Vector(Total_FD,1,.False.,1.0D0,F_InSitu)    
      print *,'    Sum of F_InSitu:', sum(F_InSitu)

      ALLOCATE(globalK_InSitu(Total_FD,Total_FD))
      print *,'    Assembling K......'
      call Assemble_Stiffness_Matrix_XFEM(1,globalK_InSitu,
     &                       Total_FD,c_Total_Num_G_P)
      print *,'    Sum of globalK: ',  sum(globalK_InSitu)

      ALLOCATE(DISP_InSitu(Total_FD))
      DISP_InSitu(1:Total_FD)= ZR
      ALLOCATE(tem_DISP_InSitu(num_FreeD_InSitu))
      
      call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &            globalK_InSitu(freeDOF_InSitu(1:num_FreeD_InSitu),
     &                           freeDOF_InSitu(1:num_FreeD_InSitu)),
     &                  F_InSitu(freeDOF_InSitu(1:num_FreeD_InSitu)),
     &                     tem_DISP_InSitu,num_FreeD_InSitu)

      DISP_InSitu(freeDOF_InSitu(1:num_FreeD_InSitu))=tem_DISP_InSitu
      print *,'    Displacements under InSitu stress obtained.'
      
      Max_Disp_x_InSitu = maxval(DISP_InSitu(1:2*num_Node:2))
      Min_Disp_x_InSitu = minval(DISP_InSitu(1:2*num_Node:2))
      Max_Disp_y_InSitu = maxval(DISP_InSitu(2:2*num_Node:2))
      Min_Disp_y_InSitu = minval(DISP_InSitu(2:2*num_Node:2))
      WRITE(*,1021) Min_Disp_x_InSitu,Max_Disp_x_InSitu
      WRITE(*,1022) Min_Disp_y_InSitu,Max_Disp_y_InSitu
      
      print *,'    Calculating InSitu stress of nodes...'
      ALLOCATE(Str_xx_InSitu(Num_Node))
      ALLOCATE(Str_yy_InSitu(Num_Node))
      ALLOCATE(Str_xy_InSitu(Num_Node))
      ALLOCATE(Str_vm_InSitu(Num_Node))
      Yes_Add_Insitu = .False.
      
      call Get_Node_Stress_XFEM_IN_OUT(Yes_Add_Insitu,1,DISP_InSitu,
     &     Str_xx_InSitu,Str_yy_InSitu,Str_xy_InSitu,Str_vm_InSitu)
     
#ifndef Silverfrost
      do i=1,Num_Node
          if (isnan(Str_xx_InSitu(i))) Str_xx_InSitu(i)= ZR 
          if (isnan(Str_yy_InSitu(i))) Str_yy_InSitu(i)= ZR 
          if (isnan(Str_xy_InSitu(i))) Str_xy_InSitu(i)= ZR 
          if (isnan(Str_vm_InSitu(i))) Str_vm_InSitu(i)= ZR        
      enddo
#endif      
      
      InSitu_x = sum(Str_xx_InSitu)/Num_Node
      InSitu_y = sum(Str_yy_InSitu)/Num_Node
      InSitu_xy = sum(Str_xy_InSitu)/Num_Node
      
      print *,"    ************************************************"
      WRITE(*,1221) -InSitu_x/1.0D6
      WRITE(*,1222) -InSitu_y/1.0D6
      print *,"    ************************************************"
      

      print *,'    Calculating InSitu stress of Gauss points...'
      ALLOCATE(InSitu_Strs_Gaus_xx(Num_Elem,Num_Gauss_P_FEM)) 
      ALLOCATE(InSitu_Strs_Gaus_yy(Num_Elem,Num_Gauss_P_FEM)) 
      ALLOCATE(InSitu_Strs_Gaus_xy(Num_Elem,Num_Gauss_P_FEM)) 
      call Cal_Gauss_Points_QUAD(4,kesi_N,yita_N,weight_N)
      G_Counter = 0
      do i_E = 1,Num_Elem
          c_D     = D(Elem_Mat(i_E),:,:)    
          
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_NN    = G_NN(:,i_E)
         
          InSitu_Strs_Gaus_xx(i_E,1:4) = sum(Str_xx_InSitu(c_NN))/FOU
          InSitu_Strs_Gaus_yy(i_E,1:4) = sum(Str_yy_InSitu(c_NN))/FOU
          InSitu_Strs_Gaus_xy(i_E,1:4) = sum(Str_xy_InSitu(c_NN))/FOU
      end do 
      
      DEALLOCATE(globalK_InSitu)
      DEALLOCATE(freeDOF_InSitu,tem_DISP_InSitu)
      
      State_InSitu = 0
      if(InSitu_x < 0.01D6 .or. InSitu_y < 0.01D6)then
          State_InSitu = 1
          InSitu_x = -InSitu_x
          InSitu_y = -InSitu_y
      endif

      return 
      end SUBROUTINE Cal_InSitu_Stress_for_Hole_and_Inclusion               
