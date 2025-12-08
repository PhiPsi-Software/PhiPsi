 
      subroutine Cal_InSitu_Stress_3D

      use Global_Float_Type
      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_DISP
      use Global_Stress
      use Global_POST
      use Global_Material
      use omp_lib
      
      implicit none
      real(kind=FT) F_InSitu(3*Num_Node)
      integer Total_Num_G_P_InSitu,num_freeDOF_InSitu,
     &        num_fixedDOF_InSitu
      real(kind=FT),ALLOCATABLE::globalK_InSitu(:,:)
      integer,ALLOCATABLE::freeDOF_InSitu(:)
      integer,ALLOCATABLE::fixedDOF_InSitu(:)
      real(kind=FT),ALLOCATABLE::tem_DISP_InSitu(:)
      real(kind=FT) Max_Disp_x_InSitu,Min_Disp_x_InSitu,
     &              Max_Disp_y_InSitu,Min_Disp_y_InSitu,
     &              Max_Disp_z_InSitu,Min_Disp_z_InSitu

      integer Total_Freedom_InSitu

      integer IS_Total_FD      
      integer i_E
      real(kind=FT) c_D(6,6),U(24)
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),
     &                 c_kesi,c_yita,c_zeta,c_Stress(6)
      integer c_NN(8)  
      integer i_G
      real(kind=FT) kesi(Num_Gauss_P_FEM_3D),
     &              yita(Num_Gauss_P_FEM_3D),
     &              zeta(Num_Gauss_P_FEM_3D),
     &              weight(Num_Gauss_P_FEM_3D)
      logical Yes_Add_Insitu
      integer isub
      real(kind=FT) Lambda
      
 1021 FORMAT(5X,'Range of InSitu displacement x:   ',F13.5,
     &                 ' m to ', F13.5,' m')  
 1022 FORMAT(5X,'Range of InSitu displacement y:   ',F13.5,
     &                 ' m to ', F13.5,' m')  
 1023 FORMAT(5X,'Range of InSitu displacement z:   ',F13.5,
     &                 ' m to ', F13.5,' m')       
 1221 FORMAT(5X,'Average Value of InSitu stress x:   ',F11.5, ' MPa')
 1222 FORMAT(5X,'Average Value of InSitu stress y:   ',F11.5, ' MPa')          
 1223 FORMAT(5X,'Average Value of InSitu stress z:   ',F11.5, ' MPa')            
      
      print *,'    ################################'
      print *,'         GEO_STATIC CALCULATION     '
      print *,'    ################################'
      print *,'    Calculating displacement under InSitu stress...'
      
      ALLOCATE(DISP_InSitu(3*Num_Node))
      DISP_InSitu = ZR
      IS_Total_FD = 3*Num_Node
      call Force_Vector_3D_no_Killed_Ele(IS_Total_FD,1,ONE,F_InSitu)
      print *,'    sum(abs(F_InSitu)):',sum(abs(F_InSitu))
      ALLOCATE(freeDOF_InSitu(IS_Total_FD))
      ALLOCATE(fixedDOF_InSitu(IS_Total_FD))
      call Boundary_Cond_3D(IS_Total_FD,1,
     &                      freeDOF_InSitu,num_freeDOF_InSitu,
     &                      fixedDOF_InSitu,num_fixedDOF_InSitu)

      if(Key_SLOE==11)then
          isub = 1
          Lambda = ONE
          call Ele_by_Ele_FEM_PCG_3D(isub,Lambda,cg_tol,max_cg,
     &                   num_freeDOF_InSitu,
     &                   freeDOF_InSitu(1:num_freeDOF_InSitu),
     &                   F_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu)),
     &                   DISP_InSitu)      
          goto 666
      endif
      
              
      ALLOCATE(globalK_InSitu(IS_Total_FD,IS_Total_FD))
      print *,'    Assembling the stiffness matrix K......'
      call Assemble_Stiffness_Matrix_FEM_3D(0,globalK_InSitu,
     &                   IS_Total_FD,Total_Num_G_P_InSitu)   
      print *,'    Sum of globalK_InSitu: ',  sum(globalK_InSitu)     

      
      DISP_InSitu(1:3*Num_Node) = ZR
      ALLOCATE(tem_DISP_InSitu(num_freeDOF_InSitu))
      call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &           globalK_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu),
     &                          freeDOF_InSitu(1:num_freeDOF_InSitu)),
     &           F_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu)),
     &           tem_DISP_InSitu,
     &           num_freeDOF_InSitu)
      
      DISP_InSitu(freeDOF_InSitu(1:num_freeDOF_InSitu))=tem_DISP_InSitu
      print *,'    Displacements under InSitu stress obtained.'
      
      
  666 continue
      
      Max_Disp_x_InSitu = maxval(DISP_InSitu(1:3*num_Node:3))
      Min_Disp_x_InSitu = minval(DISP_InSitu(1:3*num_Node:3))
      Max_Disp_y_InSitu = maxval(DISP_InSitu(2:3*num_Node:3))
      Min_Disp_y_InSitu = minval(DISP_InSitu(2:3*num_Node:3))
      Max_Disp_z_InSitu = maxval(DISP_InSitu(3:3*num_Node:3))
      Min_Disp_z_InSitu = minval(DISP_InSitu(3:3*num_Node:3))     
      
      WRITE(*,1021) Min_Disp_x_InSitu,Max_Disp_x_InSitu
      WRITE(*,1022) Min_Disp_y_InSitu,Max_Disp_y_InSitu
      WRITE(*,1023) Min_Disp_z_InSitu,Max_Disp_z_InSitu
      
      print *,'    Calculating InSitu stress of nodes...'
      ALLOCATE(Str_xx_InSitu(Num_Node))
      ALLOCATE(Str_yy_InSitu(Num_Node))
      ALLOCATE(Str_zz_InSitu(Num_Node))
      ALLOCATE(Str_xy_InSitu(Num_Node))
      ALLOCATE(Str_yz_InSitu(Num_Node))
      ALLOCATE(Str_xz_InSitu(Num_Node))
      ALLOCATE(Str_vm_InSitu(Num_Node))
      Yes_Add_Insitu = .False.
      call Get_Node_Stress_FEM_IN_OUT_3D(Yes_Add_Insitu,0,DISP_InSitu,
     &     Str_xx_InSitu,Str_yy_InSitu,Str_zz_InSitu,
     &     Str_xy_InSitu,Str_yz_InSitu,Str_xz_InSitu,Str_vm_InSitu)
      InSitu_x = sum(Str_xx_InSitu)/Num_Node
      InSitu_y = sum(Str_yy_InSitu)/Num_Node
      InSitu_z = sum(Str_zz_InSitu)/Num_Node
      
      print *,"    ****************************************************"
      WRITE(*,1221) -InSitu_x/Cn_M
      WRITE(*,1222) -InSitu_y/Cn_M
      WRITE(*,1223) -InSitu_z/Cn_M
      print *,"    ****************************************************"
  
      print *,'    Calculating InSitu stress of Gauss points...'
      if (allocated(InSitu_Strs_Gaus_xx).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_xx(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif
      if (allocated(InSitu_Strs_Gaus_yy).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_yy(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif
      if (allocated(InSitu_Strs_Gaus_zz).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_zz(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif
      if (allocated(InSitu_Strs_Gaus_xy).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_xy(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif
      if (allocated(InSitu_Strs_Gaus_yz).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_yz(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif
      if (allocated(InSitu_Strs_Gaus_xz).eqv. .false.) then
          ALLOCATE(InSitu_Strs_Gaus_xz(Num_Elem,Num_Gauss_P_FEM_3D)) 
      endif      


      call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,
     &                                kesi,yita,zeta,weight)      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_G,c_D,c_NN,
!$OMP&            c_X_NODES,c_Y_NODES,c_Z_NODES,U,
!$OMP&            c_kesi,c_yita,c_zeta,c_Stress
!$OMP&            ) SCHEDULE(static)         
      do i_E = 1,Num_Elem
          c_D     = D(Elem_Mat(i_E),1:6,1:6)     
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)    
          c_Z_NODES = G_Z_NODES(:,i_E)            
          U =
     &      [DISP_InSitu(c_NN(1)*3-2),DISP_InSitu(c_NN(1)*3-1),
     &       DISP_InSitu(c_NN(1)*3),
     &       DISP_InSitu(c_NN(2)*3-2),DISP_InSitu(c_NN(2)*3-1),
     &       DISP_InSitu(c_NN(2)*3),
     &       DISP_InSitu(c_NN(3)*3-2),DISP_InSitu(c_NN(3)*3-1),
     &       DISP_InSitu(c_NN(3)*3),
     &       DISP_InSitu(c_NN(4)*3-2),DISP_InSitu(c_NN(4)*3-1),
     &       DISP_InSitu(c_NN(4)*3),
     &       DISP_InSitu(c_NN(5)*3-2),DISP_InSitu(c_NN(5)*3-1),
     &       DISP_InSitu(c_NN(5)*3),
     &       DISP_InSitu(c_NN(6)*3-2),DISP_InSitu(c_NN(6)*3-1),
     &       DISP_InSitu(c_NN(6)*3),
     &       DISP_InSitu(c_NN(7)*3-2),DISP_InSitu(c_NN(7)*3-1),
     &       DISP_InSitu(c_NN(7)*3),
     &       DISP_InSitu(c_NN(8)*3-2),DISP_InSitu(c_NN(8)*3-1),
     &       DISP_InSitu(c_NN(8)*3)] 
          do i_G = 1,Num_Gauss_P_FEM_3D
            c_kesi = kesi(i_G)                                          
            c_yita = yita(i_G)
              c_zeta = zeta(i_G)
            call Cal_Ele_Str_N8_3D(i_E,0,1,1,
     &                               c_X_NODES,c_Y_NODES,c_Z_NODES,
     &                               c_D,c_kesi,c_yita,c_zeta,U,
     &                               c_Stress)      
              InSitu_Strs_Gaus_xx(i_E,i_G) = c_Stress(1)
            InSitu_Strs_Gaus_yy(i_E,i_G) = c_Stress(2)
            InSitu_Strs_Gaus_zz(i_E,i_G) = c_Stress(3)       
              InSitu_Strs_Gaus_xy(i_E,i_G) = c_Stress(4)
            InSitu_Strs_Gaus_yz(i_E,i_G) = c_Stress(5)
            InSitu_Strs_Gaus_xz(i_E,i_G) = c_Stress(6)                
          end do  
      end do 
!$omp end parallel do 
      
      if(allocated(globalK_InSitu)) DEALLOCATE(globalK_InSitu)     
      if(allocated(freeDOF_InSitu)) DEALLOCATE(freeDOF_InSitu)  
      if(allocated(tem_DISP_InSitu)) DEALLOCATE(tem_DISP_InSitu)  
      if(allocated(fixedDOF_InSitu)) DEALLOCATE(fixedDOF_InSitu)   
      
      State_InSitu = 0
      if(InSitu_x < 0.01D6 .or. InSitu_y < 0.01D6 .or.
     &   InSitu_z < 0.01D6)then
          State_InSitu = 1
          InSitu_x = -InSitu_x
          InSitu_y = -InSitu_y
          InSitu_z = -InSitu_z
      endif
      
      return 
      end SUBROUTINE Cal_InSitu_Stress_3D               
