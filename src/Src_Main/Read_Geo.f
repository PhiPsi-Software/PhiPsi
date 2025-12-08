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
 
      SUBROUTINE Read_Geo
          
c     ------------------------
c Read public variable module
c     ------------------------
      use Global_Float_Type
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Common
      use Global_Contact
      use Global_Field_Problem
      use Global_Inclusion
      use Global_Dynamic
      use Global_Cohesive
      use Global_Cross
      use omp_lib
      
c     ----------------------
C Variable Type Declaration
c     ----------------------
      implicit none
      LOGICAL alive
      character*200 temp_name
      integer Tool_Count_Lines
      integer temp_count,c_i,c_j
      
c     ---------------
C Temporary variable
c     ---------------
      real(kind=FT),ALLOCATABLE::Temp_DATA(:,:)
      logical Flag_Blank
      integer,ALLOCATABLE::tem1(:,:)
      integer,ALLOCATABLE::tem2(:)
      integer i,j,all_num_Outline,Out_num_Outline
      integer,ALLOCATABLE::All_Outline(:,:)
      integer,ALLOCATABLE::Temp_Outline(:,:)  
      integer N1,N2,N3,N4,NN(4)
      real(kind=FT) area,Centroid(2)
      integer c_all_num_Ele,c_Node,i_E
      integer c_all_Elem(32),Uniqued_all_Elem(32)
      integer num_Surrou
      integer i_Check
      real(kind=FT) c_L1,c_L2,c_L3,c_L4,c_L_min,c_L_max
 
c     ------------------------------
c Check the file name and path name
c     ------------------------------
      if (Filename=='*blank*')then
          print *, "    Error :: Cannot find Keyword *Filename!"
          call Warning_Message('S',Keywords_Blank) 
      end if
      if (Work_Directory == '*blank*') then
          print *, "    Error :: Cannot find Keyword *Work_Directory!"
          call Warning_Message('S',Keywords_Blank) 
      end if
      
c     ---------------------------------------------------------------
c Check if the node coordinate file exists, and if it does, read it.
c     ---------------------------------------------------------------
      print *, "    Trying to read nodal files...." 
      temp_name = trim(trim(Full_Pathname)//'.node')
      inquire(file=temp_name, exist=alive)
      if(alive.EQV..FALSE.)then
          print *, "    Error :: Can not find nodal file!!!" 
          print *, "    Check whether the following file exists or not:"
          print *, "    ",temp_name
          call Warning_Message('S',Keywords_Blank) 
      else
          ! Count the number of lines in a file
          Num_Node = Tool_Count_Lines(temp_name)
          IF(ALLOCATED(Coor)) DEALLOCATE(Coor)
          ALLOCATE(Coor(Num_Node,2))
          ALLOCATE(Temp_DATA(Num_Node,2))
          Call Tool_Read_File(temp_name,"node",Num_Node,2,Temp_DATA,
     &                        Flag_Blank)
          Coor = Temp_DATA
          DEALLOCATE(Temp_DATA)
          ! Coordinate extrema
          Max_X_Coor = maxval(Coor(:,1))
          Min_X_Coor = minval(Coor(:,1))
          Max_Y_Coor = maxval(Coor(:,2))
          Min_Y_Coor = minval(Coor(:,2))
          ! Model center coordinates. Added on 2024-06-22.
          Model_Center_x = (Max_X_Coor + Min_X_Coor)/TWO 
          Model_Center_y = (Max_Y_Coor + Min_Y_Coor)/TWO
          !Model coordinates range. Added on 2024-06-22.
          Model_X_Range = Max_X_Coor-Min_X_Coor
          Model_Y_Range = Max_Y_Coor-Min_Y_Coor
          Max_Model_Range=max(Model_X_Range,Model_Y_Range)
          Min_Model_Range=min(Model_X_Range,Model_Y_Range)
      endif
      
c     -------------------------------------------------------------------------------
c Check whether the input file is a three-dimensional problem, and if so, terminate.
c     -------------------------------------------------------------------------------
      temp_name = trim(trim(Full_Pathname)//'.bouz')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..TRUE.)then
          print *, "    Error :: *.bouz file found, please"
     &          // " check *Key_Dimension (2 or 3)!" 
          call Warning_Message('S',Keywords_Blank) 
      endif  
      temp_name = trim(trim(Full_Pathname)//'.focz')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..TRUE.)then
          print *, "    Error :: *.focz file found, please"
     &          // " check *Key_Dimension (2 or 3)!" 
          call Warning_Message('S',Keywords_Blank) 
      endif 
      
      
c     --------------------------------------------------------------
c Check if the element number file exists, and if it does, read it.
c     --------------------------------------------------------------
      print *, "    Trying to read element files..." 
      temp_name = trim(trim(Full_Pathname)//'.elem')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Error :: Can not find element file!" 
          call Warning_Message('S',Keywords_Blank) 
      else
          Num_Elem = Tool_Count_Lines(temp_name)
          IF(ALLOCATED(Elem_Node)) DEALLOCATE(Elem_Node)
          ALLOCATE(Elem_Node(Num_Elem,4))
          IF(ALLOCATED(Elem_Mat)) DEALLOCATE(Elem_Mat)
          ALLOCATE(Elem_Mat(Num_Elem)) 
          ALLOCATE(Temp_DATA(Num_Elem,5))
          Call Tool_Read_File(temp_name,"elem",Num_Elem,5,Temp_DATA,
     &                        Flag_Blank)
          Elem_Node = int(Temp_DATA(:,1:4))
          Elem_Mat  = int(Temp_DATA(:,5))
          ! num_of_Material = MaxVal(Elem_Mat)     ! Number of material types
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ---------------------------------------------------------------------------------------------
c Obtain the maximum node number difference in the element to calculate the maximum half-bandwidth
c     A first course in the finite element method, 5ed, B.4.
c     ---------------------------------------------------------------------------------------------
      Max_Diff_Elem_Num = maxval(maxval(Elem_Node(:,1:4),2)-
     &                           minval(Elem_Node(:,1:4),2))
      Max_Half_Band_Width = 2*(Max_Diff_Elem_Num + 1)
      
c     --------------------------------------------------------------------
c Check if the x boundary constraint file exists, and if it does, read it
c     --------------------------------------------------------------------
      print *, "    Trying to read boux files..." 
      temp_name = trim(trim(Full_Pathname)//'.boux')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find boux file!" 
      else
          Num_Bou_x = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Bou_x)) DEALLOCATE(Bou_x)
          ALLOCATE(Bou_x(Num_Bou_x))
          ALLOCATE(Temp_DATA(Num_Bou_x,1))
          Call Tool_Read_File(temp_name,"boux",Num_Bou_x,1,Temp_DATA,
     &                        Flag_Blank)
          Bou_x  = int(Temp_DATA(:,1))
          DEALLOCATE(Temp_DATA)
      endif  
      
c     --------------------------------------------------------------------
c Check if the y-boundary constraint file exists, and if it does, read it
c     --------------------------------------------------------------------
      print *, "    Trying to read bouy files..." 
      temp_name = trim(trim(Full_Pathname)//'.bouy')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find bouy file!" 
      else
          Num_Bou_y = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Bou_y)) DEALLOCATE(Bou_y)
          ALLOCATE(Bou_y(Num_Bou_y))
          ALLOCATE(Temp_DATA(Num_Bou_y,1))
          Call Tool_Read_File(temp_name,"bouy",Num_Bou_y,1,Temp_DATA,
     &                        Flag_Blank)
          Bou_y  = int(Temp_DATA(:,1))
          DEALLOCATE(Temp_DATA)
      endif    
      
c     -----------------------------------------------------------------------
c Check if the load file in the x-direction exists, and if it does, read it.
c     -----------------------------------------------------------------------
      print *, "    Trying to read focx files..." 
      temp_name = trim(trim(Full_Pathname)//'.focx')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find focx file!" 
      else
          Num_Foc_x = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Foc_x)) DEALLOCATE(Foc_x)
          ALLOCATE(Foc_x(Num_Foc_x,2))
          ALLOCATE(Temp_DATA(Num_Foc_x,2))
          Call Tool_Read_File(temp_name,"focx",Num_Foc_x,2,Temp_DATA,
     &                        Flag_Blank)
          Foc_x  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif  
      
c     -----------------------------------------------------------------------
c Check if the load file in the y-direction exists, and if it does, read it.
c     -----------------------------------------------------------------------
      print *, "    Trying to read focy files..." 
      temp_name = trim(trim(Full_Pathname)//'.focy')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find focy file!" 
      else
          Num_Foc_y = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Foc_y)) DEALLOCATE(Foc_y)
          ALLOCATE(Foc_y(Num_Foc_y,2))
          ALLOCATE(Temp_DATA(Num_Foc_y,2))
          Call Tool_Read_File(temp_name,"focy",Num_Foc_y,2,Temp_DATA,
     &                        Flag_Blank)
          Foc_y  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 

c     -----------------------------------------------------------------------------------------------
c Check whether the boundary condition file for non-zero displacement in the x direction exists, and
c if it does, read it.
c     -----------------------------------------------------------------------------------------------
      print *, "    Trying to read buxn files..." 
      temp_name = trim(trim(Full_Pathname)//'.buxn')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Boux_nonzero = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Bou_x_nonzero)) DEALLOCATE(Bou_x_nonzero)
          ALLOCATE(Bou_x_nonzero(Num_Boux_nonzero,2))
          ALLOCATE(Temp_DATA(Num_Boux_nonzero,2))
          Call Tool_Read_File(temp_name,"buxn",Num_Boux_nonzero,2,
     &                        Temp_DATA,Flag_Blank)
          Bou_x_nonzero  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif  

c     -----------------------------------------------------------------------------------------------
c Check whether the boundary condition file for non-zero displacement in the y direction exists, and
c if it does, read it.
c     -----------------------------------------------------------------------------------------------
      print *, "    Trying to read buyn files..." 
      temp_name = trim(trim(Full_Pathname)//'.buyn')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Bouy_nonzero = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Bou_y_nonzero)) DEALLOCATE(Bou_y_nonzero)
          ALLOCATE(Bou_y_nonzero(Num_Bouy_nonzero,2))
          ALLOCATE(Temp_DATA(Num_Bouy_nonzero,2))
          Call Tool_Read_File(temp_name,"buyn",Num_Bouy_nonzero,2,
     &                        Temp_DATA,Flag_Blank)
          Bou_y_nonzero  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     -----------------------------------------------------------------------------------------------
c Check whether the coupling file for the x-direction degree of freedom exists, and if it does, read
c it.
c     -----------------------------------------------------------------------------------------------
      print *, "    Trying to read dofx files..." 
      temp_name = trim(trim(Full_Pathname)//'.dofx ')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          temp_count = Tool_Count_Lines(temp_name) 
          ALLOCATE(Temp_DATA(temp_count,2))
          Call Tool_Read_File(temp_name,"dofx",temp_count,2,
     &                        Temp_DATA,Flag_Blank)
          num_CP_set_x = int(maxval(Temp_DATA(:,1)));
          do c_i =1,temp_count
              do c_j=1,num_CP_set_x
                  if (Temp_DATA(c_i,1) == c_j)then
                      num_nodes_CP_set_x(c_j) =num_nodes_CP_set_x(c_j)+1
                      CP_nodes_x(c_j,num_nodes_CP_set_x(c_j))=
     &                                    int(Temp_DATA(c_i,2)) 
                  endif
              enddo
          enddo
          DEALLOCATE(Temp_DATA)
      endif  
      
c     --------------------------------------------------------------------------------------------
c Check whether the y-direction degree of freedom coupling file exists; if it does, then read it.
c     --------------------------------------------------------------------------------------------
      print *, "    Trying to read dofy files..." 
      temp_name = trim(trim(Full_Pathname)//'.dofy ')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          temp_count = Tool_Count_Lines(temp_name) 
          ALLOCATE(Temp_DATA(temp_count,2))
          Call Tool_Read_File(temp_name,"dofy",temp_count,2,
     &                        Temp_DATA,Flag_Blank)
          num_CP_set_y = int(maxval(Temp_DATA(:,1)));
          do c_i =1,temp_count
              do c_j=1,num_CP_set_y
                  if (Temp_DATA(c_i,1) == c_j)then
                      num_nodes_CP_set_y(c_j) =num_nodes_CP_set_y(c_j)+1
                      CP_nodes_y(c_j,num_nodes_CP_set_y(c_j))=
     &                                    int(Temp_DATA(c_i,2)) 
                  endif
              enddo
          enddo
          DEALLOCATE(Temp_DATA)
      endif     
      
c     ----------------------------------------------------------------------------------------
c Check if the initial velocity file in the x direction exists; if it does, read it (used for
c dynamic analysis)
c     ----------------------------------------------------------------------------------------
      print *, "    Trying to read ivex files..." 
      temp_name = trim(trim(Full_Pathname)//'.ivex')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Ivex = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Ive_x)) DEALLOCATE(Ive_x)
          ALLOCATE(Ive_x(Num_Ivex,2))
          ALLOCATE(Temp_DATA(Num_Ivex,2))
          Call Tool_Read_File(temp_name,"ivex",Num_Ivex,2,Temp_DATA,
     &                        Flag_Blank)
          Ive_x  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ----------------------------------------------------------------------------------------
c Check if the initial velocity file in the y-direction exists; if it does, read it (used for
c dynamic analysis)
c     ----------------------------------------------------------------------------------------
      print *, "    Trying to read ivey files..." 
      temp_name = trim(trim(Full_Pathname)//'.ivey')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Ivey = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Ive_y)) DEALLOCATE(Ive_y)
          ALLOCATE(Ive_y(Num_Ivey,2))
          ALLOCATE(Temp_DATA(Num_Ivey,2))
          Call Tool_Read_File(temp_name,"ivey",Num_Ivey,2,Temp_DATA,
     &                        Flag_Blank)
          Ive_y  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     ----------------------------------------------------------------------------------------------
c Check whether the initial acceleration file in the x-direction exists. If it exists, read it (for
c dynamic analysis).
c     ----------------------------------------------------------------------------------------------
      print *, "    Trying to read iacx files..." 
      temp_name = trim(trim(Full_Pathname)//'.iacx')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Iacx = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Iac_x)) DEALLOCATE(Iac_x)
          ALLOCATE(Iac_x(Num_Iacx,2))
          ALLOCATE(Temp_DATA(Num_Iacx,2))
          Call Tool_Read_File(temp_name,"iacx",Num_Iacx,2,Temp_DATA,
     &                        Flag_Blank)
          Iac_x  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ----------------------------------------------------------------------------------------------
c Check whether the initial acceleration file in the Y direction exists. If it exists, read it (for
c dynamic analysis).
c     ----------------------------------------------------------------------------------------------
      print *, "    Trying to read iacy files..." 
      temp_name = trim(trim(Full_Pathname)//'.iacy')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Iacy = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Iac_y)) DEALLOCATE(Iac_y)
          ALLOCATE(Iac_y(Num_Iacy,2))
          ALLOCATE(Temp_DATA(Num_Ivey,2))
          Call Tool_Read_File(temp_name,"iacy",Num_Iacy,2,Temp_DATA,
     &                        Flag_Blank)
          Iac_y  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     --------------------------------------------------------------------------
c Check the boundary conditions of the problem domain, added on August 11, 2016
c     --------------------------------------------------------------------------
      print *, "    Trying to read fbvl file..." 
      temp_name = trim(trim(Full_Pathname)//'.fbvl')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Fd_Bou_vl = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Fd_Bou_vl)) DEALLOCATE(Fd_Bou_vl)
          ALLOCATE(Fd_Bou_vl(Num_Fd_Bou_vl,2))
          ALLOCATE(Temp_DATA(Num_Fd_Bou_vl,2))
          Call Tool_Read_File(temp_name,"fbvl",Num_Fd_Bou_vl,2,
     &                        Temp_DATA,Flag_Blank)
          Fd_Bou_vl  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
          
          ! Read the number of boundary values that are 0 and their corresponding node numbers, as well as the
          ! number of non-zero nodes and their corresponding values
          Num_Fd_Bou_fixed = 0
          Num_Fd_Bou_vl_nz = 0
          do i_Check =1,Num_Fd_Bou_vl
              if (abs(Fd_Bou_vl(i_Check,2))<=Tol_30) then
                  Num_Fd_Bou_fixed = Num_Fd_Bou_fixed +1
              else
                  Num_Fd_Bou_vl_nz  = Num_Fd_Bou_vl_nz + 1 
              endif
          enddo
          IF(ALLOCATED(Fd_Bou_fixed)) DEALLOCATE(Fd_Bou_fixed)
          ALLOCATE(Fd_Bou_fixed(Num_Fd_Bou_fixed))
          IF(ALLOCATED(Fd_Bou_vl_nz)) DEALLOCATE(Fd_Bou_vl_nz)
          ALLOCATE(Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,2))
          Num_Fd_Bou_fixed = 0
          Num_Fd_Bou_vl_nz = 0
          do i_Check =1,Num_Fd_Bou_vl
              if (abs(Fd_Bou_vl(i_Check,2))<=Tol_30) then
                  Num_Fd_Bou_fixed = Num_Fd_Bou_fixed +1
                  Fd_Bou_fixed(Num_Fd_Bou_fixed) = 
     &                                     int(Fd_Bou_vl(i_Check,1))   
              else
                  Num_Fd_Bou_vl_nz  = Num_Fd_Bou_vl_nz + 1 
                  Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,1)= Fd_Bou_vl(i_Check,1) 
                  Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,2)= Fd_Bou_vl(i_Check,2)
              endif
          enddo
      endif 
      
c     -------------------------------------------------------------------
c Check the initial value file of the field problem, added on 2016-11-08
c     -------------------------------------------------------------------
      print *, "    Trying to read .fbiv files..." 
      temp_name = trim(trim(Full_Pathname)//'.fbiv')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Fd_Ini_vl = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Fd_Ini_vl)) DEALLOCATE(Fd_Ini_vl)
          ALLOCATE(Fd_Ini_vl(Num_Fd_Ini_vl,2))
          ALLOCATE(Temp_DATA(Num_Fd_Ini_vl,2))
          Call Tool_Read_File(temp_name,"fbvl",Num_Fd_Ini_vl,2,
     &                        Temp_DATA,Flag_Blank)
          Fd_Ini_vl  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     ------------------------------------------------------------------------------
c Check the flux through the boundary normal of the field, added on August 12, 2016
c     ------------------------------------------------------------------------------
      print *, "    Trying to read fbqn files..." 
      temp_name = trim(trim(Full_Pathname)//'.fbqn')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
      else
          Num_Fd_Bou_qn = Tool_Count_Lines(temp_name) 
          IF(ALLOCATED(Fd_Bou_qn)) DEALLOCATE(Fd_Bou_qn)
          ALLOCATE( Fd_Bou_qn(Num_Fd_Bou_qn,2))
          ALLOCATE( Temp_DATA(Num_Fd_Bou_qn,2))
          Call Tool_Read_File(temp_name,"fbqn",Num_Fd_Bou_qn,2,
     &                        Temp_DATA,Flag_Blank)
          Fd_Bou_qn  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     -------------------------------------------------------------------------------------------
c Read the transient wellbore pressure variation curve data (only for analysis of No. 16 and 17)
c     -------------------------------------------------------------------------------------------
      if(Key_Analysis_Type ==16 .or. Key_Analysis_Type ==17)then
        if(Key_Gas_Production==1 .and. Key_Changing_BHP==1) then
          print *, "    Warning :: Key_Gas_Production =1" 
          print *, "    Warning :: Key_Changing_BHP =1" 
          print *, "    Trying to read bhpc file..." 
          temp_name = trim(trim(Full_Pathname)//'.bhpc')
          inquire(file=temp_name, exist=alive)  
          if(alive.EQV..FALSE.)then
              print *, "    Warning :: Can not find bhpc files!!!" 
              call Warning_Message('S',Keywords_Blank)
          else
              Num_BHP_curve_point = Tool_Count_Lines(temp_name) 
              IF(ALLOCATED(BHP_Curve)) DEALLOCATE(BHP_Curve)
              ALLOCATE( BHP_Curve(Num_BHP_curve_point,2))
              ALLOCATE( Temp_DATA(Num_BHP_curve_point,2))
              Call Tool_Read_File(temp_name,"bhpc",Num_BHP_curve_point,
     &                            2,Temp_DATA,Flag_Blank)
              BHP_Curve  = Temp_DATA(:,:)
              DEALLOCATE(Temp_DATA)
          endif
        endif
      endif
Cc     --------------------------------------------------      
C c     Check the flux in the x-direction at the problem boundary, added on August 11, 2016
Cc     --------------------------------------------------
C     print *, "    Trying to read fbqx files..." 
C     temp_name = trim(trim(Full_Pathname)//'.fbqx')
C     inquire(file=temp_name, exist=alive)  
C     if(alive.EQV..FALSE.)then
C         print *, "    Warning :: Can not find fbqx files!!!" 
C     else
C         Num_Fd_Bou_qx = Tool_Count_Lines(temp_name) 
C         ALLOCATE( Fd_Bou_qx(Num_Fd_Bou_qx,2))
C         ALLOCATE( Temp_DATA(Num_Fd_Bou_qx,2))
C         Call Tool_Read_File(temp_name,"fbqx",Num_Fd_Bou_qx,2,
C    &                        Temp_DATA,Flag_Blank)
C         Fd_Bou_qx  = Temp_DATA(:,:)
C         DEALLOCATE(Temp_DATA)
C     endif 
      
Cc     --------------------------------------------------      
C c     Check the flux in the x-direction at the problem boundary, added on August 11, 2016
Cc     --------------------------------------------------
C     print *, "    Trying to read fbqy files..." 
C     temp_name = trim(trim(Full_Pathname)//'.fbqy')
C     inquire(file=temp_name, exist=alive)  
C     if(alive.EQV..FALSE.)then
C         print *, "    Warning :: Can not find fbqy files!!!" 
C     else
C         Num_Fd_Bou_qx = Tool_Count_Lines(temp_name) 
C         ALLOCATE( Fd_Bou_qy(Num_Fd_Bou_qy,2))
C         ALLOCATE( Temp_DATA(Num_Fd_Bou_qy,2))
C         Call Tool_Read_File(temp_name,"fbqx",Num_Fd_Bou_qy,2,
C    &                        Temp_DATA,Flag_Blank)
C         Fd_Bou_qy  = Temp_DATA(:,:)
C         DEALLOCATE(Temp_DATA)
C     endif 
c     --------------------------------------------------------------------
c Determine the outer boundary of all quadrilateral elements of the model
c     --------------------------------------------------------------------
      ! The four sides of the element are stored in the temporary variable tem1
      print *, "    Finding model boundary..." 
      ALLOCATE(tem1(4*Num_Elem,2))
      ALLOCATE(tem2(4*Num_Elem))
      tem1(1:Num_Elem,1)              = Elem_Node(1:Num_Elem,1)
      tem1(1:Num_Elem,2)              = Elem_Node(1:Num_Elem,2)
      tem1(Num_Elem+1:2*Num_Elem,1)   = Elem_Node(1:Num_Elem,2)
      tem1(Num_Elem+1:2*Num_Elem,2)   = Elem_Node(1:Num_Elem,3)
      tem1(2*Num_Elem+1:3*Num_Elem,1) = Elem_Node(1:Num_Elem,3)
      tem1(2*Num_Elem+1:3*Num_Elem,2) = Elem_Node(1:Num_Elem,4)
      tem1(3*Num_Elem+1:4*Num_Elem,1) = Elem_Node(1:Num_Elem,4)
      tem1(3*Num_Elem+1:4*Num_Elem,2) = Elem_Node(1:Num_Elem,1)
      ! Sort
      call Matrix_Sort_Int(4*Num_Elem,2,tem1)
      ! Count the occurrences of each row in the matrix and store them in tem2
      !call Matrix_Count_Row_Int(4*Num_Elem,2,tem1,tem2,all_num_Outline)
      call Matrix_Count_Row_Int_m_x_2(4*Num_Elem,
     &                                tem1,tem2,all_num_Outline,0)
      ! Extract edges that appear only once, including both outer and inner edges
      ALLOCATE(All_Outline(all_num_Outline,2))
      ALLOCATE(Temp_Outline(all_num_Outline,2))
      j=0
      do i=1,4*Num_Elem
          if (tem2(i).eq.1)then
              j=j+1
              all_Outline(j,:) = tem1(i,:)
          end if
      end do
      ! Connected end to end, containing only the outer boundary
      !???????????????????????????????????????????
      !               ToDoList TDL
      ! ----------------------------------------
      ! There is a bug that needs to be fixed here; if the model has holes,
      ! Then the outer boundary obtained by Tools_Sort_by_End_to_End
      ! It's two duplicate copies, like the example of hole_regular_mesh.
      !???????????????????????????????????????????
      call Tool_Sort_by_End_to_End(all_num_Outline,all_num_Outline,
     &                             All_Outline,Temp_Outline,
     &                             Out_num_Outline)
      IF(ALLOCATED(Outline)) DEALLOCATE(Outline)
      ALLOCATE(Outline(Out_num_Outline,2))
      do i=1,Out_num_Outline
          Outline(i,1:2) = Temp_Outline(i,1:2) 
      end do
      DEALLOCATE(All_Outline)
      DEALLOCATE(Temp_Outline)
c     ------------------------------------------------------------------------------------------
c Coordinates of the four nodes of the storage unit, calculation of the element's average area,
c element side length, etc.
c     ------------------------------------------------------------------------------------------
      print *, "    Get node and element data..." 
      IF(ALLOCATED(G_NN)) DEALLOCATE(G_NN)
      ALLOCATE(G_NN(4,Num_Elem))     
      IF(ALLOCATED(G_X_NODES)) DEALLOCATE(G_X_NODES)
      ALLOCATE(G_X_NODES(4,Num_Elem))
      IF(ALLOCATED(G_Y_NODES)) DEALLOCATE(G_Y_NODES)
      ALLOCATE(G_Y_NODES(4,Num_Elem))
      IF(ALLOCATED(Elem_Area)) DEALLOCATE(Elem_Area)
      ALLOCATE(Elem_Area(Num_Elem))
      IF(ALLOCATED(Elem_Centroid)) DEALLOCATE(Elem_Centroid)
      ALLOCATE(Elem_Centroid(Num_Elem,2))
      IF(ALLOCATED(EleGaus_yes_FEM_asemd)) then
          DEALLOCATE(EleGaus_yes_FEM_asemd)
      endif
      ALLOCATE(EleGaus_yes_FEM_asemd(Num_Elem,Num_Gauss_Points))
      IF(ALLOCATED(x_max_Elements)) DEALLOCATE(x_max_Elements)
      ALLOCATE(x_max_Elements(Num_Elem)) 
      IF(ALLOCATED(x_min_Elements)) DEALLOCATE(x_min_Elements)
      ALLOCATE(x_min_Elements(Num_Elem)) 
      IF(ALLOCATED(y_max_Elements)) DEALLOCATE(y_max_Elements)
      ALLOCATE(y_max_Elements(Num_Elem)) 
      IF(ALLOCATED(y_min_Elements)) DEALLOCATE(y_min_Elements)
      ALLOCATE(y_min_Elements(Num_Elem)) 
      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points) =.false.
      G_NN(1:4,1:Num_Elem)  = 0
      Min_Ele_Edge_Length = TEN_10
      Max_Ele_Edge_Length = ZR
      do i=1,Num_Elem
          N1  = Elem_Node(i,1)                                    
          N2  = Elem_Node(i,2)                                    
          N3  = Elem_Node(i,3)                                    
          N4  = Elem_Node(i,4)                                    
          NN  = [N1,N2,N3,N4]                                           
          G_NN(1:4,i)      = NN
          G_X_NODES(1:4,i) = Coor(NN,1)
          G_Y_NODES(1:4,i) = Coor(NN,2) 
          call Tool_Area_Polygon(Coor(NN,1),Coor(NN,2),4,area)
          Elem_Area(i) = area
          call Tool_Centroid_Quad(Coor(NN,1),Coor(NN,2),Centroid)
          Elem_Centroid(i,1:2) = Centroid
          ! element side length
          c_L1= sqrt((Coor(N1,1)-Coor(N2,1))**2+
     &               (Coor(N1,2)-Coor(N2,2))**2)
          c_L2= sqrt((Coor(N2,1)-Coor(N3,1))**2+
     &               (Coor(N2,2)-Coor(N3,2))**2)
          c_L3= sqrt((Coor(N3,1)-Coor(N4,1))**2+
     &               (Coor(N3,2)-Coor(N4,2))**2)
          c_L4= sqrt((Coor(N4,1)-Coor(N1,1))**2+
     &               (Coor(N4,2)-Coor(N1,2))**2)
          c_L_min = min(c_L1,c_L2,c_L3,c_L4)
          c_L_max = max(c_L1,c_L2,c_L3,c_L4)
          if (c_L_min< Min_Ele_Edge_Length)then
              Min_Ele_Edge_Length = c_L_min
          endif
          if (c_L_max> Max_Ele_Edge_Length)then
              Max_Ele_Edge_Length = c_L_max
          endif
          ! The range of x, y, z coordinates for each element
          x_max_Elements(i) = maxval(G_X_NODES(1:4,i))
          x_min_Elements(i) = minval(G_X_NODES(1:4,i))
          y_max_Elements(i) = maxval(G_Y_NODES(1:4,i))
          y_min_Elements(i) = minval(G_Y_NODES(1:4,i))
      end do
      Max_Elem_Area = MaxVal(Elem_Area)
      Min_Elem_Area = MinVal(Elem_Area)
      Ave_Elem_Area = sum(Elem_Area)/size(Elem_Area)
      Ave_Elem_L    = sqrt(Ave_Elem_Area)
      
c     -----------------------------------------------------
c The elements around each node and the number of elements
c     -----------------------------------------------------
      IF(ALLOCATED(Node_Elements)) then
          DEALLOCATE(Node_Elements)
      endif
      ALLOCATE(Node_Elements(Num_Node,8))
      Node_Elements = 0
      IF(ALLOCATED(num_Node_Elements)) DEALLOCATE(num_Node_Elements)      
      ALLOCATE(num_Node_Elements(Num_Node))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)          
      do i=1,Num_Node
          num_Node_Elements(i) = 0
          do j=1,Num_Elem
              if (any(Elem_Node(j,:) == i))then
                  num_Node_Elements(i) = num_Node_Elements(i) +1
                  Node_Elements(i,num_Node_Elements(i)) = j
              end if
            end do
      end do      
!$omp end parallel do   
c     -------------------------------------------------------------
c The elements surrounding each element and the number of elements
c     -------------------------------------------------------------
      IF(ALLOCATED(Ele_Elements)) DEALLOCATE(Ele_Elements)
      ALLOCATE(Ele_Elements(Num_Elem,25))
      IF(ALLOCATED(num_Ele_Eles)) DEALLOCATE(num_Ele_Eles)      
      ALLOCATE(num_Ele_Eles(Num_Elem))
      Ele_Elements(1:Num_Elem,1:25) = 0
      num_Ele_Eles(1:Num_Elem) = 0
      !Loop through each element.
      do i_E =1,Num_Elem
          c_all_num_Ele      = 0
          c_all_Elem(1:32)   = 0
          !Loop through each node of the element.
          do j=1,4
              c_Node = Elem_Node(i_E,j) 
              c_all_Elem(c_all_num_Ele + 1:
     &                   c_all_num_Ele + num_Node_Elements(c_Node)) =  
     &               Node_Elements(c_Node,1:num_Node_Elements(c_Node)) 
              c_all_num_Ele = c_all_num_Ele + num_Node_Elements(c_Node)
            end do
          ! Delete duplicate unit numbers
          call Vector_Unique_Int(32,c_all_num_Ele,c_all_Elem,
     &                           Uniqued_all_Elem,num_Surrou)  
          Ele_Elements(i_E,1:num_Surrou)=Uniqued_all_Elem(1:num_Surrou)
          num_Ele_Eles(i_E) = num_Surrou
      end do
c     ------------------------------------------------------------------------------------------
c Determine dimensions for relevant public variables (mainly related to identifying enhancement
c nodes)
c     ------------------------------------------------------------------------------------------
      IF(ALLOCATED(Elem_Type)) DEALLOCATE(Elem_Type)
      ALLOCATE(Elem_Type(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Enriched_Node_Type)) DEALLOCATE(Enriched_Node_Type)
      ALLOCATE(Enriched_Node_Type(Num_Node,Max_Num_Cr))
      IF(ALLOCATED(Node_Jun_elem)) DEALLOCATE(Node_Jun_elem)
      ALLOCATE(Node_Jun_elem(Num_Node,Max_Num_Cr))
      IF(ALLOCATED(Jun_Ele_Negative_Cr_Num)) then
          DEALLOCATE(Jun_Ele_Negative_Cr_Num)
      endif
      ALLOCATE(Jun_Ele_Negative_Cr_Num(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Node_Jun_Hole)) DEALLOCATE(Node_Jun_Hole)
      ALLOCATE(Node_Jun_Hole(Num_Node,Max_Num_Cr))
      IF(ALLOCATED(Ele_Jun_Hole)) DEALLOCATE(Ele_Jun_Hole)
      ALLOCATE(Ele_Jun_Hole(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Coors_Element_Crack)) DEALLOCATE(Coors_Element_Crack)
      ALLOCATE(Coors_Element_Crack(Num_Elem,Max_Num_Cr,4))
      IF(ALLOCATED(Coors_Tip)) DEALLOCATE(Coors_Tip)
      ALLOCATE(Coors_Tip(Num_Elem,2))
      IF(ALLOCATED(Coors_Vertex)) DEALLOCATE(Coors_Vertex)
      ALLOCATE(Coors_Vertex(Num_Elem,2))
      IF(ALLOCATED(Coors_Junction)) DEALLOCATE(Coors_Junction)
      ALLOCATE(Coors_Junction(Num_Elem,Max_Num_Cr,4))
      IF(ALLOCATED(x_cr_tip_nodes)) DEALLOCATE(x_cr_tip_nodes)
      ALLOCATE(x_cr_tip_nodes(Max_Num_Cr,Num_Node))
      IF(ALLOCATED(y_cr_tip_nodes)) DEALLOCATE(y_cr_tip_nodes)
      ALLOCATE(y_cr_tip_nodes(Max_Num_Cr,Num_Node))
      IF(ALLOCATED(Ele_Num_Tip_Enriched_Node)) then 
          DEALLOCATE(Ele_Num_Tip_Enriched_Node)
      endif
      ALLOCATE(Ele_Num_Tip_Enriched_Node(Max_Num_Cr,Num_Node))
      IF(ALLOCATED(c_POS)) DEALLOCATE(c_POS)
      ALLOCATE(c_POS(Num_Node,Max_Num_Cr))
      c_POS = 0
      IF(ALLOCATED(TipEle_Adjacent_Ele)) DEALLOCATE(TipEle_Adjacent_Ele)
      ALLOCATE(TipEle_Adjacent_Ele(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(num_GP_Elem)) DEALLOCATE(num_GP_Elem)
      ALLOCATE(num_GP_Elem(Num_Elem))
      IF(ALLOCATED(Ele_GP_Start_Num)) DEALLOCATE(Ele_GP_Start_Num)
      ALLOCATE(Ele_GP_Start_Num(Num_Elem))
      ! ALLOCATE(Cracks_CalP_Elem_CalP(Max_Num_Cr, Num_Elem)) ! Starting index and number of calculation
      ! points contained in each element (relative to each crack, with local numbering for each crack)
      ! ---------Related to Holes-----------
      IF(ALLOCATED(Elem_Type_Hl)) DEALLOCATE(Elem_Type_Hl)
      ALLOCATE(Elem_Type_Hl(Num_Elem,Max_Num_Hl))
      IF(ALLOCATED(Enriched_Node_Type_Hl)) then
          DEALLOCATE(Enriched_Node_Type_Hl)
      endif
      ALLOCATE(Enriched_Node_Type_Hl(Num_Node,Max_Num_Hl))
      IF(ALLOCATED(c_POS_Hl)) DEALLOCATE(c_POS_Hl)
      ALLOCATE(c_POS_Hl(Num_Node,Max_Num_Hl))
      ! ---------Related Mixtures-----------
      IF(ALLOCATED(Elem_Type_Incl)) DEALLOCATE(Elem_Type_Incl)
      ALLOCATE(Elem_Type_Incl(Num_Elem,Max_Num_Incl))
      IF(ALLOCATED(Enriched_Node_Type_Incl)) then
          DEALLOCATE(Enriched_Node_Type_Incl)
      endif
      ALLOCATE(Enriched_Node_Type_Incl(Num_Node,Max_Num_Incl))
      IF(ALLOCATED(c_POS_Incl)) DEALLOCATE(c_POS_Incl)
      ALLOCATE(c_POS_Incl(Num_Node,Max_Num_Incl))
      ! ---------Cross-shaped Intersection Cracks Related-----
      IF(ALLOCATED(Cross_Point_Cr_num)) DEALLOCATE(Cross_Point_Cr_num)
      ALLOCATE(Cross_Point_Cr_num(Max_Num_Cross,2))
      IF(ALLOCATED(Cross_Point_Ele_num)) DEALLOCATE(Cross_Point_Ele_num)
      ALLOCATE(Cross_Point_Ele_num(Max_Num_Cross))
      ! ALLOCATE(Cross_Point_Cr_Seg_num(Max_Num_Cross,2))               ! Used for cross-shaped cracks
      IF(ALLOCATED(Cross_Point_RABCD)) DEALLOCATE(Cross_Point_RABCD)
      ALLOCATE(Cross_Point_RABCD(Max_Num_Cross,10,2))
      IF(ALLOCATED(Node_Cross_elem)) DEALLOCATE(Node_Cross_elem)
      ALLOCATE(Node_Cross_elem(Num_Node,Max_Num_Cross))
      IF(ALLOCATED(Elem_Type_Cross)) DEALLOCATE(Elem_Type_Cross)
      ALLOCATE(Elem_Type_Cross(Num_Elem,Max_Num_Cross))
      IF(ALLOCATED(Enriched_Node_Type_Cross)) then
          DEALLOCATE(Enriched_Node_Type_Cross)
      endif
      ALLOCATE(Enriched_Node_Type_Cross(Num_Node,Max_Num_Cross))
      IF(ALLOCATED(c_POS_Cross)) DEALLOCATE(c_POS_Cross)
      ALLOCATE(c_POS_Cross(Num_Node,Max_Num_Cross))
      
      ! HC_Point_Cr_num(num_HC,1:2) !The crack number and crack tip number corresponding to each HC
      ! intersection point
      ! HC_Point_Hole_num(num_HC)               !The hole number corresponding to each HC intersection point
      ! HC_Point_Ele_num(num_HC)            !The element number corresponding to each HC intersection point
      ! Elem_Type_HC(Num_Elem, num_HC); !Element type: 1--Heaviside enhanced elements additionally added
      ! at crack and hole intersections
      ! Enriched_Node_Type_HC(Num_Node, num_HC); !Type of enriched nodes: 1--Heaviside enriched nodes
      ! additionally added where cracks and holes intersect
      
      ! ---------Contact Issues Related-----------
      IF(ALLOCATED(Elem_Conta_Sta)) DEALLOCATE(Elem_Conta_Sta)
      ALLOCATE(Elem_Conta_Sta(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Elem_Coh_Sta)) DEALLOCATE(Elem_Coh_Sta)
      ALLOCATE(Elem_Coh_Sta(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Elem_Conta_Sta_Last)) then 
          DEALLOCATE(Elem_Conta_Sta_Last)
      endif
      ALLOCATE(Elem_Conta_Sta_Last(Num_Elem,Max_Num_Cr))
      IF(ALLOCATED(Ele_NumCalP)) DEALLOCATE(Ele_NumCalP)
      ALLOCATE(Ele_NumCalP(Num_Elem))
      IF(ALLOCATED(Ele_CalPNum)) DEALLOCATE(Ele_CalPNum)
      ALLOCATE(Ele_CalPNum(Num_Elem,Max_Num_Ele_CalP))
      IF(ALLOCATED(Elem_Proppant_Coor)) then
          DEALLOCATE(Elem_Proppant_Coor)
      endif
      ALLOCATE(Elem_Proppant_Coor(Num_Elem,5))
                                                            ! Among them: 1, indicates whether the element contains proppant
                                                            ! 2. Proppant Diameter
                                                            ! 3 and 4, x-coordinate and y-coordinate
                                                            ! 5. Corresponding crack segment dip angle
      ! ----------Related to Field Issues-----------
      IF(ALLOCATED(Fd_c_POS)) DEALLOCATE(Fd_c_POS)
      ALLOCATE(Fd_c_POS(Num_Node,Max_Num_Cr))
      IF(ALLOCATED(Fd_c_POS_Hl)) DEALLOCATE(Fd_c_POS_Hl)
      ALLOCATE(Fd_c_POS_Hl(Num_Node,Max_Num_Hl))
      IF(ALLOCATED(Fd_c_POS_Incl)) DEALLOCATE(Fd_c_POS_Incl)
      ALLOCATE(Fd_c_POS_Incl(Num_Node,Max_Num_Incl))
      IF(ALLOCATED(Fd_c_POS_Cross)) DEALLOCATE(Fd_c_POS_Cross)
      ALLOCATE(Fd_c_POS_Cross(Num_Node,Max_Num_Cross))
      IF(ALLOCATED(Fd_EleGaus_yes_FEM_asemd)) then
          DEALLOCATE(Fd_EleGaus_yes_FEM_asemd)
      endif
      ALLOCATE(Fd_EleGaus_yes_FEM_asemd(Num_Elem,Num_Gauss_Points))
      ! ----------Unit Destruction-------
      IF(ALLOCATED(Elem_Break)) DEALLOCATE(Elem_Break)
      ALLOCATE(Elem_Break(Num_Elem))
      IF(ALLOCATED(Elem_Ave_Gauss_Stress)) then
          DEALLOCATE(Elem_Ave_Gauss_Stress)
      endif
      ALLOCATE(Elem_Ave_Gauss_Stress(Num_Elem,3))
      IF(ALLOCATED(Elem_Ave_Gauss_S1)) DEALLOCATE(Elem_Ave_Gauss_S1)
      ALLOCATE(Elem_Ave_Gauss_S1(Num_Elem))
      Elem_Ave_Gauss_Stress(1:Num_Elem,1:3) = ZR 
      Elem_Ave_Gauss_S1(1:Num_Elem) = ZR
      Elem_Break(1:Num_Elem) = .False.
      
      !----------------------------------------------------
      ! 2D crack related (changed to Allocatable variable)
      !        2022-09-02
      !        IMPROV2022090201.
      !----------------------------------------------------
      IF(ALLOCATED(Cracks_CalP_Num)) DEALLOCATE(Cracks_CalP_Num)
      ALLOCATE(Cracks_CalP_Num(Max_Num_Cr))
      IF(ALLOCATED(Cracks_CalP_Coors)) DEALLOCATE(Cracks_CalP_Coors)
      ALLOCATE(Cracks_CalP_Coors(Max_Num_Cr,Max_Num_Cr_CalP,2))
      IF(ALLOCATED(Cracks_CalP_Orient)) DEALLOCATE(Cracks_CalP_Orient)
      ALLOCATE(Cracks_CalP_Orient(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Seg)) DEALLOCATE(Cracks_CalP_Seg)
      ALLOCATE(Cracks_CalP_Seg(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Elem)) DEALLOCATE(Cracks_CalP_Elem)
      ALLOCATE(Cracks_CalP_Elem(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Aper)) DEALLOCATE(Cracks_CalP_Aper)
      ALLOCATE(Cracks_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Pres)) DEALLOCATE(Cracks_CalP_Pres)
      ALLOCATE(Cracks_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Tractions)) then
          DEALLOCATE(Cracks_CalP_Tractions)
      endif
      ALLOCATE(Cracks_CalP_Tractions(Max_Num_Cr,Max_Num_Cr_CalP,2))
      IF(ALLOCATED(Cracks_CalP_Pgra)) DEALLOCATE(Cracks_CalP_Pgra)
      ALLOCATE(Cracks_CalP_Pgra(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Velo)) DEALLOCATE(Cracks_CalP_Velo)
      ALLOCATE(Cracks_CalP_Velo(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Quan)) DEALLOCATE(Cracks_CalP_Quan)
      ALLOCATE(Cracks_CalP_Quan(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Conc)) DEALLOCATE(Cracks_CalP_Conc)
      ALLOCATE(Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Remo_Strs)) then
          DEALLOCATE(Cracks_CalP_Remo_Strs)
      endif
      ALLOCATE(Cracks_CalP_Remo_Strs(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Contact_Strs)) then
          DEALLOCATE(Cracks_CalP_Contact_Strs)
      endif
      ALLOCATE(Cracks_CalP_Contact_Strs(Max_Num_Cr,Max_Num_Cr_CalP,2))
      IF(ALLOCATED(Cracks_CalP_Contact_Force_x)) then
          DEALLOCATE(Cracks_CalP_Contact_Force_x)
      endif
      ALLOCATE(Cracks_CalP_Contact_Force_x(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Contact_Force_y)) then
          DEALLOCATE(Cracks_CalP_Contact_Force_y)
      endif
      ALLOCATE(Cracks_CalP_Contact_Force_y(Max_Num_Cr,Max_Num_Cr_CalP))
      ! -----Added on 2016-08-25--Related to the calculation of support crack width and flow
      ! capacity-------
      IF(ALLOCATED(Cracks_CalP_wpnp)) DEALLOCATE(Cracks_CalP_wpnp)
      ALLOCATE(Cracks_CalP_wpnp(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_wpor)) DEALLOCATE(Cracks_CalP_wpor)
      ALLOCATE(Cracks_CalP_wpor(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_wdeform)) DEALLOCATE(Cracks_CalP_wdeform)
      ALLOCATE(Cracks_CalP_wdeform(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_Conductivity)) then
          DEALLOCATE(Cracks_CalP_Conductivity)
      endif
      ALLOCATE(Cracks_CalP_Conductivity(Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(Cracks_CalP_kf)) DEALLOCATE(Cracks_CalP_kf)
      ALLOCATE(Cracks_CalP_kf(Max_Num_Cr,Max_Num_Cr_CalP))
      ! real(kind=FT), ALLOCATABLE :: Cracks_CalP_Elem_CalP(:,:) ! Starting index and number of
      ! calculation points contained in each element (relative to each crack, local numbering for each
      ! crack)
      ! -----------The following variables are related to staged fracturing----------
      IF(ALLOCATED(MS_Cracks_CalP_Num)) DEALLOCATE(MS_Cracks_CalP_Num)
      ALLOCATE(MS_Cracks_CalP_Num(20,Max_Num_Cr))
      IF(ALLOCATED(MS_Cracks_CalP_Aper)) DEALLOCATE(MS_Cracks_CalP_Aper)
      ALLOCATE(MS_Cracks_CalP_Aper(20,Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(MS_Cracks_CalP_Conc)) DEALLOCATE(MS_Cracks_CalP_Conc)
      ALLOCATE(MS_Cracks_CalP_Conc(10,Max_Num_Cr,Max_Num_Cr_CalP))
      IF(ALLOCATED(MS_CalP_Propped_Aper)) then
          DEALLOCATE(MS_CalP_Propped_Aper)
      endif
      ALLOCATE(MS_CalP_Propped_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
      ! Initialize calculation point-related data
      Cracks_CalP_Pres(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
        Cracks_CalP_Tractions(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2)= ZR
      Cracks_CalP_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
      Cracks_CalP_Coors(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2)= ZR
      Cracks_CalP_Orient(1:Max_Num_Cr,1:Max_Num_Cr_CalP)   = ZR
      Cracks_CalP_Seg(1:Max_Num_Cr,1:Max_Num_Cr_CalP)      = 0
      Cracks_CalP_Elem(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = 0
      Cracks_CalP_Pgra(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
      Cracks_CalP_Velo(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
      Cracks_CalP_Quan(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
      Cracks_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = ZR
      Cracks_CalP_Remo_Strs(1:Max_Num_Cr,1:Max_Num_Cr_CalP)= ZR
      Cracks_CalP_wdeform(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  = ZR
      MS_CalP_Propped_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
      
      
      
c     ------------------------------------------------------------------------------------------
c Hydraulic fracturing related: Initialize the hydraulic fracturing status of all cracks (i.e.,
c whether they are driven by water)
c     ------------------------------------------------------------------------------------------
C     Cracks_HF_State(1:Max_Num_Cr)  = 0
C Cracks_HF_State(Inject_Crack_Num) = 1
      
      RETURN
      END SUBROUTINE Read_Geo