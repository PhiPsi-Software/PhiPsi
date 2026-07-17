!-----------------------------------------------------------
! Brief: Read the 2D geometry and problem-definition files.
!
! Notes:   Reads the node, element, boundary, focus, material, and crack
!          data files for 2D FEM/XFEM runs and populates the global model
!          modules used by downstream solvers.
!-----------------------------------------------------------

subroutine Read_Geo

!     ------------------------
!     Read public variable module
!     ------------------------
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

!     ----------------------
!     Variable Type Declaration
!     ----------------------
implicit none
logical :: alive
character(len=200) :: temp_name
integer :: Tool_Count_Lines
integer temp_count,c_i,c_j

!     ---------------
!     Temporary variable
!     ---------------
real(kind=FT),allocatable::Temp_DATA(:,:)
logical :: Flag_Blank
integer,allocatable::tem1(:,:)
integer,allocatable::tem2(:)
integer i,j,all_num_Outline,Out_num_Outline
integer,allocatable::All_Outline(:,:)
integer,allocatable::Temp_Outline(:,:)  
integer N1,N2,N3,N4,NN(4)
real(kind=FT) area,Centroid(2)
integer c_all_num_Ele,c_Node,i_E
integer c_all_Elem(32),Uniqued_all_Elem(32)
integer :: num_Surrou
integer :: i_Check
real(kind=FT) c_L1,c_L2,c_L3,c_L4,c_L_min,c_L_max
integer i_C,n_pts,i_pt

!     ------------------------------
!     Check the file name and path name
!     ------------------------------
if (Filename=='*blank*')then
    print *, "    Error :: Cannot find Keyword *Filename!"
    call Warning_Message('S',Keywords_Blank) 
end if
if (Work_Directory == '*blank*') then
    print *, "    Error :: Cannot find Keyword *Work_Directory!"
    call Warning_Message('S',Keywords_Blank) 
end if

!     ---------------------------------------------------------------
!     Check if the node coordinate file exists, and if it does, read it.
!     ---------------------------------------------------------------
print *, "    Trying to read nodal files...." 
temp_name = trim(trim(Full_Pathname)//'.node')
inquire(file=temp_name, exist=alive)
if(alive.eqv..false.)then
    print *, "    Error :: Can not find nodal file!!!" 
    print *, "    Check whether the following file exists or not:"
    print *, "    ",temp_name
    call Warning_Message('S',Keywords_Blank) 
else
    ! Count the number of lines in a file
    Num_Node = Tool_Count_Lines(temp_name)
    if (ALLOCATED(Coor)) deallocate(Coor)
    allocate(Coor(Num_Node,2))
    allocate(Temp_DATA(Num_Node,2))
    Call Tool_Read_File(temp_name,"node",Num_Node,2,Temp_DATA, Flag_Blank)
    Coor = Temp_DATA
    deallocate(Temp_DATA)
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

!     -------------------------------------------------------------------------------
!     Check whether the input file is a three-dimensional problem, and if so, terminate.
!     -------------------------------------------------------------------------------
temp_name = trim(trim(Full_Pathname)//'.bouz')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..true.)then
    print *, "    WARNING :: *.bouz file found, please check *Key_Dimension (2 or 3)!" 
    !call Warning_Message('S',Keywords_Blank) 
endif  
temp_name = trim(trim(Full_Pathname)//'.focz')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..true.)then
    print *, "    WARNING :: *.focz file found, please check *Key_Dimension (2 or 3)!" 
    !call Warning_Message('S',Keywords_Blank) 
endif 


!     --------------------------------------------------------------
!     Check if the element number file exists, and if it does, read it.
!     --------------------------------------------------------------
print *, "    Trying to read element files..." 
temp_name = trim(trim(Full_Pathname)//'.elem')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
    print *, "    Error :: Can not find element file!" 
    call Warning_Message('S',Keywords_Blank) 
else
    Num_Elem = Tool_Count_Lines(temp_name)
    if (ALLOCATED(Elem_Node)) deallocate(Elem_Node)
    allocate(Elem_Node(Num_Elem,4))
    if (ALLOCATED(Elem_Mat)) deallocate(Elem_Mat)
    allocate(Elem_Mat(Num_Elem)) 
    allocate(Temp_DATA(Num_Elem,5))
    Call Tool_Read_File(temp_name,"elem",Num_Elem,5,Temp_DATA, Flag_Blank)
    Elem_Node = int(Temp_DATA(:,1:4))
    Elem_Mat  = int(Temp_DATA(:,5))
    ! num_of_Material = MaxVal(Elem_Mat)     ! Number of material types
    deallocate(Temp_DATA)
endif  

!     ---------------------------------------------------------------------------------------------
!     Obtain the maximum node number difference in the element to calculate the maximum half-bandwidth
!     A first course in the finite element method, 5ed, B.4.
!     ---------------------------------------------------------------------------------------------
Max_Diff_Elem_Num = maxval(maxval(Elem_Node(:,1:4),2)- &
minval(Elem_Node(:,1:4),2))
Max_Half_Band_Width = 2*(Max_Diff_Elem_Num + 1)

!     --------------------------------------------------------------------
!     Check if the x boundary constraint file exists, and if it does, read it
!     --------------------------------------------------------------------
print *, "    Trying to read boux files..." 
temp_name = trim(trim(Full_Pathname)//'.boux')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
    print *, "    Warning :: Can not find boux file!" 
else
    Num_Bou_x = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Bou_x)) deallocate(Bou_x)
    allocate(Bou_x(Num_Bou_x))
    allocate(Temp_DATA(Num_Bou_x,1))
    Call Tool_Read_File(temp_name,"boux",Num_Bou_x,1,Temp_DATA, Flag_Blank)
    Bou_x  = int(Temp_DATA(:,1))
    deallocate(Temp_DATA)
endif  

!     --------------------------------------------------------------------
!     Check if the y-boundary constraint file exists, and if it does, read it
!     --------------------------------------------------------------------
print *, "    Trying to read bouy files..." 
temp_name = trim(trim(Full_Pathname)//'.bouy')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
    print *, "    Warning :: Can not find bouy file!" 
else
    Num_Bou_y = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Bou_y)) deallocate(Bou_y)
    allocate(Bou_y(Num_Bou_y))
    allocate(Temp_DATA(Num_Bou_y,1))
    Call Tool_Read_File(temp_name,"bouy",Num_Bou_y,1,Temp_DATA, Flag_Blank)
    Bou_y  = int(Temp_DATA(:,1))
    deallocate(Temp_DATA)
endif    

!     -----------------------------------------------------------------------
!     Check if the load file in the x-direction exists, and if it does, read it.
!     -----------------------------------------------------------------------
print *, "    Trying to read focx files..." 
temp_name = trim(trim(Full_Pathname)//'.focx')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
    print *, "    Warning :: Can not find focx file!" 
else
    Num_Foc_x = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Foc_x)) deallocate(Foc_x)
    allocate(Foc_x(Num_Foc_x,2))
    allocate(Temp_DATA(Num_Foc_x,2))
    Call Tool_Read_File(temp_name,"focx",Num_Foc_x,2,Temp_DATA, Flag_Blank)
    Foc_x  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif  

!     -----------------------------------------------------------------------
!     Check if the load file in the y-direction exists, and if it does, read it.
!     -----------------------------------------------------------------------
print *, "    Trying to read focy files..." 
temp_name = trim(trim(Full_Pathname)//'.focy')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
    print *, "    Warning :: Can not find focy file!" 
else
    Num_Foc_y = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Foc_y)) deallocate(Foc_y)
    allocate(Foc_y(Num_Foc_y,2))
    allocate(Temp_DATA(Num_Foc_y,2))
    Call Tool_Read_File(temp_name,"focy",Num_Foc_y,2,Temp_DATA, Flag_Blank)
    Foc_y  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     -----------------------------------------------------------------------------------------------
!     Check whether the boundary condition file for non-zero displacement in the x direction exists, and
!     if it does, read it.
!     -----------------------------------------------------------------------------------------------
print *, "    Trying to read buxn files..." 
temp_name = trim(trim(Full_Pathname)//'.buxn')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Boux_nonzero = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Bou_x_nonzero)) deallocate(Bou_x_nonzero)
    allocate(Bou_x_nonzero(Num_Boux_nonzero,2))
    allocate(Temp_DATA(Num_Boux_nonzero,2))
    Call Tool_Read_File(temp_name,"buxn",Num_Boux_nonzero,2, Temp_DATA,Flag_Blank)
    Bou_x_nonzero  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif  

!     -----------------------------------------------------------------------------------------------
!     Check whether the boundary condition file for non-zero displacement in the y direction exists, and
!     if it does, read it.
!     -----------------------------------------------------------------------------------------------
print *, "    Trying to read buyn files..." 
temp_name = trim(trim(Full_Pathname)//'.buyn')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Bouy_nonzero = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Bou_y_nonzero)) deallocate(Bou_y_nonzero)
    allocate(Bou_y_nonzero(Num_Bouy_nonzero,2))
    allocate(Temp_DATA(Num_Bouy_nonzero,2))
    Call Tool_Read_File(temp_name,"buyn",Num_Bouy_nonzero,2, Temp_DATA,Flag_Blank)
    Bou_y_nonzero  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     -----------------------------------------------------------------------------------------------
!     Check whether the coupling file for the x-direction degree of freedom exists, and if it does, read
!     it.
!     -----------------------------------------------------------------------------------------------
print *, "    Trying to read dofx files..." 
temp_name = trim(trim(Full_Pathname)//'.dofx ')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    temp_count = Tool_Count_Lines(temp_name) 
    allocate(Temp_DATA(temp_count,2))
    Call Tool_Read_File(temp_name,"dofx",temp_count,2, Temp_DATA,Flag_Blank)
    num_CP_set_x = int(maxval(Temp_DATA(:,1)));
    do c_i =1,temp_count
        do c_j=1,num_CP_set_x
            if (Temp_DATA(c_i,1) == c_j)then
                num_nodes_CP_set_x(c_j) =num_nodes_CP_set_x(c_j)+1
                CP_nodes_x(c_j,num_nodes_CP_set_x(c_j))= int(Temp_DATA(c_i,2))
            endif
        enddo
    enddo
    deallocate(Temp_DATA)
endif  

!     --------------------------------------------------------------------------------------------
!     Check whether the y-direction degree of freedom coupling file exists; if it does, then read it.
!     --------------------------------------------------------------------------------------------
print *, "    Trying to read dofy files..." 
temp_name = trim(trim(Full_Pathname)//'.dofy ')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    temp_count = Tool_Count_Lines(temp_name) 
    allocate(Temp_DATA(temp_count,2))
    Call Tool_Read_File(temp_name,"dofy",temp_count,2, Temp_DATA,Flag_Blank)
    num_CP_set_y = int(maxval(Temp_DATA(:,1)));
    do c_i =1,temp_count
        do c_j=1,num_CP_set_y
            if (Temp_DATA(c_i,1) == c_j)then
                num_nodes_CP_set_y(c_j) =num_nodes_CP_set_y(c_j)+1
                CP_nodes_y(c_j,num_nodes_CP_set_y(c_j))= int(Temp_DATA(c_i,2))
            endif
        enddo
    enddo
    deallocate(Temp_DATA)
endif     

!     ----------------------------------------------------------------------------------------
!     Check if the initial velocity file in the x direction exists; if it does, read it (used for
!     dynamic analysis)
!     ----------------------------------------------------------------------------------------
print *, "    Trying to read ivex files..." 
temp_name = trim(trim(Full_Pathname)//'.ivex')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Ivex = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Ive_x)) deallocate(Ive_x)
    allocate(Ive_x(Num_Ivex,2))
    allocate(Temp_DATA(Num_Ivex,2))
    Call Tool_Read_File(temp_name,"ivex",Num_Ivex,2,Temp_DATA, Flag_Blank)
    Ive_x  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif  

!     ----------------------------------------------------------------------------------------
!     Check if the initial velocity file in the y-direction exists; if it does, read it (used for
!     dynamic analysis)
!     ----------------------------------------------------------------------------------------
print *, "    Trying to read ivey files..." 
temp_name = trim(trim(Full_Pathname)//'.ivey')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Ivey = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Ive_y)) deallocate(Ive_y)
    allocate(Ive_y(Num_Ivey,2))
    allocate(Temp_DATA(Num_Ivey,2))
    Call Tool_Read_File(temp_name,"ivey",Num_Ivey,2,Temp_DATA, Flag_Blank)
    Ive_y  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     ----------------------------------------------------------------------------------------------
!     Check whether the initial acceleration file in the x-direction exists. If it exists, read it (for
!     dynamic analysis).
!     ----------------------------------------------------------------------------------------------
print *, "    Trying to read iacx files..." 
temp_name = trim(trim(Full_Pathname)//'.iacx')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Iacx = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Iac_x)) deallocate(Iac_x)
    allocate(Iac_x(Num_Iacx,2))
    allocate(Temp_DATA(Num_Iacx,2))
    Call Tool_Read_File(temp_name,"iacx",Num_Iacx,2,Temp_DATA, Flag_Blank)
    Iac_x  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif  

!     ----------------------------------------------------------------------------------------------
!     Check whether the initial acceleration file in the Y direction exists. If it exists, read it (for
!     dynamic analysis).
!     ----------------------------------------------------------------------------------------------
print *, "    Trying to read iacy files..." 
temp_name = trim(trim(Full_Pathname)//'.iacy')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Iacy = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Iac_y)) deallocate(Iac_y)
    allocate(Iac_y(Num_Iacy,2))
    allocate(Temp_DATA(Num_Ivey,2))
    Call Tool_Read_File(temp_name,"iacy",Num_Iacy,2,Temp_DATA, Flag_Blank)
    Iac_y  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     --------------------------------------------------------------------------
!     Check the boundary conditions of the problem domain, added on August 11, 2016
!     --------------------------------------------------------------------------
print *, "    Trying to read fbvl file..." 
temp_name = trim(trim(Full_Pathname)//'.fbvl')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Fd_Bou_vl = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Fd_Bou_vl)) deallocate(Fd_Bou_vl)
    allocate(Fd_Bou_vl(Num_Fd_Bou_vl,2))
    allocate(Temp_DATA(Num_Fd_Bou_vl,2))
    Call Tool_Read_File(temp_name,"fbvl",Num_Fd_Bou_vl,2, Temp_DATA,Flag_Blank)
    Fd_Bou_vl  = Temp_DATA(:,:)
    deallocate(Temp_DATA)

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
    if (ALLOCATED(Fd_Bou_fixed)) deallocate(Fd_Bou_fixed)
    allocate(Fd_Bou_fixed(Num_Fd_Bou_fixed))
    if (ALLOCATED(Fd_Bou_vl_nz)) deallocate(Fd_Bou_vl_nz)
    allocate(Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,2))
    Num_Fd_Bou_fixed = 0
    Num_Fd_Bou_vl_nz = 0
    do i_Check =1,Num_Fd_Bou_vl
        if (abs(Fd_Bou_vl(i_Check,2))<=Tol_30) then
            Num_Fd_Bou_fixed = Num_Fd_Bou_fixed +1
            Fd_Bou_fixed(Num_Fd_Bou_fixed) = int(Fd_Bou_vl(i_Check,1))
        else
            Num_Fd_Bou_vl_nz  = Num_Fd_Bou_vl_nz + 1 
            Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,1)= Fd_Bou_vl(i_Check,1) 
            Fd_Bou_vl_nz(Num_Fd_Bou_vl_nz,2)= Fd_Bou_vl(i_Check,2)
        endif
    enddo
endif 

!     -------------------------------------------------------------------
!     Check the initial value file of the field problem, added on 2016-11-08
!     -------------------------------------------------------------------
print *, "    Trying to read .fbiv files..." 
temp_name = trim(trim(Full_Pathname)//'.fbiv')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Fd_Ini_vl = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Fd_Ini_vl)) deallocate(Fd_Ini_vl)
    allocate(Fd_Ini_vl(Num_Fd_Ini_vl,2))
    allocate(Temp_DATA(Num_Fd_Ini_vl,2))
    Call Tool_Read_File(temp_name,"fbvl",Num_Fd_Ini_vl,2, Temp_DATA,Flag_Blank)
    Fd_Ini_vl  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     ------------------------------------------------------------------------------
!      Check the flux through the boundary normal of the field, added on August 12, 2016
!     ------------------------------------------------------------------------------
print *, "    Trying to read fbqn files..." 
temp_name = trim(trim(Full_Pathname)//'.fbqn')
inquire(file=temp_name, exist=alive)  
if(alive.eqv..false.)then
else
    Num_Fd_Bou_qn = Tool_Count_Lines(temp_name) 
    if (ALLOCATED(Fd_Bou_qn)) deallocate(Fd_Bou_qn)
    allocate( Fd_Bou_qn(Num_Fd_Bou_qn,2))
    allocate( Temp_DATA(Num_Fd_Bou_qn,2))
    Call Tool_Read_File(temp_name,"fbqn",Num_Fd_Bou_qn,2, Temp_DATA,Flag_Blank)
    Fd_Bou_qn  = Temp_DATA(:,:)
    deallocate(Temp_DATA)
endif 

!     -------------------------------------------------------------------------------------------
!     Read the transient wellbore pressure variation curve data (only for analysis of No. 16 and 17)
!     -------------------------------------------------------------------------------------------
if(Key_Analysis_Type ==16 .or. Key_Analysis_Type ==17)then
    if(Key_Gas_Production==1 .and. Key_Changing_BHP==1) then
        print *, "    Warning :: Key_Gas_Production =1" 
        print *, "    Warning :: Key_Changing_BHP =1" 
        print *, "    Trying to read bhpc file..." 
        temp_name = trim(trim(Full_Pathname)//'.bhpc')
        inquire(file=temp_name, exist=alive)  
        if(alive.eqv..false.)then
            print *, "    Warning :: Can not find bhpc files!!!" 
            call Warning_Message('S',Keywords_Blank)
        else
            Num_BHP_curve_point = Tool_Count_Lines(temp_name) 
            if (ALLOCATED(BHP_Curve)) deallocate(BHP_Curve)
            allocate( BHP_Curve(Num_BHP_curve_point,2))
            allocate( Temp_DATA(Num_BHP_curve_point,2))
            Call Tool_Read_File(temp_name,"bhpc",Num_BHP_curve_point, 2,Temp_DATA,Flag_Blank)
            BHP_Curve  = Temp_DATA(:,:)
            deallocate(Temp_DATA)
        endif
    endif
endif
!c     --------------------------------------------------      
! c     Check the flux in the x-direction at the problem boundary, added on August 11, 2016
!c     --------------------------------------------------
!     temp_name = trim(trim(Full_Pathname)//'.fbqx')
!     inquire(file=temp_name, exist=alive)  
!     if(alive.eqv..false.)then
!     else
!         Num_Fd_Bou_qx = Tool_Count_Lines(temp_name) 
!         allocate( Fd_Bou_qx(Num_Fd_Bou_qx,2))
!         allocate( Temp_DATA(Num_Fd_Bou_qx,2))
!         Call Tool_Read_File(temp_name,"fbqx",Num_Fd_Bou_qx,2,
!    &                        Temp_DATA,Flag_Blank)
!         Fd_Bou_qx  = Temp_DATA(:,:)
!         deallocate(Temp_DATA)
!     endif 

!c     --------------------------------------------------      
! c     Check the flux in the x-direction at the problem boundary, added on August 11, 2016
!c     --------------------------------------------------
!     temp_name = trim(trim(Full_Pathname)//'.fbqy')
!     inquire(file=temp_name, exist=alive)  
!     if(alive.eqv..false.)then
!     else
!         Num_Fd_Bou_qx = Tool_Count_Lines(temp_name) 
!         allocate( Fd_Bou_qy(Num_Fd_Bou_qy,2))
!         allocate( Temp_DATA(Num_Fd_Bou_qy,2))
!         Call Tool_Read_File(temp_name,"fbqx",Num_Fd_Bou_qy,2,
!    &                        Temp_DATA,Flag_Blank)
!         Fd_Bou_qy  = Temp_DATA(:,:)
!         deallocate(Temp_DATA)
!     endif 

!     --------------------------------------------------------------------
!     Determine the outer boundary of all quadrilateral elements of the model
!     --------------------------------------------------------------------
! The four sides of the element are stored in the temporary variable tem1
print *, "    Finding model boundary..." 
allocate(tem1(4*Num_Elem,2))
allocate(tem2(4*Num_Elem))
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
call Matrix_Count_Row_Int_m_x_2(4*Num_Elem, &
tem1,tem2,all_num_Outline,0)
! Extract edges that appear only once, including both outer and inner edges
allocate(All_Outline(all_num_Outline,2))
allocate(Temp_Outline(all_num_Outline,2))
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
call Tool_Sort_by_End_to_End(all_num_Outline,all_num_Outline, All_Outline,Temp_Outline, Out_num_Outline)
if (ALLOCATED(Outline)) deallocate(Outline)
allocate(Outline(Out_num_Outline,2))
do i=1,Out_num_Outline
    Outline(i,1:2) = Temp_Outline(i,1:2) 
end do
deallocate(All_Outline)
deallocate(Temp_Outline)
!     ------------------------------------------------------------------------------------------
!     Coordinates of the four nodes of the storage unit, calculation of the element's average area,
!     element side length, etc.
!     ------------------------------------------------------------------------------------------
print *, "    Get node and element data..." 
if (ALLOCATED(G_NN)) deallocate(G_NN)
allocate(G_NN(4,Num_Elem))     
if (ALLOCATED(G_X_NODES)) deallocate(G_X_NODES)
allocate(G_X_NODES(4,Num_Elem))
if (ALLOCATED(G_Y_NODES)) deallocate(G_Y_NODES)
allocate(G_Y_NODES(4,Num_Elem))
if (ALLOCATED(Elem_Area)) deallocate(Elem_Area)
allocate(Elem_Area(Num_Elem))

!2026-04-25.
if (ALLOCATED(Elem_Max_L)) deallocate(Elem_Max_L)
allocate(Elem_Max_L(Num_Elem))
if (ALLOCATED(Elem_Min_L)) deallocate(Elem_Min_L)
allocate(Elem_Min_L(Num_Elem))

if (ALLOCATED(Elem_Centroid)) deallocate(Elem_Centroid)
allocate(Elem_Centroid(Num_Elem,2))
if (ALLOCATED(EleGaus_yes_FEM_asemd)) then
    deallocate(EleGaus_yes_FEM_asemd)
endif

!2026-05-28.
if(Key_TipEnrich==5) then
    if(Num_Gauss_Points<900) then
        Num_Gauss_Points=900
    endif
end if

allocate(EleGaus_yes_FEM_asemd(Num_Elem,Num_Gauss_Points))
if (ALLOCATED(x_max_Elements)) deallocate(x_max_Elements)
allocate(x_max_Elements(Num_Elem)) 
if (ALLOCATED(x_min_Elements)) deallocate(x_min_Elements)
allocate(x_min_Elements(Num_Elem)) 
if (ALLOCATED(y_max_Elements)) deallocate(y_max_Elements)
allocate(y_max_Elements(Num_Elem)) 
if (ALLOCATED(y_min_Elements)) deallocate(y_min_Elements)
allocate(y_min_Elements(Num_Elem)) 
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
    c_L1= sqrt((Coor(N1,1)-Coor(N2,1))**2+ (Coor(N1,2)-Coor(N2,2))**2)
    c_L2= sqrt((Coor(N2,1)-Coor(N3,1))**2+ (Coor(N2,2)-Coor(N3,2))**2)
    c_L3= sqrt((Coor(N3,1)-Coor(N4,1))**2+ (Coor(N3,2)-Coor(N4,2))**2)
    c_L4= sqrt((Coor(N4,1)-Coor(N1,1))**2+ (Coor(N4,2)-Coor(N1,2))**2)
    c_L_min = min(c_L1,c_L2,c_L3,c_L4)
    c_L_max = max(c_L1,c_L2,c_L3,c_L4)
    Elem_Min_L(i) =  c_L_min
    Elem_Max_L(i) =  c_L_max


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

!     -----------------------------------------------------
!     The elements around each node and the number of elements
!     -----------------------------------------------------
if (ALLOCATED(Node_Elements)) then
    deallocate(Node_Elements)
endif
allocate(Node_Elements(Num_Node,8))
Node_Elements = 0
if (ALLOCATED(num_Node_Elements)) deallocate(num_Node_Elements)      
allocate(num_Node_Elements(Num_Node))
!$OMP PARALLEL do DEFAULT(SHARED) private(i,j)          
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
!     -------------------------------------------------------------
!     The elements surrounding each element and the number of elements
!     -------------------------------------------------------------
if (ALLOCATED(Ele_Elements)) deallocate(Ele_Elements)
allocate(Ele_Elements(Num_Elem,25))
if (ALLOCATED(num_Ele_Eles)) deallocate(num_Ele_Eles)      
allocate(num_Ele_Eles(Num_Elem))
Ele_Elements(1:Num_Elem,1:25) = 0
num_Ele_Eles(1:Num_Elem) = 0
!Loop through each element.
do i_E =1,Num_Elem
    c_all_num_Ele      = 0
    c_all_Elem(1:32)   = 0
    !Loop through each node of the element.
    do j=1,4
        c_Node = Elem_Node(i_E,j) 
        c_all_Elem(c_all_num_Ele + 1: c_all_num_Ele + num_Node_Elements(c_Node)) = &
        Node_Elements(c_Node,1:num_Node_Elements(c_Node))
        c_all_num_Ele = c_all_num_Ele + num_Node_Elements(c_Node)
    end do
    ! Delete duplicate unit numbers
    call Vector_Unique_Int(32,c_all_num_Ele,c_all_Elem, Uniqued_all_Elem,num_Surrou)
    Ele_Elements(i_E,1:num_Surrou)=Uniqued_all_Elem(1:num_Surrou)
    num_Ele_Eles(i_E) = num_Surrou
end do
!     ------------------------------------------------------------------------------------------
!     Determine dimensions for relevant public variables (mainly related to identifying enhancement
!     nodes)
!     ------------------------------------------------------------------------------------------
if (ALLOCATED(Elem_Type)) deallocate(Elem_Type)
allocate(Elem_Type(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Enriched_Node_Type)) deallocate(Enriched_Node_Type)
allocate(Enriched_Node_Type(Num_Node,Max_Num_Cr))
if (ALLOCATED(Node_Jun_elem)) deallocate(Node_Jun_elem)
allocate(Node_Jun_elem(Num_Node,Max_Num_Cr))
if (ALLOCATED(Jun_Ele_Negative_Cr_Num)) then
    deallocate(Jun_Ele_Negative_Cr_Num)
endif
allocate(Jun_Ele_Negative_Cr_Num(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Node_Jun_Hole)) deallocate(Node_Jun_Hole)
allocate(Node_Jun_Hole(Num_Node,Max_Num_Cr))
if (ALLOCATED(Ele_Jun_Hole)) deallocate(Ele_Jun_Hole)
allocate(Ele_Jun_Hole(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Coors_Element_Crack)) deallocate(Coors_Element_Crack)
allocate(Coors_Element_Crack(Num_Elem,Max_Num_Cr,4))

!2026-05-23.
if(Key_Analysis_Type==6) then
    if (ALLOCATED(Coors_Element_Crack_OLD)) DEALLOCATE(Coors_Element_Crack_OLD)
    allocate(Coors_Element_Crack_OLD(Num_Elem,Max_Num_Cr,4))
endif

if (ALLOCATED(Coors_Tip)) deallocate(Coors_Tip)
allocate(Coors_Tip(Num_Elem,2))
if (ALLOCATED(Coors_Vertex)) deallocate(Coors_Vertex)
allocate(Coors_Vertex(Num_Elem,2))
if (ALLOCATED(Coors_Junction)) deallocate(Coors_Junction)
allocate(Coors_Junction(Num_Elem,Max_Num_Cr,4))
if (ALLOCATED(x_cr_tip_nodes)) deallocate(x_cr_tip_nodes)
allocate(x_cr_tip_nodes(Max_Num_Cr,Num_Node))
if (ALLOCATED(y_cr_tip_nodes)) deallocate(y_cr_tip_nodes)
allocate(y_cr_tip_nodes(Max_Num_Cr,Num_Node))
if (ALLOCATED(Ele_Num_Tip_Enriched_Node)) then 
    deallocate(Ele_Num_Tip_Enriched_Node)
endif
allocate(Ele_Num_Tip_Enriched_Node(Max_Num_Cr,Num_Node))
if (ALLOCATED(c_POS)) deallocate(c_POS)
allocate(c_POS(Num_Node,Max_Num_Cr))
c_POS = 0
if (ALLOCATED(TipEle_Adjacent_Ele)) deallocate(TipEle_Adjacent_Ele)
allocate(TipEle_Adjacent_Ele(Num_Elem,Max_Num_Cr))
if (ALLOCATED(num_GP_Elem)) deallocate(num_GP_Elem)
allocate(num_GP_Elem(Num_Elem))
if (ALLOCATED(Ele_GP_Start_Num)) deallocate(Ele_GP_Start_Num)
allocate(Ele_GP_Start_Num(Num_Elem))
! allocate(Cracks_CalP_Elem_CalP(Max_Num_Cr, Num_Elem)) ! Starting index and number of calculation
! points contained in each element (relative to each crack, with local numbering for each crack)
! ---------Related to Holes-----------
if (ALLOCATED(Elem_Type_Hl)) deallocate(Elem_Type_Hl)
allocate(Elem_Type_Hl(Num_Elem,Max_Num_Hl))
if (ALLOCATED(Enriched_Node_Type_Hl)) then
    deallocate(Enriched_Node_Type_Hl)
endif
allocate(Enriched_Node_Type_Hl(Num_Node,Max_Num_Hl))
if (ALLOCATED(c_POS_Hl)) deallocate(c_POS_Hl)
allocate(c_POS_Hl(Num_Node,Max_Num_Hl))
! ---------Related Mixtures-----------
if (ALLOCATED(Elem_Type_Incl)) deallocate(Elem_Type_Incl)
allocate(Elem_Type_Incl(Num_Elem,Max_Num_Incl))
if (ALLOCATED(Enriched_Node_Type_Incl)) then
    deallocate(Enriched_Node_Type_Incl)
endif
allocate(Enriched_Node_Type_Incl(Num_Node,Max_Num_Incl))
if (ALLOCATED(c_POS_Incl)) deallocate(c_POS_Incl)
allocate(c_POS_Incl(Num_Node,Max_Num_Incl))
! ---------Cross-shaped Intersection Cracks Related-----
if (ALLOCATED(Cross_Point_Cr_num)) deallocate(Cross_Point_Cr_num)
allocate(Cross_Point_Cr_num(Max_Num_Cross,2))
if (ALLOCATED(Cross_Point_Ele_num)) deallocate(Cross_Point_Ele_num)
allocate(Cross_Point_Ele_num(Max_Num_Cross))
! allocate(Cross_Point_Cr_Seg_num(Max_Num_Cross,2))               ! Used for cross-shaped cracks
if (ALLOCATED(Cross_Point_RABCD)) deallocate(Cross_Point_RABCD)
allocate(Cross_Point_RABCD(Max_Num_Cross,10,2))
if (ALLOCATED(Node_Cross_elem)) deallocate(Node_Cross_elem)
allocate(Node_Cross_elem(Num_Node,Max_Num_Cross))
if (ALLOCATED(Elem_Type_Cross)) deallocate(Elem_Type_Cross)
allocate(Elem_Type_Cross(Num_Elem,Max_Num_Cross))
if (ALLOCATED(Enriched_Node_Type_Cross)) then
    deallocate(Enriched_Node_Type_Cross)
endif
allocate(Enriched_Node_Type_Cross(Num_Node,Max_Num_Cross))
if (ALLOCATED(c_POS_Cross)) deallocate(c_POS_Cross)
allocate(c_POS_Cross(Num_Node,Max_Num_Cross))

! HC_Point_Cr_num(num_HC,1:2) !The crack number and crack tip number corresponding to each HC
! intersection point
! HC_Point_Hole_num(num_HC)               !The hole number corresponding to each HC intersection point
! HC_Point_Ele_num(num_HC)            !The element number corresponding to each HC intersection point
! Elem_Type_HC(Num_Elem, num_HC); !Element type: 1--Heaviside enhanced elements additionally added
! at crack and hole intersections
! Enriched_Node_Type_HC(Num_Node, num_HC); !Type of enriched nodes: 1--Heaviside enriched nodes
! additionally added where cracks and holes intersect


!2026-06-22.
if (ALLOCATED(Arc_Crack_Coor)) deallocate(Arc_Crack_Coor)
allocate(Arc_Crack_Coor(Max_Num_Cr,Max_Num_Cr_P-1,11)) 
if (ALLOCATED(Na_Crack_Coor)) deallocate(Na_Crack_Coor)
allocate(Na_Crack_Coor(Max_Num_Cr,Max_Num_Cr_P,2) ) 
if (ALLOCATED(Edge_Disposed_Crack)) deallocate(Edge_Disposed_Crack)
allocate(Edge_Disposed_Crack(Max_Num_Cr,Max_Num_Cr_P,2)) 
if (ALLOCATED(Cracks_Coh_Ele_Type)) deallocate(Cracks_Coh_Ele_Type)
allocate(Cracks_Coh_Ele_Type(Max_Num_Cr,Max_Num_Cr_CalP-1)) 
if (ALLOCATED(Cracks_LocalNumCalP_W)) deallocate(Cracks_LocalNumCalP_W)
allocate(Cracks_LocalNumCalP_W(Max_Num_Cr*Max_Num_Cr_CalP,2))
if (ALLOCATED(Cracks_HF_Ele_L)) deallocate(Cracks_HF_Ele_L)
allocate(Cracks_HF_Ele_L(Max_Num_Cr,Max_Num_Cr_CalP-1))    
if (ALLOCATED(Cracks_GloNumCalP_W)) deallocate(Cracks_GloNumCalP_W)
allocate(Cracks_GloNumCalP_W(Max_Num_Cr,Max_Num_Cr_CalP)) 
if (ALLOCATED(Cracks_LocalNumCalP)) deallocate(Cracks_LocalNumCalP)
allocate(Cracks_LocalNumCalP(Max_Num_Cr*Max_Num_Cr_CalP,2)) 
if (ALLOCATED(Cracks_MidJunCalpNum)) deallocate(Cracks_MidJunCalpNum)
allocate(Cracks_MidJunCalpNum(Max_Num_Cr,Max_Num_Cr_CalP)) 
if (ALLOCATED(Cracks_CalP_Type)) deallocate(Cracks_CalP_Type)
allocate(Cracks_CalP_Type(Max_Num_Cr,Max_Num_Cr_CalP,2))
if (ALLOCATED(Cracks_GloNumCalP)) deallocate(Cracks_GloNumCalP)
allocate(Cracks_GloNumCalP(Max_Num_Cr,Max_Num_Cr_CalP))
! Arc Crack Coordinate Initialization
Arc_Crack_Coor(1:Max_Num_Cr,1:Max_Num_Cr_P-1,1:11)  = ZR
Cracks_Coh_Ele_Type(1:Max_Num_Cr,1:Max_Num_Cr_CalP-1) =0


! ---------Contact Issues Related-----------
if (ALLOCATED(Elem_Conta_Sta)) deallocate(Elem_Conta_Sta)
allocate(Elem_Conta_Sta(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Elem_Coh_Sta)) deallocate(Elem_Coh_Sta)
allocate(Elem_Coh_Sta(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Elem_Conta_Sta_Last)) then 
    deallocate(Elem_Conta_Sta_Last)
endif
allocate(Elem_Conta_Sta_Last(Num_Elem,Max_Num_Cr))
if (ALLOCATED(Ele_NumCalP)) deallocate(Ele_NumCalP)
allocate(Ele_NumCalP(Num_Elem))
if (ALLOCATED(Ele_CalPNum)) deallocate(Ele_CalPNum)
allocate(Ele_CalPNum(Num_Elem,Max_Num_Ele_CalP))
if (ALLOCATED(Elem_Proppant_Coor)) then
    deallocate(Elem_Proppant_Coor)
endif
allocate(Elem_Proppant_Coor(Num_Elem,5))
! Among them: 1, indicates whether the element contains proppant
! 2. Proppant Diameter
! 3 and 4, x-coordinate and y-coordinate
! 5. Corresponding crack segment dip angle
! ----------Related to Field Issues-----------
if (ALLOCATED(Fd_c_POS)) deallocate(Fd_c_POS)
allocate(Fd_c_POS(Num_Node,Max_Num_Cr))
if (ALLOCATED(Fd_c_POS_Hl)) deallocate(Fd_c_POS_Hl)
allocate(Fd_c_POS_Hl(Num_Node,Max_Num_Hl))
if (ALLOCATED(Fd_c_POS_Incl)) deallocate(Fd_c_POS_Incl)
allocate(Fd_c_POS_Incl(Num_Node,Max_Num_Incl))
if (ALLOCATED(Fd_c_POS_Cross)) deallocate(Fd_c_POS_Cross)
allocate(Fd_c_POS_Cross(Num_Node,Max_Num_Cross))
if (ALLOCATED(Fd_EleGaus_yes_FEM_asemd)) then
    deallocate(Fd_EleGaus_yes_FEM_asemd)
endif
allocate(Fd_EleGaus_yes_FEM_asemd(Num_Elem,Num_Gauss_Points))
! ----------Element Destruction-------
if (ALLOCATED(Elem_Break)) deallocate(Elem_Break)
allocate(Elem_Break(Num_Elem))
if (ALLOCATED(Elem_Ave_Gauss_Stress)) then
    deallocate(Elem_Ave_Gauss_Stress)
endif
allocate(Elem_Ave_Gauss_Stress(Num_Elem,3))
if (ALLOCATED(Elem_Ave_Gauss_S1)) deallocate(Elem_Ave_Gauss_S1)
allocate(Elem_Ave_Gauss_S1(Num_Elem))
Elem_Ave_Gauss_Stress(1:Num_Elem,1:3) = ZR 
Elem_Ave_Gauss_S1(1:Num_Elem) = ZR
Elem_Break(1:Num_Elem) = .False.

!----------------------------------------------------
! 2D crack related (changed to Allocatable variable)
!        2022-09-02
!        IMPROV2022090201.
!----------------------------------------------------
if (ALLOCATED(Cracks_CalP_Num)) deallocate(Cracks_CalP_Num)
allocate(Cracks_CalP_Num(Max_Num_Cr))
if (ALLOCATED(Cracks_CalP_Coors)) deallocate(Cracks_CalP_Coors)
allocate(Cracks_CalP_Coors(Max_Num_Cr,Max_Num_Cr_CalP,2))
if (ALLOCATED(Cracks_CalP_Orient)) deallocate(Cracks_CalP_Orient)
allocate(Cracks_CalP_Orient(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Seg)) deallocate(Cracks_CalP_Seg)
allocate(Cracks_CalP_Seg(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Elem)) deallocate(Cracks_CalP_Elem)
allocate(Cracks_CalP_Elem(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Aper)) deallocate(Cracks_CalP_Aper)
allocate(Cracks_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Pres)) deallocate(Cracks_CalP_Pres)
allocate(Cracks_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Tractions)) then
    deallocate(Cracks_CalP_Tractions)
endif
allocate(Cracks_CalP_Tractions(Max_Num_Cr,Max_Num_Cr_CalP,2))
if (ALLOCATED(Cracks_CalP_Pgra)) deallocate(Cracks_CalP_Pgra)
allocate(Cracks_CalP_Pgra(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Velo)) deallocate(Cracks_CalP_Velo)
allocate(Cracks_CalP_Velo(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Quan)) deallocate(Cracks_CalP_Quan)
allocate(Cracks_CalP_Quan(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Conc)) deallocate(Cracks_CalP_Conc)
allocate(Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Remo_Strs)) then
    deallocate(Cracks_CalP_Remo_Strs)
endif
allocate(Cracks_CalP_Remo_Strs(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Contact_Strs)) then
    deallocate(Cracks_CalP_Contact_Strs)
endif
allocate(Cracks_CalP_Contact_Strs(Max_Num_Cr,Max_Num_Cr_CalP,2))
if (ALLOCATED(Cracks_CalP_Contact_Force_x)) then
    deallocate(Cracks_CalP_Contact_Force_x)
endif
allocate(Cracks_CalP_Contact_Force_x(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Contact_Force_y)) then
    deallocate(Cracks_CalP_Contact_Force_y)
endif
allocate(Cracks_CalP_Contact_Force_y(Max_Num_Cr,Max_Num_Cr_CalP))
! -----Added on 2016-08-25--Related to the calculation of support crack width and flow
! capacity-------
if (ALLOCATED(Cracks_CalP_wpnp)) deallocate(Cracks_CalP_wpnp)
allocate(Cracks_CalP_wpnp(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_wpor)) deallocate(Cracks_CalP_wpor)
allocate(Cracks_CalP_wpor(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_wdeform)) deallocate(Cracks_CalP_wdeform)
allocate(Cracks_CalP_wdeform(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_Conductivity)) then
    deallocate(Cracks_CalP_Conductivity)
endif
allocate(Cracks_CalP_Conductivity(Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(Cracks_CalP_kf)) deallocate(Cracks_CalP_kf)
allocate(Cracks_CalP_kf(Max_Num_Cr,Max_Num_Cr_CalP))
! real(kind=FT), allocatable :: Cracks_CalP_Elem_CalP(:,:) ! Starting index and number of
! calculation points contained in each element (relative to each crack, local numbering for each
! crack)
! -----------The following variables are related to staged fracturing----------
if (ALLOCATED(MS_Cracks_CalP_Num)) deallocate(MS_Cracks_CalP_Num)
allocate(MS_Cracks_CalP_Num(20,Max_Num_Cr))
if (ALLOCATED(MS_Cracks_CalP_Aper)) deallocate(MS_Cracks_CalP_Aper)
allocate(MS_Cracks_CalP_Aper(20,Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(MS_Cracks_CalP_Conc)) deallocate(MS_Cracks_CalP_Conc)
allocate(MS_Cracks_CalP_Conc(10,Max_Num_Cr,Max_Num_Cr_CalP))
if (ALLOCATED(MS_CalP_Propped_Aper)) then
    deallocate(MS_CalP_Propped_Aper)
endif
allocate(MS_CalP_Propped_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
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


! Get initial length of cracks. 2026-06-29. 
Initial_Cracks_Length(Max_Num_Cr) = ZR
do i_C = 1, num_Crack
    n_pts = Each_Cr_Poi_Num(i_C)
    if (n_pts < 2) cycle
    ! Compute total crack length as polyline segment sum
    do i_pt = 1, n_pts - 1
        Initial_Cracks_Length(i_C) = Initial_Cracks_Length(i_C) + &
        sqrt((Crack_Coor(i_C, i_pt+1, 1) - Crack_Coor(i_C, i_pt, 1))**2 + &
             (Crack_Coor(i_C, i_pt+1, 2) - Crack_Coor(i_C, i_pt, 2))**2)
    end do
enddo



!     ------------------------------------------------------------------------------------------
!     Hydraulic fracturing related: Initialize the hydraulic fracturing status of all cracks (i.e.,
!     whether they are driven by water)
!     ------------------------------------------------------------------------------------------
!     Cracks_HF_State(1:Max_Num_Cr)  = 0
!     Cracks_HF_State(Inject_Crack_Num) = 1  ! Crack containing injection point

return
END subroutine Read_Geo
