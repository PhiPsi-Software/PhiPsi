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
 
SUBROUTINE D3_Check_and_Adjust_Cracks
! Adjust the coordinates of the fracture surface according to the nodal symbol distances.
! Purpose: To prevent the stiffness matrix from being singular.
!Firstly written on 2022-08-01. NEWFTU2022080102.



!----------------------------------
! Read the public variable module.
!----------------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_XFEM_Elements
use omp_lib
!
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane

!----------------------------
! Variable type declaration.
!----------------------------
implicit none
integer i_C,i_Node,i_E,i_Try
integer c_Crack_Type
real(kind=FT) c_Singularity,Singularity_Factor
logical Yes_Singularity
real(kind=FT) delta_L,tem_L
real(kind=FT) Tool_Function_Odd_Positive_Even_Negative
integer num_Try_1,num_Try_2
real(kind=FT) c_Crack_Center(3)
integer c_OUT_Elem,c_NN(8)
integer c_Node
integer c_Num_Eles
integer c_Logical_Positive_Num,c_Logical_Negative_Num
logical c_Logical_Positive,c_Logical_Negative
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT),ALLOCATABLE::kesi_Enr(:),yita_Enr(:),zeta_Enr(:),weight_Enr(:)
integer Tem_Num_Gauss_3D,Num_Gauss_Cubes
integer i_G
real(kind=FT) N(3,24),Global_coor_Gauss(3)
real(kind=FT) c_PER_Node_to_FS(3)
real(kind=FT) Check_Ball_R
real(kind=FT) c_Signed_Dis_v2
real(kind=FT) ,ALLOCATABLE::c_Distance(:)
logical Yes_Found_Min_Signed_Dis
logical c_Yes_Node_PER_in_FS
real(kind=FT) c_n_Vector(3)
real(kind=FT) c_Sign
integer Ele_Num_Cache

!-----------------------------
! Formatted output statement.
!-----------------------------
101 FORMAT(5X,'Attention: possible Heaviside singularity found in crack',I5)
102 FORMAT(5X,'Attention: center adjustment succesfully done after ',I5,' tries')
103 FORMAT(5X,'Attention: center adjustment failed after ',I5,' tries')



!------------------------
! Basic data parameters.
!------------------------
Singularity_Factor = 20.0D0
num_Try_1          = 50
num_Try_2          = 50

!---------------------------------
! Get the crack coordinate range.
!---------------------------------
call D3_Get_Cracks_Coor_Ranges

!---------------------------------------------------------------
! Point data, used only for Key_Check_and_Adjust_Cracks_3D = 5.
! NEWFTU2022092901.
!---------------------------------------------------------------
if(Key_Check_and_Adjust_Cracks_3D ==5) then
    ! Depending on the points system:
    if(allocated(kesi_Enr)) deallocate(kesi_Enr)
    if(allocated(yita_Enr)) deallocate(yita_Enr)
    if(allocated(zeta_Enr)) deallocate(zeta_Enr)
    if(allocated(weight_Enr)) deallocate(weight_Enr)
    select case(Key_Integral_Sol)
    ! Fixed points
    case(2)
        allocate(kesi_Enr(Num_Gau_Points_3D))
        allocate(yita_Enr(Num_Gau_Points_3D))
        allocate(zeta_Enr(Num_Gau_Points_3D))
        allocate(weight_Enr(Num_Gau_Points_3D))
        call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
        Tem_Num_Gauss_3D = Num_Gau_Points_3D
    ! Fixed block integration. NEWFTU2022072901.
    case(3)
        Num_Gauss_Cubes = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
        allocate(kesi_Enr(Num_Gauss_Cubes))
        allocate(yita_Enr(Num_Gauss_Cubes))
        allocate(zeta_Enr(Num_Gauss_Cubes))
        allocate(weight_Enr(Num_Gauss_Cubes))
        call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,&
                                              kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
        Tem_Num_Gauss_3D = Num_Gauss_Cubes
    ! Integration by Parts (Temporarily Abandoned). 2022-07-27.
    case(4)
    ! If it is a block integration scheme, skip this step.
    print *, ' To be done!'
    end select
endif

!----------------------
! Each fracture cycle.
!----------------------
do i_C=1,num_Crack

    !*************************************************
    ! If the crack has already been inspected before,
    ! proceed to the next loop.
    !*************************************************
    if(Cracks_Checked_and_Adjusted(i_C)==1)then
        cycle
    endif

    !******************************
    ! Mark this crack as detected.
    !******************************
    Cracks_Checked_and_Adjusted(i_C)=1

    !*********************************
    ! Obtain the type of crack shape.
    !*********************************
    if(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>Tol_20)then
        c_Crack_Type = 1
    elseif(sum(abs(Crack3D_Cir_Coor(i_C,1:7)))>Tol_20)then
        c_Crack_Type = 2
    elseif(sum(abs(Crack3D_Ellip_Coor(i_C,1:8)))>Tol_20)then
        c_Crack_Type = 3
    endif



    select case(Key_Check_and_Adjust_Cracks_3D)

    !**************************************************************************
    !                                                                        *
    !                                                                        *
    !                                                                        *
    ! Plan 1: Continuously try fine-tuning, divided into (1)centroid .       *
    !coordinate fine-tuning and (2)normal direction fine-tuning (To be done) *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !**************************************************************************
    case(1)
        !------------------------------------------------------------
        !-                                                         -
        !- Step 1: Try adjusting the crack coordinates several     -
        !- times (along the normal direction of the crack surface, -
        !- Adjust back and forth around the centroid.              -
        !-                                                         -
        !------------------------------------------------------------
        DO i_Try =1,num_Try_1
            !///////////////////////////////////////////////////////////////////////////////
            ! Determine whether the Heaviside enrichment function distance for the crack is
            ! singular.
            !///////////////////////////////////////////////////////////////////////////////
            call  D3_Check_Crack_Heaviside_Singularity(i_C,Singularity_Factor,Yes_Singularity)

            !////////////////////////////////////////////////////////////////////////////////////////
            ! If there are anomalies, adjust the crack coordinates and re-mesh.
            ! Note: Adjust around the original crack face centroid coordinates to the left or right.
            ! Ref: My PhiPsi Development Notebook V1-P60-(1).
            !////////////////////////////////////////////////////////////////////////////////////////
            if(Yes_Singularity)then
                if(i_Try ==1) write(*,101)  i_C
                !======================================================
                ! Adjust the initial crack along the normal direction.
                !======================================================
                ! Tool_Function_Odd_Positive_Even_Negative: Odd numbers are 1, even numbers are -1
                tem_L = 0.25D0/dble(num_Try_1)
                if(i_Try==1) then
                    delta_L = Ave_Elem_L*(tem_L+1.0D-3)*DBLE(i_Try)*Tool_Function_Odd_Positive_Even_Negative(i_Try)
                else
                    delta_L = TWO*Ave_Elem_L*(tem_L+1.0D-3)*DBLE(i_Try)*Tool_Function_Odd_Positive_Even_Negative(i_Try)
                endif
                ! Make sure to return to the original position in the final step
                if(i_Try==num_Try_1) then
                    delta_L = Ave_Elem_L*(tem_L+1.0D-3)*DBLE(i_Try-1)*Tool_Function_Odd_Positive_Even_Negative(i_Try)
                endif
                ! Adjust along the normal direction.
                call D3_Offset_Initial_Crack_Along_Normal(i_C,c_Crack_Type,delta_L)


                ! Set the mark of the crack Cracks_Initial_Meshed to 0.
                Cracks_Initial_Meshed(i_C)=0

                ! Reinitialize the discretization for this crack. Since the other cracks have already been
                ! discretized,
                ! So here, i_C is actually discretized again.
                CALL D3_Mesh_Initial_Crack

            endif

            !//////////////////////////////////////////////////////////
            ! If there are no more singularities, exit the i_Try loop.
            !//////////////////////////////////////////////////////////
            if(Yes_Singularity .eqv. .False.)then
                if(i_Try >1) write(*,102) i_Try
                exit
            endif
            ! If there are still anomalies after multiple attempts, a prompt will be displayed.
            if(i_Try==50) write(*,103) i_Try
        enddo

    !**************************************************************************
    !                                                                        *
    !                                                                        *
    !                                                                        *
    ! Option 2: For Heaviside singular cracks, adjust the crack centroid to  *
    !       the element center and normalize the normal vector.              *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !**************************************************************************
    case(2)
        ! Determine whether the Heaviside enrichment function distance for the crack is singular.
        call  D3_Check_Crack_Heaviside_Singularity(i_C,Singularity_Factor,Yes_Singularity)

        ! If there are anomalies, adjust the crack coordinates and orientation, and then regrid.
        if(Yes_Singularity)then
            ! Adjust the centroid to the center of the element.
            call D3_Offset_Initial_Crack_to_Ele_Center(i_C,c_Crack_Type )

            ! Adjust the normal vector to the nearest standard cardinal direction (0, 45, 90, 135, 180, 225,
            ! 270, etc.).
            !call D3_Offset_Initial_Crack_to_Regular_n_Vector(i_C,c_Crack_Type)

            ! Set the mark of the crack Cracks_Initial_Meshed to 0.
            Cracks_Initial_Meshed(i_C)=0

            ! Reinitialize the discretization of this crack. Since the other cracks have already been
            ! discretized,
            ! So here, i_C is actually discretized again.
            CALL D3_Mesh_Initial_Crack
        endif


    !*************************************************************************
    !                                                                        *
    !                                                                        *
    !                                                                        *
    ! Option 3: Regardless of whether it is singular, adjust the crack       *
    !        centroid to the element center and normalize the normal vector. *
    !       NEWFTU2022080201.                                                *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !*************************************************************************
    case(3)
        ! Adjust the centroid to the center of the element.
        call D3_Offset_Initial_Crack_to_Ele_Center(i_C,c_Crack_Type)

        ! Adjust the normal vector to the nearest standard direction (0, 45, 90, 135, 180, 225, 270, etc.).
        call D3_Offset_Initial_Crack_to_Regular_n_Vector(i_C,c_Crack_Type)

        ! Set the mark of the crack Cracks_Initial_Meshed to 0.
        Cracks_Initial_Meshed(i_C)=0

        ! Reinitialize the discretization of this crack. Since the other cracks have already been
        ! discretized,
        ! So here, i_C is actually discretized again.
        CALL D3_Mesh_Initial_Crack


    !**************************************************************************
    !                                                                        *
    !                                                                        *
    !                                                                        *
    ! Scheme 4: Normalize the crack normal vector regardless of whether      *
    !           it is singular. 2022-08-05.                                  *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !**************************************************************************
    case(4)
        ! Adjust the normal vector to the nearest standard cardinal direction (0, 45, 90, 135, 180, 225,
        ! 270, etc.).
        call D3_Offset_Initial_Crack_to_Regular_n_Vector(i_C,c_Crack_Type)
        ! Set the mark of the crack Cracks_Initial_Meshed to 0.
        Cracks_Initial_Meshed(i_C)=0
        ! Reinitialize the discretization of this crack. Since the other cracks have already been
        ! discretized,
        ! So here, i_C is actually discretized again.
        CALL D3_Mesh_Initial_Crack



    !*****************************************************************************
    !                                                                        *
    !                                                                        *
    !                                                                        *
    ! Plan 5: Adjust only the centroids of crack surfaces that are very close*
    !       to the element boundaries. NEWFTU2022092901. 2022-09-29.            *
    !                                                                        *
    !                                                                        *
    !                                                                        *
    !*****************************************************************************
    case(5)
        ! Obtain the centroid coordinates of the crack.
        call D3_Get_Crack_Center(i_C,c_Crack_Center)
#ifndef Silverfrost
        if (isnan(c_Crack_Center(1))) then
            print *, '    ERROR2024112901 :: Coodinates of crack center is NaN in D3_Check_and_Adjust_Cracks.f90! '
            print *, '                       Crack number:',i_C
            call Warning_Message('S',Keywords_Blank)
        end if
#endif
        ! Calculate the element number where the centroid of the crack surface is located.
        Ele_Num_Cache = 1
        call Cal_Ele_Num_by_Coors_3D(c_Crack_Center(1),c_Crack_Center(2),c_Crack_Center(3),Ele_Num_Cache,c_OUT_Elem)

        c_Logical_Positive_Num = 0
        c_Logical_Negative_Num = 0
        c_Logical_Positive = .False.
        c_Logical_Negative = .False.

        ! Loop through each node of the element.
        c_NN  = G_NN(:,c_OUT_Elem)
        c_X_NODES = G_X_NODES(1:8,c_OUT_Elem)
        c_Y_NODES = G_Y_NODES(1:8,c_OUT_Elem)
        c_Z_NODES = G_Z_NODES(1:8,c_OUT_Elem)

        if(allocated(c_Distance)) deallocate(c_Distance)
        allocate(c_Distance(Tem_Num_Gauss_3D))
        c_Distance(1:Tem_Num_Gauss_3D) = ZR

        do i_G = 1,Tem_Num_Gauss_3D
          call Cal_N_3D(kesi_Enr(i_G),yita_Enr(i_G),zeta_Enr(i_G),N)
          !Global coordinates of the gauss point.
          Global_coor_Gauss(1) = DOT_PRODUCT(N(1,1:24:3),c_X_NODES(1:8))
          Global_coor_Gauss(2) = DOT_PRODUCT(N(1,1:24:3),c_Y_NODES(1:8))
          Global_coor_Gauss(3) = DOT_PRODUCT(N(1,1:24:3),c_Z_NODES(1:8))

          !********************************************
          ! Calculate the signed distance from the 
          ! current Gauss point to the fracture plane.
          !********************************************
          Check_Ball_R  = 3.0D0*Ave_Elem_L

          ! If out-of-plane extension.
          if(Key_InPlane_Growth== 0 ) then
              call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,i_C,Check_Ball_R,c_Distance(i_G),   &
                            c_Signed_Dis_v2,                                                             &
                            c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis,c_n_Vector)
          ! For in-plane extension, there is no need to find the minimum symbol distance, resulting in low
          ! resource consumption. 2022-09-08. NEWFTU2022090801.
          elseif(Key_InPlane_Growth== 1) then
              call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Global_coor_Gauss,i_C,&
                                       Check_Ball_R,c_Distance(i_G),                        &
                                       c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),          &
                                       Yes_Found_Min_Signed_Dis,c_n_Vector)
              c_Signed_Dis_v2 = c_Distance(i_G)
          endif

          call Cal_Sign(c_Distance(i_G),c_Sign)

          if(c_Sign>=  (ONE -Tol_6)) then
              c_Logical_Positive = .True.
              c_Logical_Positive_Num = c_Logical_Positive_Num + 1
          endif
          if(c_Sign<= (-ONE +Tol_6)) then
              c_Logical_Negative = .True.
              c_Logical_Negative_Num = c_Logical_Negative_Num + 1
          endif

          ! If the current fracture surface has both + and - Gauss points relative to the element where the
          ! fracture centroid is located, and the number exceeds Must_Gauss_Number_3D,
          ! Then exit the loop.
          if(c_Logical_Positive .and. c_Logical_Negative .and. &
                c_Logical_Positive_Num>=Must_Gauss_Number_3D .and. &
                c_Logical_Negative_Num>=Must_Gauss_Number_3D) then
              goto 1000
              !exit
          endif
        enddo

        ! Adjust the centroid to the center of the element.
        call D3_Offset_Initial_Crack_to_Ele_Center(i_C,c_Crack_Type)
        print *,'    Center of crack ',i_C,' has been adjusted!'


        ! Set the mark of the crack Cracks_Initial_Meshed to 0.
        Cracks_Initial_Meshed(i_C)=0

        ! Reinitialize the discretization of this crack. Since the other cracks have already been
        ! discretized,
        ! So here, i_C is actually discretized again.
        CALL D3_Mesh_Initial_Crack
    end select

    1000 continue
enddo



RETURN
END SUBROUTINE D3_Check_and_Adjust_Cracks
