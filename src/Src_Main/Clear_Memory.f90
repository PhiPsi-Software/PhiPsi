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
 
SUBROUTINE Clear_Memory
! Clear the allocatable array variables in memory.
! 2022-08-04.
      
!.............................
! Read public variable module
!.............................
use Global_Float_Type
use Global_Common   
use Global_Filename
use Global_Model
use Global_XFEM_Elements
use Global_Elem_Area_Vol
use Global_Crack_3D
use Global_HF
use Global_Inclusion
use Global_Contact
use Global_Material

!...........................................
! Clear the memory of the related variables
!...........................................
if(allocated(material_type))  deallocate(material_type)
if(allocated(Coor))           deallocate(Coor)
if(allocated(Elem_Node))      deallocate(Elem_Node)

if(allocated(Elem_Mat))       deallocate(Elem_Mat)
if(allocated(Bou_x))          deallocate(Bou_x)
if(allocated(Bou_y))          deallocate(Bou_y)
if(allocated(Bou_z))          deallocate(Bou_z)

if(allocated(Foc_x))          deallocate(Foc_x)
if(allocated(Foc_y))          deallocate(Foc_y)
if(allocated(Foc_z))          deallocate(Foc_z)

if(allocated(G_NN))           deallocate(G_NN)
if(allocated(G_X_NODES))      deallocate(G_X_NODES)
if(allocated(G_Y_NODES))      deallocate(G_Y_NODES)
if(allocated(G_Z_NODES))      deallocate(G_Z_NODES)

if(allocated(Elem_Vol))               deallocate(Elem_Vol)
if(allocated(Elem_Centroid))          deallocate(Elem_Centroid)
if(allocated(EleGaus_yes_FEM_asemd))  deallocate(EleGaus_yes_FEM_asemd)
if(allocated(x_max_Elements))         deallocate(x_max_Elements)
if(allocated(x_min_Elements))         deallocate(x_min_Elements)
if(allocated(y_max_Elements))         deallocate(y_max_Elements)
if(allocated(y_min_Elements))         deallocate(y_min_Elements)
if(allocated(z_max_Elements))         deallocate(z_max_Elements)
if(allocated(z_min_Elements))         deallocate(z_min_Elements)

if(allocated(Elements_Gauss_Num))              deallocate(Elements_Gauss_Num)
if(allocated(Elems_Block_Bou))                 deallocate(Elems_Block_Bou)
if(allocated(Outline))                         deallocate(Outline)
if(allocated(OutArea))                         deallocate(OutArea)
if(allocated(Surface_Nodes))                   deallocate(Surface_Nodes)
if(allocated(Node_Elements))                   deallocate(Node_Elements)
if(allocated(num_Node_Elements))               deallocate(num_Node_Elements)
if(allocated(Element_Edges))                   deallocate(Element_Edges)
if(allocated(c_POS_3D))                        deallocate(c_POS_3D)
if(allocated(num_GP_Elem))                     deallocate(num_GP_Elem)
if(allocated(Ele_GP_Start_Num))                deallocate(Ele_GP_Start_Num)
if(allocated(Elem_num_Related_Cracks))         deallocate(Elem_num_Related_Cracks)
if(allocated(Elem_Conta_Sta))                  deallocate(Elem_Conta_Sta)
if(allocated(Flag_FreeDOF))                    deallocate(Flag_FreeDOF)
if(allocated(Location_FreeDOF))                deallocate(Location_FreeDOF)

if(allocated(D))                    deallocate(D)
if(allocated(D_Comp))               deallocate(D_Comp)
if(allocated(S))                    deallocate(S)
if(allocated(St))                   deallocate(St)
if(allocated(Sc))                   deallocate(Sc)
if(allocated(T_Alpha))              deallocate(T_Alpha)
if(allocated(KIc))                  deallocate(KIc)
if(allocated(E))                    deallocate(E)
if(allocated(v))                    deallocate(v)
if(allocated(density))              deallocate(density)
if(allocated(D_for_cylindrical))    deallocate(D_for_cylindrical)

if(allocated(size_local_Old))               deallocate(size_local_Old)
if(allocated(Dis_Node_to_FS))               deallocate(Dis_Node_to_FS)
if(allocated(Elem_Type_3D))                    deallocate(Elem_Type_3D)
if(allocated(Enriched_Node_Type_3D))           deallocate(Enriched_Node_Type_3D)
if(allocated(Enriched_Node_Crack_n_Vector_3D)) deallocate(Enriched_Node_Crack_n_Vector_3D)
if(allocated(Ele_Num_Tip_Enriched_Node_3D))    deallocate(Ele_Num_Tip_Enriched_Node_3D)      
RETURN
END SUBROUTINE Clear_Memory
