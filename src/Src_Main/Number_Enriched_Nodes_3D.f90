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
 
SUBROUTINE Number_Enriched_Nodes_3D(isub)
!     Number the enriched nodes(3D).

use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Filename

implicit none

integer,intent(in)::isub
integer i_C,i_N,i_H,i_Incl
integer i,j,max_Cr_Node
character(200) c_File_name_1
character(5) temp

integer el_max_dof(Num_Elem),i_E
integer c_NN(8)   
integer num_t,num_h,num_j,cnt,tem,i_f
integer max_DOF_element
integer tem_MDOF_3D


print *,'    Numbering enriched nodes...'

n_h_Node =0
n_t_Node =0
n_j_Node =0
n_hl_Node=0
!n_Incl_Node    = 0
Total_FD       = 0
Usual_Freedom  = 0
Enrich_Freedom = 0

Usual_Freedom  = 3*Num_Node


!------------------------------------
! Standard Tip Enhancement (4 items)
!------------------------------------
if (Key_TipEnrich==1) then      
      Num_F_Functions = 4
!-----------------------------
! Keep only the first item F1
!-----------------------------
elseif (Key_TipEnrich==2) then       
      Num_F_Functions = 1                 
endif
              
! Initialize c_POS
IF(ALLOCATED(c_POS_3D)) DEALLOCATE(c_POS_3D)
ALLOCATE(c_POS_3D(Num_Node,num_Crack))
c_POS_3D(1:Num_Node,1:num_Crack) = 0

!Loop through each crack.
! Strong data dependency, not suitable for OpenMP parallelization. 2022-08-22.
do i_C = 1,num_Crack
      !Loop through each node.
      do i_N = 1,Num_Node
          if (Enriched_Node_Type_3D(i_N,i_C).eq.2) then
              c_POS_3D(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*Num_F_Functions + n_j_Node) + 1
              n_h_Node = n_h_Node + 1
          elseif (Enriched_Node_Type_3D(i_N,i_C) .eq.1)then
              c_POS_3D(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*Num_F_Functions+n_j_Node)+1
              n_t_Node = n_t_Node + 1               
          elseif (Enriched_Node_Type_3D(i_N,i_C) .eq.3)then
              c_POS_3D(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*Num_F_Functions + n_j_Node) + 1
              n_j_Node = n_j_Node + 1
          end if
      end do
end do     

! Maximum enhancement node number (of the crack)
if(num_Crack >=1)then
      max_Cr_Node = maxval(c_POS_3D(1:Num_Node,1:num_Crack))
else
      max_Cr_Node = Num_Node
endif

!Total degrees of freedom. 
Total_FD=3*(Num_Node+n_h_Node + n_t_Node*Num_F_Functions + n_j_Node)

Enrich_Freedom = Total_FD - Usual_Freedom

print *,'    Number of tip enriched nodes: ',n_t_Node
print *,'    Number of Heaviside enriched nodes: ',n_h_Node 
print *,'    Number of Junction enriched nodes: ',n_j_Node 


!--------------------------------------------------------------------------------------------------
! Loop through each element and check whether it exceeds the total number of structural degrees of
! freedom mdof_3d of the element.
! 2022-05-01. NEWFTU2022050101.
!--------------------------------------------------------------------------------------------------
el_max_dof(1:Num_Elem) = 0
! OpenMP Parallelization. 2022-08-22. IMPROV2022082201.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_NN,i_C,num_t,num_h,  &
!$OMP     num_j,tem,cnt,i_N,i_f)        
do i_E = 1,Num_Elem
      c_NN    = G_NN(:,i_E)
      do i_C=1,num_Crack
          !If the element has no enriched nodes, then:
          if(maxval(Enriched_Node_Type_3D(c_NN,i_C)).eq.0 .and. &
             minval(Enriched_Node_Type_3D(c_NN,i_C)).eq.0) then
              if (i_C .eq. 1)then
                  el_max_dof(i_E) =24
              end if
          !If the element has enriched nodes, then:      
          else
              !Get the number of the tip enriched nodes of the element.
              num_t = count(Enriched_Node_Type_3D(c_NN,i_C).eq.1)
              !Get the number of the Heaviside enriched nodes of the element.
              num_h = count(Enriched_Node_Type_3D(c_NN,i_C).eq.2)
              !Get the number of the tip junction nodes of the element.     
              num_j = count(Enriched_Node_Type_3D(c_NN,i_C).eq.3)
              
              tem = 3*(num_h*1 + num_t*Num_F_Functions + num_j*1) 
              cnt = 0
              
              do i_N = 1,8
                if ( Enriched_Node_Type_3D(c_NN(i_N),i_C)==2 )then
                    cnt = cnt + 1
                elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)==1)then
                    do i_f=1,Num_F_Functions
                          cnt = cnt + 1
                    end do               
                elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)==3)then
                       cnt = cnt + 1
                end if
              end do   
              if (i_C .eq. 1)then
                  ! Includes FEM degrees of freedom
                  el_max_dof(i_E)      = el_max_dof(i_E)+24+tem    
              else
                  el_max_dof(i_E)      = el_max_dof(i_E)+tem     
              end if
          end if
      enddo
enddo
!$OMP END PARALLEL DO


!IMPROV2023011101. 2023-01-11.
201 FORMAT(5X,'Element number of max DOF: ',I7)  
if(isub==1)then
      MDOF_3D = maxval(el_max_dof(1:Num_Elem))
      max_DOF_element =  maxloc(el_max_dof(1:Num_Elem),1)
      print *, '    MDOF_3D: ' ,MDOF_3D
      write(*,201) max_DOF_element 
      Old_MDOF_3D = MDOF_3D
else
      Old_MDOF_3D = MDOF_3D
      tem_MDOF_3D = maxval(el_max_dof(1:Num_Elem))
      if(tem_MDOF_3D > Old_MDOF_3D)then
          ! Update MDOF_3D
          MDOF_3D = tem_MDOF_3D
          max_DOF_element =  maxloc(el_max_dof(1:Num_Elem),1)
          print *, '    Updated MDOF_3D: ' ,MDOF_3D
          write(*,201) max_DOF_element  
      else
          print *, '    MDOF_3D: ' ,MDOF_3D
      endif
endif

RETURN
END SUBROUTINE Number_Enriched_Nodes_3D
