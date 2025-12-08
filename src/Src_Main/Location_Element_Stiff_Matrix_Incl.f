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
 
      SUBROUTINE Location_Element_Stiff_Matrix_Incl(i_E,i_Incl,
     &                                       POS_c_Incl,
     &                                       Location_ESM_C_Crack,
     &                                       num_Loc_ESM_C_Crack,
     &                                       Location_ESM_C_Cr_NoFEM,
     &                                       num_Loc_ESM_C_Cr_NoFEM)
c     The position of the element stiffness matrix in the total stiffness.
      use Global_Float_Type
      use Global_Crack
      use Global_Model
      use Global_Common
      use Global_Inclusion
   
      implicit none
      
      integer,intent(in)::i_E,i_Incl,POS_c_Incl(Num_Node)
      integer,intent(out)::Location_ESM_C_Crack(80)
      integer,intent(out)::Location_ESM_C_Cr_NoFEM(60)
      integer,intent(out)::num_Loc_ESM_C_Crack,num_Loc_ESM_C_Cr_NoFEM
      integer::Location_FEM(8)     
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4),i_N,cnt      
      integer Enriched_Node_c_Incl(Num_Node)
      integer num_Incl
      integer:: Location_XFEM(60)
      integer tem
      
      Enriched_Node_c_Incl = Enriched_Node_Type_Incl(:,i_Incl)
      
      c_NN    = G_NN(:,i_E)
      c_X_NODES = G_X_NODES(:,i_E)
      c_Y_NODES = G_Y_NODES(:,i_E)
      ! Vector initialization
      Location_FEM(1:8)=0
      Location_XFEM(1:60)=0
      Location_ESM_C_Crack(1:80)=0
      Location_ESM_C_Cr_NoFEM(1:60) = 0
      
      num_Loc_ESM_C_Crack = 0
      num_Loc_ESM_C_Cr_NoFEM = 0
      
      if (i_Incl .eq. 1)then
          do i_N =1,4
              Location_FEM(2*i_N-1) = 2*c_NN(i_N) - 1
              Location_FEM(2*i_N)   = 2*c_NN(i_N)
            end do
      end if
      
      !If the element has no enriched nodes, then:
      if(maxval(Enriched_Node_c_Incl(c_NN)).eq.0 .and.
     &   minval(Enriched_Node_c_Incl(c_NN)).eq.0) then
          if (i_Incl.eq. 1)then
              Location_ESM_C_Crack(1:8) = Location_FEM
              num_Loc_ESM_C_Crack =8
              num_Loc_ESM_C_Cr_NoFEM = 0
          end if
      else
      
          ! Get the number of Inclusion enhancement nodes in the current cell.
          num_Incl = count(Enriched_Node_c_Incl(c_NN).eq.1)
          tem = 2*num_Incl  
          cnt = 0
          do i_N = 1,4
              if (Enriched_Node_c_Incl(c_NN(i_N)) .eq. 1)then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Incl(c_NN(i_N)) - 1
            Location_XFEM(2*cnt  ) = 2*POS_c_Incl(c_NN(i_N))
               end if
          end do   
          !-----------------------
          if (i_Incl .eq. 1)then
              ! Includes FEM degrees of freedom
              Location_ESM_C_Crack(1:8)     = Location_FEM
              Location_ESM_C_Crack(9:8+tem) = Location_XFEM(1:tem)
              num_Loc_ESM_C_Crack           = 8+tem    
              ! Does not include FEM degrees of freedom
              Location_ESM_C_Cr_NoFEM(1:tem)= Location_XFEM(1:tem)
              num_Loc_ESM_C_Cr_NoFEM        = tem
          else
              ! Includes FEM degrees of freedom
              Location_ESM_C_Crack(1:tem)   = Location_XFEM(1:tem)
              num_Loc_ESM_C_Crack           = tem 
              ! Does not include FEM degrees of freedom
              Location_ESM_C_Cr_NoFEM(1:tem)= Location_XFEM(1:tem)
              num_Loc_ESM_C_Cr_NoFEM        = tem              
          end if          
      end if
      
      RETURN
      END SUBROUTINE Location_Element_Stiff_Matrix_Incl
