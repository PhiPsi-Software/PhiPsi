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
 
      SUBROUTINE Force_Vector_3D_no_Killed_Ele(c_Total_FD,
     &                                         isub,Lambda,globalF)
      ! Load vector (excluding active and inactive elements, and not considering the external load
      ! correction for initial stress from Zienkiewicz_7ed_P216; used for initial stress field
      ! calculation).

c     ----------------------------
c     Read Public Variable Module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      
      implicit none

      integer,intent(in)::isub,c_Total_FD
      real(kind=FT),intent(in)::Lambda
      real(kind=FT),intent(out)::globalF(c_Total_FD)
      
      real(kind=FT) h,density
      integer i,i_E,i_N,cur_Node,mat_num,c_mat_type,c_node
      
      print *,'    Constructing global force vector...'        
      
      globalF(1:c_Total_FD) = ZR
      
      do i = 1,Num_Foc_x
            cur_Node = int(Foc_x(i,1))
          globalF(3*cur_Node-2) = Lambda*Foc_x(i,2)                   
      end do

      do i = 1,Num_Foc_y
            cur_Node = int(Foc_y(i,1))
            globalF(3*cur_Node-1) =   Lambda*Foc_y(i,2)                  
      end do

      do i = 1,Num_Foc_z
            cur_Node = int(Foc_z(i,1))
            globalF(3*cur_Node) =   Lambda*Foc_z(i,2)                  
      end do     
      
      if(Key_Gravity==1) then
          do i_E = 1,Num_Elem
                mat_num    = Elem_Mat(i_E)
                c_mat_type = Material_Type(mat_num)
              if (c_mat_type ==1) then
                  density = Material_Para(mat_num,3)
              else
            density = Material_Para(mat_num,3); 
              end if
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node-2)  = globalF(3*c_node-2)-
     &                          g_X_Y_Z(1)*density*Elem_Vol(i_E)/EIG
              end do              
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node-1)  = globalF(3*c_node-1)-
     &                          g_X_Y_Z(2)*density*Elem_Vol(i_E)/EIG
              end do
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node)  = globalF(3*c_node)-
     &                          g_X_Y_Z(3)*density*Elem_Vol(i_E)/EIG
              end do
            end do
      end if
      

      
      
  199 continue 
  
  
      RETURN
      END SUBROUTINE Force_Vector_3D_no_Killed_Ele
