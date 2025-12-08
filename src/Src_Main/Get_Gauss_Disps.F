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
 
      SUBROUTINE Get_Gauss_Disps(isub,c_DISP,Total_Num_G_P)
c     Calculate Gauss point displacement
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use Global_DISP
      use Global_Inclusion
      use Global_Cross
#ifndef Silverfrost
      use omp_lib
#endif      
      implicit none
      !include 'omp_lib.h'
      integer,intent(in)::isub,Total_Num_G_P
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      integer i_E,i_G,c_NN(4),c_Num_Gauss_Point,G_Counter
      real(kind=FT)    kesi_Enr(Num_Gauss_Points),
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)    
      real(kind=FT)    kesi_N_Enr(Num_Gauss_P_FEM),
     &                 yita_N_Enr(Num_Gauss_P_FEM),
     &                 weight_N_Enr(Num_Gauss_P_FEM)    
      real(kind=FT) kesi(900),yita(900)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),P_Disp(2)
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64)      
     
      print *,'    Calculating displacements of Gauss points...'
      ! Standard 64 Gauss points
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,
     &                           weight_Enr)
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      ! Subdivision of quadrilaterals to calculate Gauss points (Gauss points are unrelated to cracks,
      ! consisting of several regular quadrilaterals)
      elseif (Key_Integral_Sol  == 3)then
          call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                           kesi_Enr,yita_Enr,
     &                           weight_Enr)
          Num_Gauss_Points = Num_Sub_Quads*4
          Num_Gauss_P_Inc  = Num_Sub_Quads*4
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      endif
     
     
      G_Counter = 0
      
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Crack-free static analysis (the simplest case)
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (Yes_XFEM.eqv..False.) then 
          !.............................
          ! OpenMP multi-core computing
          !.............................
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_G,
!$OMP&         c_NN,c_X_NODES,c_Y_NODES,kesi,yita,P_Disp,G_Counter) 
          do i_E = 1,Num_Elem
              !c_NN    = G_NN(i_E,:)
              !c_X_NODES = G_X_NODES(i_E,:)
              !c_Y_NODES = G_Y_NODES(i_E,:)
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)
              kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
              yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
              ! Loop over each Gauss point.
              do i_G = 1,Num_Gauss_P_FEM
                  !G_Counter =  G_Counter +1 
                  G_Counter =  (i_E-1)*Num_Gauss_P_FEM +i_G 
                  call Cal_Any_Point_Disp_KesiYita(i_E,kesi(i_G),
     &                                             yita(i_G),i_G,
     &                                             c_DISP,P_Disp)
                  DISP_x_Gauss(G_Counter) = P_Disp(1)
                  DISP_y_Gauss(G_Counter) = P_Disp(2)                  
              end do 
          enddo   
!$omp end parallel do
c     %%%%%%%%%
c       Others
c     %%%%%%%%%
      else
          do i_E = 1,Num_Elem
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)              
              ! ------------------------------
              ! Point Scheme 1: Triangulation
              ! ------------------------------
              if(Key_Integral_Sol.eq.1)then
                  !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)	
                  !c_Num_Gauss_Point  = size(kesi,2);	
              ! ----------------------------------------------------------------------------------
              ! Integration scheme 2 or 3: a fixed number of integration points, Num_Gauss_Points
              ! ----------------------------------------------------------------------------------
              elseif(Key_Integral_Sol.eq.2 .or.
     &               Key_Integral_Sol.eq.3)then

                  !if the current element are enriched element, then 8x8 gauss points is suggested:
                  if((num_Crack/= 0) .and. 
     &                maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)
     &                                                     then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  ! If it is a Hole enhanced node, then 8x8 Gauss points are suggested:
                  elseif(num_Hole/= 0 .and.
     &                (maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole))
     &                              .ne.0))then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  ! If it is a Cross enhanced node, then 8x8 Gauss points are suggested:
                  elseif(num_Cross/= 0 .and.
     &               (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross))
     &                              .ne.0))then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  ! If it is a hybrid reinforcement node, it should be divided into two situations. For reinforcement
                  ! elements that include material interfaces, use at least 400 Gauss points.
                  ! For a general unit, 64 integration points are used.
                  elseif(num_Inclusion/= 0 .and.
     &                (maxval(Enriched_Node_Type_Incl
     &                                    (c_NN,1:num_Inclusion)).ne.0))
     &                                                             then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  !if the current element are not enriched element, then 2x2 gauss points:
                  else 
                      kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
                      yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
                      c_Num_Gauss_Point = Num_Gauss_P_FEM
                  end if 
                            
                  ! Loop over each Gauss point.
                  do i_G = 1,c_Num_Gauss_Point
                      G_Counter =  G_Counter +1 
                      call Cal_Any_Point_Disp_KesiYita(i_E,kesi(i_G),
     &                                                 Yita(i_G),i_G,
     &                                                 c_DISP,P_Disp)
                      DISP_x_Gauss(G_Counter) = P_Disp(1)
                      DISP_y_Gauss(G_Counter) = P_Disp(2)
                  end do             
              end if            
          end do   
      end if
     
      RETURN
      END SUBROUTINE Get_Gauss_Disps
