 
      subroutine Cal_Gauss_Coors(isub,Total_Num_G_P,
     &                           c_G_CoorX,c_G_CoorY)
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use Global_Inclusion
      use Global_Cross
      use omp_lib
      
      implicit none
      integer,intent(in)::isub,Total_Num_G_P
      real(kind=FT),intent(out)::c_G_CoorX(Total_Num_G_P),
     &                           c_G_CoorY(Total_Num_G_P)
      
      integer i_E,i_G,c_NN(4),c_Num_Gauss_Point,Gauss_Counter
      real(kind=FT) kesi_Enr(Num_Gauss_Points),
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)    
      real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM),
     &                 yita_N_Enr(Num_Gauss_P_FEM),
     &                 weight_N_Enr(Num_Gauss_P_FEM)    
      real(kind=FT) kesi(900),yita(900)
      real(kind=FT) Out_x,Out_y    
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64) 
      integer G_Counter
     
      print *,'    Calculating coordinates of Gauss points...'
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,
     &                             kesi_Enr,
     &                             yita_Enr,
     &                             weight_Enr)
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                             weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,
     &                             kesi_N_Enr,
     &                             yita_N_Enr,
     &                             weight_N_Enr)
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
     
     
      Gauss_Counter = 0
      if (Yes_XFEM.eqv..False.) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_G,
!$OMP&        c_NN,c_X_NODES,c_Y_NODES,kesi,yita,G_Counter,Out_x,Out_y) 
          do i_E = 1,Num_Elem
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)                
              kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
              yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
              do i_G = 1,Num_Gauss_P_FEM
                  G_Counter =  (i_E-1)*Num_Gauss_P_FEM +i_G 
                  call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G),
     &                                      c_X_NODES,c_Y_NODES,
     &                                      Out_x,Out_y)
                  c_G_CoorX(G_Counter) = Out_x
                  c_G_CoorY(G_Counter) = Out_y
              end do 
          enddo
!$omp end parallel do
      elseif(Yes_XFEM.eqv..True.) then
          do i_E = 1,Num_Elem
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)                 
              if(Key_Integral_Sol.eq.1)then
              elseif(Key_Integral_Sol.eq.2 .or. 
     &               Key_Integral_Sol.eq.3)then
                  if (num_Crack/= 0 .and. 
     &                (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).
     &                                                       ne.0))then
                      kesi(1:Num_Gauss_Points)    = kesi_Enr
                      yita(1:Num_Gauss_Points)    = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  elseif(num_Hole/= 0 .and.
     &                (maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).
     &                                                     ne.0))then
                      kesi(1:Num_Gauss_Points)    = kesi_Enr
                      yita(1:Num_Gauss_Points)    = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  elseif(num_Cross/= 0 .and.
     &               (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross))
     &                              .ne.0))then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  elseif(num_Inclusion/= 0 .and.
     &                (maxval(Enriched_Node_Type_Incl
     &                                    (c_NN,1:num_Inclusion)).ne.0))
     &                                                             then
                      kesi(1:Num_Gauss_Points)   = kesi_Enr
                      yita(1:Num_Gauss_Points)   = yita_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  else 
                      kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
                      yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
                      c_Num_Gauss_Point = Num_Gauss_P_FEM
                  end if 
                  do i_G = 1,c_Num_Gauss_Point
                      Gauss_Counter =  Gauss_Counter +1 
                      call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G),
     &                                          C_X_NODES,C_Y_NODES,
     &                                          Out_x,Out_y)
                      c_G_CoorX(Gauss_Counter) = Out_x
                      c_G_CoorY(Gauss_Counter) = Out_y

                  end do             
              end if            
          end do   
      end if
      
      return 
      end SUBROUTINE Cal_Gauss_Coors   
