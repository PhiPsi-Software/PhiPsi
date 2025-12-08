 
      subroutine Cal_B_Matrix_Hl(kesi,yita,
     &                        i_H,i_E,i_G,
     &                        c_NN,c_X_NODES,c_Y_NODES,
     &                        tem_B,num_tem_B)  

      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      
      implicit none
      
      integer,intent(in)::i_H,i_E,i_G
      real(kind=FT),intent(in)::c_X_NODES(4),c_Y_NODES(4)
      integer,intent(in)::c_NN(4)
      real(kind=FT),intent(in)::kesi,yita
      real(kind=FT),intent(out)::tem_B(3,80)
      integer,intent(out)::num_tem_B        
      real(kind=FT) detJ, J(2,2), Inverse_J(2,2)
      real(kind=FT) N(2,8),dNdkesi(4,2),dNdx(4,2)
      integer mat_num,c_mat_type
      real(kind=FT) B_FEM(3,8),B_XFEM(3,60),BI_enr(3,2)
      integer num_B_XFEM
      integer i_N
      real(kind=FT) Hl_gp
      real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r
      real(kind=FT) Hl_i,c_G_x,c_G_y,Tool_Function_2Point_Dis
      real(kind=FT) c_Dis_G,c_Dis_N
      integer c_NODE
      real(kind=FT) x0_Hole,y0_Hole,R0_Hole
      real(kind=FT) a_Hole,b_Hole,theta_Hole,Tol  
      integer Return_Statu_G,Return_Statu_N
      
      call Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                           c_X_NODES,c_Y_NODES,
     &                           detJ,J,N,dNdkesi)    
     
      call Matrix_Inverse_2x2(J,Inverse_J)

      dNdx = MATMUL(dNdkesi,Inverse_J)
      
      c_G_x = DOT_PRODUCT(N(1,1:7:2),c_X_NODES(1:4))
      c_G_y = DOT_PRODUCT(N(1,1:7:2),c_Y_NODES(1:4))      
      
      
      mat_num    = Elem_Mat(i_E)
      
      c_mat_type = Material_Type(mat_num)     
      
      B_FEM(1:3,1:8) = ZR

      if (EleGaus_yes_FEM_asemd(i_E,i_G) .eqv. .False.) then
          B_FEM(1,1:8:2)   =  dNdx(:,1)
            B_FEM(2,2:8:2)   =  dNdx(:,2)
            B_FEM(3,1:8:2)   =  dNdx(:,2)
            B_FEM(3,2:8:2)   =  dNdx(:,1)
          EleGaus_yes_FEM_asemd(i_E,i_G)=.True.
      end if

      if(maxval(Enriched_Node_Type_Hl(c_NN,i_H)).eq.0 .and.
     &   minval(Enriched_Node_Type_Hl(c_NN,i_H)).eq.0) then
          if (i_H.eq.1) then
              tem_B(1:3,1:8) = B_FEM
              num_tem_B = 8
          else
              num_tem_B = 0
          end if
      elseif(maxval(Enriched_Node_Type_Hl(c_NN,i_H)).gt.0 .or.
     &       minval(Enriched_Node_Type_Hl(c_NN,i_H)).gt.0) then 
          B_XFEM(1:3,1:60) = ZR
          num_B_XFEM = 0
          do i_N = 1,4
              if (Enriched_Node_Type_Hl(c_NN(i_N),i_H).eq.1)then
                if (Elem_Type_Hl(i_E,i_H).eq.1) then
                      c_NODE  = Elem_Node(i_E,i_N)     
                      if (num_Circ_Hole>=1)then
                          c_Hole_x = Hole_Coor(i_H,1)
                          c_Hole_y = Hole_Coor(i_H,2)
                          c_Hole_r = Hole_Coor(i_H,3)
                          c_Dis_G=Tool_Function_2Point_Dis(
     &                          [c_G_x,c_G_y],[c_Hole_x,c_Hole_y])
                          c_Dis_N=Tool_Function_2Point_Dis(
     &                           Coor(c_NODE,1:2),[c_Hole_x,c_Hole_y])
                          if (c_Dis_G<c_Hole_r)then
                              if(Key_Hole_Value==0)then
                                  Hl_gp = ZR
                              elseif(Key_Hole_Value==-1)then
                                  Hl_gp = -ONE
                              endif
                          else
                              Hl_gp = ONE
                          endif
                          
                          if (c_Dis_N<c_Hole_r)then
                              if(Key_Hole_Value==0)then
                                  Hl_i = ZR
                              elseif(Key_Hole_Value==-1)then
                                  Hl_i = -ONE
                              endif
                          else
                              Hl_i = ONE
                          endif
                      endif
                      if (num_Ellip_Hole>=1)then
                          x0_Hole =  Ellip_Hole_Coor(i_H,1)
                          y0_Hole =  Ellip_Hole_Coor(i_H,2)
                          a_Hole  =  Ellip_Hole_Coor(i_H,3) 
                          b_Hole  =  Ellip_Hole_Coor(i_H,4) 
                          theta_Hole =  Ellip_Hole_Coor(i_H,5)

                          Tol = 1.0D-10
                          call Tool_Yes_Point_in_Oblique_Ellipse(
     &                         [c_G_x,c_G_y],
     &                         x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,
     &                         Return_Statu_G,Tol)
                          call Tool_Yes_Point_in_Oblique_Ellipse(
     &                         Coor(c_NODE,1:2),
     &                         x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,
     &                         Return_Statu_N,Tol)
     
                          if (Return_Statu_G==1)then
                              if(Key_Hole_Value==0)then
                                  Hl_gp = ZR
                              elseif(Key_Hole_Value==-1)then
                                  Hl_gp = -ONE
                              endif
                          else
                              Hl_gp = ONE
                          endif
                          
                          if (Return_Statu_N==1)then
                              if(Key_Hole_Value==0)then
                                  Hl_i = ZR
                              elseif(Key_Hole_Value==-1)then
                                  Hl_i = -ONE
                              endif
                          else
                              Hl_i = ONE
                          endif
                      endif
                      
            BI_enr(1,:) = [dNdx(i_N,1)*(Hl_gp-Hl_i), ZR]
                      BI_enr(2,:) = [ZR,dNdx(i_N,2)*(Hl_gp-Hl_i)]
                      BI_enr(3,:) = [dNdx(i_N,2)*(Hl_gp-Hl_i), 
     &                               dNdx(i_N,1)*(Hl_gp-Hl_i)]
                  else
                      BI_enr(1:3,1:2)=ZR
                  end if
                  B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
                  num_B_XFEM = num_B_XFEM + 2                  
              end if
          end do
          if (i_H.eq.1) then
              tem_B(1:3,1:8)            = B_FEM
              tem_B(1:3,9:8+num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
              num_tem_B = 8 + num_B_XFEM
          else
              tem_B(1:3,1:num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
              num_tem_B = num_B_XFEM
          end if
      end if   

      RETURN
      end subroutine Cal_B_Matrix_Hl
