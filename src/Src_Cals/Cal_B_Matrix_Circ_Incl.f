 
      subroutine Cal_B_Matrix_Circ_Incl(kesi,yita,
     &                        i_Incl,i_E,i_G,
     &                        c_NN,c_X_NODES,c_Y_NODES,
     &                        tem_B,num_tem_B)  

      use Global_Float_Type
      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use Global_Inclusion
      
      implicit none
      
      integer,intent(in)::i_Incl,i_E,i_G
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
      real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r
      real(kind=FT) c_G_x,c_G_y,Tool_Function_2Point_Dis
      real(kind=FT) c_Dis_N
      real(kind=FT) Incl_gp
      integer c_NODE
      real(kind=FT) c_N(4)
      real(kind=FT) Zeta_Node(4),Zeta_Node_abs(4)
      real(kind=FT) Zm,Za,Zmn,Ex,Ey
      
      tem_B(1:3,1:80) = ZR
       
      call Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                           c_X_NODES,c_Y_NODES,
     &                           detJ,J,N,dNdkesi)    
     
      call Matrix_Inverse_2x2(J,Inverse_J)

      dNdx = MATMUL(dNdkesi,Inverse_J)
      
      call Cal_Coor_by_KesiYita(kesi,yita,c_X_NODES,c_Y_NODES,
     &                                    c_G_x,c_G_y) 
      
      mat_num    = Elem_Mat(i_E)
      
      c_mat_type = Material_Type(mat_num)     
      
      B_FEM(1:3,1:8) = ZR
      
      if (EleGaus_yes_FEM_asemd(i_E,i_G) .eqv. .False.) then
          B_FEM(1,1:8:2)   =  dNdx(:,1)
            B_FEM(2,2:8:2)   =  dNdx(:,2)
            B_FEM(3,1:8:2)   =  dNdx(:,2)
            B_FEM(3,2:8:2)   =  dNdx(:,1)
          EleGaus_yes_FEM_asemd(i_E,i_G)= .True.
      end if
      
      if(maxval(Enriched_Node_Type_Incl(c_NN,i_Incl)).eq.0 .and.
     &   minval(Enriched_Node_Type_Incl(c_NN,i_Incl)).eq.0) then
          if (i_Incl.eq.1) then

              tem_B(1:3,1:8) = B_FEM
              num_tem_B = 8
          else
              num_tem_B = 0
          end if
      elseif(maxval(Enriched_Node_Type_Incl(c_NN,i_Incl)).gt.0 .or.
     &       minval(Enriched_Node_Type_Incl(c_NN,i_Incl)).gt.0) then 
          B_XFEM(1:3,1:60) = ZR
          num_B_XFEM = 0
          Zeta_Node(1:4) = ZR
          do i_N = 1,4 
              c_NODE  = Elem_Node(i_E,i_N)  
              c_Incl_x = Circ_Inclu_Coor(i_Incl,1)
              c_Incl_y = Circ_Inclu_Coor(i_Incl,2)
              c_Incl_r = Circ_Inclu_Coor(i_Incl,3)
              c_Dis_N=Tool_Function_2Point_Dis(Coor(c_NODE,1:2),
     &                          [c_Incl_x,c_Incl_y])
              Zeta_Node(i_N) = c_Dis_N-c_Incl_r
              if(abs(Zeta_Node(i_N))<1.0D-8)then
                  Zeta_Node(i_N) = ZR
              endif
          enddo
          Zeta_Node_abs = abs(Zeta_Node)
          c_N(1) = N(1,1)
          c_N(2) = N(1,3)
          c_N(3) = N(1,5)
          c_N(4) = N(1,7)
          Zm = Zeta_Node(1)*c_N(1) + Zeta_Node(2)*c_N(2)+
     &         Zeta_Node(3)*c_N(3) + Zeta_Node(4)*c_N(4)
          Za = Zeta_Node_abs(1)*c_N(1)+Zeta_Node_abs(2)*c_N(2)+
     &         Zeta_Node_abs(3)*c_N(3)+Zeta_Node_abs(4)*c_N(4)
          Zmn  = Zm/abs(Zm)
          do i_N = 1,4 
              if (Enriched_Node_Type_Incl(c_NN(i_N),i_Incl).eq.1)then
                if (Elem_Type_Incl(i_E,i_Incl).eq.1) then
                      Incl_gp = Za-abs(Zm)                 
                      Ex =     dNdx(1,1)*Zeta_Node_abs(1)+
     &                         dNdx(2,1)*Zeta_Node_abs(2)+
     &                         dNdx(3,1)*Zeta_Node_abs(3)+
     &                         dNdx(4,1)*Zeta_Node_abs(4)-
     &                    Zmn*(dNdx(1,1)*Zeta_Node(1)+
     &                         dNdx(2,1)*Zeta_Node(2)+
     &                         dNdx(3,1)*Zeta_Node(3)+
     &                         dNdx(4,1)*Zeta_Node(4))
                      Ey =     dNdx(1,2)*Zeta_Node_abs(1)+
     &                         dNdx(2,2)*Zeta_Node_abs(2)+
     &                         dNdx(3,2)*Zeta_Node_abs(3)+
     &                         dNdx(4,2)*Zeta_Node_abs(4)-
     &                    Zmn*(dNdx(1,2)*Zeta_Node(1)+
     &                         dNdx(2,2)*Zeta_Node(2)+
     &                         dNdx(3,2)*Zeta_Node(3)+
     &                         dNdx(4,2)*Zeta_Node(4))
            BI_enr(1,:)=[dNdx(i_N,1)*Incl_gp+c_N(i_N)*Ex,ZR]
                      BI_enr(2,:)=[ZR,dNdx(i_N,2)*Incl_gp+c_N(i_N)*Ey]
                      BI_enr(3,:)=[dNdx(i_N,2)*Incl_gp+c_N(i_N)*Ey, 
     &                             dNdx(i_N,1)*Incl_gp+c_N(i_N)*Ex]
                  else
                      BI_enr(1:3,1:2)=ZR
                  end if
                  B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
                  num_B_XFEM = num_B_XFEM + 2                  
              end if
          end do
          
          if (i_Incl.eq.1) then
            tem_B(1:3,1:8)            = B_FEM
            tem_B(1:3,9:8+num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
            num_tem_B = 8 + num_B_XFEM
          else
            tem_B(1:3,1:num_B_XFEM)   = B_XFEM(1:3,1:num_B_XFEM)
            num_tem_B = num_B_XFEM
          end if
      end if   

      RETURN
      end SUBROUTINE Cal_B_Matrix_Circ_Incl