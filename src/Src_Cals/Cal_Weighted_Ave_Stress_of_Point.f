 
      subroutine Cal_Weighted_Ave_Stress_of_Point(c_Key_Ave_Stress,
     &                       Point,Search_R,a_Weight,l_Ordina,
     &                       WA_S_xx,WA_S_yy,WA_S_xy)
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Stress
      
      implicit none
      integer,intent(in)::c_Key_Ave_Stress
      real(kind=FT),intent(in)::Point(2),Search_R,a_Weight,
     &                             l_Ordina
      real(kind=FT),intent(out)::WA_S_xx,WA_S_yy,WA_S_xy
      integer Tip_Elem,num_Surround_Ele,c_E
      integer i_E,i_G,c_num_Gauss,c_Global_G
      real(kind=FT) x_GP,y_GP,Sxx_GP,Syy_GP,Sxy_GP,c_L
      real(kind=FT)  All_weight,weight,tem1
      integer num_G
      
      call Cal_Ele_Num_by_Coors(Point(1),Point(2),Tip_Elem)
      
      num_Surround_Ele = num_Ele_Eles(Tip_Elem)
      
      All_weight  = ZR
      WA_S_xx     = ZR
      WA_S_yy     = ZR
      WA_S_xy     = ZR
      tem1 = ONE/((TWO*pi)**(1.5D0))/(l_Ordina**3)
      num_G = 0
      do i_E = 1,num_Surround_Ele
          c_E = Ele_Elements(Tip_Elem,i_E)
          do i_G = 1,num_GP_Elem(c_E)
              c_Global_G = Ele_GP_Start_Num(c_E) + i_G -1
              x_GP = Gauss_CoorX(c_Global_G)
              y_GP = Gauss_CoorY(c_Global_G)
              Sxx_GP = Stress_xx_Gauss(c_Global_G)
              Syy_GP = Stress_yy_Gauss(c_Global_G)
              Sxy_GP = Stress_xy_Gauss(c_Global_G)
              c_L = sqrt((x_GP-Point(1))**2 + (y_GP-Point(2))**2)
              if (c_L <= Search_R) then
                  if(c_Key_Ave_Stress==1)then
                      weight = (ONE-((c_L/Search_R)**a_Weight))**3
                  elseif(c_Key_Ave_Stress==2)then
                      weight = tem1*exp(-c_L**2/(TWO*l_Ordina**2))
                  endif
                  All_weight = All_weight + weight
                  num_G = num_G + 1
                  WA_S_xx     = WA_S_xx + weight * Sxx_GP
                  WA_S_yy     = WA_S_yy + weight * Syy_GP
                  WA_S_xy     = WA_S_xy + weight * Sxy_GP
              endif
          enddo
      enddo
      WA_S_xx = WA_S_xx /All_weight
      WA_S_yy = WA_S_yy /All_weight
      WA_S_xy = WA_S_xy /All_weight
      RETURN
      END SUBROUTINE Cal_Weighted_Ave_Stress_of_Point