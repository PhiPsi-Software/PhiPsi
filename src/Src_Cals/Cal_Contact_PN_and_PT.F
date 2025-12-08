 
      subroutine Cal_Contact_PN_and_PT(iter,ifra,Counter_Iter,i_NR_P,
     &                         c_Total_Freedom,c_num_freeDOF,
     &                         Kn,Kn_Gauss,Kt_Gauss,fric_mu,c_DISP,
     &                         delta_u_a,
     &                         CT_State_Gauss,c_Elem_Conta_Sta,
     &                         PC_Gauss_x,PC_Gauss_y)
      
      
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Material
      use Global_Field_Problem
      
      implicit none
      integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
      integer, intent(in)::c_Total_Freedom,c_num_freeDOF
      real(kind=FT),intent(in)::c_DISP(c_Total_Freedom)
      real(kind=FT),intent(in)::delta_u_a(c_Total_Freedom)
      real(kind=FT),intent(in)::Kn,fric_mu
      real(kind=FT),intent(inout)::
     &        Kn_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT),intent(inout)::
     &        Kt_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      integer,intent(inout)::
     &        CT_State_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)             
      integer,intent(inout):: 
     &        c_Elem_Conta_Sta(Num_Elem,num_Crack)
      real(kind=FT),intent(inout)::
     &              PC_Gauss_x(num_Crack,Max_Num_Cr_CalP-1,2),
     &              PC_Gauss_y(num_Crack,Max_Num_Cr_CalP-1,2)
      integer i_E,i_C
      integer Num_Div_Points,i_CT_Elem
      real(kind=FT) Aper_Tol
      real(kind=FT) x1,y1,x2,y2,ori_CT_elem
      integer i_CT_GAUSS,num_CT_Gauss
      real(kind=FT) kesi_CT(2),Weight_CT(2)
      real(kind=FT) CT_GAUSS_x,CT_GAUSS_y
      integer CT_GAUSS_Elem
      real(kind=FT) c_Aperture,c_Tangent_disp,ele_Sign
      real(kind=FT) c_PN,c_PT
      integer c_Adj_Ele
      real(kind=FT) c_Kt
      real(kind=FT) delta_Aperture,delta_Tangent_disp
      integer Old_Elem_Conta_Sta(Num_Elem,num_Crack)
      real(kind=FT) tem_Aper_Tol,Prop_C1,Prop_C2,P_closure,Prop_C0
      real(kind=FT) c_Kn,Curve_Slope
      integer c_Yes_Proppant
      real(kind=FT) CT_Length
 
      
      if(i_NR_P>=2)then
          Old_Elem_Conta_Sta = c_Elem_Conta_Sta
          CT_State_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2)  = 0
          c_Elem_Conta_Sta(1:Num_Elem,1:num_Crack) = 0
      endif
      PC_Gauss_x(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = ZR
      PC_Gauss_y(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = ZR
      
      if(Conta_Integ_Point==2)then
          num_CT_Gauss = 2
          kesi_CT(1)  = -0.577350269189626D0
          kesi_CT(2)  =  0.577350269189626D0
          Weight_CT(1)=  ONE
          Weight_CT(2)=  ONE
      elseif(Conta_Integ_Point==1)then
          num_CT_Gauss = 1
          kesi_CT(1)  =  0.0D0
          Weight_CT(1)=  TWO
      endif
      c_Yes_Proppant = 0
      
      do i_C=1,num_Crack
          Num_Div_Points = Cracks_CalP_Num(i_C)
          do i_CT_Elem = 1,Num_Div_Points - 1
            x1= Cracks_CalP_Coors(i_C,i_CT_Elem,1)
            y1= Cracks_CalP_Coors(i_C,i_CT_Elem,2)
            x2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,1)
            y2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,2)
              ori_CT_elem = atan2(y2-y1,x2-x1)
              CT_Length = Cracks_HF_Ele_L(i_C,i_CT_Elem)

              do i_CT_GAUSS  = 1,num_CT_Gauss
                  call Cal_HF_Coor_by_Kesi(kesi_CT(i_CT_GAUSS),
     &                                     x1,y1,x2,y2,
     &                                     CT_GAUSS_x,CT_GAUSS_y)
                  call Cal_Ele_Num_by_Coors(CT_GAUSS_x,CT_GAUSS_y,
     &                                      CT_GAUSS_Elem)

                  call Cal_Point_Aperture_and_Tangent_disp(i_C,
     &                             [CT_GAUSS_x,CT_GAUSS_y],
     &                              c_DISP,ori_CT_elem,c_Aperture,
     &                              c_Tangent_disp)   
                  
                  
                  if(Key_Propped_Width==0)then
                      Aper_Tol = ZR
                  endif
                  
#ifndef github
                  if(Key_Propped_Width==1)then
                      c_Yes_Proppant = 1
                      tem_Aper_Tol = (Cracks_CalP_wpnp(i_C,i_CT_Elem) +
     &                        Cracks_CalP_wpnp(i_C,i_CT_Elem+1))/TWO
                      Aper_Tol = tem_Aper_Tol
                      call Cal_C1_C2_for_Width_Propped_Frac(Aper_Tol,
     &                                               Prop_C1,Prop_C2)
                  endif
#endif
#ifndef github
                  if(Key_Analysis_Type==17 .and. 
     &               Key_Proppant_Active==1)  then
                      c_Yes_Proppant=1
                      tem_Aper_Tol = (Cracks_CalP_wpnp(i_C,i_CT_Elem) +
     &                        Cracks_CalP_wpnp(i_C,i_CT_Elem+1))/TWO
                      Aper_Tol = tem_Aper_Tol
                      call Cal_C1_C2_for_Width_Propped_Frac(Aper_Tol,
     &                                               Prop_C1,Prop_C2)
                  endif
#endif
                  if(c_Aperture <= Aper_Tol)then
                      c_Kt = Kt_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) 
                      c_Kn = Kn_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) 
                      ele_Sign =  sign(ONE,cos(ori_CT_elem))
                      if(c_Yes_Proppant==0)then
                          c_PN = Kn*(c_Aperture-Aper_Tol)   
                          c_PT = c_Kt*c_Tangent_disp*ele_Sign
                      elseif(c_Yes_Proppant==1)then
#ifndef github
                          if(Aper_Tol > 0.01D-3) then
                              Prop_C0 = c_Aperture - Aper_Tol
                              call Cal_Closure_P_Given_C1_C2_Wp_Secant(
     &                             Prop_C0,Prop_C1,Prop_C2,P_closure,
     &                             Curve_Slope)
                              
                              c_PN = -P_closure
                              c_PT = ZR
                              
                              Cracks_CalP_wdeform(i_C,i_CT_Elem) = 
     &                          1.04D0*tem_Aper_Tol*
     &                          (P_closure**(TWO/THR))*
     &                          ((ONE-Prop_Poss)/Prop_Elas)**(TWO/THR)
                              
                              Cracks_CalP_wdeform(i_C,i_CT_Elem+1) = 
     &                           Cracks_CalP_wdeform(i_C,i_CT_Elem)                         
                              c_Kn  = Curve_Slope
                              
                              
                              c_Kt  = ZR
                              
                          elseif(Aper_Tol <= 0.01D-3)then
                              c_PN = Kn*(c_Aperture)   
                              c_PT = c_Kt*c_Tangent_disp*ele_Sign
                          endif
#endif
                      endif
                      if(i_NR_P>=2)then
                          CT_State_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) = 1
                          if(old_Elem_Conta_Sta(CT_GAUSS_Elem,i_C)/=2)
     &                                                            then
                              c_Elem_Conta_Sta(CT_GAUSS_Elem,i_C) = 1
                          else
                              c_Elem_Conta_Sta(CT_GAUSS_Elem,i_C) = 2
                          endif
                      endif
                      if((Key_Propped_Width==0) .or.
     &                  (Key_Propped_Width==1 .and. 
     &                   Aper_Tol <= 0.01D-3)) then
                        if(abs(c_PT)-fric_mu*abs(c_PN) >ZR)then
                          call Cal_Point_Aperture_and_Tangent_disp(i_C,
     &                             [CT_GAUSS_x,CT_GAUSS_y],
     &                              delta_u_a,ori_CT_elem,
     &                              delta_Aperture,delta_Tangent_disp)    
                          c_PT = fric_mu*abs(c_PN)*
     &                       delta_Tangent_disp/abs(delta_Tangent_disp)
                          c_Kt  = c_PT/c_Tangent_disp
                          
                          
                          
                          if(i_NR_P>=2)then
                              CT_State_Gauss(i_C,i_CT_Elem,i_CT_GAUSS)=2
                              c_Elem_Conta_Sta(CT_GAUSS_Elem,i_C) = 2
                          endif
                        endif
                      endif
                      PC_Gauss_x(i_C,i_CT_Elem,i_CT_GAUSS) =
     &                       c_PT*cos(ori_CT_elem)-c_PN*sin(ori_CT_elem)
                      PC_Gauss_y(i_C,i_CT_Elem,i_CT_GAUSS) =
     &                       c_PT*sin(ori_CT_elem)+c_PN*cos(ori_CT_elem)
     
                      Kn_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) = c_Kn
                      Kt_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) = c_Kt
                  endif
              enddo
          enddo
      enddo
      
      do i_C=1,num_Crack
          Num_Div_Points = Cracks_CalP_Num(i_C)
          do i_CT_Elem = 1,Num_Div_Points - 1
              do i_CT_GAUSS  = 1,num_CT_Gauss         
                Cracks_CalP_Contact_Force_x(i_C,i_CT_Elem) = 
     &                         sum(PC_Gauss_x(i_C,i_CT_Elem,1:2))
                Cracks_CalP_Contact_Force_y(i_C,i_CT_Elem) = 
     &                         sum(PC_Gauss_y(i_C,i_CT_Elem,1:2))
              enddo
          enddo
      enddo
 
      
      if(i_NR_P>=2)then
        do i_E = 1,Num_Elem
          do i_C=1,num_Crack
              c_Adj_Ele = TipEle_Adjacent_Ele(i_E,i_C) 
              if (c_Adj_Ele/= 0) then
                  if(c_Elem_Conta_Sta(c_Adj_Ele,i_C) /= 0 )then 
                     c_Elem_Conta_Sta(i_E,i_C)=
     &                    c_Elem_Conta_Sta(c_Adj_Ele,i_C) 
                  endif
              end if
          enddo
        end do
      endif
      
      return 
      end SUBROUTINE Cal_Contact_PN_and_PT         
