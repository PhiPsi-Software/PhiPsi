 
      SUBROUTINE Cal_Contact_Contact_State_Gauss(
     &                iter,ifra,Counter_Iter,i_NR_P,
     &                c_DISP,Yes_Contact,c_Elem_Conta_Sta,
     &                CT_State_Gauss) 



      use Global_Float_Type      
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Common
      use Global_Contact
      use Global_HF
      use Global_Elem_Area_Vol
      
      implicit none
      integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      logical,intent(out)::Yes_Contact
      integer,intent(out)::CT_State_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      integer,intent(out):: c_Elem_Conta_Sta(Num_Elem,num_Crack)
      
      integer i_E,i_C
      integer Num_Div_Points,i_CT_Elem
      real(kind=FT) Aper_Tol
      real(kind=FT) x1,y1,x2,y2,ori_n(2),ori_CT_elem
      integer i_CT_GAUSS,num_CT_Gauss
      real(kind=FT) kesi_CT(2),Weight_CT(2)
      real(kind=FT) CT_GAUSS_x,CT_GAUSS_y
      integer CT_GAUSS_Elem
      real(kind=FT) c_Aperture
      integer c_Adj_Ele
      
      Aper_Tol = 1.0D-5*Ave_Elem_L_Enrich
      
      CT_State_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2)  = 0
      c_Elem_Conta_Sta(1:Num_Elem,1:num_Crack) = 0
      
      Yes_Contact = .False.
      
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
      
      do i_C=1,num_Crack
          Num_Div_Points = Cracks_CalP_Num(i_C)
          do i_CT_Elem = 1,Num_Div_Points - 1
            x1= Cracks_CalP_Coors(i_C,i_CT_Elem,1)
            y1= Cracks_CalP_Coors(i_C,i_CT_Elem,2)
            x2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,1)
            y2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,2)
              ori_CT_elem = atan2(y2-y1,x2-x1)
              ori_n = [cos(pi/TWO + ori_CT_elem),
     &                 sin(pi/TWO + ori_CT_elem)]
              do i_CT_GAUSS  = 1,num_CT_Gauss
                  call Cal_HF_Coor_by_Kesi(kesi_CT(i_CT_GAUSS),
     &                                     x1,y1,x2,y2,
     &                                     CT_GAUSS_x,CT_GAUSS_y)
                  call Cal_Ele_Num_by_Coors(CT_GAUSS_x,CT_GAUSS_y,
     &                                      CT_GAUSS_Elem)

                  call Cal_Point_Aperture(i_C,[CT_GAUSS_x,CT_GAUSS_y],
     &                              c_DISP,ori_CT_elem,c_Aperture)

                  if ((Key_Analysis_Type == 3) .or.
     &                (Key_Analysis_Type == 4) .or.
     &                (Key_Analysis_Type == 5))  then
                      if (c_Aperture<=Aper_Tol) then     
                          CT_State_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) = 1
                          c_Elem_Conta_Sta(CT_GAUSS_Elem,i_C) = 1
                          Yes_Contact = .True.
                      end if
                  else   
                      if (c_Aperture<Aper_Tol) then 
                          CT_State_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) = 1
                          c_Elem_Conta_Sta(CT_GAUSS_Elem,i_C) = 1
                          Yes_Contact = .True.
                      end if                     
                  end if
              enddo
          enddo
      enddo
      
      
      do i_E = 1,Num_Elem
          do i_C=1,num_Crack
              c_Adj_Ele = TipEle_Adjacent_Ele(i_E,i_C) 
              if (c_Adj_Ele/= 0) then
                  if(c_Elem_Conta_Sta(c_Adj_Ele,i_C) == 1 )then 
                      c_Elem_Conta_Sta(i_E,i_C) = 1 
                  endif
              end if
          enddo
      end do
      
      
      
      RETURN
      END SUBROUTINE Cal_Contact_Contact_State_Gauss
