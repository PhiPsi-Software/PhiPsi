 
      subroutine Cal_Crack_Aperture(isub,c_DISP)
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Common
      use Global_HF
      
      implicit none
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      integer,intent(in)::isub
      integer i_C,i_P
      integer Num_CalP
      real(kind=FT) CalP_Points(Max_Num_Cr_CalP,2)   
      real(kind=FT) P_Aperture,Cr_Omega,L_Aperture,R_Aperture
      real(kind=FT) L_Coor(2),R_Coor(2),L_Omega,R_Omega
      real(kind=FT) Cr_p1(2),Cr_p2(2),L_Length,R_Length
      real(kind=FT) c_Line_AB(2,2),delta_L,shorted_Line_AB(2,2)
      Cracks_CalP_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
      
      do i_C = 1,num_Crack
          Num_CalP = Cracks_CalP_Num(i_C)
          CalP_Points(1:Num_CalP,1:2) = 
     &               Cracks_CalP_Coors(i_C,1:Num_CalP,1:2)
          
          do i_P =1,Num_CalP
              if(i_P==1 .or. i_P==Num_CalP) then
                  Cr_Omega = Cracks_CalP_Orient(i_C,i_P)
                  
                  call Cal_Point_Aperture(i_C,CalP_Points(i_P,1:2),
     &                                    c_DISP,Cr_Omega,P_Aperture)
     
                  Cracks_CalP_Aper(i_C,i_P) = P_Aperture
              else
                  Cr_p1 = [Cracks_CalP_Coors(i_C,i_P-1,1),
     &                     Cracks_CalP_Coors(i_C,i_P-1,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,i_P,1),
     &                     Cracks_CalP_Coors(i_C,i_P,2)]
                  L_Length = sqrt((Cr_p2(2)-Cr_p1(2))**2 + 
     &                            (Cr_p2(1)-Cr_p1(1))**2)
                  delta_L = ZPZ1*Ave_Elem_L_Enrich
                  if (delta_L >= L_Length) then
                      delta_L = L_Length*HLF
                  endif
                  c_Line_AB(1,1:2) = Cr_p1; c_Line_AB(2,1:2) = Cr_p2
                  call Tool_Shorten_or_Extend_Line
     &                              (c_Line_AB, -delta_L, 'B',
     &                              shorted_Line_AB,L_Coor)
                  L_Omega= atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
                  Cr_p1 = [Cracks_CalP_Coors(i_C,i_P,1),
     &                     Cracks_CalP_Coors(i_C,i_P,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,i_P+1,1),
     &                     Cracks_CalP_Coors(i_C,i_P+1,2)]
                  R_Length = sqrt((Cr_p2(2)-Cr_p1(2))**2 + 
     &                            (Cr_p2(1)-Cr_p1(1))**2)
                  delta_L = ZPZ1*Ave_Elem_L_Enrich
                  if (delta_L >= R_Length) then
                      delta_L = R_Length*HLF
                  endif
                  c_Line_AB(1,1:2) = Cr_p1; c_Line_AB(2,1:2) = Cr_p2
                  call Tool_Shorten_or_Extend_Line
     &                              (c_Line_AB, -delta_L, 'A',
     &                              shorted_Line_AB,R_Coor)
                  R_Omega= atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
                  call Cal_Point_Aperture(i_C,L_Coor,
     &                                    c_DISP,L_Omega,L_Aperture)
                  call Cal_Point_Aperture(i_C,R_Coor,
     &                                    c_DISP,R_Omega,R_Aperture)
                  Cracks_CalP_Aper(i_C,i_P) = HLF*(L_Aperture+
     &                                               R_Aperture)
              endif
          end do
      end do
      
      if(Key_Min_Aperture==1 .and. Key_Analysis_Type==3)then
          do i_C = 1,num_Crack
              do i_P =1,Cracks_CalP_Num(i_C)
                  if(Cracks_CalP_Aper(i_C,i_P) <= Min_Aperture)then
                      Cracks_CalP_Aper(i_C,i_P) = Min_Aperture
                  endif
              enddo
          enddo
      endif
      
      return 
      end SUBROUTINE Cal_Crack_Aperture                
