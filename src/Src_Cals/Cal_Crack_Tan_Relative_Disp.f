 
      subroutine Cal_Crack_Tan_Relative_Disp(isub,c_DISP,
     &                                       Cracks_CalP_Tan_Aper)
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Common
      use Global_HF
      
      implicit none
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      real(kind=FT),intent(out)::Cracks_CalP_Tan_Aper(num_Crack,
     &                                                Max_Num_Cr_CalP)
      integer,intent(in)::isub
      integer i_C,i_P
      integer Num_CalP
      real(kind=FT) CalP_Points(Max_Num_Cr_CalP,2)   
      real(kind=FT) P_Aperture,Cr_Omega
      real(kind=FT) c_Tangent_disp
      
      Cracks_CalP_Tan_Aper(1:num_Crack,1:Max_Num_Cr_CalP) = ZR
      
      do i_C = 1,num_Crack
          Num_CalP = Cracks_CalP_Num(i_C)
          
          CalP_Points(1:Num_CalP,1:2)  =
     &         Cracks_CalP_Coors(i_C,1:Num_CalP,1:2)  
          
          do i_P =1,Num_CalP
              Cr_Omega = Cracks_CalP_Orient(i_C,i_P)
              call Cal_Point_Aperture(i_C,CalP_Points(i_P,1:2),
     &                                c_DISP,Cr_Omega,P_Aperture)
                  call Cal_Point_Aperture_and_Tangent_disp(i_C,
     &                             CalP_Points(i_P,1:2),
     &                              c_DISP,Cr_Omega,P_Aperture,
     &                              c_Tangent_disp)  
              Cracks_CalP_Tan_Aper(i_C,i_P) = c_Tangent_disp
          enddo
          

      end do
      
      return
      end SUBROUTINE Cal_Crack_Tan_Relative_Disp               
