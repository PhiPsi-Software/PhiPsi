 
      subroutine Tool_HF_Pres_Transform(iter,CalP_Pres,
     &                                  temp_Pres)
      use Global_Float_Type      
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      
      integer,intent(in)::iter
      real(kind=FT),intent(in)::CalP_Pres(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::temp_Pres(Max_Num_Cr,
     &                                        Max_Num_Cr_CalP)
      
      integer i_C,i_CalP,num_Count
      
      temp_Pres = ZR
      
      num_Count = 0
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then         
              do i_CalP=1,Cracks_CalP_Num(i_C)
                  num_Count = num_Count + 1
                  temp_Pres(i_C,i_CalP) = CalP_Pres(num_Count)
              end do
          end if
      end do

      return 
      end SUBROUTINE Tool_HF_Pres_Transform            
