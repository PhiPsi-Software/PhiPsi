 
      subroutine Tool_Length_of_Crack(i_Crack,Length_Crack)
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      
      implicit none
      integer,intent(in)::i_Crack
      real(kind=FT),intent(out)::Length_Crack
      integer c_C,i_S
      real(kind=FT) crack_p1(2),crack_p2(2),delta_L
      
      if(i_Crack>num_Crack)then
          print *,'     Error:: wrong crack number for '
     &            // 'Tool_Length_of_Crack.f!'
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      Length_Crack = ZR
      c_C = i_Crack   
      
      do i_S = 1,Each_Cr_Poi_Num(c_C)-1
          crack_p1 = [Crack_Coor(c_C,i_S,1),Crack_Coor(c_C,i_S,2)]
          crack_p2 = [Crack_Coor(c_C,i_S+1,1),Crack_Coor(c_C,i_S+1,2)]
          delta_L  =  sqrt((crack_p2(2)-crack_p1(2))**2 +
     &                     (crack_p2(1)-crack_p1(1))**2)
          Length_Crack = Length_Crack +  delta_L
      end do
      
      return 
      end SUBROUTINE Tool_Length_of_Crack                          
