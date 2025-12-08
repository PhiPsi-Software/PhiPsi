 
      subroutine Tool_Smooth_Data(Data_Points,num_points,Smooth_method,
     &                            Moving_Average_Method_n,logical_Close)

      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points,Smooth_method
      integer,intent(in)::Moving_Average_Method_n
      logical,intent(in)::logical_Close
      real(kind=FT),intent(inout)::Data_Points(num_points)
      real(kind=FT) in_Data(num_points),out_Data(num_points)
      integer n
      integer i,j,k
      integer twice_n1
      real(kind=FT),ALLOCATABLE::tem_Data(:)
      real(kind=FT) c_average
      
      in_Data = Data_Points
      n       = Moving_Average_Method_n
      
      
      twice_n1 = 2*n+1
      
      if(twice_n1 >= num_points) then
        print *,'    Error :: n is too big or num_points is too small!'
        print *,'             in Tool_Smooth_Data.f!'
        print *,'             n:',n
        print *,'             2*n+1:',2*n+1
        print *,'             num_points:',num_points
        call Warning_Message('S',Keywords_Blank)   
      endif
      
      if (Smooth_method==1) then
          ALLOCATE(tem_Data(-num_points:2*num_points))
          if (logical_Close .eqv. .False.) then
              tem_Data(-num_points:0) = ZR
              tem_Data(1:num_points)  = in_Data(1:num_points)
              tem_Data((num_points+1):2*num_points) = ZR
          elseif (logical_Close .eqv. .True.) then
              tem_Data(-num_points+1:0)  = in_Data(1:num_points)
              tem_Data(1:num_points)     = in_Data(1:num_points)
              tem_Data((num_points+1):2*num_points) =  
     &                                   in_Data(1:num_points)
          endif
          do i=1,num_points
             c_average   = sum(tem_Data((i-n):(i+n)))/dble(twice_n1)
             out_Data(i) = c_average 
          enddo
          DEALLOCATE(tem_Data)
      endif
      
      Data_Points = out_Data
      
      return 
      end SUBROUTINE Tool_Smooth_Data             
