 
      subroutine Tool_Denoise_Data(Data_Points,num_points,
     &                            Denoise_method, n_Sigma,
     &                            logical_Close)   
     

      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points,Denoise_method
      integer,intent(in)::n_Sigma
      logical,intent(in)::logical_Close
      real(kind=FT),intent(inout)::Data_Points(num_points)
      real(kind=FT) in_Data(num_points)
      integer i
      real(kind=FT) c_average,SD,c_x,low_limit,up_limit
      integer c_count
          
      in_Data = Data_Points

      
      if (Denoise_method==1) then
          c_average = sum(in_Data)/dble(num_points)
          call Tool_Standard_Deviation(in_Data,num_points,SD)
          
          
          c_count = 0 
          do i=1,num_points
             c_x = in_Data(i)
             low_limit = c_average - n_Sigma*SD
             up_limit  = c_average + n_Sigma*SD 
             if((c_x < low_limit) .or. (c_x > up_limit)) then
                 if (c_x > c_average ) then 
                     Data_Points(i) = c_average + SD
                 else
                     Data_Points(i) = c_average - SD
                 endif
                 c_count =  c_count +1
             endif
          enddo  
          
      endif
      
      return 
      end SUBROUTINE Tool_Denoise_Data             
