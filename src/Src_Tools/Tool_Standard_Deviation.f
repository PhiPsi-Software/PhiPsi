 
      subroutine Tool_Standard_Deviation(Data_Points,num_points,SD)

      use Global_Float_Type
      implicit none
      integer,intent(in)::num_points
      real(kind=FT),intent(in)::Data_Points(num_points)
      real(kind=FT),intent(out)::SD
      integer i
      real(kind=FT)  c_average,tem 
      
      c_average = sum(Data_Points)/dble(num_points)
      
      tem  = ZR
      
      do i=1,num_points
          tem = tem +(Data_Points(i)-c_average)**2
      enddo
      
      tem = tem/dble(num_points)
      
      if(tem<Tol_11)then
          tem=Tol_11
      endif
      
      SD = sqrt(tem)
      
      return 
      end SUBROUTINE Tool_Standard_Deviation          
