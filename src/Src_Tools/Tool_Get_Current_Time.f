 
      subroutine Tool_Get_Current_Time(c_data,date_time,second)
      use Global_Float_Type 
      implicit none
      
      integer date_time(8)
      integer(kind=LIT) second
      character*10 b(3),c_data
      
      call date_and_time(b(1), b(2), b(3), date_time)
      
      second = date_time(7)          +
     &         date_time(6)*60       +
     &         date_time(5)*60*60    +
     &         date_time(3)*60*60*24 +
     &         date_time(2)*60*60*24*30
      
      
      c_data= b(1)

      return 
      end SUBROUTINE Tool_Get_Current_Time                          
