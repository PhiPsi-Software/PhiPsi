 
      function tool_get_mac_address()  result(line) 
          character(80) :: line 
          call system('net config workstation >mac.dat ' ) 
          open (158,file='mac.dat') 
          read (158,*) 
          read (158,*) 
          read (158,*) 
          read (158,*) 
          read (158,*) 
          read (158,*) 
          read (158,*) line 
          close (158,status='delete') 

      end function                       
