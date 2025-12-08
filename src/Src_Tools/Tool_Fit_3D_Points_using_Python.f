 
      subroutine Tool_Fit_3D_Points_using_Python(Key_Twice,
     &                                    num_Point,
     &                                    In_Points,Out_Points)

      use Global_Float_Type     
      use Global_Model
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point,Key_Twice
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      real(kind=FT) In_Points_1(num_Point,3),In_Points_2(num_Point,3)
      real(kind=FT) Out_Points_1(num_Point,3),Out_Points_2(num_Point,3)   
      real(kind=FT) New_Out_Points_2(num_Point,3)   
      integer mid_num
      
      if (Key_Twice==0) then
          call Fit_3D_Points_using_Python(num_Point,
     &                                    In_Points,Out_Points)  
          Out_Points(1,1:3)         = In_Points(1,1:3)
          Out_Points(num_Point,1:3) = In_Points(num_Point,1:3)     
      endif

      if (Key_Twice==1) then
          In_Points_1(1:num_Point,1:3)=In_Points(1:num_Point,1:3)
          call Fit_3D_Points_using_Python(num_Point,
     &                                    In_Points_1,Out_Points_1)  
     
          mid_num = int(num_Point/2)
          In_Points_2(1:mid_num,1:3)=In_Points(mid_num:1:-1,1:3)
          In_Points_2(mid_num+1:num_Point,1:3)=
     &          In_Points(num_Point:(mid_num-1):-1,1:3)
          call Fit_3D_Points_using_Python(num_Point,
     &                                    In_Points_2,Out_Points_2)  
          
          New_Out_Points_2(1:mid_num,1:3)=Out_Points_2(mid_num:1:-1,1:3)
          New_Out_Points_2(mid_num+1:num_Point,1:3)=
     &                       Out_Points_2(num_Point:(mid_num-1):-1,1:3)
     
          Out_Points = (Out_Points_1 + New_Out_Points_2)/TWO
      endif      

      
      return 
      end SUBROUTINE Tool_Fit_3D_Points_using_Python   
      
      
      subroutine Fit_3D_Points_using_Python(num_Point,
     &                                      In_Points,Out_Points)      
      use Global_Float_Type     
      use Global_Model
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      
      character(200) c_File_name_1,pathname
      character(200) c_File_name_2
      integer i_P,istat
      logical alive



      open(unit=123, iostat=istat, file=c_File_name_1, status='old')
      if (istat == 0) close(123, status='delete')
      open(unit=124, iostat=istat, file=c_File_name_2, status='old')
      if (istat == 0) close(124, status='delete')

      pathname =  trim(PhiPsi_Current_Directory)//string_connector  
     &          //'Python_Tools'//string_connector 
      c_File_name_1   =  trim(pathname)//'Tem_Input_Points.txt'
      open(101,file=c_File_name_1,status='unknown')  
      do i_P=1,num_Point
          write(101,'(3E20.12)') In_Points(i_P,1:3)
      end do
      close(101)
      
      if(Key_Smooth_Front==3) then
        call system 
     &     ('cd .\Python_Tools\ && Smooth_3D_Points_CSAPS_v1.0.py')
      elseif(Key_Smooth_Front==4) then
        call system 
     &     ('cd .\Python_Tools\ && Smooth_3D_Points_Scipy_v1.0.py')
      endif
      
      
      c_File_name_2   =  trim(pathname)//'Tem_Output_Points.txt'
      inquire(file=c_File_name_2, exist=alive)
      if(alive.eqv..False.) then
        print *,'    Error :: can not find Tem_Output_Points.txt!'
        print *,'           in Tool_Fit_3D_Points_using_Python_CSAPS.f!'
        call Warning_Message('S',Keywords_Blank)
      endif
      
      open(112,file=c_File_name_2,status='old')
      Read(112,*,IOSTAT=istat)
      close(112)
      
      if ( istat /= 0 ) then
        print *,'    Error :: blank Tem_Output_Points.txt!'
        print *,'           in Tool_Fit_3D_Points_using_Python_CSAPS.f!'
        call Warning_Message('S',Keywords_Blank)
      else
          open(122,file=c_File_name_2,status='old')  
          do i_P=1,num_Point
              read(122,*) Out_Points(i_P,1:3)
          end do
          close(122)
      end if
      
      return 
      end SUBROUTINE Fit_3D_Points_using_Python         
