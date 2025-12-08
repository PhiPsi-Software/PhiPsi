 
      subroutine Tool_Get_Nearest_Point_3D(num_Point,Points,A,
     &                                     Point_Num) 
 
      use Global_Float_Type
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::Points(num_Point,3),A(3)       
      integer,intent(out)::Point_Num
      real(kind=FT) Dis_Vector(num_Point)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      integer i_P
      
      if(num_Point<=1)then
        print *, '    Error :: num_Point<=1!'
        print *, '             in Tool_Get_Nearest_Point_3D.f!'
        call Warning_Message('S',Keywords_Blank)      
      endif
      
      
      do i_P=1,num_Point
          Dis_Vector(i_P) = Tool_Function_2Point_Dis_3D(
     &                              A,Points(i_P,1:3))
      enddo
      
      Point_Num = MINLOC(Dis_Vector(1:num_Point),1)
      
      return 
      end SUBROUTINE Tool_Get_Nearest_Point_3D                     
