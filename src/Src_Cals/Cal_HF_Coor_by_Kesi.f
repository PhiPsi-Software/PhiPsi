 
      subroutine Cal_HF_Coor_by_Kesi(kesi,x1,y1,x2,y2,
     &                                Out_x,Out_y)
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::kesi,x1,y1,x2,y2
      real(kind=FT),intent(out)::Out_x,Out_y
      real(kind=FT) x0,y0
      real(kind=FT) angle,length
      
      x0     = HLF*(x1+x2)
      y0     = HLF*(y1+y2)
      angle  = atan2(y2-y1,x2-x1)
      length = sqrt((x2-x1)**2 + (y2-y1)**2)
      
      Out_x  = x0  + HLF*length*kesi*cos(angle)
      Out_y  = y0  + HLF*length*kesi*sin(angle)
     
      return 
      end SUBROUTINE Cal_HF_Coor_by_Kesi                         
