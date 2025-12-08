 
      SUBROUTINE Vector_Sort_Dou(n,Vector)
      use Global_Float_Type 
      implicit none
      integer,intent(in)::n
      real(kind=FT),intent(inout)::Vector(n)
      integer i,j
      real(kind=FT) a

      do j=2, n
          a=Vector(j)
          do i=j-1,1,-1
              if (Vector(i)<=a) goto 10
              Vector(i+1)=Vector(i)
          end do
          i=0
   10     Vector(i+1)=a
      end do
      
      return

      END SUBROUTINE Vector_Sort_Dou