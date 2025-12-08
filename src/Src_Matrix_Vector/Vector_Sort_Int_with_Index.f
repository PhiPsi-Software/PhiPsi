 
      SUBROUTINE Vector_Sort_Int_with_Index(n,Vector,V_Index)   

      use Global_Float_Type     
      implicit none
      integer,intent(in)::n
      integer,intent(inout):: Vector(n)
      integer,intent(out):: V_Index(n)
      integer Vector_Old(n)
      integer i,j,a
      logical c_Yes_In
      
      Vector_Old = Vector
      
      do j=2, n
          a=Vector(j)
          do i=j-1,1,-1
              if (Vector(i).le.a) goto 10
              Vector(i+1)=Vector(i)
              
              
          end do
          i=0
          
   10     Vector(i+1)=a
   
      end do
      
      
      do i=1,n
          call Vector_Location_Int(n,Vector_Old,Vector(i),
     &                             V_Index(i),c_Yes_In)   
      enddo
      
      
      return
      END SUBROUTINE Vector_Sort_Int_with_Index
    


