 
      subroutine Tool_Intersection_of_AB_and_Brick_Ele(A,B,c_El,
     &                 Yes_Inter,InterSection_P) 

 
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Model
      
      implicit none
      real(kind=FT),intent(in)::A(3),B(3)
      integer,intent(in)::c_El
      real(kind=FT),intent(out)::InterSection_P(3)
      logical,intent(out)::Yes_Inter
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
      integer Triangle(12,3)
      real(kind=FT) Point1(3),Point2(3),Point3(3)
      logical c_Yes_Inter
      real(kind=FT) c_InterSection_P(3)
      integer i_Tri
      
      Yes_Inter = .False.
      InterSection_P(1:3)=ZR 
      
      c_X_NODES = G_X_NODES(1:8,c_El)
      c_Y_NODES = G_Y_NODES(1:8,c_El)  
      c_Z_NODES = G_Z_NODES(1:8,c_El)  
      
      Triangle(1,1:3) = [1,2,5]
      Triangle(2,1:3) = [2,6,5]
      Triangle(3,1:3) = [2,3,6]
      Triangle(4,1:3) = [3,7,6]
      Triangle(5,1:3) = [3,7,8]
      Triangle(6,1:3) = [3,8,4]
      Triangle(7,1:3) = [4,8,5]
      Triangle(8,1:3) = [4,5,1]
      Triangle(9,1:3) = [5,6,8]
      Triangle(10,1:3) = [6,7,8]
      Triangle(11,1:3) = [1,2,4]
      Triangle(12,1:3) = [2,3,4]    
      
      do i_Tri = 1,12
          Point1 =  [c_X_NODES(Triangle(i_Tri,1)),
     &               c_Y_NODES(Triangle(i_Tri,1)),
     &               c_Z_NODES(Triangle(i_Tri,1))]
          Point2 =  [c_X_NODES(Triangle(i_Tri,2)),
     &               c_Y_NODES(Triangle(i_Tri,2)),
     &               c_Z_NODES(Triangle(i_Tri,2))]
          Point3 =  [c_X_NODES(Triangle(i_Tri,3)),
     &               c_Y_NODES(Triangle(i_Tri,3)),
     &               c_Z_NODES(Triangle(i_Tri,3))] 
          c_Yes_Inter = .False.
          call Tool_Intersection_of_AB_and_Triangle_3D(A,B,
     &                 Point1,Point2,Point3,
     &                 c_Yes_Inter,c_InterSection_P) 
          if (c_Yes_Inter .eqv. .True.) then
              Yes_Inter = .True.
              InterSection_P =c_InterSection_P
              return
          endif
      enddo
      
      return 
      end SUBROUTINE Tool_Intersection_of_AB_and_Brick_Ele                         
