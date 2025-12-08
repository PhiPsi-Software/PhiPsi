 
      subroutine Tool_Intersection_of_AB_and_Triangle_3D(A,B,
     &                 Point1,Point2,Point3,
     &                 Yes_Inter,InterSection_P) 
 
      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::A(3),B(3),Point1(3),Point2(3),Point3(3)
      real(kind=FT),intent(out)::InterSection_P(3)
      logical,intent(out)::Yes_Inter
      real(kind=FT) Tool_Function_2Point_Dis_3D
      real(kind=FT) t0(3),u(3),v(3),n(3),vectorNorm3d
      real(kind=FT) dir(3),w0(3),value_a,value_b,value_pos
      real(kind=FT) uu,uv,vv,w(3),wu,wv,value_D
      real(kind=FT) value_t,value_s
      logical Yes_on_AB
      real(kind=FT) DIS_AB
      
      Yes_Inter = .False.
      InterSection_P(1:3)=ZR 

      DIS_AB = Tool_Function_2Point_Dis_3D(A,B)
      if (DIS_AB <=Tol_11)then
          return
      endif
      
      
      t0  = Point1(1:3)
      u   = Point2(1:3) - t0
      v   = Point3(1:3) - t0
 
      call Vector_Cross_Product_3(u,v,n)   
      call Vector_Norm2(3,n,vectorNorm3d)   
      if (vectorNorm3d < Tol_11) then
           return
      endif
 
      dir = B - A
 
      w0 = A - t0
      value_a = -DOT_PRODUCT(n, w0)
      value_b =  DOT_PRODUCT(n, dir)
  
      if (abs(value_b) < Tol_11)then
          return    
      endif
  
      value_pos = value_a / value_b
 
      InterSection_P = A + value_pos * dir
 
  
      uu  = DOT_PRODUCT(u, u);
      uv  = DOT_PRODUCT(u, v);
      vv  = DOT_PRODUCT(v, v);
 
      w   = InterSection_P - t0
      wu  = DOT_PRODUCT(w, u)
      wv  = DOT_PRODUCT(w, v)
  
      value_D = uv**2 - uu * vv
    
      value_s = (uv * wv - vv * wu) / value_D
      
      
      if (value_s < ZR .or. value_s > ONE) then
         return
      endif
    
      value_t = (uv * wu - uu * wv) / value_D

      if (value_t < -Tol_11  .or. (value_s + value_t) > ONE)then
         return
      endif
      
      
      
      call Tool_Yes_Point_on_Line_Segment_3D(A,B,InterSection_P,
     &                                       Yes_on_AB)
      if(Yes_on_AB .eqv. .False.) then
          return 
      endif
      
      Yes_Inter = .True.
      

      
      return 
      end SUBROUTINE Tool_Intersection_of_AB_and_Triangle_3D                         
