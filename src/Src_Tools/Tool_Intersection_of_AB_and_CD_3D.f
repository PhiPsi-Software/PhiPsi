 
      subroutine Tool_Intersection_of_AB_and_CD_3D(A,B,C,D,
     &                 Yes_Inter,InterSection_P) 

 
      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::A(3),B(3),C(3),D(3)
      real(kind=FT),intent(out)::InterSection_P(3)
      logical,intent(out)::Yes_Inter
      
      real(kind=FT) A1,A2,A3,B1,B2,B3
      real(kind=FT) C1,C2,C3,D1,D2,D3
      real(kind=FT) t,s
      real(kind=FT) K(2,2),F(2),U(2),det_K
      real(kind=FT) check_value
      logical Yes_on_AB,Yes_on_CD
      
      real(kind=FT) max_x_AB,max_x_CD,min_x_AB,min_x_CD
      real(kind=FT) max_y_AB,max_y_CD,min_y_AB,min_y_CD
      real(kind=FT) max_z_AB,max_z_CD,min_z_AB,min_z_CD
      
      logical Logical_Yes_x,Logical_Yes_y,Logical_Yes_z
      
      Yes_Inter = .False.
      InterSection_P(1:3)=ZR 
      

      max_x_AB= max(A(1),B(1))
      min_x_AB= min(A(1),B(1))
      max_y_AB= max(A(2),B(2))
      min_y_AB= min(A(2),B(2))
      max_z_AB= max(A(3),B(3))
      min_z_AB= min(A(3),B(3))
      
      max_x_CD= max(C(1),D(1))
      min_x_CD= min(C(1),D(1))
      max_y_CD= max(C(2),D(2))
      min_y_CD= min(C(2),D(2))
      max_z_CD= max(C(3),D(3))
      min_z_CD= min(C(3),D(3))     
      
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_x_AB,max_x_AB],
     &                            [min_x_CD,max_x_CD],Logical_Yes_x) 
      if(Logical_Yes_x .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_y_AB,max_y_AB],
     &                            [min_y_CD,max_y_CD],Logical_Yes_y) 
      if(Logical_Yes_y .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_z_AB,max_z_AB],
     &                            [min_z_CD,max_z_CD],Logical_Yes_z) 
      if(Logical_Yes_z .eqv. .False.)then
          return
      endif    
      
      a1 = A(1);a2 = A(2);a3 = A(3)
      b1 = B(1);b2 = B(2);b3 = B(3)
      c1 = C(1);c2 = C(2);c3 = C(3)
      d1 = D(1);d2 = D(2);d3 = D(3)


      K(1,1) = b1-a1
      K(1,2) = c1-d1
      K(2,1) = b2-a2
      K(2,2) = c2-d2  
      F(1)   = c1-a1
      F(2)   = c2-a2
      det_K  = K(1,1)*K(2,2) - K(1,2)*K(2,1) 
      if (det_K /= ZR) then
          call Matrix_Solve_LSOE(-1,0,2,K,F,U,2)
          t = U(1)
          s = U(2)
          check_value = abs(a3+t*(b3-a3)-c3-s*(d3-c3))
          if(check_value<=Tol_11) then
              InterSection_P(1) = a1+t*(b1-a1)
              InterSection_P(2) = a2+t*(b2-a2)
              InterSection_P(3) = a3+t*(b3-a3)
 
              
              call Tool_Yes_Point_on_Line_Segment_3D(A,B,
     &                     InterSection_P,Yes_on_AB)
              call Tool_Yes_Point_on_Line_Segment_3D(C,D,
     &                     InterSection_P,Yes_on_CD)
              if(Yes_on_AB .and. Yes_on_CD)then
                  
                  Yes_Inter = .True.
                  return
              endif
          endif
      endif
      
      K(1,1) = b2-a2
      K(1,2) = c2-d2
      K(2,1) = b3-a3
      K(2,2) = c3-d3  
      F(1)   = c2-a2
      F(2)   = c3-a3
      det_K  = K(1,1)*K(2,2) - K(1,2)*K(2,1) 
      if (det_K /= ZR) then
          call Matrix_Solve_LSOE(-1,0,2,K,F,U,2)
          t = U(1)
          s = U(2)
          check_value = abs(a1+t*(b1-a1)-c1-s*(d1-c1))
          if(check_value<=Tol_11) then
              InterSection_P(1) = a1+t*(b1-a1)
              InterSection_P(2) = a2+t*(b2-a2)
              InterSection_P(3) = a3+t*(b3-a3)

              
              call Tool_Yes_Point_on_Line_Segment_3D(A,B,
     &                     InterSection_P,Yes_on_AB)
              call Tool_Yes_Point_on_Line_Segment_3D(C,D,
     &                     InterSection_P,Yes_on_CD)
              if(Yes_on_AB .and. Yes_on_CD)then
                  
                  Yes_Inter = .True.
                  return
              endif
          endif
      endif    
      
      K(1,1) = b1-a1
      K(1,2) = c1-d1
      K(2,1) = b3-a3
      K(2,2) = c3-d3  
      F(1)   = c1-a1
      F(2)   = c3-a3
      det_K  = K(1,1)*K(2,2) - K(1,2)*K(2,1) 
      if (det_K /= ZR) then
          call Matrix_Solve_LSOE(-1,0,2,K,F,U,2)
          t = U(1)
          s = U(2)
          check_value = abs(a2+t*(b2-a2)-c2-s*(d2-c2))
          if(check_value<=Tol_11) then
              InterSection_P(1) = a1+t*(b1-a1)
              InterSection_P(2) = a2+t*(b2-a2)
              InterSection_P(3) = a3+t*(b3-a3)

              
              call Tool_Yes_Point_on_Line_Segment_3D(A,B,
     &                     InterSection_P,Yes_on_AB)
              call Tool_Yes_Point_on_Line_Segment_3D(C,D,
     &                     InterSection_P,Yes_on_CD)
              if(Yes_on_AB .and. Yes_on_CD)then
                  
                  Yes_Inter = .True.
                  return
              endif
          endif
      endif        
  
      return 
      end SUBROUTINE Tool_Intersection_of_AB_and_CD_3D                        
