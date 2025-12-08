 
      subroutine Tool_Dis_Point_to_3D_Tri(Point,
     &                            Tri_P1,Tri_P2,Tri_P3,
     &                            Distance,PER,Yes_PER_in,Yes_PER_on)

      use Global_Float_Type    
      use Global_Common
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                          Tri_P1(3),Tri_P2(3),Tri_P3(3)
      real(kind=FT),intent(out):: Distance,PER(3)
      logical,intent(out):: Yes_PER_in,Yes_PER_on
      real(kind=FT) P1P2(3),P1P3(3),P1P0(3),Np(3),cosa,P0P01(3)
      real(kind=FT) abs_P1P0,abs_Np,abs_P0P01,P01(3)
      real(kind=FT) Distance_Sign
      real(kind=FT) B_i(3),C_i(3),X_i(3)
      real(kind=FT) c_DET
      integer c_Orient
      real(kind=FT) v0(3),v1(3),v2(3)
      real(kind=FT) dot00,dot01,dot02,dot11,dot12,inverDeno
      real(kind=FT) u,v
      logical c_Yes_on
      
      Yes_PER_in = .False.
      Yes_PER_on = .False.
      Distance   =  ZR
      
      P1P2 = Tri_P2 - Tri_P1
      P1P3 = Tri_P3 - Tri_P1
      Np(1) = P1P2(2) * P1P3(3) - P1P2(3) * P1P3(2)
      Np(2) = P1P2(3) * P1P3(1) - P1P2(1) * P1P3(3)
      Np(3) = P1P2(1) * P1P3(2) - P1P2(2) * P1P3(1)
      
      P1P0     = Point - Tri_P1
      
      abs_P1P0 = sqrt(P1P0(1)**2 + P1P0(2)**2 + P1P0(3)**2)
      
      if(abs_P1P0<=Tol_10) abs_P1P0=Tol_10
      
      
      abs_Np   = sqrt(Np(1)**2   + Np(2)**2   + Np(3)**2)
      
      cosa     = dot_product(P1P0,Np)/abs_P1P0/abs_Np 
      
      abs_P0P01 = abs_P1P0 * cosa
      
      P0P01 = -abs_P0P01*Np/abs_Np
      
      P01 = Point + P0P01
      PER = P01

      
      B_i = Tri_P2-Tri_P1
      C_i = Tri_P3-Tri_P1
      X_i = Point-Tri_P1
      
      
      
      
      c_DET = 
     &      B_i(1)*C_i(2)*X_i(3) +
     &      B_i(2)*C_i(3)*X_i(1) +
     &      C_i(1)*X_i(2)*B_i(3) -
     &      B_i(1)*C_i(3)*X_i(2) -
     &      B_i(2)*C_i(1)*X_i(3) -
     &      B_i(3)*C_i(2)*X_i(1)

     
     
      if(c_DET>ZR)then
          c_Orient =  1
      elseif(c_DET<ZR)then
          c_Orient = -1
      else
          c_Orient =  0
      endif
      
      Distance_Sign = dble(c_Orient)
      
      Distance = Distance_Sign*abs(abs_P0P01)
      
#ifndef Silverfrost
      if (isnan(Distance)) then
           print *, '    Error-2023061502 :: Distance is NAN in '//
     & 'Tool_Dis_Point_to_3D_Tri.f90!'
          call Warning_Message('S',Keywords_Blank)
      endif
#endif      
      
      
      v0 = Tri_P3 - Tri_P1
      v1 = Tri_P2 - Tri_P1
      v2 = PER    - Tri_P1
      
      dot00 = dot_product(v0,v0)
      dot01 = dot_product(v0,v1)
      dot02 = dot_product(v0,v2)
      dot11 = dot_product(v1,v1)
      dot12 = dot_product(v1,v2)   
      
      inverDeno = ONE / (dot00 * dot11 - dot01 * dot01)
      u = (dot11 * dot02 - dot01 * dot12) * inverDeno
      v = (dot00 * dot12 - dot01 * dot02) * inverDeno
      
      if((u >= ZR) .and. (v >= ZR) .and. (u + v <= ONE))then
          Yes_PER_in = .True.
          return
      endif
      
      call Tool_Yes_Point_on_Line_Segment_3D(Tri_P1,Tri_P2,PER,c_Yes_on)
      if (c_Yes_on) then
          Yes_PER_on = .True.
          return
      endif
      call Tool_Yes_Point_on_Line_Segment_3D(Tri_P2,Tri_P3,PER,c_Yes_on)
      if (c_Yes_on) then
          Yes_PER_on = .True.
          return
      endif      
      call Tool_Yes_Point_on_Line_Segment_3D(Tri_P1,Tri_P3,PER,c_Yes_on)
      if (c_Yes_on) then
          Yes_PER_on = .True.
          return
      endif  
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_3D_Tri                   
