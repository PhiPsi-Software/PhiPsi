 
      subroutine Cal_Point_Aperture(c_Crack,Point,
     &                              c_DISP,Cr_Omega,Aperture)
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::Point(2),c_DISP(Total_FD),Cr_Omega
      integer,intent(in)::c_Crack
      real(kind=FT),intent(out)::Aperture
      real(kind=FT) Up_Point(2),Low_Point(2),Cr_p1(2),Cr_p2(2),
     &                 Up_Disp(2),Low_Disp(2),delta_L,
     &                 Delta_Disp(2),cos_angle,x1,y1,x2,y2
      real(kind=FT) Vector_a_b(2),Vector_ap_bp(2),sign_Aper,
     &                 Magni_Aper,Vector_n_Up(2),Vector_n_Low(2),
     &                 Up_Disp_Pro(2),Low_Disp_Pro(2),Disp_Normal(2)
      real(kind=FT) c_Kesi,c_Yita
      real(kind=FT) c_Kesi_Up,c_Yita_Up
      real(kind=FT) c_Kesi_Low,c_Yita_Low
      integer c_Elem_Up,c_Elem_Low
      integer Elem_Before_Offset
      real(kind=FT) Elem_Before_Off_Center(2)
      
      delta_L = Delta_Factor_Aper*Ave_Elem_L_Enrich
      
      Up_Point(1) = Point(1)  -  delta_L*sin(Cr_Omega)
      Up_Point(2) = Point(2)  +  delta_L*cos(Cr_Omega)
      
      
      Low_Point(1) = Point(1)  +  delta_L*sin(Cr_Omega)
      Low_Point(2) = Point(2)  -  delta_L*cos(Cr_Omega) 
      
      
      call Cal_Ele_Num_by_Coors(Point(1),Point(2),Elem_Before_Offset)
      Elem_Before_Off_Center(1:2)=Elem_Centroid(Elem_Before_Offset,1:2)
      
      
      call Cal_Ele_Num_by_Coors(Up_Point(1),Up_Point(2),c_Elem_Up)
      if (c_Elem_Up<=0) then
          
          Up_Point(1) = Elem_Before_Off_Center(1) -delta_L*sin(Cr_Omega)
          Up_Point(2) = Elem_Before_Off_Center(2) +delta_L*cos(Cr_Omega)
          call Cal_Ele_Num_by_Coors(Up_Point(1),Up_Point(2),c_Elem_Up)
          
          if (c_Elem_Up<=0) then
              Aperture = -TEN_15
              print *, '    Warning:  illegal element number found'
     &               //' in Cal_Point_Aperture.f!'
              return 
          endif
      endif
      
      
      call Cal_KesiYita_by_Coor(Up_Point,c_Elem_Up,c_Kesi_Up,c_Yita_Up)
      
      call Cal_Any_Point_Disp_KesiYita(c_Elem_Up,c_Kesi_Up,
     &                                 c_Yita_Up,0,c_DISP,Up_Disp)
      
      call Cal_Ele_Num_by_Coors(Low_Point(1),Low_Point(2),c_Elem_Low)
      if (c_Elem_Low<=0) then
          Low_Point(1) =Elem_Before_Off_Center(1) +delta_L*sin(Cr_Omega)
          Low_Point(2) =Elem_Before_Off_Center(2) -delta_L*cos(Cr_Omega) 
          call Cal_Ele_Num_by_Coors(Low_Point(1),Low_Point(2),
     &                              c_Elem_Low)
          
          if (c_Elem_Low<=0) then
              Aperture = -TEN_15
              print *, '    Warning:  illegal element number found'
     &               //' in Cal_Point_Aperture.f!'
              return 
          endif
      endif
      
      call Cal_KesiYita_by_Coor(Low_Point,c_Elem_Low,
     &                          c_Kesi_Low,c_Yita_Low)
      call Cal_Any_Point_Disp_KesiYita(c_Elem_Low,c_Kesi_Low,c_Yita_Low,
     &                                 0,c_DISP,Low_Disp)
      
      Vector_n_Up  = [-sin(Cr_Omega),cos(Cr_Omega)]
      Vector_n_Low = -Vector_n_Up
      Up_Disp_Pro(1:2) = dot_product(Up_Disp,Vector_n_Up)*Vector_n_Up
      Low_Disp_Pro(1:2)= dot_product(Low_Disp,Vector_n_Low)*Vector_n_Low
      Disp_Normal = Up_Disp_Pro - Low_Disp_Pro
      Magni_Aper = sqrt(Disp_Normal(1)**2+Disp_Normal(2)**2)
      x1 = Disp_Normal(1) 
      y1 = Disp_Normal(2) 
      x2 = Vector_n_Up(1) 
      y2 = Vector_n_Up(2) 
      cos_angle = (x1*x2+y1*y2)/(sqrt(x1**2+y1**2)*sqrt(x2**2+y2**2))
      if (cos_angle > 0) then
          sign_Aper = ONE
      else
          sign_Aper =-ONE
      end if

      Aperture = sign_Aper*Magni_Aper
      
      return 
      end SUBROUTINE Cal_Point_Aperture               
