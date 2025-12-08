 
      subroutine Cal_F_dFdx_dFdy(r,theta,omega,
     &                                 c_mat_type,
     &                                 F,dFdx,dFdy)
     
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Common
      implicit none
      
      real(kind=FT),intent(in)::r,theta,omega
      integer,intent(in)::c_mat_type
      real(kind=FT),intent(out)::dFdx(4),dFdy(4),F(4)
      
      real(kind=FT) fac,st2,ct2,s3t2,c3t2,st,ct,
     &                 r2,c_theta,dF1dx1,dF1dx2,dF2dx1,dF2dx2,
     &                 dF3dx1,dF3dx2,dF4dx1,dF4dx2,
     &                 dx1dx,dx2dx,dx1dy,dx2dy
     
      c_theta = theta
      if (r.ne.ZR) then
          r2 = sqrt(r)
      else
          r2 = sqrt(Ave_Elem_Area)*0.1D-4
          c_theta = ZR
      end if
      
      if (c_mat_type ==1 ) then
            fac = HLF/r2
            st2 = sin(c_theta/TWO)
            ct2 = cos(c_theta/TWO)
            s3t2= sin(1.5D0*c_theta)
            c3t2= cos(1.5D0*c_theta)
            st  = sin(c_theta)
            ct  = cos(c_theta)

          
            F(1) = r2*st2
            F(2) = r2*ct2
            F(3) = r2*st2*st
            F(4) = r2*ct2*st

            dF1dx1 = -fac*st2
            dF1dx2 =  fac*ct2
            dF2dx1 =  dF1dx2
            dF2dx2 = -dF1dx1
            dF3dx1 = -fac*s3t2*st
            dF3dx2 =  fac*(st2 + s3t2*ct)
            dF4dx1 = -fac*c3t2*st
            dF4dx2 =  fac*(ct2 + c3t2*ct)

            dx1dx =  cos(omega)
            dx2dx = -sin(omega)
            dx1dy =  sin(omega)
            dx2dy =  cos(omega)

            dFdx(1) = dF1dx1*dx1dx + dF1dx2*dx2dx
            dFdy(1) = dF1dx1*dx1dy + dF1dx2*dx2dy
            dFdx(2) = dF2dx1*dx1dx + dF2dx2*dx2dx
            dFdy(2) = dF2dx1*dx1dy + dF2dx2*dx2dy
            dFdx(3) = dF3dx1*dx1dx + dF3dx2*dx2dx
            dFdy(3) = dF3dx1*dx1dy + dF3dx2*dx2dy
            dFdx(4) = dF4dx1*dx1dx + dF4dx2*dx2dx
            dFdy(4) = dF4dx1*dx1dy + dF4dx2*dx2dy
      elseif (c_mat_type ==2 .or. c_mat_type ==3) then
      
      end if
      
      
      if(Key_TipEnrich==4)then
            st2 = sin(c_theta/TWO)
            ct2 = cos(c_theta/TWO)  
          

            F(1) = r*st2
            dF1dx1 = st2
            dF1dx2 = r/TWO*ct2
          
          
            dx1dx =  cos(omega)
            dx2dx = -sin(omega)
            dx1dy =  sin(omega)
            dx2dy =  cos(omega)
            dFdx(1) = dF1dx1*dx1dx + dF1dx2*dx2dx
            dFdy(1) = dF1dx1*dx1dy + dF1dx2*dx2dy
            dFdx(2) = ZR
            dFdy(2) = ZR
            dFdx(3) = ZR
            dFdy(3) = ZR
            dFdx(4) = ZR
            dFdy(4) = ZR
      endif
      
      return 
      end SUBROUTINE Cal_F_dFdx_dFdy                 
