 
      subroutine Cal_F_dFdx_dFdy_dFdz_3D(r,theta,Matrix_T,
     &                                 c_mat_type,
     &                                 F,dFdx,dFdy,dFdz)
     
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Common
      implicit none
      
      real(kind=FT),intent(in)::r,theta,Matrix_T(3,3)
      integer,intent(in)::c_mat_type
      real(kind=FT),intent(out)::dFdx(4),dFdy(4),dFdz(4),F(4)
      
      real(kind=FT) fac,st2,ct2,s3t2,c3t2,st,ct,
     &                 r2,c_theta,dF1dx1,dF1dx2,dF2dx1,dF2dx2,
     &                 dF3dx1,dF3dx2,dF4dx1,dF4dx2,
     &                 dx1dx,dx2dx,dx1dy,dx2dy
      real(kind=FT) dF1dx3,dF2dx3,dF3dx3,dF4dx3,dx1dz,dx2dz,
     &              dx3dx,dx3dy,dx3dz
     
      c_theta = theta
      if (r.ne.ZR) then
          r2 = sqrt(r)
      else
          r2 = sqrt(Ave_Elem_Area)*0.1D-4
            c_theta = ZR
      end if
      
      dFdx(1:4) = ZR
      dFdy(1:4) = ZR
      dFdz(1:4) = ZR
      F(1:4)    = ZR
      
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
            dF1dx3 = ZR
            dF2dx1 =  dF1dx2
            dF2dx2 = -dF1dx1
            dF2dx3 = ZR
            dF3dx1 = -fac*s3t2*st
            dF3dx2 =  fac*(st2 + s3t2*ct)
            dF3dx3 = ZR
            dF4dx1 = -fac*c3t2*st
            dF4dx2 =  fac*(ct2 + c3t2*ct)
            dF4dx3 = ZR

            dx1dx =  Matrix_T(1,1)
            dx1dy =  Matrix_T(1,2)
            dx1dz =  Matrix_T(1,3)          
            dx2dx =  Matrix_T(2,1)   
            dx2dy =  Matrix_T(2,2)
            dx2dz =  Matrix_T(2,3)
            dx3dx =  Matrix_T(3,1)   
            dx3dy =  Matrix_T(3,2)
            dx3dz =  Matrix_T(3,3)
          
            dFdx(1) = dF1dx1*dx1dx + dF1dx2*dx2dx + dF1dx3*dx3dx 
            dFdx(2) = dF2dx1*dx1dx + dF2dx2*dx2dx + dF2dx3*dx3dx 
            dFdx(3) = dF3dx1*dx1dx + dF3dx2*dx2dx + dF3dx3*dx3dx 
            dFdx(4) = dF4dx1*dx1dx + dF4dx2*dx2dx + dF4dx3*dx3dx 
          
            dFdy(1) = dF1dx1*dx1dy + dF1dx2*dx2dy + dF1dx3*dx3dy 
            dFdy(2) = dF2dx1*dx1dy + dF2dx2*dx2dy + dF2dx3*dx3dy 
            dFdy(3) = dF3dx1*dx1dy + dF3dx2*dx2dy + dF3dx3*dx3dy 
            dFdy(4) = dF4dx1*dx1dy + dF4dx2*dx2dy + dF4dx3*dx3dy 
          
            dFdz(1) = dF1dx1*dx1dz + dF1dx2*dx2dz + dF1dx3*dx3dz 
            dFdz(2) = dF2dx1*dx1dz + dF2dx2*dx2dz + dF2dx3*dx3dz 
            dFdz(3) = dF3dx1*dx1dz + dF3dx2*dx2dz + dF3dx3*dx3dz 
            dFdz(4) = dF4dx1*dx1dz + dF4dx2*dx2dz + dF4dx3*dx3dz           
      elseif (c_mat_type ==2 .or. c_mat_type ==3) then
      
      end if
      
      
      
      return 
      end SUBROUTINE Cal_F_dFdx_dFdy_dFdz_3D              
