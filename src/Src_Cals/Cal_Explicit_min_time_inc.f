 
      subroutine Cal_Explicit_min_time_inc
      

      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Material
      use Global_Common
      use Global_Model
      use Global_Dynamic
      
      implicit none
      
      integer i_mat
      real(kind=FT)  c_E, c_v,c_density,c(500),c_max 
      
      c(1:500) = ZR
      do i_mat = 1,num_of_Material
          c_E  = Material_Para(i_mat,1)
          c_v  = Material_Para(i_mat,2)
          c_density = Material_Para(i_mat,3)    
          c(i_mat) = sqrt(c_E/(ONE-c_v**2)/c_density)
      enddo
      c_max = maxval(c) 
      
      Explicit_time_inc = 0.9D0*Min_Ele_Edge_Length/c_max
      
      
      return 
      end SUBROUTINE Cal_Explicit_min_time_inc               
