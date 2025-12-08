!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
Module CMD_Progress
!http://fcode.cn/code_gen-118-1.html      
Implicit None
private
Logical , parameter , public :: CMD_PROGRESS_ABSOLUTE = .true.
Type , public :: CLS_CMD_Progress
Integer , private :: N , lens , i
Character :: M = "*" , O = "."
Character(len=64) :: Prefix
Contains
Procedure :: Set
Procedure :: Put
End Type CLS_CMD_Progress

contains

Subroutine Set( this , N , L )
  Class( CLS_CMD_Progress ) :: this
  Integer , Intent( IN ) :: N , L
  this % N    = N
  this % lens = L
  this % i = 0
  this % Prefix = " Progress: "
End Subroutine Set

Subroutine Put( this , K , bAbsol )
  Class( CLS_CMD_Progress ) :: this
  Integer , Intent( IN ) :: K
  Logical , optional :: bAbsol
  Character(len=1) :: br
  integer :: jm
  this % i = this % i + K
  if ( present( bAbsol ) ) then
      if ( bAbsol ) this % i = K
  end if
  if ( this % i > this % n ) this % i = this % n    
  jm = Nint( real( this%i * this%lens ) / real( this%N ) )
  if ( this%i < this%n ) then
    br = char(13)
  else
    br = char(10)
  end if
write( * , '(5a,f6.2,2a)',advance="no") trim(this%Prefix) , '[' , &
repeat(this%M , jm ) , repeat( this%O , this%lens-jm ) , '] ' , this%i*100.0/this%N , "%" , br
End Subroutine Put

End Module CMD_Progress

