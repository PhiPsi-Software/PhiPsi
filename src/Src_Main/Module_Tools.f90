!-----------------------------------------------------------
! Brief: Provide a small ASCII progress-bar type CLS_CMD_Progress
!        and its Set / Put methods for printing percentage progress
!        to the console (with carriage-return overwriting).
!
! Notes:   The bar draws a [***...] / percentage line; the Put
!          method supports absolute or relative advancement.
!-----------------------------------------------------------

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

