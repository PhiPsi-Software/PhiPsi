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
#ifndef Silverfrost 
module fhash_fnv
  use iso_fortran_env, only: int32, int64
  use iso_c_binding, only: c_char
  implicit none

  private
  public :: fnv_1a, hash_string

  !> Starting seed
  integer(int64), parameter :: FNV_OFFSET_32 = 2166136261_int64

  !> Hashing prime
  integer(int64), parameter :: FNV_PRIME_32 = 16777619_int64

  !> Generic interface to perform hashing
  !>
  !> Usage:
  !>```fortran
  !> fnv_1a([seed],input)
  !>```
  !> where `input` is any of the supported types
  interface fnv_1a
    module procedure fnv_1a_char_scalar
    module procedure fnv_1a_char_scalar_seed
    module procedure fnv_1a_int32_scalar
    module procedure fnv_1a_int32_scalar_seed
    module procedure fnv_1a_int32_1d
    module procedure fnv_1a_int32_1d_seed
    module procedure fnv_1a_int64_scalar
    module procedure fnv_1a_int64_scalar_seed
    module procedure fnv_1a_int64_1d
    module procedure fnv_1a_int64_1d_seed
  end interface fnv_1a

contains


  !> Hash a single default kind character variable
  pure function fnv_1a_char_scalar(input) result(hash)
    character(*), intent(in) :: input
    integer(int64) :: hash

    hash = fnv_1a(FNV_OFFSET_32,input)

  end function fnv_1a_char_scalar


  !> Hash a character(*) string of default kind
  pure function fnv_1a_char_scalar_seed(seed, input) result(hash)
    integer(int64), intent(in) :: seed
    character(*), intent(in) :: input
    integer(int64) :: hash

    integer :: i
    integer(int64) :: item

    hash = seed

    do i=1,len(input)
      item = transfer([iachar(input(i:i),int32),0_int32],item)
      hash = ieor(hash,item) * fnv_prime_32
    end do

  end function fnv_1a_char_scalar_seed


  !> Hash a single 32bit integer
  pure function fnv_1a_int32_scalar(input) result(hash)
    integer(int32), intent(in) :: input
    integer(int64) :: hash

    hash =  fnv_1a(FNV_OFFSET_32,input)

  end function fnv_1a_int32_scalar


  !> Hash a single 32bit integer with a starting seed
  pure function fnv_1a_int32_scalar_seed(seed,input) result(hash)
    integer(int64), intent(in) :: seed
    integer(int32), intent(in) :: input
    integer(int64) :: hash

    character(len=4,kind=c_char) :: chars

    chars = transfer(input,chars)

    hash =  fnv_1a(seed,chars)

  end function fnv_1a_int32_scalar_seed


  !> Hash a 1D array of 32bit integers
  pure function fnv_1a_int32_1d(input) result(hash)
    integer(int32), intent(in) :: input(:)
    integer(int64) :: hash

    hash =  fnv_1a(FNV_OFFSET_32,input)

  end function fnv_1a_int32_1d


  !> Hash a 1D array of 32bit integers with a starting seed
  pure function fnv_1a_int32_1d_seed(seed,input) result(hash)
    integer(int64), intent(in) :: seed
    integer(int32), intent(in) :: input(:)
    integer(int64) :: hash

    integer :: i

    hash = seed
    do i=1,size(input)
      hash = fnv_1a(hash,input(i))
    end do

  end function fnv_1a_int32_1d_seed


  !> Hash a single 64bit integer
  pure function fnv_1a_int64_scalar(input) result(hash)
    integer(int64), intent(in) :: input
    integer(int64) :: hash

    hash =  fnv_1a(FNV_OFFSET_32,input)

  end function fnv_1a_int64_scalar


  !> Hash a single 64bit integer with a starting seed
  pure function fnv_1a_int64_scalar_seed(seed,input) result(hash)
    integer(int64), intent(in) :: seed
    integer(int64), intent(in) :: input
    integer(int64) :: hash

    character(len=8,kind=c_char) :: chars

    chars = transfer(input,chars)

    hash =  fnv_1a(seed,chars)

  end function fnv_1a_int64_scalar_seed


  !> Hash a 1D array of 64bit integers
  pure function fnv_1a_int64_1d(input) result(hash)
    integer(int64), intent(in) :: input(:)
    integer(int64) :: hash

    hash =  fnv_1a(FNV_OFFSET_32,input)

  end function fnv_1a_int64_1d


  !> Hash a 1D array of 64bit integers with a starting seed
  pure function fnv_1a_int64_1d_seed(seed,input) result(hash)
    integer(int64), intent(in) :: seed
    integer(int64), intent(in) :: input(:)
    integer(int64) :: hash

    integer :: i

    hash = seed
    do i=1,size(input)
      hash = fnv_1a(hash,input(i))
    end do

  end function fnv_1a_int64_1d_seed


  !> Help fcn to convert hash to hex representation
  function hash_string(hash_value) result(str)
    integer(int64), intent(in) :: hash_value
    character(:), allocatable :: str

    allocate(character(len=10) :: str)
    write(str,'(Z0)') int(hash_value,int32)

  end function hash_string


end module fhash_fnv



!////////////////////////////////////////////
!> Implements an abstract type for hash keys
!////////////////////////////////////////////
module fhash_key_base
  use iso_fortran_env, only: int32, int64
  implicit none

  private
  public fhash_key_t

  !> Abstract base type for defining hash keys
  type, abstract :: fhash_key_t
  contains
    procedure(hash_proc), deferred :: hash
    procedure(equality_proc), deferred :: equals
    generic, public :: operator(==) => equals
  end type fhash_key_t

  abstract interface

    pure function equality_proc(key1,key2) result(keys_equal)
      import
      class(fhash_key_t), intent(in) :: key1
      class(fhash_key_t), intent(in) :: key2
      logical :: keys_equal
    end function equality_proc

    pure function hash_proc(key) result(hash)
      import
      class(fhash_key_t), intent(in) :: key
      integer(int64) :: hash
    end function hash_proc

  end interface

end module fhash_key_base


!////////////////////////////////////////////////////////
!> Implements a concrete type for scalar int32 hash keys
!////////////////////////////////////////////////////////
module fhash_key_char
  use iso_fortran_env, only: int32, int64
  use fhash_key_base, only: fhash_key_t
  use fhash_fnv, only: fnv_1a
  implicit none

  private
  public fhash_key_char_t
  public fhash_key

  !> Hash table key container
  type, extends(fhash_key_t) :: fhash_key_char_t
    private
    character(:), allocatable :: value
  contains
    procedure, pass :: hash => key_hash_char
    procedure, pass :: equals => key_equal_char
  end type fhash_key_char_t

  interface fhash_key
    module procedure :: key_from_char
  end interface fhash_key

  contains


  !> Check if two keys are equal
  pure function key_equal_char(key1,key2) result(keys_equal)
    class(fhash_key_char_t), intent(in) :: key1
    class(fhash_key_t), intent(in) :: key2
    logical :: keys_equal

    keys_equal = .false.

    select type(k2=>key2)
    type is (fhash_key_char_t)
      if (allocated(key1%value) .and. allocated(k2%value)) then
        if (key1%value == k2%value) then
          keys_equal = .true.
          return
        end if
      end if
    end select

  end function key_equal_char


  !> Generate hash of key
  pure function key_hash_char(key) result(hash)
    class(fhash_key_char_t), intent(in) :: key
    integer(int64) :: hash

    hash = fnv_1a(key%value)

  end function key_hash_char


  !> Create new key container from a scalar int32
  function key_from_char(source) result(key)
    character(*), intent(in) :: source
    type(fhash_key_char_t) :: key

    key%value = source

  end function key_from_char


end module fhash_key_char

!////////////////////////////////////////////////////////
!> Implements a concrete type for scalar int32 hash keys
!////////////////////////////////////////////////////////
module fhash_key_int32
  use iso_fortran_env, only: int32, int64
  use fhash_key_base, only:fhash_key_t
  use fhash_fnv, only: fnv_1a
  implicit none

  private
  public fhash_key_int32_t
  public fhash_key

  !> Hash table key container
  type, extends(fhash_key_t) :: fhash_key_int32_t
    private
    integer(int32) :: value
  contains
    procedure, pass :: hash => key_hash_int32
    procedure, pass :: equals => key_equal_int32
  end type fhash_key_int32_t

  interface fhash_key
    module procedure :: key_from_int32
  end interface fhash_key

contains


  !> Check if two keys are equal
  pure function key_equal_int32(key1,key2) result(keys_equal)
    class(fhash_key_int32_t), intent(in) :: key1
    class(fhash_key_t), intent(in) :: key2
    logical :: keys_equal

    keys_equal = .false.

    select type(k2=>key2)
    type is (fhash_key_int32_t)
      if (key1%value == k2%value) then
        keys_equal = .true.
        return
      end if
    end select

  end function key_equal_int32


  !> Generate hash of key
  pure function key_hash_int32(key) result(hash)
      class(fhash_key_int32_t), intent(in) :: key
      integer(int64) :: hash

      hash = fnv_1a(key%value)

  end function key_hash_int32


  !> Create new key container from a scalar int32
  function key_from_int32(source) result(key)
    integer(int32), intent(in) :: source
    type(fhash_key_int32_t) :: key

    key%value = source

  end function key_from_int32


end module fhash_key_int32




!////////////////////////////////////////////////////////
!> Implements a concrete type for scalar int64 hash keys
!////////////////////////////////////////////////////////
module fhash_key_int64
  use iso_fortran_env, only: int64
  use fhash_key_base, only:fhash_key_t
  use fhash_fnv, only: fnv_1a
  implicit none

  private
  public fhash_key_int64_t
  public fhash_key

  !> Hash table key container
  type, extends(fhash_key_t) :: fhash_key_int64_t
    private
    integer(int64) :: value
  contains
    procedure, pass :: hash => key_hash_int64
    procedure, pass :: equals => key_equal_int64
  end type fhash_key_int64_t

  interface fhash_key
    module procedure :: key_from_int64
  end interface fhash_key

contains


  !> Check if two keys are equal
  pure function key_equal_int64(key1,key2) result(keys_equal)
    class(fhash_key_int64_t), intent(in) :: key1
    class(fhash_key_t), intent(in) :: key2
    logical :: keys_equal

    keys_equal = .false.

    select type(k2=>key2)
    type is (fhash_key_int64_t)
      if (key1%value == k2%value) then
        keys_equal = .true.
        return
      end if
    end select

  end function key_equal_int64


  !> Generate hash of key
  pure function key_hash_int64(key) result(hash)
      class(fhash_key_int64_t), intent(in) :: key
      integer(int64) :: hash

      hash = fnv_1a(key%value)

  end function key_hash_int64


  !> Create new key container from a scalar int64
  function key_from_int64(source) result(key)
    integer(int64), intent(in) :: source
    type(fhash_key_int64_t) :: key

    key%value = source

  end function key_from_int64


end module fhash_key_int64

!//////////////////////////////////////////////////////////
!> Implements a concrete type for 1D int32 array hash keys
!//////////////////////////////////////////////////////////
module fhash_key_int32_1d
  use iso_fortran_env, only: int32, int64
  use fhash_key_base, only: fhash_key_t
  use fhash_fnv, only: fnv_1a
  implicit none

  private
  public fhash_key_int32_1d_t
  public fhash_key

  !> Hash table key container
  type, extends(fhash_key_t) :: fhash_key_int32_1d_t
    private
    integer(int32), allocatable :: value(:)
  contains
    procedure, pass :: hash => key_hash_int32_1d
    procedure, pass :: equals => key_equal_int32_1d
  end type fhash_key_int32_1d_t

  interface fhash_key
    module procedure :: key_from_int32_1d
  end interface fhash_key

contains


  !> Check if two keys are equal
  pure function key_equal_int32_1d(key1,key2) result(keys_equal)
    class(fhash_key_int32_1d_t), intent(in) :: key1
    class(fhash_key_t), intent(in) :: key2
    logical :: keys_equal

    keys_equal = .false.

    select type(k2=>key2)
    type is (fhash_key_int32_1d_t)
      if (.not.(allocated(key1%value) .and. allocated(k2%value))) then
        return
      end if
      if (size(key1%value) /= size(k2%value)) then
        return
      end if
      if (all(key1%value == k2%value)) then
        keys_equal = .true.
        return
      end if
    end select

  end function key_equal_int32_1d


  !> Generate hash of key
  pure function key_hash_int32_1d(key) result(hash)
      class(fhash_key_int32_1d_t), intent(in) :: key
      integer(int64) :: hash

      hash = fnv_1a(key%value)

  end function key_hash_int32_1d


  !> Create new key container from a scalar int32
  function key_from_int32_1d(source) result(key)
    integer(int32), intent(in) :: source(:)
    type(fhash_key_int32_1d_t) :: key

    key%value = source

  end function key_from_int32_1d


end module fhash_key_int32_1d



!//////////////////////////////////////////////////////////
!> Implements a concrete type for 1D int64 array hash keys
!//////////////////////////////////////////////////////////
module fhash_key_int64_1d
  use iso_fortran_env, only: int64
  use fhash_key_base, only: fhash_key_t
  use fhash_fnv, only: fnv_1a
  implicit none

  private
  public fhash_key_int64_1d_t
  public fhash_key

  !> Hash table key container
  type, extends(fhash_key_t) :: fhash_key_int64_1d_t
    private
    integer(int64), allocatable :: value(:)
  contains
    procedure, pass :: hash => key_hash_int64_1d
    procedure, pass :: equals => key_equal_int64_1d
  end type fhash_key_int64_1d_t

  interface fhash_key
    module procedure :: key_from_int64_1d
  end interface fhash_key

contains


  !> Check if two keys are equal
  pure function key_equal_int64_1d(key1,key2) result(keys_equal)
    class(fhash_key_int64_1d_t), intent(in) :: key1
    class(fhash_key_t), intent(in) :: key2
    logical :: keys_equal

    keys_equal = .false.

    select type(k2=>key2)
    type is (fhash_key_int64_1d_t)
      if (.not.(allocated(key1%value) .and. allocated(k2%value))) then
        return
      end if
      if (size(key1%value) /= size(k2%value)) then
        return
      end if
      if (all(key1%value == k2%value)) then
        keys_equal = .true.
        return
      end if
    end select

  end function key_equal_int64_1d


  !> Generate hash of key
  pure function key_hash_int64_1d(key) result(hash)
      class(fhash_key_int64_1d_t), intent(in) :: key
      integer(int64) :: hash

      hash = fnv_1a(key%value)

  end function key_hash_int64_1d


  !> Create new key container from a scalar int64
  function key_from_int64_1d(source) result(key)
    integer(int64), intent(in) :: source(:)
    type(fhash_key_int64_1d_t) :: key

    key%value = source

  end function key_from_int64_1d


end module fhash_key_int64_1d



!/////////////////////////////////////////
!> Implements simple container type
!>  for polymorphic scalars and 1D arrays
!/////////////////////////////////////////
module fhash_data_container
  use iso_fortran_env, only: sp=>real32, dp=>real64, int32, int64
  implicit none

  private
  public fhash_container_t
  public fhash_container

  !> Generic container for scalar and 1D data
  type fhash_container_t

    class(*), allocatable :: scalar_data
    class(*), pointer :: scalar_ptr => NULL()

  contains

    procedure :: allocated => fhash_container_allocated
    procedure :: get => fhash_container_get_scalar
    procedure :: get_ptr => fhash_container_get_scalar_ptr

  end type fhash_container_t

  !> Create a fhash_container object from a polymorphic value
  interface fhash_container
    module procedure fhash_container_scalar
  end interface fhash_container

contains

  !> Helper to initialise a polymorphic data container with scalar
  function fhash_container_scalar(value,pointer) result(container)

    !> Value to store
    class(*), intent(in), target :: value

    !> If .true., store pointer to value instead of copying
    logical, intent(in), optional :: pointer

    type(fhash_container_t) :: container

    if (present(pointer)) then
      if (pointer) then
        container%scalar_ptr => value
      else
        if (allocated(container%scalar_data)) deallocate(container%scalar_data)
        allocate(container%scalar_data, source = value)
      end if
    else
      if (allocated(container%scalar_data)) deallocate(container%scalar_data)
      allocate(container%scalar_data, source = value)
    end if

  end function fhash_container_scalar


  !> Helper to determine if container contains anything
  function fhash_container_allocated(container) result(alloc)
    class(fhash_container_t), intent(in) :: container
    logical :: alloc

    alloc = allocated(container%scalar_data) .OR. &
                associated(container%scalar_ptr)

  end function fhash_container_allocated


  !> Helper to return container value as intrinsic type
  subroutine fhash_container_get_scalar(container,i32,i64,r32,r64,char,bool,raw,match,type_string)
    class(fhash_container_t), intent(in), target :: container
    integer(int32), intent(out), optional :: i32
    integer(int64), intent(out), optional :: i64
    real(sp), intent(out), optional :: r32
    real(dp), intent(out), optional :: r64
    character(:), allocatable, intent(out), optional :: char
    logical, intent(out), optional :: bool
    class(*), allocatable, intent(out), optional :: raw
    logical, intent(out), optional :: match
    character(:), allocatable, intent(out), optional :: type_string

    class(*), pointer :: data

    if (present(match)) match = .false.

    if (.not.container%allocated()) return

    if (allocated(container%scalar_data)) then
      data => container%scalar_data
    else
      data => container%scalar_ptr
    end if

    if (present(raw)) then
      if (present(match)) match = .true.
      raw = data
    end if

    select type(d=>data)
    type is(integer(int32))
      if (present(type_string)) type_string = 'integer32'
      if (present(i32)) then
        if (present(match)) match = .true.
        i32 = d
        return
      end if

    type is (integer(int64))
      if (present(type_string)) type_string = 'integer64'
      if (present(i64)) then
        if (present(match)) match = .true.
        i64 = d
        return
      end if

    type is (real(sp))
      if (present(type_string)) type_string = 'real32'
      if (present(r32)) then
        if (present(match)) match = .true.
        r32 = d
        return
      end if

    type is (real(dp))
      if (present(type_string)) type_string = 'real64'
      if (present(r64)) then
        if (present(match)) match = .true.
        r64 = d
        return
      end if

    type is (character(*))
      if (present(type_string)) type_string = 'character*'
      if (present(char)) then
        if (present(match)) match = .true.
        char = d
        return
      end if

    type is (logical)
      if (present(type_string)) type_string = 'logical'
      if (present(bool)) then
        if (present(match)) match = .true.
        bool = d
        return
      end if

    class default
      if (present(type_string)) type_string = 'unknown'

    end select

  end subroutine fhash_container_get_scalar


  !> Helper to return pointer to container value as intrinsic type
  subroutine fhash_container_get_scalar_ptr(container,i32,i64,r32,r64,char,bool,raw,match,type_string)
    class(fhash_container_t), intent(in), target :: container
    integer(int32), pointer, intent(out), optional :: i32
    integer(int64), pointer, intent(out), optional :: i64
    real(sp), pointer, intent(out), optional :: r32
    real(dp), pointer, intent(out), optional :: r64
    character(:), pointer, intent(out), optional :: char
    logical, pointer, intent(out), optional :: bool
    class(*), pointer, intent(out), optional :: raw
    logical, intent(out), optional :: match
    character(:), allocatable, intent(out), optional :: type_string

    class(*), pointer :: data

    if (present(match)) match = .false.

    if (.not.container%allocated()) return

    if (allocated(container%scalar_data)) then
      data => container%scalar_data
    else
      data => container%scalar_ptr
    end if

    if (present(raw)) then
      if (present(match)) match = .true.
      raw => data
    end if

    select type(d=>data)
    type is(integer(int32))
      if (present(i32)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'integer32'
        i32 => d
        return
      end if

    type is (integer(int64))
      if (present(i64)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'integer64'
        i64 => d
        return
      end if

    type is (real(sp))
      if (present(r32)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'real32'
        r32 => d
        return
      end if

    type is (real(dp))
      if (present(r64)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'real64'
        r64 => d
        return
      end if

    type is (character(*))
      if (present(char)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'character*'
        char => d
        return
      end if

    type is (logical)
      if (present(bool)) then
        if (present(match)) match = .true.
        if (present(type_string)) type_string = 'logical'
        bool => d
        return
      end if

    class default
      if (present(type_string)) type_string = 'unknown'

    end select

  end subroutine fhash_container_get_scalar_ptr


end module fhash_data_container



!///////////////////////////////////////////////////////////////////////
!> Implements singly-linked list (sll) node with generic data container
!///////////////////////////////////////////////////////////////////////
module fhash_sll
  use iso_fortran_env, only: int32, int64
  use fhash_key_base, only: fhash_key_t
  use fhash_data_container, only: fhash_container_t
  implicit none

  !> Node type for hash table singly linked list
  type fhash_node_t

    class(fhash_key_t), allocatable :: key
    type(fhash_container_t) :: value
    type(fhash_node_t), pointer :: next => NULL()

  end type fhash_node_t

contains

  !> Append node to SLL
  recursive subroutine sll_push_node(node,key,value,pointer)

    !> Node to which to add data
    type(fhash_node_t), intent(inout) :: node

    !> Key to add
    class(fhash_key_t), intent(in) :: key

    !> Value to add
    class(*), intent(in), target :: value

    !> Store only a point if .true.
    logical, intent(in), optional :: pointer


    if (allocated(node%key)) then

      if (node%key == key) then

        call sll_node_set(node,value,pointer)
        return

      end if

      if (.not.associated(node%next)) then
        allocate(node%next)
      end if

      call sll_push_node(node%next,key,value,pointer)

    else

      node%key = key
      call sll_node_set(node,value,pointer)

    end if

  end subroutine sll_push_node


  !> Set container value in node
  !>
  subroutine sll_node_set(node,value,pointer)

    !> Node to which to add data
    type(fhash_node_t), intent(inout) :: node

    !> Value to set
    class(*), intent(in), target :: value

    !> Store only a pointer if .true.
    logical, intent(in), optional :: pointer

    if (present(pointer)) then
      if (pointer) then
        node%value%scalar_ptr => value
        return
      end if
    end if

    if (allocated(node%value%scalar_data)) deallocate(node%value%scalar_data)
    allocate(node%value%scalar_data, source = value)

  end subroutine sll_node_set


  !> Search for a node with a specific key.
  !> Returns a pointer to the 'data' component of the corresponding node.
  !> Pointer is not associated if node cannot be found
  recursive subroutine sll_find_in(node,key,data,found)

    !> Node to search in
    type(fhash_node_t), intent(in), target :: node

    !> Key to look for
    class(fhash_key_t) :: key

    !> Pointer to value container if found.
    !> (Unassociated if the key is not found in node)
    type(fhash_container_t), pointer, intent(out) :: data

    logical, intent(out), optional :: found

    data => NULL()

    if (present(found)) found = .false.

    if (.not.allocated(node%key)) then

      return

    else if (node%key == key) then

      if (present(found)) found = .true.
      data => node%value
      return

    else if (associated(node%next)) then

      call sll_find_in(node%next,key,data,found)

    end if

  end subroutine sll_find_in


  !> Search for a node with a specific key and remove
  recursive subroutine sll_remove(node,key,found,parent_node)

    !> Node to remove from
    type(fhash_node_t), intent(inout) :: node

    !> Key to remove
    class(fhash_key_t) :: key

    !> Indicates if the key was found in node and removed
    logical, optional, intent(out) :: found

    !> Used internally
    type(fhash_node_t), intent(inout), optional :: parent_node

    type(fhash_node_t), pointer :: next_temp

    if (present(found)) then
      found = .false.
    end if

    if (.not.allocated(node%key)) then

      return

    else if (node%key == key) then

      if (present(found)) then
        found = .true.
      end if

      if (.not.present(parent_node)) then
        ! This is the top-level node
        if (associated(node%next)) then
          ! Replace with next
          next_temp => node%next
          node = next_temp
          deallocate(next_temp)
          return
        else
          ! No children, just deallocate
          deallocate(node%key)
          return
        end if

      else
        ! Not top-level node
        if (associated(node%next)) then
          ! Join previous with next
          next_temp => node%next
          deallocate(parent_node%next)
          parent_node%next => next_temp
          return
        else
          ! No children, just deallocate
          deallocate(node%key)
          deallocate(parent_node%next)
          return
        end if
      end if

    else if (associated(node%next)) then
      ! Look further down
      call sll_remove(node%next,key,found,node)

    end if

  end subroutine sll_remove


  !> Deallocate node components and those of its children
  recursive subroutine sll_clean(node)

    !> Node to search in
    type(fhash_node_t), intent(inout) :: node

    if (associated(node%next)) then

      call sll_clean(node%next)
      deallocate(node%next)

    end if

  end subroutine sll_clean


  !> Determine depth of SLL
  function node_depth(node) result(depth)

    !> Node to check depth
    type(fhash_node_t), intent(in), target :: node

    integer :: depth

    type(fhash_node_t), pointer :: current

    if (.not.allocated(node%key)) then

      depth = 0
      return

    else

      depth = 1
      current => node
      do while(associated(current%next))
        depth = depth + 1
        current => current%next
      end do

    end if

  end function node_depth


end module fhash_sll


!//////////////////////////////////////
!                                    !
!                                    !
!        module fhash_tbl            !
!                                    !
!                                    !
!//////////////////////////////////////
module fhash_tbl
  use iso_fortran_env, only: int32, int64, sp=>real32, dp=>real64
  use fhash_data_container, only: fhash_container
  use fhash_sll
  implicit none

  private
  public fhash_tbl_t

  !> This condition should be unreachable by the public interface
  integer, parameter, public :: FHASH_INTERNAL_ERROR = -4

  !> Error flag for operating on an unallocated table
  integer, parameter, public :: FHASH_EMPTY_TABLE = -3

   !> Error flag for when retrieved data-type does not
  !>  match that expected by the invoked getter function
  !>  (`get_int32`,`get_int63`,`get_float`,'get_double`,`get_char`)
  integer, parameter, public :: FHASH_FOUND_WRONG_TYPE = -2

  !> Error flag for when specified key is not found in the hash table
  integer, parameter, public :: FHASH_KEY_NOT_FOUND = -1

  !> Default allocation size
  integer, parameter :: FHASH_DEFAULT_ALLOCATION = 127

  type fhash_tbl_t

    type(fhash_node_t), allocatable :: buckets(:)

  contains

    procedure :: allocate => fhash_tbl_allocate
    procedure :: unset => fhash_tbl_unset
    procedure :: check_key => fhash_tbl_check_key
    procedure :: stats => fhash_tbl_stats

    procedure :: fhash_tbl_set_scalar
    generic :: set => fhash_tbl_set_scalar

    procedure :: fhash_tbl_set_scalar_ptr
    generic :: set_ptr => fhash_tbl_set_scalar_ptr

    procedure :: fhash_tbl_get_int32, fhash_tbl_get_int64
    procedure :: fhash_tbl_get_float, fhash_tbl_get_double
    procedure :: fhash_tbl_get_char,fhash_tbl_get_logical
    procedure :: fhash_tbl_get_data,fhash_tbl_get_raw

    generic :: get => fhash_tbl_get_int32, fhash_tbl_get_int64
    generic :: get => fhash_tbl_get_float, fhash_tbl_get_double
    generic :: get => fhash_tbl_get_char, fhash_tbl_get_logical
    generic :: get => fhash_tbl_get_data
    generic :: get_raw => fhash_tbl_get_raw

    procedure :: fhash_tbl_get_int32_ptr, fhash_tbl_get_int64_ptr
    procedure :: fhash_tbl_get_float_ptr, fhash_tbl_get_double_ptr
    procedure :: fhash_tbl_get_char_ptr,fhash_tbl_get_logical_ptr
    procedure :: fhash_tbl_get_raw_ptr

    generic :: get_ptr => fhash_tbl_get_int32_ptr, fhash_tbl_get_int64_ptr
    generic :: get_ptr => fhash_tbl_get_float_ptr, fhash_tbl_get_double_ptr
    generic :: get_ptr => fhash_tbl_get_char_ptr, fhash_tbl_get_logical_ptr
    generic :: get_raw_ptr => fhash_tbl_get_raw_ptr

    final :: fhash_tbl_cleanup

  end type fhash_tbl_t

  contains

    !> Allocate hash table
    subroutine fhash_tbl_allocate(tbl,size)

    !> Table object to allocate
    class(fhash_tbl_t), intent(inout) :: tbl

    !> Number of buckets in hash table
    !> If ommited, `tbl` is allocated with `FHASH_DEFAULT_ALLOCATION`
    integer, intent(in), optional :: size

    if (present(size)) then
        allocate(tbl%buckets(size))
    else
        allocate(tbl%buckets(FHASH_DEFAULT_ALLOCATION))
    end if

    end subroutine fhash_tbl_allocate


    !> Finalizer for fhash_tbl_t
    subroutine fhash_tbl_cleanup(tbl)

    !> Table object to allocate
    type(fhash_tbl_t), intent(inout) :: tbl

    integer :: i

    if (.not.allocated(tbl%buckets)) return

    do i=1,size(tbl%buckets)

        call sll_clean(tbl%buckets(i))

    end do

    end subroutine fhash_tbl_cleanup


    !> Unset a value in the table
    !>
    subroutine fhash_tbl_unset(tbl,key,stat)

    !> Hash table object
    class(fhash_tbl_t), intent(inout) :: tbl

    !> Key to remove
    class(fhash_key_t), intent(in) :: key

    !> Status flag. Zero if successful.
    !> Unsuccessful: FHASH_EMPTY_TABLE | `FHASH_KEY_NOT_FOUND`
    integer, intent(out), optional :: stat

    integer :: index
    logical :: found

    if (present(stat)) stat = 0

    if (.not.allocated(tbl%buckets)) then
        if (present(stat)) stat = FHASH_EMPTY_TABLE
        return
    end if

    index = modulo(key%hash(),size(tbl%buckets)) + 1
    call sll_remove(tbl%buckets(index),key,found)

    if (present(stat)) stat = merge(0,FHASH_KEY_NOT_FOUND,found)

    end subroutine fhash_tbl_unset


    !> Check if key exists in table
    subroutine fhash_tbl_check_key(tbl,key,stat)

    !> Hash table object
    class(fhash_tbl_t), intent(in) :: tbl

    !> Key to retrieve
    class(fhash_key_t), intent(in) :: key

    !> Status flag. Zero if key is found.
    !> Unsuccessful: `FHASH_EMPTY_TABLE` | `FHASH_KEY_NOT_FOUND`
    integer, intent(out) :: stat

    integer :: index
    logical :: found
    type(fhash_container_t), pointer :: data

    if (.not.allocated(tbl%buckets)) then
        stat = FHASH_EMPTY_TABLE
        return
    end if

    stat = 0

    index = modulo(key%hash(),size(tbl%buckets)) + 1

    call sll_find_in(tbl%buckets(index),key,data,found)

    stat = merge(0,FHASH_KEY_NOT_FOUND,found)

    return

    end subroutine fhash_tbl_check_key


    !> Get stats about the hash table
    subroutine fhash_tbl_stats(tbl,num_buckets,num_items,num_collisions,max_depth)

    !> Hash table object
    class(fhash_tbl_t), intent(in) :: tbl

    !> Number of buckets allocated in table
    integer, intent(out), optional :: num_buckets

    !> Number of key-value pairs stored in table
    integer, intent(out), optional :: num_items

    !> Number of hash collisions
    integer, intent(out), optional :: num_collisions

    !> Maximum depth of bucket in table
    integer, intent(out), optional :: max_depth

    integer :: i, depth

    ! Initialise stats
    if (present(num_items)) num_items = 0
    if (present(num_collisions)) num_collisions = 0
    if (present(max_depth)) max_depth = 0
    if (present(num_buckets)) num_buckets = 0

    if (.not.allocated(tbl%buckets)) return

    if (present(num_buckets)) then
        num_buckets = size(tbl%buckets)
    end if

    do i=1,size(tbl%buckets)

        depth = node_depth(tbl%buckets(i))

        if (present(num_items)) num_items = num_items + depth

        if (present(num_collisions)) num_collisions = num_collisions + &
                                              merge(depth-1,0,depth > 1)

        if (present(max_depth)) max_depth = max(max_depth,depth)

    end do

    end subroutine fhash_tbl_stats


    !> Set/update a polymorphic scalar value in the table
    !>
    !> `tbl` is allocated with default size if not already allocated
    subroutine fhash_tbl_set_scalar(tbl,key,value,pointer)

    !> Hash table object
    class(fhash_tbl_t), intent(inout) :: tbl

    !> Key to set/update
    class(fhash_key_t), intent(in) :: key

    !> Value for key
    class(*), intent(in), target :: value

    !> If .true., store a pointer to value instead of copying
    logical, intent(in), optional :: pointer

    integer :: index

    if (.not.allocated(tbl%buckets)) call fhash_tbl_allocate(tbl)

    index = modulo(key%hash(),size(tbl%buckets)) + 1

    call sll_push_node(tbl%buckets(index),key,value,pointer)

    end subroutine fhash_tbl_set_scalar


    !> Get wrapper routine for generic 'set_ptr'
    !>
    !> `tbl` is allocated with default size if not already allocated
    subroutine fhash_tbl_set_scalar_ptr(tbl,key,value)

    !> Hash table object
    class(fhash_tbl_t), intent(inout) :: tbl

    !> Key to set/update
    class(fhash_key_t), intent(in) :: key

    !> Value for key
    class(*), intent(in), target :: value

    call fhash_tbl_set_scalar(tbl,key,value,pointer=.true.)

    end subroutine fhash_tbl_set_scalar_ptr



    !> Retrieve data container from the hash table
    subroutine fhash_tbl_get_data(tbl,key,data,stat)

    !> Hash table object
    class(fhash_tbl_t), intent(in) :: tbl

    !> Key to retrieve
    class(fhash_key_t), intent(in) :: key

    !> Copy of value retrieved for key
    type(fhash_container_t), pointer :: data

    !> Status flag. Zero if successful.
    !> Unsuccessful: `FHASH_EMPTY_TABLE` | `FHASH_KEY_NOT_FOUND`
    integer, intent(out), optional :: stat

    integer :: index
    logical :: found

    if (.not.allocated(tbl%buckets)) then
        if (present(stat)) stat = FHASH_EMPTY_TABLE
        return
    end if

    if (present(stat)) stat = 0

    index = modulo(key%hash(),size(tbl%buckets)) + 1

    call sll_find_in(tbl%buckets(index),key,data,found)

    if (.not.found) then

        if (present(stat)) stat = FHASH_KEY_NOT_FOUND
        return

    end if

    end subroutine fhash_tbl_get_data



    !> Get wrapper to retrieve a scalar intrinsic type value
    subroutine fhash_tbl_get_intrinsic_scalar(tbl,key,i32,i64,r32,r64,char,raw,bool,stat)

    !> Hash table object
    class(fhash_tbl_t), intent(in) :: tbl

    !> Key to retrieve
    class(fhash_key_t), intent(in) :: key

    !> Value to retrieve
    integer(int32), intent(out), optional :: i32
    integer(int64), intent(out), optional :: i64
    real(sp), intent(out), optional :: r32
    real(dp), intent(out), optional :: r64
    character(:), allocatable, intent(out), optional :: char
    logical, intent(out), optional :: bool
    class(*), allocatable, intent(out), optional :: raw

    !> Status flag. Zero if successful.
    !> Unsuccessful: `FHASH_EMPTY_TABLE` |
    !>  `FHASH_FOUND_WRONG_TYPE` | `FHASH_KEY_NOT_FOUND`
    integer, intent(out), optional :: stat

    logical :: type_match
    integer :: local_stat
    type(fhash_container_t), pointer :: data

    character(:), allocatable :: char_temp

    if (present(stat)) stat = 0

    call fhash_tbl_get_data(tbl,key,data,local_stat)

    if (local_stat /= 0) then
        if (present(stat)) stat = local_stat
        return
    end if

    if (present(char)) then

        call data%get(i32,i64,r32,r64,char_temp,bool,raw,type_match)

        if (type_match) char = char_temp

    else

        call data%get(i32,i64,r32,r64,bool=bool,raw=raw,match=type_match)

    end if

    if (.not.type_match) then
        if (present(stat)) stat = FHASH_FOUND_WRONG_TYPE
        return
    end if

    end subroutine fhash_tbl_get_intrinsic_scalar


    !> Get wrapper to retrieve a scalar intrinsic type pointer
    subroutine fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,i32,i64,r32,r64,char,bool,raw,stat)

    !> Hash table object
    class(fhash_tbl_t), intent(in) :: tbl

    !> Key to retrieve
    class(fhash_key_t), intent(in) :: key

    !> Value to retrieve
    integer(int32), pointer, intent(out), optional :: i32
    integer(int64), pointer, intent(out), optional :: i64
    real(sp), pointer, intent(out), optional :: r32
    real(dp), pointer, intent(out), optional :: r64
    character(:), pointer, intent(out), optional :: char
    logical, pointer, intent(out), optional :: bool
    class(*), pointer, intent(out), optional :: raw

    !> Status flag. Zero if successful.
    !> Unsuccessful: `FHASH_EMPTY_TABLE` |
    !>  `FHASH_FOUND_WRONG_TYPE` | `FHASH_KEY_NOT_FOUND`
    integer, intent(out), optional :: stat

    logical :: type_match
    integer :: local_stat
    type(fhash_container_t), pointer :: data

    character(:), pointer :: char_temp

    if (present(stat)) stat = 0

    call fhash_tbl_get_data(tbl,key,data,local_stat)

    if (local_stat /= 0) then
        if (present(stat)) stat = local_stat
        return
    end if

    if (present(char)) then

        call data%get_ptr(i32,i64,r32,r64,char_temp,bool,raw,type_match)

        if (type_match) char => char_temp

    else

        call data%get_ptr(i32,i64,r32,r64,bool=bool,raw=raw,match=type_match)

    end if

    if (.not.type_match) then
        if (present(stat)) stat = FHASH_FOUND_WRONG_TYPE
        return
    end if

    end subroutine fhash_tbl_get_intrinsic_scalar_ptr


    !> Get wrapper to directly retrieve a scalar int32 value
    subroutine fhash_tbl_get_int32(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    integer(int32), intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,i32=value,stat=stat)

    end subroutine fhash_tbl_get_int32


    !> Get wrapper to directly retrieve a scalar int64 value
    subroutine fhash_tbl_get_int64(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    integer(int64), intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,i64=value,stat=stat)

    end subroutine fhash_tbl_get_int64


    !> Get wrapper to directly retrieve a scalar float value
    subroutine fhash_tbl_get_float(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    real(sp), intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,r32=value,stat=stat)

    end subroutine fhash_tbl_get_float


    !> Get wrapper to directly retrieve a scalar double value
    subroutine fhash_tbl_get_double(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    real(dp), intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,r64=value,stat=stat)

    end subroutine fhash_tbl_get_double


    !> Get wrapper to directly retrieve a scalar character value
    subroutine fhash_tbl_get_char(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    character(:), allocatable, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,char=value,stat=stat)

    end subroutine fhash_tbl_get_char


    !> Get wrapper to directly retrieve a scalar logical value
    subroutine fhash_tbl_get_logical(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    logical, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,bool=value,stat=stat)

    end subroutine fhash_tbl_get_logical


    !> Get wrapper to directly retrieve underlying polymorhpic scalar value
    subroutine fhash_tbl_get_raw(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    class(*), allocatable, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar(tbl,key,raw=value,stat=stat)

    end subroutine fhash_tbl_get_raw


    !> Get wrapper to directly retrieve a scalar int32 value
    subroutine fhash_tbl_get_int32_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    integer(int32), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,i32=value,stat=stat)

    end subroutine fhash_tbl_get_int32_ptr


    !> Get wrapper to directly retrieve a scalar int64 value
    subroutine fhash_tbl_get_int64_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    integer(int64), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,i64=value,stat=stat)

    end subroutine fhash_tbl_get_int64_ptr


    !> Get wrapper to directly retrieve a scalar float value
    subroutine fhash_tbl_get_float_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    real(sp), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,r32=value,stat=stat)

    end subroutine fhash_tbl_get_float_ptr


    !> Get wrapper to directly retrieve a scalar double value
    subroutine fhash_tbl_get_double_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    real(dp), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,r64=value,stat=stat)

    end subroutine fhash_tbl_get_double_ptr


    !> Get wrapper to directly retrieve a scalar character value
    subroutine fhash_tbl_get_char_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    character(:), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,char=value,stat=stat)

    end subroutine fhash_tbl_get_char_ptr


    !> Get wrapper to directly retrieve a scalar logical value
    subroutine fhash_tbl_get_logical_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    logical, pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,bool=value,stat=stat)

    end subroutine fhash_tbl_get_logical_ptr


    !> Get wrapper to directly retrieve underlying polymorhpic scalar value
    subroutine fhash_tbl_get_raw_ptr(tbl,key,value,stat)
    class(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: key
    class(*), pointer, intent(out) :: value
    integer, intent(out), optional :: stat

    call fhash_tbl_get_intrinsic_scalar_ptr(tbl,key,raw=value,stat=stat)

    end subroutine fhash_tbl_get_raw_ptr


end module fhash_tbl



!//////////////////////////////////////
!                                    !
!                                    !
!           module fhash             !
!                                    !
!                                    !
!//////////////////////////////////////
module fhash
use fhash_tbl, only: fhash_tbl_t
use fhash_key_base, only: fhash_key_t
use fhash_key_char, only: fhash_key_char_t, fhash_key
use fhash_key_int32, only: fhash_key_int32_t, fhash_key
use fhash_key_int64, only: fhash_key_int64_t, fhash_key
use fhash_key_int32_1d, only: fhash_key_int32_1d_t, fhash_key
use fhash_key_int64_1d, only: fhash_key_int64_1d_t, fhash_key
implicit none


end module fhash
#endif


#ifdef Silverfrost
module fhash
end module fhash
#endif
