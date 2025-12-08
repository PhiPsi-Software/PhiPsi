 
module fpm_environment
#ifndef Silverfrost
    use,intrinsic :: iso_fortran_env, only : stdin=>input_unit,   &
                                           & stdout=>output_unit, &
                                           & stderr=>error_unit
    implicit none
    private
    public :: get_os_type
    public :: os_is_unix
    public :: get_env
    public :: get_command_arguments_quoted
    public :: separator

    integer, parameter, public :: OS_UNKNOWN = 0
    integer, parameter, public :: OS_LINUX   = 1
    integer, parameter, public :: OS_MACOS   = 2
    integer, parameter, public :: OS_WINDOWS = 3
    integer, parameter, public :: OS_CYGWIN  = 4
    integer, parameter, public :: OS_SOLARIS = 5
    integer, parameter, public :: OS_FREEBSD = 6
    integer, parameter, public :: OS_OPENBSD = 7
contains
    integer function get_os_type() result(r)
        character(len=32) :: val
        integer           :: length, rc
        logical           :: file_exists
        logical, save     :: first_run = .true.
        integer, save     :: ret = OS_UNKNOWN
        !$omp threadprivate(ret, first_run)

        if (.not. first_run) then
            r = ret
            return
        end if

        first_run = .false.
        r = OS_UNKNOWN

        call get_environment_variable('OSTYPE', val, length, rc)

        if (rc == 0 .and. length > 0) then
            if (index(val, 'linux') > 0) then
                r = OS_LINUX
                ret = r
                return
            end if

            if (index(val, 'darwin') > 0) then
                r = OS_MACOS
                ret = r
                return
            end if

            if (index(val, 'win') > 0 .or. index(val, 'msys') > 0) then
                r = OS_WINDOWS
                ret = r
                return
            end if

            if (index(val, 'cygwin') > 0) then
                r = OS_CYGWIN
                ret = r
                return
            end if

            if (index(val, 'SunOS') > 0 .or. index(val, 'solaris') > 0) then
                r = OS_SOLARIS
                ret = r
                return
            end if

            if (index(val, 'FreeBSD') > 0 .or. index(val, 'freebsd') > 0) then
                r = OS_FREEBSD
                ret = r
                return
            end if

            if (index(val, 'OpenBSD') > 0 .or. index(val, 'openbsd') > 0) then
                r = OS_OPENBSD
                ret = r
                return
            end if
        end if

        call get_environment_variable('OS', val, length, rc)

        if (rc == 0 .and. length > 0 .and. index(val, 'Windows_NT') > 0) then
            r = OS_WINDOWS
            ret = r
            return
        end if

        inquire (file='/etc/os-release', exist=file_exists)

        if (file_exists) then
            r = OS_LINUX
            ret = r
            return
        end if

        inquire (file='/usr/bin/sw_vers', exist=file_exists)

        if (file_exists) then
            r = OS_MACOS
            ret = r
            return
        end if

        inquire (file='/bin/freebsd-version', exist=file_exists)

        if (file_exists) then
            r = OS_FREEBSD
            ret = r
            return
        end if
    end function get_os_type

    logical function os_is_unix(os)
        integer, intent(in), optional :: os
        integer :: build_os
        if (present(os)) then
            build_os = os
        else
            build_os = get_os_type()
        end if
        os_is_unix = build_os /= OS_WINDOWS
    end function os_is_unix

    function get_env(NAME,DEFAULT) result(VALUE)
    implicit none
    character(len=*),intent(in)          :: NAME
    character(len=*),intent(in),optional :: DEFAULT
    character(len=:),allocatable         :: VALUE
    integer                              :: howbig
    integer                              :: stat
    integer                              :: length
        length=0
        if(NAME/='')then
           call get_environment_variable(NAME, length=howbig,status=stat,trim_name=.true.)
           select case (stat)
           case (1)
               VALUE=''
           case (2)
               VALUE=''
           case default
               allocate(character(len=max(howbig,1)) :: VALUE)
               call get_environment_variable(NAME,VALUE,status=stat,trim_name=.true.)
               if(stat/=0)VALUE=''
           end select
        else
           VALUE=''
        endif
        if(VALUE==''.and.present(DEFAULT))VALUE=DEFAULT
     end function get_env

    function get_command_arguments_quoted() result(args)
    character(len=:),allocatable :: args
    character(len=:),allocatable :: arg
    character(len=1)             :: quote
    integer                      :: ilength, istatus, i
    ilength=0
    args=''
        quote=merge('"',"'",separator()=='\')
        do i=2,command_argument_count()
            call get_command_argument(number=i,length=ilength,status=istatus)
            if(istatus /= 0) then
                write(stderr,'(*(g0,1x))')'<ERROR>*get_command_arguments_stack* error obtaining argument ',i
                exit
            else
                if(allocated(arg))deallocate(arg)
                allocate(character(len=ilength) :: arg)
                call get_command_argument(number=i,value=arg,length=ilength,status=istatus)
                if(istatus /= 0) then
                    write(stderr,'(*(g0,1x))')'<ERROR>*get_command_arguments_stack* error obtaining argument ',i
                    exit
                elseif(ilength>0)then
                    if(index(arg//' ','-')/=1)then
                        args=args//quote//arg//quote//' '
                    elseif(index(arg,' ')/=0)then
                        args=args//quote//arg//quote//' '
                    else
                        args=args//arg//' '
                    endif
                else
                    args=args//repeat(quote,2)//' '
                endif
             endif
         enddo
    end function get_command_arguments_quoted

function separator() result(sep)

implicit none
character(len=:),allocatable :: arg0
integer                      :: arg0_length
integer                      :: istat
logical                      :: existing
character(len=1)             :: sep
character(len=4096)          :: name
character(len=:),allocatable :: fname
   arg0_length=0
   name=' '
   call get_command_argument(0,length=arg0_length,status=istat)
   if(allocated(arg0))deallocate(arg0)
   allocate(character(len=arg0_length) :: arg0)
   call get_command_argument(0,arg0,status=istat)
   if(index(arg0,'\')/=0)then
      sep='\'
   elseif(index(arg0,'/')/=0)then
      sep='/'
   else
      existing=.false.
      name=' '
      inquire(file=arg0,iostat=istat,exist=existing,name=name)
      if(index(name,'\')/=0)then
         sep='\'
      elseif(index(name,'/')/=0)then
         sep='/'
      else
         fname='.\'//arg0
         inquire(file=fname,iostat=istat,exist=existing)
         if(existing)then
            sep='\'
         else
            fname='./'//arg0
            inquire(file=fname,iostat=istat,exist=existing)
            if(existing)then
               sep='/'
            else
               sep=merge('\','/',index(get_env('PATH'),'\')/=0)
            endif
         endif
      endif
   endif
end function separator
#endif
end module fpm_environment
