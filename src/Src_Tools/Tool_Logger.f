 
      module Module_Tool_Logger
      implicit none
      private
      public :: log_msg
      public :: log_startup
      public :: log_shutdown
      public :: log_get_unit
      public :: log_delimiter
      public :: log_get_delimiter
      public :: log_isinitialized
      public :: log_set_stoponerror
      public :: log_configure
      public :: log_cget
      public :: log_reset
      public :: log_inactivate_file
      public :: log_activate_file
      public :: log_inactivate_screen
      public :: log_activate_screen
      public :: log_init
      integer :: log_fileunit = 6
      integer :: log_stdout = -1
      logical :: activate_screen = .true.
      logical :: activate_file = .true.
      logical :: log_timestamp = .false.
      logical :: logger_initialized = .false.
      integer , parameter , public :: LOG_LEVEL_DELIMITER_LENGTH = 50
      character (len=LOG_LEVEL_DELIMITER_LENGTH) :: 
     &            log_level_string_volume = "==============="
      character (len=LOG_LEVEL_DELIMITER_LENGTH) :: 
     &            log_level_string_chapter = "---------------"
      logical, save :: logger_stoponerror = .true.
      integer , parameter , public :: LOG_LEVEL_VOLUME = 1
      integer , parameter , public :: LOG_LEVEL_CHAPTER = 2
      integer , parameter , public :: LOG_LEVEL_SECTION = 3
      integer , parameter , public :: LOG_LEVEL_SUBSECTION = 4
      interface log_configure
      module procedure log_configure_logical
      module procedure log_configure_integer
      module procedure log_configure_character
      end interface log_configure
      interface log_cget
      module procedure log_cget_logical
      module procedure log_cget_integer
      module procedure log_cget_character
      end interface log_cget

      contains
      subroutine log_startup (log_file, append )
      character(len=*), intent(in) :: log_file
      logical, intent(in), optional :: append
      logical :: append_real
      if ( present ( append ) ) then
           append_real = append
      else
           append_real = .true.
      endif
      if ( logger_initialized ) then
           call log_error
     &       ( "Logger is allready initialized in log_startup." )
      else
           log_fileunit = log_get_freeunit()
           if ( append_real ) then
              open (log_fileunit, FILE= log_file , 
     &              ACTION='WRITE', STATUS='UNKNOWN', 
     &              POSITION ='APPEND')
           else
              open (log_fileunit, FILE= log_file ,
     &              ACTION='WRITE', STATUS='UNKNOWN')
           endif
           logger_initialized = .true.
      endif
      end subroutine log_startup
      subroutine log_shutdown ()
      close ( log_fileunit )
      logger_initialized = .false.
      end subroutine log_shutdown
      subroutine log_reset ()
      activate_screen = .true.
      activate_file = .true.
      log_timestamp = .false.
      end subroutine log_reset
      subroutine log_msg( msg )
      character(len=*), intent(in) :: msg
      character(len=40)            :: date_string
      character(len=40)            :: time_string
      character(len=40)            :: stamp
      if ( log_timestamp ) then
           call date_and_time( date = date_string, time = time_string )
           write( stamp, '(11a)' ) 
     &           date_string(1:4), '-', 
     &           date_string(5:6), '-', date_string(7:8), ' ',
     &           time_string(1:2), ':', time_string(3:4),
     &            ':', time_string(5:6)
      else
           stamp = ' '
      endif
      if (activate_screen) then
           if ( log_timestamp ) then
              call log_write ( log_stdout, trim(stamp) // ' ' // msg )
           else
              call log_write ( log_stdout, msg )
           endif
      endif
      if (activate_file) then
           if ( log_timestamp ) then
              call log_write ( log_fileunit, trim(stamp) // ' ' // msg )
           else
              call log_write ( log_fileunit, msg )
           endif
      endif
      end subroutine log_msg
      subroutine log_write( unit, msg )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: msg
      character(len=500)           :: filename
      if (  unit == -1 ) then
           write ( *, '(a)' ) trim(msg)
      else
           write ( unit, '(a)' ) trim(msg)
           inquire( unit, name = filename )
           close( unit )
           open( unit, FILE=filename, ACTION='WRITE', STATUS='UNKNOWN', 
     &           POSITION ='APPEND')
      endif
      end subroutine log_write
      subroutine log_delimiter( level )
      integer , intent(in), optional :: level
      character(len=40)              :: msg
      integer                        :: used_level
      if (present(level)) then
           used_level = level
      else
           used_level = LOG_LEVEL_VOLUME
      endif
      call log_get_delimiter( used_level , msg )
      call log_msg( msg )
      end subroutine log_delimiter
      integer function log_get_freeunit ( )
      integer :: iunit
      integer :: ios
      logical :: lopen
      logical :: unit_found
      iunit = 0
      unit_found = .false.
      log_get_freeunit = 0
      do iunit = 1, 100
           if ( iunit /= 5 .and. iunit /= 6 .and. iunit /= 9 ) then
              inquire ( UNIT = iunit, opened = lopen, iostat = ios )
              if ( ios == 0 ) then
                 if ( .not. lopen ) then
                    log_get_freeunit = iunit
                    unit_found = .true.
                    exit
                 endif
              endif
           endif
      enddo
      if (.NOT.unit_found) then
           write(*,*) "Logging: No free logical unit for log file"
      endif
      end function log_get_freeunit
      function log_isinitialized ( ) result ( isinitialized )
      implicit none
      logical :: isinitialized
      isinitialized = logger_initialized
      end function log_isinitialized
      subroutine log_error ( message )
      implicit none
      character (len=*), intent(in) :: message
      write ( 6, "(A)" ) "Error in m_logger."
      write ( 6 , "(A)" ) message
      call log_error_stop ( )
      end subroutine log_error
      subroutine log_error_stop ( )
      if ( logger_stoponerror ) then
           stop
      endif
      end subroutine log_error_stop
      subroutine log_set_stoponerror ( stoponerror )
      logical , intent(in) :: stoponerror
      logger_stoponerror = stoponerror
      end subroutine log_set_stoponerror
      subroutine log_configure_logical ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      logical, intent(in) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "timestamp" )
           log_timestamp = value
      case ( "writeonstdout" )
           activate_screen = value
      case ( "writeonlogfile" )
           activate_file = value
      case ( "stoponerror" )
           logger_stoponerror = value
      case default
           write (message,"(A,A,A,l5,A)") "Unknown option ", option, 
     &   " for value ", value, " in log_configure_logical"
           call log_error ( message )
      end select
      end subroutine log_configure_logical
      subroutine log_configure_integer ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      integer, intent(in) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "logfileunit" )
           log_fileunit = value
      case default
           write (message,"(A,A,A,I5,A)") "Unknown option ", option,
     &     " for value ", value, " in log_configure_integer"
           call log_error ( message )
      end select
      end subroutine log_configure_integer
      subroutine log_configure_character ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      character ( len = * ) , intent(in) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "level_string_volume" )
           log_level_string_volume = value
      case ( "level_string_chapter" )
           log_level_string_chapter = value
      case ( "level_string_section" )
           log_level_string_chapter = value
      case ( "level_string_subsection" )
           log_level_string_chapter = value
      case default
           write (message,"(A,A,A,A,A)") "Unknown option ", option, 
     &     " for value ", value, " in log_configure_character"
           call log_error ( message )
      end select
      end subroutine log_configure_character
      subroutine log_cget_logical ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      logical, intent(out) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "timestamp" )
           value = log_timestamp
      case ( "writeonstdout" )
           value = activate_screen
      case ( "writeonlogfile" )
           value = activate_file
      case ( "stoponerror" )
           value = logger_stoponerror
      case default
           write (message,"(A,l5,A)") "Unknown option ", option, 
     & " in log_cget_logical"
           call log_error ( message )
      end select
      end subroutine log_cget_logical
      subroutine log_cget_integer ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      integer, intent(out) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "logfileunit" )
           value = log_fileunit
      case default
           write (message,"(A,I5,A)") "Unknown option ", option,
     &     " in log_cget_integer"
           call log_error ( message )
      end select
      end subroutine log_cget_integer
      subroutine log_cget_character ( option , value )
      implicit none
      character ( len = * ) , intent(in) :: option
      character ( len = * ) , intent(out) :: value
      character ( len = 500 ) :: message
      select case ( option )
      case ( "level_string_volume" )
           value = log_level_string_volume 
      case ( "level_string_chapter" )
           value = log_level_string_chapter 
      case ( "level_string_section" )
           value = log_level_string_chapter
      case ( "level_string_subsection" )
           value = log_level_string_chapter
      case default
           write (message,"(A,A,A)") "Unknown option ", option,
     &           " in log_cget_character"
           call log_error ( message )
      end select
      end subroutine log_cget_character
      subroutine log_inactivate_file ()
      activate_file = .false.
      end subroutine log_inactivate_file
      subroutine log_activate_file ()
      activate_file = .true.
      end subroutine log_activate_file
      subroutine log_inactivate_screen ()
      activate_screen = .false.
      end subroutine log_inactivate_screen
      subroutine log_activate_screen ()
      activate_screen = .true.
      end subroutine log_activate_screen
      subroutine log_init (log_file, append )
      character(len=*), intent(in) :: log_file
      logical, intent(in), optional :: append
      call log_startup ( log_file, append)
      end subroutine log_init
      subroutine log_get_delimiter ( level , msg )
      implicit none
      integer , intent(in) :: level
      character(len=*), intent(out) :: msg
      select case (level)
      case (LOG_LEVEL_VOLUME)
           write(msg,*) "==============="
      case (LOG_LEVEL_CHAPTER)
           write(msg,*) "---------------"
      case (LOG_LEVEL_SECTION)
           write(msg,*) "***************"
      case (LOG_LEVEL_SUBSECTION)
           write(msg,*) "+++++++++++++++"
      case default
           write(*,*) "Bad value for the message level:" , level
           write(*,*)
           stop
      end select
      end subroutine log_get_delimiter
      function log_get_unit () result ( logger_unit )
      integer :: logger_unit
      logger_unit = log_fileunit
      end function log_get_unit
      end module Module_Tool_Logger

