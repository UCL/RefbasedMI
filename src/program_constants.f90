!#####################################################################
module program_constants
   ! Programming constants used throughout the NORM  program.
   ! Unlike most modules, everything here is public.
   implicit none
   public
   ! Define compiler-specific KIND numbers for integers,
   ! single and double-precision reals to help ensure consistency of
   ! performance across platforms:
   integer, parameter :: our_int = selected_int_kind(9), &
        our_sgle = selected_real_kind(6,37), &
        our_dble = selected_real_kind(15,307)
   ! Define UNIT numbers for Fortran I/O:
   integer, parameter :: ctrl_file_handle = 11, &
        data_file_handle = 12, names_file_handle = 13, &
        out_file_handle = 13
   ! Define maximum lengths for various types of character strings:
   integer, parameter :: file_name_length = 256, &
        var_name_length = 8, case_id_length = 12
   ! Define the maximum line widths for various types of files:
   integer, parameter :: ctrl_line_width = 80, &
        data_line_width = 2048, names_line_width = 80, &
        out_line_width = 80
   ! Common integer values returned by all functions to indicate
   ! success or failure:
   integer(kind=our_int), parameter :: RETURN_SUCCESS = 0, &
        RETURN_FAIL = -1
   ! Module names for error messaging
   character(len=*), parameter :: modname_list(6) = (/ &
        ! list all modules used by this program,
        ! except for program_constants and error_handler
        ! make sure these are same length
        "dynalloc            ", &
        "quick_sort          ", &
        "matrix_methods      ", &
        "random_generator    ", &
        "tabulate            ", &
        "norm_engine         " &
        /)
end module program_constants
!#####################################################################
