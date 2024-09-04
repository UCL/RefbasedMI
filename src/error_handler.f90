!#####################################################################
module error_handler
   ! Generic error message handler for both console and non-console
   ! applications. The routines in this module do not halt program
   ! execution; they merely store error messages for subsequent
   ! retrieval.
   use program_constants
   implicit none
   private ! by default
   ! declare public types
   public :: error_type
   ! declare public subroutines and functions
   public :: err_reset, err_handle, err_msg_present, err_get_msgs, &
        err_get_codes
   ! Parameters private to this module
   integer, parameter :: &
        ! max width of any single error message line
        err_msg_width=72
   !##################################################################
   type :: err_msg_line
      ! Private type for holding one line in the linked list
      sequence
      integer(our_int) :: err_type = 0
      !    err_type = 1 : a comment line
      !               2 : OCCURED IN: (whichsub) in MOD (whichmod)
      !               3 : Observation (iobs)
      !               4 : Variable (ivar)
      !               5 : Iteration (iiter)
      !               6 : Iteration (iiter), Cycle (icycle)
      character(len=err_msg_width) :: line = "" ! text message line
      character(len=err_msg_width) :: comment = "" ! text comment
      integer(our_int) :: icomment = 0             ! integer code
      character(len=err_msg_width) :: whichmod = "" ! module name
      integer(our_int) :: iwhichmod = 0             ! integer code
      character(len=err_msg_width) :: whichsub = "" ! procedure name
      integer(our_int) :: iwhichsub = 0             ! integer code
      integer(our_int) :: iobs = 0, ivar = 0, iiter=0, icycle=0
      type(err_msg_line), pointer :: next=>null()
   end type err_msg_line
   !##################################################################
   type :: error_type
      ! Public type for holding a linked list of messages
      sequence
      private ! contents of this type are private to this module
      logical :: msg_present=.false.
      type(err_msg_line), pointer :: head=>null(), tail=>null()
   end type error_type
   !##################################################################
contains
   !##################################################################
   integer(our_int) function get_modname_icode( whichmod ) &
         result(answer)
      implicit none
      character(len=*), intent(in) :: whichmod
	  character(len=*), parameter :: &
		   nonredundant_modnames(6) = (/ &
		   "dynalloc                                                                ", &
		   "matrix_methods                                                          ", &
		   "norm_engine                                                             ", &
		   "quick_sort                                                              ", &
		   "random_generator                                                        ", &
		   "                                    " // &
		   "                                    " /)
      integer(our_int) :: i
      answer = 0
      if( whichmod == "" ) return
      do i = 1, size( nonredundant_modnames )
         if( whichmod == nonredundant_modnames(i) ) then
            answer = i
            exit
         end if
      end do
    end function get_modname_icode
   !##################################################################
   integer(our_int) function get_subname_icode( whichsub ) &
         result(answer)
      implicit none
      character(len=*), intent(in) :: whichsub
	  character(len=*), parameter :: nonredundant_subnames(78) = (/ &
		   "allocate_workspace_items                                                ", &
		   "cholesky_saxpy                                                          ", &
		   "dbl1_alloc                                                              ", &
		   "dbl1_dealloc                                                            ", &
		   "dbl2_alloc                                                              ", &
		   "dbl2_dealloc                                                            ", &
		   "dbl3_alloc                                                              ", &
		   "dbl3dealloc                                                             ", &
		   "dbl4_alloc                                                              ", &
		   "dbl4dealloc                                                             ", &
		   "estimate_worst_frac                                                     ", &
		   "fill_sort_strings_and_sort                                              ", &
		   "find_max_rel_diff                                                       ", &
		   "find_means_and_variances                                                ", &
		   "house                                                                   ", &
		   "householder_ols                                                         ", &
		   "householder_qr                                                          ", &
		   "householder_qr2                                                         ", &
		   "householder_qr2_pivot                                                   ", &
		   "householder_qr_pivot                                                    ", &
		   "initgn                                                                  ", &
		   "int1_alloc                                                              ", &
		   "int1_dealloc                                                            ", &
		   "int2_alloc                                                              ", &
		   "int2_dealloc                                                            ", &
		   "int3_alloc                                                              ", &
		   "int3_dealloc                                                            ", &
		   "invert_lower                                                            ", &
		   "log1_alloc                                                              ", &
		   "log1_dealloc                                                            ", &
		   "log2_alloc                                                              ", &
		   "log2_dealloc                                                            ", &
		   "log3_alloc                                                              ", &
		   "log3_dealloc                                                            ", &
		   "make_xsort_and_ysort                                                    ", &
		   "make_xtxinv                                                             ", &
		   "matmul_boundcheck                                                       ", &
		   "mltmod                                                                  ", &
		   "nullify_workspace_type                                                  ", &
		   "premult_lower_by_transpose                                              ", &
		   "qsort_double                                                            ", &
		   "qsort_integer                                                           ", &
		   "qsort_strings                                                           ", &
		   "ran_genchi                                                              ", &
		   "ran_gengam                                                              ", &
		   "ran_genint                                                              ", &
		   "ran_genunf                                                              ", &
		   "ran_phrasd                                                              ", &
		   "ran_setall                                                              ", &
		   "ran_sexp                                                                ", &
		   "ran_sgamma                                                              ", &
		   "ran_snorm                                                               ", &
		   "ran_timeseed                                                            ", &
		   "regress_univariate                                                      ", &
		   "row_house                                                               ", &
		   "row_house_vec                                                           ", &
		   "run_estep                                                               ", &
		   "run_istep                                                               ", &
		   "run_mstep                                                               ", &
		   "run_norm_engine_em                                                      ", &
		   "run_norm_engine_impute_mean                                             ", &
		   "run_norm_engine_impute_random                                           ", &
		   "run_norm_engine_logpost                                                 ", &
		   "run_norm_engine_mcmc                                                    ", &
		   "run_pstep                                                               ", &
		   "sort_cases_by_missingness                                               ", &
		   "str1_alloc                                                              ", &
		   "str1_dealloc                                                            ", &
		   "str2_alloc                                                              ", &
		   "str2_dealloc                                                            ", &
		   "str3_alloc                                                              ", &
		   "str3_dealloc                                                            ", &
		   "sweep_forward                                                           ", &
		   "sweep_reverse                                                           ", &
		   "update_rates                                                            ", &
		   "update_theta                                                            ", &
		   "update_theta_and_oldtheta                                               ", &
		   "                                    " // &  ! Concatenate additional literal strings
		   "                                    " /)
      integer(our_int) :: i
      answer = 0
      if( whichsub == "" ) return
      do i = 1, size( nonredundant_subnames )
         if( whichsub == nonredundant_subnames(i) ) then
            answer = i
            exit
         end if
      end do
    end function get_subname_icode
   !##################################################################
   integer(our_int) function get_comment_icode( comment ) &
         result(answer)
      implicit none
      character(len=*), intent(in) :: comment
	  character(len=*), parameter :: nonredundant_comments(90) = (/ &
		   "A zero-length string is invalid                                         ", &
		   "Argument beta has incorect length                                       ", &
		   "Argument beta has incorrect shape                                       ", &
		   "Argument beta_series has incorrect shape                                ", &
		   "Argument criterion cannot be negative                                   ", &
		   "Argument imp_list has incorrect shape                                   ", &
		   "Argument impute_every cannot be negative                                ", &
		   "Argument impute_every cannot exceed max_iter                            ", &
		   "Argument iwkp has incorect length                                       ", &
		   "Argument loglik has incorrect size                                      ", &
		   "Argument logpost has incorrect size                                     ", &
		   "Argument max_iter cannot be negative                                    ", &
		   "Argument max_iter must be positive                                      ", &
		   "Argument mis has incorrect shape                                        ", &
		   "Argument mis_int has incorrect size                                     ", &
		   "Argument multicycle must be positive                                    ", &
		   "Argument n_in_patt has incorrect size                                   ", &
		   "Argument nobs has incorrect size                                        ", &
		   "Argument omit has incorect length                                       ", &
		   "Argument piv has incorrect length                                       ", &
		   "Argument prior_sscp has incorrect dimensions                            ", &
		   "Argument r has incorect dimensions                                      ", &
		   "Argument rate_beta has incorrect shape                                  ", &
		   "Argument rate_sigma has incorrect shape                                 ", &
		   "Argument sigma has incorrect shape                                      ", &
		   "Argument sigma_series has incorrect shape                               ", &
		   "Argument which_patt has incorrect size                                  ", &
		   "Argument wkn has incorect length                                        ", &
		   "Argument wkn1 or wkn2 has incorect length                               ", &
		   "Argument wknp has incorect dimensions                                   ", &
		   "Argument wkp has incorect length                                        ", &
		   "Argument wkp1 or wkp2 has incorect length                               ", &
		   "Argument wkpp has incorect dimensions                                   ", &
		   "Argument worst_linear_coef has incorrect size                           ", &
		   "Argument worst_series has incorrect size                                ", &
		   "Argument x has more columns than rows                                   ", &
		   "Argument ybar has incorrect size                                        ", &
		   "Argument yimp has incorrect shape                                       ", &
		   "Argument ysdv has incorrect size                                        ", &
		   "Arguments a and b not conformable                                       ", &
		   "Arguments a and c not conformable                                       ", &
		   "Arguments b and c not conformable                                       ", &
		   "Arguments x and q have different dimensions                             ", &
		   "Arguments x and y have incompatible dimensions                          ", &
		   "Attempted division by zero                                              ", &
		   "Attempted division by zero; procedure aborted                           ", &
		   "Attempted logarithm of non-positive number                              ", &
		   "Cannot estimate variance; fewer than 2 cases                            ", &
		   "Computation of loglikelihood aborted                                    ", &
		   "Could not initialize random seeds                                       ", &
		   "Cov. matrix became non-positive definite                                ", &
		   "Cov. matrix became singular or negative definite                        ", &
		   "Degrees of freedom not positive                                         ", &
		   "Dimensions of matrix arguments not conformable                          ", &
		   "Dimensions of x and y not conformable                                   ", &
		   "EM algorithm aborted                                                    ", &
		   "Eigen power method failed to converge                                   ", &
		   "Error: stack size is too small                                          ", &
		   "Fewer than two observed cases available                                 ", &
		   "Finite-differencing procedure strayed outside                           ", &
		   "First element of v is not 1.D0                                          ", &
		   "Imputation procedure aborted                                            ", &
		   "Incorrect dimensions for argument beta_start                            ", &
		   "Incorrect dimensions for argument sigma_start                           ", &
		   "Incorrect dimensions for array argument beta                            ", &
		   "Incorrect dimensions for array argument sigma                           ", &
		   "Incorrect size for argument worst_linear_coef                           ", &
		   "Lower bound exceeds upper bound                                         ", &
		   "MCMC procedure aborted                                                  ", &
		   "Matrix apparently singular                                              ", &
		   "Matrix not positive definite                                            ", &
		   "Non-square matrix encountered; square expected                          ", &
		   "Operation failed                                                        ", &
		   "Parameters out of order                                                 ", &
		   "Pivot out of bounds                                                     ", &
		   "Posterior distribution may not be proper                                ", &
		   "Predictor (X) matrix does not have full rank                            ", &
		   "Prior type not recognized                                               ", &
		   "Prior_df cannot be negative for ridge prior                             ", &
		   "Random generator seeds have not been set                                ", &
		   "Shape or scale parameter not positive                                   ", &
		   "Unable to allocate memory for object                                    ", &
		   "Unable to deallocate memory for object                                  ", &
		   "Unable to regress response on x                                         ", &
		   "Workspace wkn has incorrect length                                      ", &
		   "Workspace wkp has incorrect length                                      ", &
		   "Workspace wkp1 or wkp2 has incorrect length                             ", &
		   "Zero variance; observed values are identical                            ", &
		   "parameter space; solution at or near boundary                           ", &
		   "                                    " // &  ! Concatenate additional literal strings
		   "                                    " /)
      integer(our_int) :: i
      answer = 0
      if( comment == "" ) return
      do i = 1, size( nonredundant_comments )
         if( comment == nonredundant_comments(i) ) then
            answer = i
            exit
         end if
      end do
    end function get_comment_icode
   !##################################################################
   subroutine insert_msg_line(new, head, tail)
      ! Private: Inserts a message line at the end of the current
      ! list. The new line must already be allocated.
      implicit none
      ! declare arguments
      type(err_msg_line), pointer :: new, head, tail
      ! begin
      if( associated(head) ) then
         ! linked list not empty
         tail%next => new
         nullify(new%next)
         tail => new
      else
         ! list currently empty
         head => new
         tail => new
         nullify(tail%next)
      end if
   end subroutine insert_msg_line
   !##################################################################
   subroutine append_msg_line(line, err, &
        err_type, comment, whichsub, whichmod, &
        iobs, ivar, iiter, icycle )
      ! Private: Allocates and appends a new text message line into
      ! the error_type, along with extra information if desired
      implicit none
      ! declare arguments
      character(len=*), intent(in) :: line
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer(our_int), optional, intent(in) :: err_type
      character (len=*), optional, intent(in) :: &
           comment, whichsub, whichmod
      integer(our_int), optional, intent(in) :: &
           iobs, ivar, iiter, icycle
      ! declare local variables
      type(err_msg_line), pointer :: new => null()
      ! begin
      allocate(new)
      new%line = line
      if( present( err_type ) ) new%err_type = err_type
      if( present( comment ) ) then
         new%comment = comment
         new%icomment = get_comment_icode( comment )
      end if
      if( present( whichsub ) ) then
         new%whichsub = whichsub
         new%iwhichsub = get_subname_icode( whichsub )
      end if
      if( present( whichmod ) ) then
         new%whichmod = whichmod
         new%iwhichmod = get_modname_icode( whichmod )
      end if
      if( present( iobs ) ) new%iobs = iobs
      if( present( ivar ) ) new%ivar = ivar
      if( present( iiter ) ) new%iiter = iiter
      if( present( icycle ) ) new%icycle = icycle
      call insert_msg_line(new, err%head, err%tail)
   end subroutine append_msg_line
   !##################################################################
   subroutine remove_msg_line(head, tail, first)
      ! Private: Removes the first line from the list, but does not
      ! deallocate it.
      implicit none
      ! declare arguments
      type(err_msg_line), pointer :: head, tail, first
      ! begin
      if ( associated(head) ) then
         if ( associated(head%next) ) then
            ! list has two or more lines
            first => head
            head => head%next
         else
            ! list has one line
            first => head
            nullify(head,tail)
         end if
      else
         ! list is empty
         nullify(first)
      end if
   end subroutine remove_msg_line
   !##################################################################
   recursive subroutine kill_err_list(head, tail)
      ! Private: Recursively remove the first line from the list and
      ! deallocate it until the list is empty
      implicit none
      ! declare arguments
      type(err_msg_line), pointer :: head, tail
      ! declare local variables
      type(err_msg_line), pointer :: first
      ! begin
      if( associated(head) ) then
         call remove_msg_line(head, tail, first)
         deallocate(first)
         call kill_err_list(head, tail)
      end if
   end subroutine kill_err_list
   !##################################################################
   subroutine err_reset(err)
      ! Public: Wipes out any error messages present
      implicit none
      ! declare arguments
      type(error_type), intent(inout) :: err
      ! begin
      call kill_err_list(err%head, err%tail)
      err%msg_present = .false.
   end subroutine err_reset
   !##################################################################
   subroutine err_handle(err, err_type, comment, whichsub, whichmod, &
        iobs, ivar, iiter, icycle )
      ! Public: Stores various types of message lines in the
      ! error handler
      !    err_type = 1 : just store a comment line
      !               2 : OCCURED IN: (whichsub) in MOD (whichmod)
      !               3 : Observation (iobs)
      !               4 : Variable (ivar)
      !               5 : Iteration (iiter)
      !               6 : Iteration (iiter), Cycle (icycle)
      implicit none
      ! declare required arguments
      type(error_type), intent(inout) :: err
      integer(our_int), intent(in) :: err_type
      ! declare optional arguments
      character (len=*), optional, intent(in) :: &
           comment, whichsub, whichmod
      integer(our_int), optional, intent(in) :: &
           iobs, ivar, iiter, icycle
      ! local variables
      character(len=12) :: sInt, sInt2
      character(len=err_msg_width) :: msg_line
      integer(our_int) :: posn
      ! begin
      err%msg_present = .true.
      ! to enable string-based messaging, un-comment the lines
      ! where we write to sInt and sInt2
      sInt  = "???"
      sInt2 = "???"
      select case(err_type)
         case(1)
            ! just store the comment
            if( present(comment) ) then
               msg_line = comment
            else
               msg_line = "???"
            end if
            call append_msg_line(msg_line, err, &
                 err_type = err_type, comment = comment )
         case(2)
            ! OCCURRED IN:
            msg_line = "OCCURRED IN:"
            posn = index( msg_line, "  ") + 1
            if( present(whichsub) ) then
               msg_line(posn:) = adjustl( whichsub )
            else
               msg_line(posn:) = "???"
            end if
            posn = index( msg_line, "  ") + 1
            if( present(whichmod) ) then
               msg_line(posn:) = "in MOD " // adjustl( whichmod )
            else
               msg_line(posn:) = "in MOD ???"
            end if
            call append_msg_line(msg_line, err, &
                 err_type = err_type, whichsub = whichsub, whichmod = whichmod )
         case(3)
            ! Observation (iobs)
            if( present( iobs ) ) then
               ! write(sInt,"(I12)") iobs
               sInt = adjustl(sInt)
            else
               sInt = "???"
            end if
            msg_line = "Observation " // sInt
            call append_msg_line(msg_line, err, &
                 err_type = err_type, iobs = iobs )
         case(4)
            ! Variable (ivar)
            if( present( ivar ) ) then
               ! write(sInt,"(I12)") ivar
               sInt = adjustl(sInt)
            else
               sInt = "???"
            end if
            msg_line = "Variable " // sInt
            call append_msg_line(msg_line, err, &
                 err_type = err_type, ivar = ivar )
         case(5)
            ! Iteration (iiter)
            if( present( iiter ) ) then
               ! write(sInt,"(I12)") iiter
               sInt = adjustl(sInt)
            else
               sInt = "???"
            end if
            msg_line = "Iteration " // sInt
            call append_msg_line(msg_line, err, &
                 err_type = err_type, iiter = iiter )
         case(6)
            ! Iteration (iiter), Cycle (icycle)
            if( present( iiter ) ) then
               ! write(sInt,"(I12)") iiter
               sInt = adjustl(sInt)
            else
               sInt = "???"
            end if
            if( present( icycle ) ) then
               ! write(sInt2,"(I12)") icycle
               sInt2 = adjustl(sInt2)
            else
               sInt2 = "???"
            end if
            msg_line = "Iteration " // trim(sInt) // ", " // &
                 trim(sInt2)
            call append_msg_line(msg_line, err, &
                 err_type = err_type, iiter = iiter, icycle = icycle )
         case default
            call append_msg_line("Unknown error type", err)
      end select
   end subroutine err_handle
   !##################################################################
   logical function err_msg_present(err)
      ! Public: Queries the error_type to see if a message is present
      implicit none
      type(error_type), intent(inout) :: err
      err_msg_present = err%msg_present
   end function err_msg_present
   !##################################################################
   subroutine err_get_msgs(err, msg_string, platform)
      ! Public: Retrieves all stored messages as a single character
      ! string, with message lines separated by platform-appropriate
      ! ASCII carriage control characters.
      ! Values for platform may be "UNIX", "MAC", "PC" or "S+"
      implicit none
      ! required arguments
      type(error_type), intent(inout) :: err
      character(len=*), intent(out) :: msg_string
      ! optional arguments
      character(len=*), intent(in), optional :: platform
      ! local variables
      character(len=4) :: plat
      integer :: posn
      logical :: first_time
      type(err_msg_line), pointer :: cur_line=>null()
      ! determine platform
      if( present(platform) ) then
         plat = platform
      else
         plat = "PC"
      end if
      ! clean out msg_string
      msg_string = ""
      ! step through the linked list, appending the lines
      first_time = .true.
      cur_line => err%head
      do
         if( .not.associated(cur_line) ) exit
         posn = len_trim(msg_string)
         if( (posn+3) >= len(msg_string) ) exit ! out of space
         posn = posn + 1
         if( .not.first_time) then
            select case(plat)
               case("UNIX")
                  ! Separate lines with LF
                  msg_string(posn:) = achar(10)
                  posn = posn + 1
               case("MAC")
                  ! Separate lines with CR
                  msg_string(posn:) = achar(13)
                  posn = posn + 1
               case("S+")
                  ! Separate lines with LF
                  msg_string(posn:) = char(10)
                  posn = posn + 1
               case default
                  msg_string(posn:) = char(13) // char(10)
                  posn = posn + 2
            end select
         end if
         msg_string(posn:) = trim(cur_line%line)
         first_time = .false.
         cur_line => cur_line%next
      end do
   end subroutine err_get_msgs
   !##################################################################
   subroutine err_get_codes(err, msg_codes, nlines)
      ! Public: Retrieves stored messages as a matrix of integer
      ! codes. The correspondence between integer codes and
      ! character strings is determined by the errcodes program,
      ! which scans the module source files and compiles lists
      ! of modules, procedures and comments.
      !
      !    msg_codes(n,p) should have at least 8 columns:
      !
      !    1. err_type = 1 : (comment)
      !                  2 : OCCURED IN: (whichsub) in MOD (whichmod)
      !                  3 : Observation (iobs)
      !                  4 : Variable (ivar)
      !                  5 : Iteration (iiter)
      !                  6 : Iteration (iiter), Cycle (icycle)
      !    2. iwhichmod : integer code for whichmod (icodes_1.h)
      !    3. iwhichsub : integer code for whichsub (icodes_2.h)
      !    4. icomment  : integer code for comment  (icodes_3.h)
      !    5. iobs
      !    6. ivar
      !    7. iiter
      !    8. icycle
      !
      implicit none
      ! required arguments
      type(error_type), intent(inout) :: err
      integer(our_int), intent(out) :: msg_codes(:,:)
      integer(our_int), intent(out) :: nlines ! no. of message lines
      ! local variables
      type(err_msg_line), pointer :: cur_line=>null()
      ! clean out msg_codes
      msg_codes(:,:) = 0
      ! step through the linked list, appending the lines
      nlines = 0
      cur_line => err%head
      do
         if( .not.associated(cur_line) ) exit
         if( nlines > size( msg_codes, 2 ) ) exit ! out of space
         nlines = nlines + 1
         msg_codes( nlines, 1 ) = cur_line%err_type
         msg_codes( nlines, 2 ) = cur_line%iwhichmod
         msg_codes( nlines, 3 ) = cur_line%iwhichsub
         msg_codes( nlines, 4 ) = cur_line%icomment
         msg_codes( nlines, 5 ) = cur_line%iobs
         msg_codes( nlines, 6 ) = cur_line%ivar
         msg_codes( nlines, 7 ) = cur_line%iiter
         msg_codes( nlines, 8 ) = cur_line%icycle
         cur_line => cur_line%next
      end do
    end subroutine err_get_codes
   !##################################################################
end module error_handler
!#####################################################################
