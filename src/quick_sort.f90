!######################################################################
module quick_sort
   ! Routines for performing quick sorting of arrays
   ! 3/16: modified to include sorting of integer arrays
   use program_constants
   use error_handler
   implicit none
   private ! by default
   public qsort  ! generic function
   interface qsort
      module procedure qsort_strings
      module procedure qsort_double
      module procedure qsort_integer
   end interface
   ! Module parameters
   integer, parameter :: MIN_PART=7
   character(len=*), parameter :: modname = "quick_sort"
   !##################################################################
contains
   !##################################################################
   integer(kind=our_int) function qsort_strings( sarr1, iarr, idim,&
        ilen, qopt, qcase, err) result(answer)
      !	Performs a quicksort on the character array sarr1. 
      !	If qopt ==.true. then it returns the sorted array.
      ! If qopt == .false. then it does not change sarr1.
      !
      !	An integer array, iarr, containing ordinal numbers, 
      !	is sorted along with darr, providing an array of 
      !	the original indices of the values in the array.
      !
      !	An additional option, qcase, gives the option of 
      !	ignoring letter case by specifying .true. .
      !
      !	argument declarations
      integer, intent(in) :: idim, ilen
      character(len=ilen), dimension(idim), intent(inout) :: sarr1
      integer, dimension(idim), intent(out) :: iarr
      logical, intent(in) :: qopt, qcase
      type(error_type), intent(inout) :: err
      !	local declarations
      integer :: i, j, k, ir, ia, l, istk_ptr, stk_size
      integer, dimension(:), allocatable :: istk
      character(len=ilen) :: sa
      character(len=ilen), dimension(idim) :: sarr
      character(len=*), parameter :: subname = "qsort_strings"
      ! begin
      answer = RETURN_FAIL
      istk_ptr = 0
      l = 1
      !	Allocate the stack
      stk_size = (abs(idim)*2)/MIN_PART
      allocate( istk(stk_size) )
      !	Fill the string and index arrays
      do i=1,idim
         if(qcase) then
            sarr(i) = tolower(sarr1(i), ilen)
         else
            sarr(i) = sarr1(i)
         endif
         iarr(i) = i
      end do
      !	Start sorting
      ir = idim
      sort: do
         if(ir-l < MIN_PART) then
            do j=(l+1),ir
               sa = sarr(j)
               ia = iarr(j)
               do i=(j-1), l, -1
                  !if(sarr(i) <= sa) exit
                  if(lle(sarr(i), sa)) exit
                  sarr(i+1) = sarr(i)
                  iarr(i+1) = iarr(i)
               end do
               sarr(i+1) = sa
               iarr(i+1) = ia
            end do
            if(istk_ptr == 0) exit sort
            ir = istk(istk_ptr)
            l = istk(istk_ptr-1)
            istk_ptr = istk_ptr - 2
         else
            ! Use k as the index of the next
            ! partitioning element
            k = (l + ir)/2
            call sswap( sarr, k, l+1, ilen )
            call iswap( iarr, k, l+1 )
            !if(sarr(l) > sarr(ir)) then
            if(lgt(sarr(l) , sarr(ir))) then
               call sswap( sarr, l, ir, ilen )
               call iswap( iarr, l, ir )
            endif
            !if(sarr(l+1) > sarr(ir)) then
            if(lgt(sarr(l+1) , sarr(ir))) then
               call sswap( sarr, l+1, ir, ilen )
               call iswap( iarr, l+1, ir )
            endif
            !if(sarr(l) > sarr(l+1)) then
            if(lgt(sarr(l) , sarr(l+1))) then
               call sswap( sarr, l, l+1, ilen )
               call iswap( iarr, l, l+1 )
            endif
            i= l + 1
            j= ir
            sa = sarr(l+1)
            ia = iarr(l+1)
            !	Begin innermost loop
            innermost: do
               countup: do
                  i = i + 1
                  if( lge(sarr(i),sa) ) exit countup
               end do countup
               countdown: do
                  j = j - 1
                  if(lle(sarr(j) , sa)) exit countdown
               end do countdown
               if(j < i) exit innermost
               call sswap(sarr, i, j, ilen)
               call iswap(iarr, i, j)
            end do innermost
            !
            sarr(l+1) = sarr(j)
            iarr(l+1) = iarr(j)
            sarr(j) = sa
            iarr(j) = ia
            istk_ptr = istk_ptr + 2
            !	Deal with stack bounds errors
            if( istk_ptr > stk_size ) goto 700
            !
            if((ir-i+1) >= (j-l)) then
               istk(istk_ptr) = ir
               istk(istk_ptr-1) = i
               ir = j - 1
            else
               istk(istk_ptr) = j-1
               istk(istk_ptr-1) = l
               l = i
            end if
         end if
      end do sort
      !	Substitute strings if requested
      if(qopt) then
         if( qcase ) then
            ! use strings as they were before changing case
            do i=1,idim
               sarr(i) = sarr1( iarr(i) )
            end do
         end if
         do i=1,idim
            sarr1(i) = sarr(i)
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
700   call err_handle(err, 1, &
           comment = "Error: stack size is too small")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   if( allocated( istk ) ) deallocate( istk )
   end function qsort_strings
   !##################################################################
   integer(kind=our_int) function qsort_double( darr1, iarr, idim, &
        qopt, err) result(answer)
      !	Performs a quicksort on the double precision array darr. 
      !	If qopt ==.true. then it returns the sorted array.
      ! If qopt == .false. then it does not change darr1.
      !
      !	An integer array, iarr, containing ordinal numbers, 
      !	is sorted along with darr, providing an array of 
      !	the original indices of the values in the array.
      !
      !	argument declarations
      integer, intent(in) :: idim
      real(kind=our_dble), dimension(idim), intent(inout) :: darr1
      integer, dimension(idim), intent(out) :: iarr
      logical, intent(in) :: qopt
      type(error_type), intent(inout) :: err
      !	local declarations
      integer :: i, j, k, ir, ia, l, istk_ptr, stk_size
      integer, dimension(:), allocatable :: istk
      real(kind=our_dble) :: da
      real(kind=our_dble), dimension(idim) :: darr
      character(len=*), parameter :: subname = "qsort_double"
      ! begin
      answer = RETURN_FAIL
      l = 1
      istk_ptr = 0
      !	Allocate the stack
      stk_size = (abs(idim)*2)/MIN_PART
      allocate( istk(stk_size) )
      !	Fill the arrays
      do i=1,idim
         darr(i) = darr1(i)
         iarr(i) = i
      end do
      !	Start sorting
      ir = idim
      sort: do
         if(ir-l < MIN_PART) then
            do j=(l+1),ir
               da = darr(j)
               ia = iarr(j)
               do i=(j-1), l, -1
                  if(darr(i) <= da) exit
                  darr(i+1) = darr(i)
                  iarr(i+1) = iarr(i)
               end do
               darr(i+1) = da
               iarr(i+1) = ia
            end do
            if(istk_ptr == 0) exit sort
            ir = istk(istk_ptr)
            l = istk(istk_ptr-1)
            istk_ptr = istk_ptr - 2
         else
            !	Use k as the index of the next partitioning element
            k = (l + ir)/2
            call dswap( darr, k, l+1 )
            call iswap( iarr, k, l+1 )
            if(darr(l) > darr(ir)) then
               call dswap( darr, l, ir )
               call iswap( iarr, l, ir )
            endif
            if(darr(l+1) > darr(ir)) then
               call dswap( darr, l+1, ir )
               call iswap( iarr, l+1, ir )
            endif
            if(darr(l) > darr(l+1)) then
               call dswap( darr, l, l+1 )
               call iswap( iarr, l, l+1 )
            endif
            i= l + 1
            j= ir
            da = darr(l+1)
            ia = iarr(l+1)
            !	Begin innermost loop
            innermost: do
               countup: do
                  i = i + 1
                  if( darr(i) >= da ) exit countup
               end do countup
               countdown: do
                  j = j - 1
                  if( darr(j) <= da ) exit countdown
               end do countdown
               if(j < i) exit innermost
               call dswap(darr, i, j)
               call iswap(iarr, i, j)
            end do innermost
            !
            darr(l+1) = darr(j)
            iarr(l+1) = iarr(j)
            darr(j) = da
            iarr(j) = ia
            istk_ptr = istk_ptr + 2
            !	Deal with stack bounds errors
            if( istk_ptr > stk_size ) goto 700
            ! 
            if((ir-i+1) >= (j-l)) then
               istk(istk_ptr) = ir
               istk(istk_ptr-1) = i
               ir = j - 1
            else
               istk(istk_ptr) = j-1
               istk(istk_ptr-1) = l
               l = i
            end if
         end if
      end do sort
      !	Substitute values if requested
      if(qopt) then
         do i=1,idim
            darr1(i) = darr(i)
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
700   call err_handle(err, 1, &
           comment = "Error: stack size is too small")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   if( allocated( istk ) ) deallocate( istk )
   end function qsort_double
   !##################################################################
   integer(kind=our_int) function qsort_integer( darr1, iarr, idim, &
        qopt, err) result(answer)
      !	Performs a quicksort on the integer array darr. 
      !	If qopt ==.true. then it returns the sorted array.
      ! If qopt == .false. then it does not change darr1.
      !
      !	An integer array, iarr, containing ordinal numbers, 
      !	is sorted along with darr, providing an array of 
      !	the original indices of the values in the array.
      !
      !	argument declarations
      integer, intent(in) :: idim
      integer(kind=our_int), dimension(idim), intent(inout) :: darr1
      integer, dimension(idim), intent(out) :: iarr
      logical, intent(in) :: qopt
      type(error_type), intent(inout) :: err
      !	local declarations
      integer :: i, j, k, ir, ia, l, istk_ptr, stk_size
      integer, dimension(:), allocatable :: istk
      integer(kind=our_int) :: da
      integer(kind=our_int), dimension(idim) :: darr
      character(len=*), parameter :: subname = "qsort_integer"
      ! begin
      answer = RETURN_FAIL
      l = 1
      istk_ptr = 0
      !	Allocate the stack
      stk_size = (abs(idim)*2)/MIN_PART
      allocate( istk(stk_size) )
      !	Fill the arrays
      do i=1,idim
         darr(i) = darr1(i)
         iarr(i) = i
      end do
      !	Start sorting
      ir = idim
      sort: do
         if(ir-l < MIN_PART) then
            do j=(l+1),ir
               da = darr(j)
               ia = iarr(j)
               do i=(j-1), l, -1
                  if(darr(i) <= da) exit
                  darr(i+1) = darr(i)
                  iarr(i+1) = iarr(i)
               end do
               darr(i+1) = da
               iarr(i+1) = ia
            end do
            if(istk_ptr == 0) exit sort
            ir = istk(istk_ptr)
            l = istk(istk_ptr-1)
            istk_ptr = istk_ptr - 2
         else
            !	Use k as the index of the next partitioning element
            k = (l + ir)/2
            call iswap( darr, k, l+1 )
            call iswap( iarr, k, l+1 )
            if(darr(l) > darr(ir)) then
               call iswap( darr, l, ir )
               call iswap( iarr, l, ir )
            endif
            if(darr(l+1) > darr(ir)) then
               call iswap( darr, l+1, ir )
               call iswap( iarr, l+1, ir )
            endif
            if(darr(l) > darr(l+1)) then
               call iswap( darr, l, l+1 )
               call iswap( iarr, l, l+1 )
            endif
            i= l + 1
            j= ir
            da = darr(l+1)
            ia = iarr(l+1)
            !	Begin innermost loop
            innermost: do
               countup: do
                  i = i + 1
                  if( darr(i) >= da ) exit countup
               end do countup
               countdown: do
                  j = j - 1
                  if( darr(j) <= da ) exit countdown
               end do countdown
               if(j < i) exit innermost
               call iswap(darr, i, j)
               call iswap(iarr, i, j)
            end do innermost
            !
            darr(l+1) = darr(j)
            iarr(l+1) = iarr(j)
            darr(j) = da
            iarr(j) = ia
            istk_ptr = istk_ptr + 2
            !	Deal with stack bounds errors
            if( istk_ptr > stk_size ) goto 700
            ! 
            if((ir-i+1) >= (j-l)) then
               istk(istk_ptr) = ir
               istk(istk_ptr-1) = i
               ir = j - 1
            else
               istk(istk_ptr) = j-1
               istk(istk_ptr-1) = l
               l = i
            end if
         end if
      end do sort
      !	Substitute values if requested
      if(qopt) then
         do i=1,idim
            darr1(i) = darr(i)
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
700   call err_handle(err, 1, &
           comment = "Error: stack size is too small")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   if( allocated( istk ) ) deallocate( istk )
    end function qsort_integer
   !##################################################################
   subroutine sswap( sarr, i, j, ilen )
      !	Swaps elements i and j in string array sarr
      integer, intent(in) :: i,j,ilen
      character(len=ilen), dimension(:), intent(inout) :: sarr
      character(len=ilen) :: stemp
      stemp = sarr(i)
      sarr(i) = sarr(j)
      sarr(j) = stemp
   end subroutine sswap
   !##################################################################
   subroutine iswap( iarr, i, j )
      !	Swaps elements i and j in integer array iarr
      integer, intent(in) :: i,j
      integer, dimension(:), intent(inout) :: iarr
      integer :: itemp
      itemp = iarr(i)
      iarr(i) = iarr(j)
      iarr(j) = itemp
   end subroutine iswap
   !##################################################################
   subroutine dswap( darr, i, j )
      !	Swaps elements i and j in integer array iarr
      integer, intent(in) :: i,j
      real(kind=our_dble), dimension(:), intent(inout) :: darr
      real(kind=our_dble) :: dtemp
      dtemp = darr(i)
      darr(i) = darr(j)
      darr(j) = dtemp
   end subroutine dswap
   !##################################################################
   function tolower( sstr1, ilen ) result(sstr2)
      !	Returns the character string in all lower case
      integer, intent(in) :: ilen
      character(len=ilen), intent(in):: sstr1
      character(len=ilen) :: sstr2
      integer :: i, j
      ! begin
      do i=1,ilen
         !	Get the ASCII value
         j = iachar(sstr1(i:i))
         !	Is it upper case? (65-90)
         if( (j >64) .and. (j<91) ) then
            !  If so, find corresponding lower case value
            j = j + 32
         end if
         sstr2(i:i) = achar(j)
      end do
   end function tolower
   !##################################################################
!   function toupper( sstr1, ilen ) result(sstr2)
!      !	Returns the character string in all upper case
!      integer, intent(in) :: ilen
!      character(len=ilen), intent(in):: sstr1
!      character(len=ilen) :: sstr2
!      integer :: i, j
!      ! begin
!      do i=1,ilen
!         !	Get the ASCII value
!         j = iachar(sstr1(i:i))
!         !	Is it lower case? (97-122)
!         if( (j >96) .and. (j<123) ) then
!            !  If so, find corresponding upper case value
!            j = j - 32
!         end if
!         sstr2(i:i) = achar(j)
!      end do
!   end function toupper
   !##################################################################
end module quick_sort
!######################################################################
