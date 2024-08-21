!#####################################################################
module dynalloc
   ! Routines for allocating and deallocating pointers to arrays
   use error_handler
   use program_constants
   implicit none
   private ! by default
   public :: dyn_alloc, dyn_dealloc
   ! For allocating pointers to arrays
   interface dyn_alloc
      module procedure int1_alloc
      module procedure int2_alloc
      module procedure int3_alloc
      module procedure dbl1_alloc
      module procedure dbl2_alloc
      module procedure dbl3_alloc
      module procedure dbl4_alloc
      module procedure str1_alloc
      module procedure str2_alloc
      module procedure str3_alloc
      module procedure log1_alloc
      module procedure log2_alloc
      module procedure log3_alloc
   end interface
   ! For deallocating pointers to arrays
   interface dyn_dealloc
      module procedure int1_dealloc
      module procedure int2_dealloc
      module procedure int3_dealloc
      module procedure dbl1_dealloc
      module procedure dbl2_dealloc
      module procedure dbl3_dealloc
      module procedure dbl4_dealloc
      module procedure str1_dealloc
      module procedure str2_dealloc
      module procedure str3_dealloc
      module procedure log1_dealloc
      module procedure log2_dealloc
      module procedure log3_dealloc
   end interface
   ! parameters private to this module
   character(len=*), parameter :: modname = "dynalloc"
   !##################################################################
contains
   !##################################################################
   integer(our_int) function int1_alloc(intArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates an integer array of rank 1
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "int1_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( intArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int1_alloc
   !##################################################################
   integer(our_int) function int2_alloc(intArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates an integer array of rank 2
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "int2_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int2_alloc
   !##################################################################
   integer(our_int) function int3_alloc(intArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates an integer array of rank 3
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "int3_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int3_alloc
   !##################################################################
   integer(our_int) function dbl1_alloc(dblArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates a double-precision real array of rank 1
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "dbl1_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( dblArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl1_alloc
   !##################################################################
   integer(our_int) function dbl2_alloc(dblArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates a double-precision real array of rank 2
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "dbl2_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl2_alloc
   !##################################################################
   integer(our_int) function dbl3_alloc(dblArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates a double-precision real array of rank 3
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "dbl3_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl3_alloc
   !##################################################################
   integer(our_int) function dbl4_alloc(dblArray, dim1, dim2, dim3, &
        dim4, err, lbound1, lbound2, lbound3, lbound4) result(answer)
      ! Allocates a double-precision real array of rank 4
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
           lbound4
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4
      character(len=*), parameter :: subname = "dbl4_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl4_alloc
   !##################################################################
   integer(our_int) function str1_alloc(strArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates a character-string array of rank 1
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "str1_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( strArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str1_alloc
   !##################################################################
   integer(our_int) function str2_alloc(strArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates a character-string array of rank 2
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "str2_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str2_alloc
   !##################################################################
   integer(our_int) function str3_alloc(strArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates a character-string array of rank 3
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "str3_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str3_alloc
   !##################################################################
   integer(our_int) function log1_alloc(logArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates a logical array of rank 1
      implicit none
      ! declare required arguments
      logical, pointer :: logArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "log1_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( logArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log1_alloc
   !##################################################################
   integer(our_int) function log2_alloc(logArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates a logical array of rank 2
      implicit none
      ! declare required arguments
      logical, pointer :: logArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "log2_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( logArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log2_alloc
   !##################################################################
   integer(our_int) function log3_alloc(logArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates a logical array of rank 3
      implicit none
      ! declare required arguments
      logical, pointer :: logArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "log3_alloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( logArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
810   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log3_alloc
   !##################################################################
   integer(kind=our_int) function int1_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 1
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int1_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int1_dealloc
   !##################################################################
   integer(kind=our_int) function int2_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 2
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int2_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int2_dealloc
   !##################################################################
   integer(kind=our_int) function int3_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 3
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int3_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function int3_dealloc
   !##################################################################
   integer(kind=our_int) function dbl1_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 1
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl1_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl1_dealloc
   !##################################################################
   integer(kind=our_int) function dbl2_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 2
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl2_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl2_dealloc
   !##################################################################
   integer(kind=our_int) function dbl3_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 3
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl3dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl3_dealloc
   !##################################################################
   integer(kind=our_int) function dbl4_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 4
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl4dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function dbl4_dealloc
   !##################################################################
   integer(kind=our_int) function str1_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 1
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str1_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str1_dealloc
   !##################################################################
   integer(kind=our_int) function str2_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 2
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str2_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str2_dealloc
   !##################################################################
   integer(kind=our_int) function str3_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 3
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str3_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function str3_dealloc
   !##################################################################
   integer(kind=our_int) function log1_dealloc(logArray, err) &
        result(answer)
      ! Deallocates a logical array of rank 1
      implicit none
      ! declare arguments
      logical, pointer :: logArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "log1_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log1_dealloc
   !##################################################################
   integer(kind=our_int) function log2_dealloc(logArray, err) &
        result(answer)
      ! Deallocates a logical array of rank 2
      implicit none
      ! declare arguments
      logical, pointer :: logArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "log2_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log2_dealloc
   !##################################################################
   integer(kind=our_int) function log3_dealloc(logArray, err) &
        result(answer)
      ! Deallocates a logical array of rank 3
      implicit none
      ! declare arguments
      logical, pointer :: logArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "log3_dealloc"
      ! begin
      answer = RETURN_FAIL
      status = 0
      if( associated(logArray) ) deallocate(logArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function log3_dealloc
   !##################################################################
end module dynalloc
!#####################################################################
