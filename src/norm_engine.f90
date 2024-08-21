!#####################################################################
module norm_engine
   ! This module contains the computational routines for NORM 2.
   use error_handler
   use program_constants
   use dynalloc
   use quick_sort
   use matrix_methods
   use random_generator
   implicit none
   private ! by default
   public :: run_norm_engine_em, run_norm_engine_mcmc, &
        run_norm_engine_impute_random, run_norm_engine_impute_mean, &
        run_norm_engine_logpost
   ! parameters private to this module
   character(len=*), parameter :: modname = "norm_engine"
   ! private data types
   type :: workspace_type
      !   n    = number of cases or rows in dataset
      !   r    = number of response variables (i.e. columns of y)
      !   ntot = n*r = total size of y matrix
      !   p    = number of columns of the x-matrix
      !   case_order(n) = array giving order of cases
      !       after they have been sorted by missingness pattern
      !   npatt = number of distinct missingness patterns
      !   mvcode = missing-value code
      !   mis(npatt,r) = logical array storing the missingness
      !       patterns (T=missing, F=observed)
      !   first_case_in_patt(npatt) = integer array indicating
      !       the locations of the first case (after sorting) 
      !       in each missingness pattern
      !   last_case_in_patt(npatt) = integer array indicating
      !       the locations of the last case (after sorting) 
      !       in each missingness pattern
      !   n_in_patt(npatt) = integer array indicating
      !       the number of cases in each missingness pattern
      !   which_patt(n) = integer array indicating the missingness
      !       pattern to which each case belongs
      !   ysort(n,r) = sorted version of y-matrix
      !   yimp(n,r) = imputed version of y-matrix in original order
      !   xsort(n,p) = sorted version of x-matrix
      !   nobs(r) = no. of observed cases for each response variable
      !   ybar(r) = means of the observed responses
      !   ysdv(r) = stdevs of the observed responses
      !   beta(p,r) = regression coefficients
      !   sigma(r,r) = cov matrix for y
      !   oldbeta(p,r) = regression coefficients
      !   oldsigma(r,r) = cov matrix for y
      !   oldoldbeta(p,r) = regression coefficients
      !   oldoldsigma(r,r) = cov matrix for y
      !   sigma_ridge(r,r) = diagonal estimate of sigma for ridging
      !   ridge = logical flag if ridge prior is invoked
      !   prior_df = prior degrees of freedom (a.k.a. xi)
      !   prior_sscp(r,r) = prior sums of squares and cross products
      !      (a.k.a. lambda-inverse)
      !   wkp(p) = p x 1 workspace
      !   wkppA(p,p) = p x p workspace
      !   wkprA(p,p) = p x r workspace
      !   wkrrA(p,p) = r x r workspace
      !   wkrrB(p,p) = r x r workspace
      !   wkrprpA(r*p,r*p) = r*p x r*p workspace
      !   wkrpA(r*p) = r*p workspace
      !   wkrpB(r*p) = r*p workspace
      !   wknpA(n,p) = n x p workspace
      !   wknpB(n,p) = n x p workspace
      !   wknA(n) = n x 1 workspace
      !   wknB(n) = n x 1 workspace
      !   wknC(n) = n x 1 workspace
      !   wkpA(p) = p x 1 workspace
      !   wkpB(p) = p x 1 workspace
      !   wkpC(p) = p x 1 workspace
      !   iwkp(p) = p x 1 workspace (integer)
      !   lwkp(p) = p x 1 workspace (logical)
      !   xtxinv(p,p) = workspace for holding X^TX inverse
      !   xtxinvfac(p,p) = workspace for holding cholesky factor 
      !        of X^TX inverse
      !   wkrA(r) = workspace of length r
      !   epsteps(r,r) = r x r workspace for holding sum of
      !        ( Yi - t(beta)*xi) * t(yi - t(beta)*xi )
      !   yhat(n,r) = n x r workspace for holding X * beta
      !   eps(n,r) =  n x r workspace for holding y - yhat
      !   xty(p,r) =  p x p workspace for accumulating t(x) * y
      !   yty(r,r) =  r x r workspace for accumulating t(y) * y
      !   sigswp(r,r) = workspace for sweeping sigma
      !   swept(r) = logical array indicating which variables
      !       have been swept on
      !   loglik = current value of observed-data loglikelihood
      !   logpri = current value of log-prior density
      !   nparam = number of nonredundant model parameterr
      !          = p*r + r*(r+1)/2
      !   theta(nparam) = vectorized parameter
      !   oldtheta(nparam) = parameter from last iteration
      !   worst_weights(nparam) = weights for computing worst linear
      !          function of parameters
      !   reldiff = max. relative difference between oldtheta and theta
      ! These components are used to estimate the largest fraction
      ! of missing information, using the power method of Fraley et al.
      ! (2007)
      !      em_store_incr = iteration increment at which we update the
      !          stored value of theta (real)
      !      em_store_counter = integer
      !      next_em_store_iter = EM iteration for next stored theta
      !      em_store_iter = EM iteration for currently stored theta
      !      em_store_theta(nparam) = stored value of theta, 
      !         from approximately halfway along in the sequence
      !         of EM iterates
      !      next_em_store_theta(nparam) = next stored value of theta, 
      !         held in reserve until needed
      !      em_thetahat(nparam) = convergent value of theta from em
      !      em_worst_ok = did worst fraction procedure work? (T or F)
      !      uvec(nparam) = workspace
      !      uvec_new(nparam) = workspace
      !      vvec(nparam) = workspace
      !      vvec_new(nparam) = workspace
      !      worst_frac = worst fraction of missing information
      !      worst_linear_coef(nparam) = unit vector from worst-fraction
      !          procedure
      !  The following code cycles through the observations
      !  within patterns.
      ! do patt = 1, work%npatt
      !    do kase = work%first_case_in_patt(patt), &
      !         work%last_case_in_patt(patt)
      !       i = work%case_order(kase)  ! actual case number
      !       !  Now the yi vector can  be found both in ysort(kase,:)
      !       !  and in y(i,:)
      !       !  and the xi vector is in xsort(kase,:) and in x(i,:)
      !    end do
      ! end do
      !
      integer(kind=our_int) :: ntot=0, n=0, r=0, p=0, npatt=0, &
           nparam = 0
      integer(kind=our_int), pointer :: case_order(:) => null()
      real(kind=our_dble) :: mvcode = huge(0.D0)
      logical, pointer :: mis(:,:) => null(), lwkp(:) => null()
      integer(kind=our_int), pointer :: &
           first_case_in_patt(:) => null(), &
           last_case_in_patt(:) => null(), &
           n_in_patt(:) => null(), which_patt(:) => null(), &
           nobs(:) => null(), iwkp(:) => null()
      real(kind=our_dble), pointer :: ysort(:,:) => null(), &
           yimp(:,:) => null(), &
           xsort(:,:) => null(), ybar(:) => null(), &
           ysdv(:) => null(), &
           beta(:,:) => null(), sigma(:,:) => null(), &
           oldbeta(:,:) => null(), oldsigma(:,:) => null(), &
           oldoldbeta(:,:) => null(), oldoldsigma(:,:) => null(), &
           ratebeta(:,:) => null(), ratesigma(:,:) => null(), &
           sigma_ridge(:,:) => null(), &
           prior_sscp(:,:) => null(), &
           wkp(:) => null(), &
           wkppA(:,:) => null(), wkprA(:,:) => null(), &
           wkrrA(:,:) => null(), wkrrB(:,:) => null(), &
           wkrprpA(:,:) => null(), &
           wkrpA(:) => null(), wkrpB(:) => null(), &
           wknpA(:,:) => null(), wknpB(:,:) => null(), &
           wknA(:) => null(), wknB(:) => null(), wknC(:) => null(), &
           wkpA(:) => null(), wkpB(:) => null(), wkpC(:) => null(), &
           xtxinv(:,:) => null(), xtxinvfac(:,:) => null(), &
           wkrA(:) => null(), epsteps(:,:) => null(), &
           yhat(:,:) => null(), eps(:,:) => null(), &
           xty(:,:) => null(), yty(:,:) => null(), &
           sigswp(:,:) => null(), &
           theta(:) => null(), oldtheta(:) => null(), &
           worst_weights(:) => null()
      character(len=12) :: prior_type = ""
      logical, pointer :: swept(:) => null()
      real(kind=our_dble) :: loglik = 0.D0, logpri = 0.D0, &
           reldiff = 0.D0, prior_df = 0.D0
      real(kind=our_dble) :: em_store_incr = 0.D0
      integer(kind=our_int) :: em_store_counter = 0, &
           em_store_iter = 0, next_em_store_iter = 0
      real(kind=our_dble), pointer :: em_store_theta(:) => null(), &
           next_em_store_theta(:) => null(), em_thetahat(:) => null()
      logical :: em_worst_ok = .false.
      real(kind=our_dble), pointer :: uvec(:) => null(), &
           uvec_new(:) => null(), vvec(:) => null(), &
           vvec_new(:) => null()
      real(kind=our_dble) :: worst_frac = 0.D0
      real(kind=our_dble), pointer :: worst_linear_coef(:) => null()
   end type workspace_type
   !##################################################################
contains
   !##################################################################
   integer(kind=our_int) function nullify_workspace_type( work, &
        err) result(answer)
      ! Returns a workspace_type to its unititialized (null) state
      implicit none
      ! declare arguments
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = &
           "nullify_workspace_type"
      ! begin
      answer = RETURN_FAIL
      work%ntot = 0
      work%n = 0
      work%r = 0
      work%p = 0
      work%npatt = 0
      work%nparam = 0
      if( dyn_dealloc(work%case_order, err) == RETURN_FAIL ) goto 800
      work%mvcode = huge(0.D0)
      if( dyn_dealloc(work%mis, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%first_case_in_patt, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%last_case_in_patt, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%n_in_patt, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%which_patt, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%nobs, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%ysort, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%yimp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%xsort, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%ybar, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%ysdv, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%beta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%sigma, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%oldbeta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%oldsigma, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%oldoldbeta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%oldoldsigma, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%ratebeta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%ratesigma, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%sigma_ridge, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%prior_sscp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkppA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkprA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkrrA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkrrB, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkrprpA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkrpA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wknpA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wknpB, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wknA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wknB, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wknC, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkpA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkpB, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkpC, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%iwkp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%lwkp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%xtxinv, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%xtxinvfac, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%wkrA, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%epsteps, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%yhat, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%eps, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%xty, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%yty, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%sigswp, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%swept, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%theta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%oldtheta, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%worst_weights, err) == RETURN_FAIL ) goto 800
      work%prior_type = ""
      work%loglik = 0.D0
      work%logpri = 0.D0
      work%reldiff = 0.D0
      work%prior_df = 0.D0
      work%em_store_incr = 0.D0
      work%em_store_counter = 0
      work%em_store_iter = 0
      work%next_em_store_iter = 0
      if( dyn_dealloc(work%em_store_theta, err) &
           == RETURN_FAIL ) goto 800 
      if( dyn_dealloc(work%em_thetahat, err) &
           == RETURN_FAIL ) goto 800
      work%em_worst_ok = .false.
      if( dyn_dealloc(work%uvec, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%uvec_new, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%vvec, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(work%vvec_new, err) == RETURN_FAIL ) goto 800
      work%worst_frac = 0.D0
      if( dyn_dealloc(work%worst_linear_coef, err) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function nullify_workspace_type
   !##################################################################
   integer(kind=our_int) function run_norm_engine_em( &
        ! required inputs
        x, y, mvcode, &
        startval_present_int, beta_start, sigma_start, &
        prior_type, prior_df, prior_sscp, &
        max_iter, criterion, estimate_worst_int, &
        ! required outputs
        err, &
        ! optional outputs
        iter, converged_int, reldiff, &
        loglik, logpost, &
        beta, sigma, yimp, &
        npatt, mis_int, n_in_patt, &
        nobs, which_patt, &
        ybar, ysdv, &
        rate_beta, rate_sigma, &
        em_worst_ok_int, worst_frac, worst_linear_coef &
        ) result(answer)
      ! EM algorithm for ML estimation
      ! Required input arguments:
      !    x = n x p matrix of predictors (completely observed)
      !    y = n x r matrix of responses (missing values allowed)
      !    mvcode = missing-value code for y
      !    startval_present_int = integer flag indicating whether user
      !        is supplying starting values for parameters
      !    beta_start = starting value for beta
      !        (ignored if .not.startval_present, in which case
      !        default starting values are computed from the data)
      !    sigma_start = starting value for sigma
      !        (ignored if .not.startval_present, in which case
      !        default starting values are computed from the data)
      !    prior_type = character string: must be "uniform", 
      !        "jeffreys", "ridge" or "invwish"
      !           "uniform" maximizes the loglikelihood; other
      !           choices will maximize a log-posterior density
      !           "jeffreys" is the standard noniformative prior
      !           "invwish" is an inverted wishart. If this option
      !             is selected, then prior_df and prior_sscp must
      !             be provided.
      !           "ridge" is a data-determined prior.  If this is
      !            selected, then prior_df must be
      !             provided.
      !    max_iter = maximum number of iterations
      !        If this is not given, it defaults to the
      !        module parameter em_maxits_default
      !    criterion = relative convergence criterion
      !    estimate_worst_int = upon convergence, estimate the worst
      !        fraction of missing information? (T=1, F=0)
      ! Outputs (all are optional except err):
      !    err = object of type(error) to hold error message
      !    iter = number of iterations performed
      !    converged = T or F
      !    reldiff = maximum relative difference in elements of
      !       the parameter vector at the last two iterations
      !    loglik  = array of size max_iter giving loglikelihood at
      !              each iteration
      !    logpost = array of size max_iter giving log-posterior
      !              density at each iteration
      !    beta = p x r matrix of regression coefficients
      !    sigma = r x r covariance matrix
      !    yimp = imputed version of y before the last iteration
      !        (conditional-mean imputation)
      !    npatt = number of distinct missingness patterns
      !    mis_int = logical matrix of missingness patterns
      !    n_in_patt = number of observations per pattern
      !    nobs = numbet of observed values per variable
      !    which_patt = missing-data patterns to which the cases belong
      !    ybar = means of the observed values for each response variable
      !    ysdv = stdevs of the observed values for each response variable
      !    rate_beta = elementwise rates of convergence for
      !     beta, estimated from the last two iterations of EM
      !    rate_sigma = elementwise rates of convergence for
      !     sigma, estimated from the last two iterations of EM
      !    em_worst_ok_int = was procedure for estimating the worst
      !     fraction of missing information successful? (T/F)
      !    worst_frac = estimated worst fraction of missing info
      !    worst_linear_coef(nparam) = estimated unit vector
      !      corresponding to largest eigenvalue of rate matrix
      implicit none
      ! declare required inputs
      real(kind=our_dble), intent(in) :: x(:,:), y(:,:), mvcode
      integer(kind=our_int), intent(in) :: startval_present_int
      real(kind=our_dble), intent(in) :: beta_start(:,:), &
           sigma_start(:,:)
      character(len=*), intent(in) :: prior_type
      real(kind=our_dble), intent(in) :: prior_df, prior_sscp(:,:)
      integer(kind=our_int), intent(in) :: max_iter
      real(kind=our_dble), intent(in) :: criterion
      integer(kind=our_int), intent(in) :: estimate_worst_int
      ! declare required outputs
      type(error_type), intent(inout) :: err
      ! declare optional outputs
      integer(kind=our_int), optional, intent(out) :: iter
      integer(kind=our_int), optional, intent(out) :: converged_int
      real(kind=our_dble), optional, intent(out) :: reldiff
      real(kind=our_dble), optional, intent(out) :: loglik(:), logpost(:)
      real(kind=our_dble), optional, intent(out) :: &
           beta(:,:), sigma(:,:), yimp(:,:)
      integer(kind=our_int), optional, intent(out) :: npatt
      integer(kind=our_int), optional, intent(out) :: mis_int(:,:)
      integer(kind=our_int), optional, intent(out) :: n_in_patt(:)
      integer(kind=our_int), optional, intent(out) :: &
           nobs(:), which_patt(:)
      real(kind=our_dble), optional, intent(out) :: ybar(:), ysdv(:)
      real(kind=our_dble), optional, intent(out) :: rate_beta(:,:), &
           rate_sigma(:,:)
      integer(kind=our_int), optional, intent(out) :: em_worst_ok_int
      real(kind=our_dble), optional, intent(out) :: worst_frac
      real(kind=our_dble), optional, intent(out) :: worst_linear_coef(:)
      ! declare locals
      integer(kind=our_int) :: iter_local, i, j
      logical :: startval_present, estimate_worst, &
           converged_local, aborted
      real(kind=our_dble), pointer :: loglik_local(:)=>null(), &
           logpost_local(:)=>null()
      type(workspace_type) :: work
      integer(kind=our_int) :: ijunk
      character(len=12) :: sInt
      character(len=*), parameter :: subname = "run_norm_engine_em"
      ! begin
      answer = RETURN_FAIL
      sInt = "???"
      ! to enable text-based messaging, un-comment the lines where we
      ! write to sInt
      ! check args and prepare for EM
      if( size(y,1) /= size(x,1) ) goto 200
      !
      if( adjustl(prior_type) == "uniform" ) then
         work%prior_type = "uniform"
      else if( adjustl( prior_type ) == "jeffreys" ) then
         work%prior_type = "jeffreys"
      else if( adjustl( prior_type ) == "ridge" ) then
         work%prior_type = "ridge"
         if( prior_df < 0.D0 ) goto 223
      else if( adjustl( prior_type ) == "invwish" ) then
         work%prior_type = "invwish"
         if( ( size( prior_sscp, 1) /= size(y,2) ) .or. &
              ( size( prior_sscp, 2) /= size(y,2) ) ) goto 224
      else
         goto 228
      end if
      !
      if( max_iter < 0 ) goto 230
      !
      if( criterion < 0.D0 ) goto 240
      !
      if( startval_present_int /= 0 ) then
         startval_present = .true.
      else
         startval_present = .false.
      end if
      if( startval_present ) then
         if( ( size(beta_start,1) /= size(x,2) ) .or. &
              ( size(beta_start,2) /= size(y,2) ) ) goto 250
         if( ( size(sigma_start, 1) /= size(y,2) ) .or. &
             ( size(sigma_start, 2) /= size(y,2) ) ) goto 260
      end if
      if( estimate_worst_int /= 0 ) then
         estimate_worst = .true.
      else
         estimate_worst = .false.
      end if
      ! set up workspace
      work%n = size(y,1)
      work%r = size(y,2)
      work%ntot = work%n * work%r
      work%p = size(x,2)
      work%mvcode = mvcode
      if( sort_cases_by_missingness(y, work, err) == RETURN_FAIL ) goto 800
      if( allocate_workspace_items( work, err) == RETURN_FAIL ) goto 800
      if( make_xsort_and_ysort(x, y, work ) == RETURN_FAIL ) goto 800
      if( find_means_and_variances( work, err) == RETURN_FAIL ) goto 800
      if( regress_univariate( work, startval_present, &
           beta_start, sigma_start, err ) &
           == RETURN_FAIL ) goto 800
      if( make_xtxinv( work, err) == RETURN_FAIL ) goto 800
      if( work%prior_type == "uniform" ) then
         work%prior_df = real( -(work%r + 1), our_dble)
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "jeffreys" ) then
         work%prior_df = 0.D0
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "ridge" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_df * work%sigma_ridge(:,:)
      else if( work%prior_type == "invwish" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_sscp(:,:)         
      end if
      ! prepare for main iteration
      work%yimp(:,:) = y(:,:)  ! need to do this once at the
      ! beginning, because the e-step changes the missing values 
      ! but not the observed values
      iter_local = 0
      converged_local = .false.
      if( dyn_alloc( loglik_local, max_iter, err ) &
           == RETURN_FAIL ) goto 800
      loglik_local(:) = 0.D0
      if( dyn_alloc( logpost_local, max_iter, err ) &
           == RETURN_FAIL ) goto 800
      logpost_local(:) = 0.D0
      aborted = .false.
      ! for estimating worst fraction of missing information
      work%em_store_incr = 1.D0
      work%next_em_store_iter = 0
      work%em_store_iter = 0
      work%em_store_counter = 0
      ijunk = update_theta( work )
      work%em_store_theta(:) = work%theta(:) ! starting values
      work%next_em_store_theta(:) = work%theta(:)
      do
         if( max_iter == 0 ) exit
         iter_local = iter_local + 1
         ! write(sInt,"(I12)") iter_local
         sInt = adjustl(sInt)
         work%oldoldbeta(:,:) = work%oldbeta(:,:)
         work%oldoldsigma(:,:) = work%oldsigma(:,:)
         work%oldbeta(:,:) = work%beta(:,:)
         work%oldsigma(:,:) = work%sigma(:,:)
         if( run_estep( work, err ) == RETURN_FAIL) then
            aborted = .true.
            exit
         end if
         loglik_local( iter_local )  = work%loglik
         logpost_local( iter_local ) = work%loglik + work%logpri
         if( run_mstep( work, err ) == RETURN_FAIL) then
            aborted = .true.
            exit
         end if
         if( update_theta_and_oldtheta( work ) &
              == RETURN_FAIL ) then
            aborted = .true.
            exit
         end if
         ! for estimating worst frac miss inf
         work%em_store_counter = work%em_store_counter + 1
         if( real( work%em_store_counter, kind=our_dble) >= &
              work%em_store_incr ) then
            work%em_store_iter = work%next_em_store_iter
            work%em_store_theta(:) = work%next_em_store_theta(:)
            work%next_em_store_iter = iter_local
            work%next_em_store_theta(:) = work%theta(:)
            work%em_store_counter = 0
            work%em_store_incr = work%em_store_incr * 1.7D0
         end if
         if( iter_local > 1 ) then
            if( update_rates( work ) == RETURN_FAIL ) then
               aborted = .true.
               exit
            end if
         end if
         if( find_max_rel_diff( work ) == RETURN_FAIL) then
            aborted = .true.
            exit
         end if
         if( work%reldiff <= criterion ) &
              converged_local = .true.
         if( converged_local .or. &
              (iter_local >= max_iter) ) exit
      end do
      if( aborted ) then
         call err_handle(err, 1, &
              comment = "EM algorithm aborted" )
         call err_handle(err, 5, iiter = iter_local )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
         iter_local = iter_local - 1
      end if
      ! report the desired results, even if EM was aborted
      ! or didn't run at all
      if( present( iter ) ) iter = iter_local
      if( present( converged_int ) ) then
         if( converged_local ) then
            converged_int = 1
         else
            converged_int = 0
         end if
      end if
      if( present( reldiff ) ) reldiff = work%reldiff
      if( present( loglik ) ) then
         if( size( loglik ) /= max_iter) goto 300
         loglik(:) = loglik_local(:)  
      end if
      if( present( logpost ) ) then
         if( size( logpost ) /= max_iter ) goto 310
         logpost(:) = logpost_local(:)  
      end if
      if( present( beta ) ) then
         if( size( beta, 1) /= work%p ) goto 320
         if( size( beta, 2) /= work%r ) goto 320
         beta(:,:) = work%beta(:,:)
      end if
      if( present( sigma ) ) then
         if( size( sigma, 1 ) /= work%r ) goto 330
         if( size( sigma, 2 ) /= work%r ) goto 330
         sigma(:,:) = work%sigma(:,:)
      end if
      if( present( yimp ) ) then
         if( size( yimp, 1 ) /= work%n ) goto 340
         if( size( yimp, 2 ) /= work%r ) goto 340
         yimp(:,:) = work%yimp(:,:)
      end if
      if( present( npatt ) ) then
         npatt = work%npatt
      end if
      if( present( mis_int ) ) then
         if( size( mis_int, 1 ) /= work%n ) goto 350
         if( size( mis_int, 2 ) /= work%r ) goto 350
         mis_int(:,:) = 0
         do i = 1, work%npatt
            do j = 1, work%r
               if( work%mis(i,j) ) mis_int(i,j) = 1
            end do
         end do
      end if
      if( present( n_in_patt ) ) then
         if( size( n_in_patt ) /= work%n ) goto 360
         n_in_patt(:) = 0
         n_in_patt( 1:work%npatt ) = work%n_in_patt(:)
      end if
      if( present( nobs ) ) then
         if( size( nobs ) /= work%r ) goto 370
         nobs(:) = work%nobs(:)
      end if
      if( present( which_patt ) ) then
         if( size( which_patt ) /= work%n ) goto 380
         which_patt(:) = work%which_patt(:)
      end if
      if( present( ybar ) ) then
         if( size( ybar ) /= work%r ) goto 390
         ybar(:) = work%ybar(:)
      end if
      if( present( ysdv ) ) then
         if( size( ysdv ) /= work%r ) goto 400
         ysdv(:) = work%ysdv(:)
      end if
      if( present( rate_beta ) ) then
         if( size( rate_beta, 1 ) /= work%p ) goto 410
         if( size( rate_beta, 2 ) /= work%r ) goto 410
         if( ( iter_local > 2 ) .or. &
              ( ( iter_local == 2) .and. (.not. aborted ) ) ) then
            ! only report the rates if at least two iterations were
            ! successfully completed
            rate_beta(:,:) = work%ratebeta(:,:)
         end if
      end if
      if( present( rate_sigma ) ) then
         if( size( rate_sigma, 1 ) /= work%r ) goto 420
         if( size( rate_sigma, 2 ) /= work%r ) goto 420
         if( ( iter_local > 2 ) .or. &
              ( ( iter_local == 2) .and. (.not. aborted ) ) ) then
            ! only report the rates if at least two iterations were
            ! successfully completed
            rate_sigma(:,:) = work%ratesigma(:,:)
         end if
      end if
      !!! estimate worst fraction of missing information 
      work%em_worst_ok = .false.
      if ( all( work%nobs == work%n ) ) then
         ! no missing information
         work%em_worst_ok = .true.
         work%worst_frac = 0.D0
         work%worst_linear_coef(:) = 0.D0
      else
         ! use power method
         if( converged_local .and. estimate_worst ) then
            ijunk = estimate_worst_frac( work, err )
         end if
      end if
      if( present( em_worst_ok_int ) ) then
         if( work%em_worst_ok ) then
            em_worst_ok_int = 1
         else
            em_worst_ok_int = 0
         end if
      end if
      if( present( worst_frac ) ) then
         if( work%em_worst_ok ) then
            worst_frac = work%worst_frac
         else
            worst_frac = 0.D0
         end if
      end if
      if( present( worst_linear_coef ) ) then
         if( size( worst_linear_coef ) /= work%nparam ) goto 430
         if( work%em_worst_ok ) then
            worst_linear_coef(:) = work%worst_linear_coef(:)
         else
            worst_linear_coef(:) = 0.D0
         end if
      end if
      ! normal exit even if EM was aborted
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Dimensions of x and y not conformable" )
      goto 800
223   call err_handle(err, 1, &
           comment = "Prior_df cannot be negative for ridge prior")
      goto 800
224   call err_handle(err, 1, &
           comment = "Argument prior_sscp has incorrect dimensions")
      goto 800
228   call err_handle(err, 1, &
           comment = "Prior type not recognized")
      goto 800
230   call err_handle(err, 1, &
           comment = "Argument max_iter cannot be negative")
      goto 800
240   call err_handle(err, 1, &
           comment = "Argument criterion cannot be negative")
      goto 800
250   call err_handle(err, 1, &
           comment = "Incorrect dimensions for argument beta_start")
      goto 800
260   call err_handle(err, 1, &
           comment = "Incorrect dimensions for argument sigma_start")
      goto 800
300   call err_handle(err, 1, &
           comment = "Argument loglik has incorrect size")
      goto 800
310   call err_handle(err, 1, &
           comment = "Argument logpost has incorrect size")
      goto 800
320   call err_handle(err, 1, &
           comment = "Argument beta has incorrect shape")
      goto 800
330   call err_handle(err, 1, &
           comment = "Argument sigma has incorrect shape")
      goto 800
340   call err_handle(err, 1, &
           comment = "Argument yimp has incorrect shape")
      goto 800
350   call err_handle(err, 1, &
           comment = "Argument mis has incorrect shape")
      goto 800
360   call err_handle(err, 1, &
           comment = "Argument n_in_patt has incorrect size")
      goto 800
370   call err_handle(err, 1, &
           comment = "Argument nobs has incorrect size")
      goto 800
380   call err_handle(err, 1, &
           comment = "Argument which_patt has incorrect size")
      goto 800
390   call err_handle(err, 1, &
           comment = "Argument ybar has incorrect size")
      goto 800
400   call err_handle(err, 1, &
           comment = "Argument ysdv has incorrect size")
      goto 800
410   call err_handle(err, 1, &
           comment = "Argument rate_beta has incorrect shape")
      goto 800
420   call err_handle(err, 1, &
           comment = "Argument rate_sigma has incorrect shape")
      goto 800
430   call err_handle(err, 1, &
           comment = "Argument worst_linear_coef has incorrect size")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   ijunk = nullify_workspace_type(work, err)
   end function run_norm_engine_em
   !##################################################################
   integer(kind=our_int) function run_norm_engine_mcmc( &
        ! required inputs
        x, y, mvcode, &
        rand, &
        beta_start, sigma_start, &
        prior_type, prior_df, prior_sscp, &
        max_iter, multicycle, impute_every, &
        save_all_series_int, save_worst_series_int, worst_linear_coef, &
        ! required outputs
        err, &
        ! optional outputs
        iter, beta, sigma, loglik, logpost, yimp, &
        npatt, mis_int, n_in_patt, nobs, which_patt, ybar, ysdv, &
        beta_series, sigma_series, worst_series, &
        nimp, imp_list) &
        result(answer)
      ! MCMC algorithm
      ! Note: before calling this function, the random generator
      ! must be seeded.
      ! Required inputs:
      !    x = n x p matrix of predictors (completely observed)
      !    y = n x r matrix of responses (missing values allowed)
      !    mvcode = missing-value code for y
      !    rand = random generator object
      !    beta_start = starting value for beta
      !    sigma_start = starting value for sigma
      !    prior_type = character string: must be "uniform", 
      !        "jeffreys", "ridge" or "invwish"
      !           "uniform" maximizes the loglikelihood; other
      !           choices will maximize a log-posterior density
      !           "jeffreys" is the standard noniformative prior
      !           "invwish" is an inverted wishart. If this option
      !             is selected, then prior_df and prior_sscp must
      !             be provided.
      !           "ridge" is a data-determined prior.  If this is
      !            selected, then prior_df must be
      !             provided.
      !         If no prior_type is given, it defaults to
      !          a uniform prior.
      !    max_iter = number of MCMC iterations to perform
      !    multicycle = number of times to cycle the MCMC within
      !           each iteration (default = 1 )
      !    impute_every = no. of MCMC cycles between imputed datasets
      !           If multiple imputations are not desired, 
      !           then this argument should be set to zero.
      !    save_all_series = should full param series be saved? (T or F)
      !    save_worst_series = should worst linear function
      !           series be saved? (T or F)
      !    worst_linear_coef = weights determining the 
      !           worst linear function
      ! Outputs (all are optional except err) 
      !    err = object of type(error) to hold error message
      !    iter = number of iterations actually performed
      !          (will be equal to max_iter unless the algorithm
      !          was aborted)
      !    beta = p x r matrix of regression coefficients,
      !           last simulated value
      !           (if present upon input, used as starting value)
      !    sigma = r x r covariance matrix
      !           last simulated value
      !           (if present upon input, used as starting value)
      !    loglik  = array of size iter giving loglikelihood at
      !              each iteration
      !    logpost = array of size iter giving log-posterior
      !              density at each iteration
      !    yimp = imputed version of y, last simulated value
      !    npatt = number of distinct missingness patterns
      !    mis = logical matrix of missingness patterns
      !    n_in_patt = number of observations per pattern
      !    nobs = numbet of observed values per variable
      !    which_patt = missing-data patterns to which the cases belong
      !    ybar = means of the observed values for each response variable
      !    ysdv = stdevs of the observed values for each response variable
      !    beta_series =  array of size p x r x iter
      !           returning all the simulated values of beta
      !    sigma_series =  array of size r x r x iter
      !           returning all the simulated values of sigma
      !    worst_series = output array of length iter
      !           returning all the simulated values of the
      !           worst linear function
      !    nimp = number of imputations saved in imp_list
      !         (determined by max_iter and impute_every)
      !    imp_list = n x r x nimp array holding all imputed datasets
      implicit none
      ! declare required inputs
      real(kind=our_dble), intent(in) :: x(:,:), y(:,:), mvcode
      type(random_gendata), intent(inout) :: rand
      real(kind=our_dble), intent(in) :: beta_start(:,:), &
           sigma_start(:,:)
      character(len=*), intent(in) :: prior_type
      real(kind=our_dble), intent(in) :: prior_df, &
           prior_sscp(:,:)
      integer(kind=our_int), intent(in) :: max_iter
      integer(kind=our_int), intent(in) :: multicycle
      integer(kind=our_int), intent(in) :: impute_every
      integer(kind=our_int), intent(in) :: save_all_series_int, &
           save_worst_series_int
      real(kind=our_dble), intent(in) :: worst_linear_coef(:)
      ! declare required outputs
      type(error_type), intent(inout) :: err
      ! declare optional outputs
      integer(kind=our_int), optional, intent(out) :: iter
      real(kind=our_dble), optional, intent(out) :: beta(:,:), sigma(:,:), &
           yimp(:,:), loglik(:), logpost(:)
      integer, optional, intent(out) :: npatt
      integer(kind=our_int), optional, intent(out) :: mis_int(:,:)
      integer(kind=our_int), optional, intent(out) :: n_in_patt(:)
      integer(kind=our_int), optional, intent(out) :: &
           nobs(:), which_patt(:)
      real(kind=our_dble), optional, intent(out) :: ybar(:), ysdv(:)
      real(kind=our_dble), optional, intent(out) :: beta_series(:,:,:), &
           sigma_series(:,:,:), worst_series(:)
      integer(kind=our_int), optional, intent(out) :: nimp
      real(kind=our_dble), optional, intent(out) :: imp_list(:,:,:)
      ! declare locals
      integer(kind=our_int) :: iter_local, iter_complete
      integer(kind=our_int) :: ijunk, nimps, impno, r, p, nparam, &
           which_cycle
      real (kind=our_dble), pointer :: loglik_local(:)=>null(), &
           logpost_local(:)=>null(), &
           beta_series_local(:,:,:)=>null(), &
           sigma_series_local(:,:,:)=>null(), &
           worst_series_local(:)=>null(), &
           imp_list_local(:,:,:)=>null()
      type(workspace_type) :: work
      character(len=12) :: sInt, sInt2
      logical :: aborted, save_all_series, save_worst_series
      real(kind=our_dble) :: num, denom
      character(len=*), parameter :: subname = "run_norm_engine_mcmc"
      integer(our_int) :: i,j,k
      ! begin
      answer = RETURN_FAIL
      sInt = "???"
      sInt2 = "???"
      save_all_series = .false.
      if( save_all_series_int /= 0 ) save_all_series = .true.
      save_worst_series = .false.
      if( save_worst_series_int /= 0 ) save_worst_series = .true.
      ! check args
      if( size(y,1) /= size(x,1) ) goto 200
      !
      if( ( size(beta_start,1) /= size(x,2) ) .or. &
           ( size(beta_start,2) /= size(y,2) ) ) goto 250
      !
      if( ( size(sigma_start, 1) /= size(y,2) ) .or. &
           ( size(sigma_start, 2) /= size(y,2) ) ) goto 260
      !
      if( adjustl(prior_type) == "uniform" ) then
         work%prior_type = "uniform"
      else if( adjustl( prior_type ) == "jeffreys" ) then
         work%prior_type = "jeffreys"
      else if( adjustl( prior_type ) == "ridge" ) then
         work%prior_type = "ridge"
         if( prior_df < 0.D0 ) goto 223
      else if( adjustl( prior_type ) == "invwish" ) then
         work%prior_type = "invwish"
         if( ( size( prior_sscp, 1) /= size(y,2) ) .or. &
              ( size( prior_sscp, 2) /= size(y,2) ) ) goto 224
      else
         goto 228
      end if
      !
      if( max_iter < 1 ) goto 230
      if( multicycle < 1 ) goto 235
      if( impute_every < 0 ) goto 240
      if( impute_every > max_iter ) goto 245
      !
      p = size(x,2)
      r = size(y,2)
      nparam = p*r + r*(r+1)/2
      !
      if( save_worst_series ) then
         if( size( worst_linear_coef ) /= nparam ) goto 270
      end if
      ! set up workspace
      work%n = size(y,1)
      work%r = size(y,2)
      work%ntot = work%n * work%r
      work%p = size(x,2)
      work%mvcode = mvcode
      if( sort_cases_by_missingness(y, work, err) == RETURN_FAIL ) goto 800
      if( allocate_workspace_items( work, err) == RETURN_FAIL ) goto 800
      if( make_xsort_and_ysort(x, y, work ) == RETURN_FAIL ) goto 800
      if( find_means_and_variances( work, err) == RETURN_FAIL ) goto 800
      if( regress_univariate( work, .true., &
           beta_start, sigma_start, err ) &
           == RETURN_FAIL ) goto 800
      if( make_xtxinv( work, err) == RETURN_FAIL ) goto 800
      if( work%prior_type == "uniform" ) then
         work%prior_df = real( -(work%r + 1), our_dble)
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "jeffreys" ) then
         work%prior_df = 0.D0
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "ridge" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_df * work%sigma_ridge(:,:)
      else if( work%prior_type == "invwish" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_sscp(:,:)         
      end if
      work%yimp(:,:) = y(:,:)  ! need to do this once at the
      !  beginning, because the e/i-step changes the missing values 
      ! but not the observed values
      !!!!!!!
      ! allocate space for series and imputations
      if( present( beta_series ) .and. save_all_series ) then
         if( dyn_alloc( beta_series_local, work%p, work%r, &
              max_iter, err ) == RETURN_FAIL ) goto 800
         beta_series_local(:,:,:) = 0.D0
      end if
      if( present( sigma_series ) .and. save_all_series ) then
         if( dyn_alloc( sigma_series_local, work%r, work%r, &
              max_iter, err ) == RETURN_FAIL ) goto 800
         sigma_series_local(:,:,:) = 0.D0
      end if
      if( present( worst_series  ) .and. save_worst_series ) then
         if( dyn_alloc( worst_series_local, &
              max_iter, err ) == RETURN_FAIL ) goto 800
         worst_series_local(:) = 0.D0
         ! store normalized coefficients
         denom = sqrt( sum( worst_linear_coef(:)**2 ) )
         if( denom /= 0.D0 ) then
            work%worst_weights(:) = worst_linear_coef(:) / denom
         else
            work%worst_weights(:) = 0.D0
         end if
      end if
      if( present( imp_list ) .and. ( impute_every > 0 ) ) then
         nimps = max_iter / impute_every  ! integer division
         if( dyn_alloc( imp_list_local, work%n, work%r, nimps, err ) &
              == RETURN_FAIL ) goto 800
         imp_list_local(:,:,:) = 0.D0
      end if
      ! set starting values
      work%beta(:,:) = beta_start(:,:)
      work%sigma(:,:) = sigma_start(:,:)
      ! store staring values in oldtheta for all iterations
      if( update_theta( work ) == RETURN_FAIL) goto 800
      work%oldtheta(:) = work%theta(:)
      ! prepare for mcmc iteration
      if( dyn_alloc( loglik_local, max_iter, err ) &
           == RETURN_FAIL ) goto 800
      loglik_local(:) = 0.D0
      if( dyn_alloc( logpost_local, max_iter, err ) &
           == RETURN_FAIL ) goto 800
      logpost_local(:) = 0.D0
      aborted = .false.
      impno = 0
      iter_complete = 0
      do iter_local = 1, max_iter
         ! write(sInt,"(I12)") iter_local
         sInt = adjustl(sInt)
         do which_cycle = 1, multicycle
            ! write(sInt2,"(I12)") which_cycle
            sInt2 = adjustl(sInt2)
            if( run_istep( work, rand, err ) == RETURN_FAIL) then
               aborted = .true.
               exit
            end if
            if( which_cycle == 1 ) then
               loglik_local( iter_local )  = work%loglik
               logpost_local( iter_local ) = work%loglik + work%logpri
            end if
            if( run_pstep( work, rand, err) == RETURN_FAIL) then
               aborted = .true.
               call err_handle(err, 1, &
                    comment = "Posterior distribution may not be proper")
               exit
            end if
            if( which_cycle == multicycle ) then
               ! save the relevant results
               if( present( beta_series ) .and. save_all_series ) then
                  beta_series_local(:,:, iter_local) = work%beta(:,:)
               end if
               if( present( sigma_series ) .and. save_all_series ) then
                  sigma_series_local(:,:, iter_local) = work%sigma(:,:)
               end if
               if( present( worst_series ) .and. save_worst_series ) then
                  if( update_theta( work ) == RETURN_FAIL) then
                     aborted = .true.
                     exit
                  end if
                  num = sum( work%worst_weights(:) * work%theta(:) )
                  denom = sqrt( sum( work%theta(:)**2 ) )
                  if( denom /= 0.D0 ) then
                     worst_series_local( iter_local ) = num / denom
                  else
                     aborted = .true.
                     call err_handle(err, 1, &
                          comment = "Attempted division by zero" )
                      exit
                  end if
               end if
               if( ( impute_every > 0 ) .and. present( imp_list ) ) then
                  if( mod( iter_local, impute_every ) == 0 ) then
                     impno = impno + 1
                     imp_list_local(:,:,impno) = work%yimp(:,:)
                  end if
               end if
            end if
         end do
         if( aborted ) exit
         iter_complete = iter_local
      end do
      if( aborted ) then
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 6, iiter = iter_local, icycle = which_cycle )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      ! report the desired results, even if MCMC was aborted
      ! or didn't run at all
      if( present( iter ) ) iter = iter_complete
      if( present( loglik ) ) then
         if( size( loglik ) /= max_iter ) goto 330
         loglik(:) = loglik_local(:)
      end if
      if( present( logpost ) ) then
         if( size( logpost ) /= max_iter ) goto 340
         logpost(:) = logpost_local(:)
      end if
      if( present( beta ) ) then
         if( size( beta, 1 ) /= work%p ) goto 300
         if( size( beta, 2 ) /= work%r ) goto 300
         beta(:,:) = work%beta(:,:)
      end if
      if( present( sigma ) ) then
         if( size( sigma, 1 ) /= work%r ) goto 310
         if( size( sigma, 2 ) /= work%r ) goto 310
         sigma(:,:) = work%sigma(:,:)
      end if
      if( present( yimp ) ) then
         if( size( yimp, 1 ) /= work%n ) goto 320
         if( size( yimp, 2 ) /= work%r ) goto 320
         yimp(:,:) = work%yimp(:,:)
      end if
      if( present( npatt ) ) then
         npatt = work%npatt
      end if
      if( present( mis_int ) ) then
         if( size( mis_int, 1 ) /= work%n ) goto 350
         if( size( mis_int, 2 ) /= work%r ) goto 350
         mis_int(:,:) = 0
         do i = 1, work%npatt
            do j = 1, work%r
               if( work%mis(i,j) ) mis_int(i,j) = 1
            end do
         end do
      end if
      if( present( n_in_patt ) ) then
         if( size( n_in_patt ) /= work%n ) goto 360
         n_in_patt(:) = 0
         n_in_patt( 1:work%npatt ) = work%n_in_patt( 1:work%npatt )
      end if
      if( present( nobs ) ) then
         if( size( nobs ) /= work%r ) goto 370
         nobs(:) = work%nobs(:)
      end if
      if( present( which_patt ) ) then
         if( size( which_patt ) /= work%n ) goto 380
         which_patt(:) = work%which_patt(:)
      end if
      if( present( ybar ) ) then
         if( size( ybar ) /= work%r ) goto 390
         ybar(:) = work%ybar(:)
      end if
      if( present( ysdv ) ) then
         if( size( ysdv ) /= work%r ) goto 400
         ysdv(:) = work%ysdv(:)
      end if
      if( present( beta_series ) .and. save_all_series ) then
         if( size( beta_series, 1 ) /= work%p ) goto 410
         if( size( beta_series, 2 ) /= work%r ) goto 410
         if( size( beta_series, 3 ) < iter_complete ) goto 410
         beta_series(:,:,:) = 0.D0
         do k = 1, iter_complete
            do j = 1, work%r
               do i = 1, work%p
                  beta_series(i,j,k) = beta_series_local(i,j,k)
               end do
            end do
         end do
      end if
      if( present( sigma_series ) .and. save_all_series ) then
         if( size( sigma_series, 1 ) /= work%r ) goto 420
         if( size( sigma_series, 2 ) /= work%r ) goto 420
         if( size( sigma_series, 3 ) < iter_complete ) goto 420
         sigma_series(:,:,:) = 0.D0
         do k = 1, iter_complete
            do j = 1, work%r
               do i = 1, work%r
                  sigma_series(i,j,k) = sigma_series_local(i,j,k)
               end do
            end do
         end do
      end if
      if( present( worst_series ) .and. save_worst_series ) then
         if( size( worst_series ) /= max_iter ) goto 430
         worst_series(:) = worst_series_local(:)
      end if
      if( present( nimp ) ) nimp = impno
      if( present( imp_list ) .and. ( impno > 0 ) ) then
         if( size( imp_list, 1 ) /= work%n ) goto 440
         if( size( imp_list, 2 ) /= work%r ) goto 440
         if( size( imp_list, 3 ) < impno ) goto 440
         imp_list(:,:,:) = 0.D0
         do k = 1, impno
            do j = 1, work%r
               do i = 1, work%n
                  imp_list(i,j,k) = imp_list_local(i,j,k)
               end do
            end do
         end do
      end if
      ! normal exit even if MCMC was aborted
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Dimensions of x and y not conformable" )
      goto 800
223   call err_handle(err, 1, &
           comment = "Prior_df cannot be negative for ridge prior")
      goto 800
224   call err_handle(err, 1, &
           comment = "Argument prior_sscp has incorrect dimensions")
      goto 800
228   call err_handle(err, 1, &
           comment = "Prior type not recognized")
      goto 800
230   call err_handle(err, 1, &
           comment = "Argument max_iter must be positive")
      goto 800
235   call err_handle(err, 1, &
           comment = "Argument multicycle must be positive")
      goto 800
240   call err_handle(err, 1, &
           comment = "Argument impute_every cannot be negative")
      goto 800
245   call err_handle(err, 1, &
           comment = "Argument impute_every cannot exceed max_iter")
      goto 800
250   call err_handle(err, 1, &
           comment = "Incorrect dimensions for argument beta_start")
      goto 800
260   call err_handle(err, 1, &
           comment = "Incorrect dimensions for argument sigma_start")
      goto 800
270   call err_handle(err, 1, &
           comment = "Incorrect size for argument worst_linear_coef")
      goto 800
300   call err_handle(err, 1, &
           comment = "Argument beta has incorrect shape")
      goto 800
310   call err_handle(err, 1, &
           comment = "Argument sigma has incorrect shape")
      goto 800
320   call err_handle(err, 1, &
           comment = "Argument yimp has incorrect shape")
      goto 800
330   call err_handle(err, 1, &
           comment = "Argument loglik has incorrect size")
      goto 800
340   call err_handle(err, 1, &
           comment = "Argument logpost has incorrect size")
      goto 800
350   call err_handle(err, 1, &
           comment = "Argument mis_int has incorrect size")
      goto 800
360   call err_handle(err, 1, &
           comment = "Argument n_in_patt has incorrect size")
      goto 800
370   call err_handle(err, 1, &
           comment = "Argument nobs has incorrect size")
      goto 800
380   call err_handle(err, 1, &
           comment = "Argument which_patt has incorrect size")
      goto 800
390   call err_handle(err, 1, &
           comment = "Argument ybar has incorrect size")
      goto 800
400   call err_handle(err, 1, &
           comment = "Argument ysdv has incorrect size")
      goto 800
410   call err_handle(err, 1, &
           comment = "Argument beta_series has incorrect shape")
      goto 800
420   call err_handle(err, 1, &
           comment = "Argument sigma_series has incorrect shape")
      goto 800
430   call err_handle(err, 1, &
           comment = "Argument worst_series has incorrect size")
      goto 800
440   call err_handle(err, 1, &
           comment = "Argument imp_list has incorrect shape")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   ijunk = nullify_workspace_type(work, err)
      ijunk = dyn_dealloc( loglik_local, err )
      ijunk = dyn_dealloc( logpost_local, err )
      ijunk = dyn_dealloc( beta_series_local, err )
      ijunk = dyn_dealloc( sigma_series_local, err )
      ijunk = dyn_dealloc( worst_series_local, err )
      ijunk = dyn_dealloc( imp_list_local, err )
   end function run_norm_engine_mcmc
   !##################################################################
   integer(kind=our_int) function run_norm_engine_impute_random(x, &
        y, mvcode, beta, sigma, yimp, rand, err) result(answer)
      ! Impute data from fixed parameters
      ! Required inputs:
      !    x = n x p matrix of predictors (completely observed)
      !    y = n x r matrix of responses (missing values allowed)
      !    mvcode = missing-value code for y
      !    iseed1 = first integer seed for random number generator
      !    iseed2 = second integer seed for random number generator
      !    beta = p x r matrix of regression coefficients,
      !    sigma = r x r covariance matrix
      ! Outputs:
      !    yimp = imputed version of y
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: x(:,:), y(:,:), mvcode
      real(kind=our_dble), intent(in) :: beta(:,:), sigma(:,:)
      ! declare outputs
      real(kind=our_dble), pointer :: yimp(:,:)
      type(random_gendata), intent(inout) :: rand
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: ijunk
      type(workspace_type) :: work
      character(len=*), parameter :: &
           subname = "run_norm_engine_impute_random"
      ! begin
      answer = RETURN_FAIL
      ! check args
      if( size(y,1) /= size(x,1) ) goto 200
      if( ( size(beta,1) /= size(x,2) ) .or. &
           ( size(beta,2) /= size(y,2) ) ) goto 250
      if( ( size(sigma,1) /= size(y,2) ) .or. &
           ( size(sigma,2) /= size(y,2) ) ) goto 260
      ! set up workspace
      work%n = size(y,1)
      work%r = size(y,2)
      work%ntot = work%n * work%r
      work%p = size(x,2)
      work%mvcode = mvcode
      if( sort_cases_by_missingness(y, work, err) == RETURN_FAIL ) goto 800
      if( allocate_workspace_items( work, err) == RETURN_FAIL ) goto 800
      if( make_xsort_and_ysort(x, y, work ) == RETURN_FAIL ) goto 800
      if( find_means_and_variances( work, err) == RETURN_FAIL ) goto 800
      if( make_xtxinv( work, err) == RETURN_FAIL ) goto 800
      work%yimp(:,:) = y(:,:)  ! need to do this once at the
      !  beginning, because the e/i-step changes the missing values 
      ! but not the observed values
      work%beta(:,:) = beta(:,:)
      work%sigma(:,:) = sigma(:,:)
      ! impute
      if( run_istep( work, rand, err ) == RETURN_FAIL) goto 700
      ! allocate space for results
      if( dyn_alloc( yimp, work%n, work%r, err ) == RETURN_FAIL ) goto 800
      yimp(:,:) = work%yimp(:,:)
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Dimensions of x and y not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
250   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument beta")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
260   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument sigma")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
700   call err_handle(err, 1, &
           comment = "Imputation procedure aborted")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ! final cleanup
999   ijunk = nullify_workspace_type(work, err)
   end function run_norm_engine_impute_random
   !##################################################################
   integer(kind=our_int) function run_norm_engine_impute_mean( x, &
        y, mvcode, beta, sigma, yimp, err) result(answer)
      ! Impute data from fixed parameters (conditional mean imputation)
      ! Required inputs:
      !    x = n x p matrix of predictors (completely observed)
      !    y = n x r matrix of responses (missing values allowed)
      !    mvcode = missing-value code for y
      !    beta = p x r matrix of regression coefficients,
      !    sigma = r x r covariance matrix
      ! Outputs:
      !    yimp = imputed version of y
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: x(:,:), y(:,:), mvcode
      real(kind=our_dble), intent(in) :: beta(:,:), sigma(:,:)
      ! declare outputs
      real(kind=our_dble), pointer :: yimp(:,:)
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: ijunk
      type(workspace_type) :: work
      character(len=*), parameter :: &
           subname = "run_norm_engine_impute_mean"
      ! begin
      answer = RETURN_FAIL
      ! check args
      if( size(y,1) /= size(x,1) ) goto 200
      if( ( size(beta,1) /= size(x,2) ) .or. &
           ( size(beta,2) /= size(y,2) ) ) goto 250
      if( ( size(sigma,1) /= size(y,2) ) .or. &
           ( size(sigma,2) /= size(y,2) ) ) goto 260
      ! set up workspace
      work%n = size(y,1)
      work%r = size(y,2)
      work%ntot = work%n * work%r
      work%p = size(x,2)
      work%mvcode = mvcode
      if( sort_cases_by_missingness(y, work, err) == RETURN_FAIL ) goto 800
      if( allocate_workspace_items( work, err) == RETURN_FAIL ) goto 800
      if( make_xsort_and_ysort(x, y, work ) == RETURN_FAIL ) goto 800
      if( find_means_and_variances( work, err) == RETURN_FAIL ) goto 800
      if( make_xtxinv( work, err) == RETURN_FAIL ) goto 800
      work%yimp(:,:) = y(:,:)  ! need to do this once at the
      !  beginning, because the e/i-step changes the missing values 
      ! but not the observed values
      work%beta(:,:) = beta(:,:)
      work%sigma(:,:) = sigma(:,:)
      ! impute
      if( run_estep( work, err ) == RETURN_FAIL) goto 700
      ! allocate space for results
      if( dyn_alloc( yimp, work%n, work%r, err ) == RETURN_FAIL ) goto 800
      yimp(:,:) = work%yimp(:,:)
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Dimensions of x and y not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
250   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument beta")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
260   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument sigma")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
700   call err_handle(err, 1, &
           comment = "Imputation procedure aborted")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ! final cleanup
999   ijunk = nullify_workspace_type(work, err)
   end function run_norm_engine_impute_mean
!##################################################################
   integer(kind=our_int) function run_norm_engine_logpost( x, &
        y, mvcode, beta, sigma, logpost, err, &
        prior_type, prior_df, prior_sscp) result(answer)
      ! Evaluates the log-posterior density function
      ! Required inputs:
      !    x = n x p matrix of predictors (completely observed)
      !    y = n x r matrix of responses (missing values allowed)
      !    mvcode = missing-value code for y
      !    beta = p x r matrix of regression coefficients,
      !    sigma = r x r covariance matrix
      !    prior_type = character string: must be "uniform", 
      !        "jeffreys", "ridge" or "invwish"
      !           "uniform" computes the loglikelihood; other
      !           choices will compute a log-posterior density
      !           "jeffreys" is the standard noniformative prior
      !           "invwish" is an inverted wishart. If this option
      !             is selected, then prior_df and prior_sscp must
      !             be provided.
      !           "ridge" is a data-determined prior.  If this is
      !            selected, then prior_df must be
      !             provided.
      !         If no prior_type is given, it defaults to
      !          a uniform prior.
      ! Outputs:
      !    logpost = answer
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: x(:,:), y(:,:), mvcode
      real(kind=our_dble), intent(in) :: beta(:,:), sigma(:,:)
      ! declare outputs
      real(kind=our_dble), intent(out) :: logpost
      type(error_type), intent(inout) :: err
      character(len=*), intent(in) :: prior_type
      real(kind=our_dble), intent(in) :: prior_df, prior_sscp(:,:)
      ! declare locals
      type(workspace_type) :: work
      character(len=*), parameter :: subname = "run_norm_engine_logpost"
      ! begin
      answer = RETURN_FAIL
      ! check args
      if( size(y,1) /= size(x,1) ) goto 200
      !
      if( adjustl(prior_type) == "uniform" ) then
         work%prior_type = "uniform"
      else if( adjustl( prior_type ) == "jeffreys" ) then
         work%prior_type = "jeffreys"
      else if( adjustl( prior_type ) == "ridge" ) then
         work%prior_type = "ridge"
         if( prior_df < 0.D0 ) goto 223
      else if( adjustl( prior_type ) == "invwish" ) then
         work%prior_type = "invwish"
         if( ( size( prior_sscp, 1) /= size(y,2) ) .or. &
              ( size( prior_sscp, 2) /= size(y,2) ) ) goto 224
      else
         goto 228
      end if
      !
      if( ( size(beta,1) /= size(x,2) ) .or. &
           ( size(beta,2) /= size(y,2) ) ) goto 250
      if( ( size(sigma,1) /= size(y,2) ) .or. &
           ( size(sigma,2) /= size(y,2) ) ) goto 260
      ! set up workspace
      work%n = size(y,1)
      work%r = size(y,2)
      work%ntot = work%n * work%r
      work%p = size(x,2)
      work%mvcode = mvcode
      if( sort_cases_by_missingness(y, work, err) == RETURN_FAIL ) goto 800
      if( allocate_workspace_items( work, err) == RETURN_FAIL ) goto 800
      if( make_xsort_and_ysort(x, y, work ) == RETURN_FAIL ) goto 800
      if( find_means_and_variances( work, err) == RETURN_FAIL ) goto 800
      if( regress_univariate( work, .true., beta, sigma, err ) &
           == RETURN_FAIL ) goto 800
      if( make_xtxinv( work, err) == RETURN_FAIL ) goto 800
      if( work%prior_type == "uniform" ) then
         work%prior_df = real( -(work%r + 1), our_dble)
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "jeffreys" ) then
         work%prior_df = 0.D0
         work%prior_sscp(:,:) = 0.D0
      else if( work%prior_type == "ridge" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_df * work%sigma_ridge(:,:)
      else if( work%prior_type == "invwish" ) then
         work%prior_df = prior_df
         work%prior_sscp(:,:) = prior_sscp(:,:)         
      end if
      work%yimp(:,:) = y(:,:)  ! need to do this once at the
      !  beginning, because the e/i-step changes the missing values 
      ! but not the observed values
      work%beta(:,:) = beta(:,:)
      work%sigma(:,:) = sigma(:,:)
      ! impute
      if( run_estep( work, err ) == RETURN_FAIL) goto 700
      logpost = work%loglik + work%logpri
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Dimensions of x and y not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
223   call err_handle(err, 1, &
           comment = "Prior_df cannot be negative for ridge prior")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
224   call err_handle(err, 1, &
           comment = "Argument prior_sscp has incorrect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
228   call err_handle(err, 1, &
           comment = "Prior type not recognized")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
250   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument beta")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
260   call err_handle(err, 1, &
           comment = "Incorrect dimensions for array argument sigma")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
700   call err_handle(err, 1, &
           comment = "Computation of loglikelihood aborted")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
      999 continue
   end function run_norm_engine_logpost
   !##################################################################
   integer(kind=our_int) function sort_cases_by_missingness(y, work, &
        err) result(answer)
      implicit none
      ! declare args
      real(kind=our_dble), intent(in) :: y(:,:)
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      ! Notice that sort_strings is an array of len=1 strings of
      ! size ntot=n*r. But we will pass it to a subroutine
      ! "fill sort_strings_and_sort" which treats it as len=r
      ! and dimension n.
      ! The reason we do this is that allocation of deferred-length 
      ! character strings is a Fortran-2003 feature, not yet supported
      ! by some compilers including gfortran.
      character(len=*), parameter :: subname = &
           "sort_cases_by_missingness"
      character(len=1), allocatable :: sort_strings(:)
      integer(kind=our_int) :: status
      ! begin
      answer = RETURN_FAIL
      allocate( sort_strings(work%ntot), stat=status)
      if(status/=0) goto 700
      if( fill_sort_strings_and_sort(work%r, work%n, sort_strings, &
           y, work%mvcode, work, err) == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      deallocate(sort_strings, stat=status)
      if(status/=0) goto 750
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
750   call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function sort_cases_by_missingness
   !##################################################################
   integer(kind=our_int) function fill_sort_strings_and_sort( &
        r, n, sort_strings, y, mvcode, work, err) result(answer)
      implicit none
      ! declare args
      integer(kind=our_int), intent(in) :: r, n
      character(len=r), intent(inout) :: sort_strings(n)
      real(kind=our_dble), intent(in) :: y(:,:), mvcode
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=r) :: tmp_string
      integer(kind=our_int) :: i, j, patt, kase
      character(len=*), parameter :: &
           subname = "fill_sort_strings_and_sort"
      ! begin
      answer = RETURN_FAIL
      ! fill sort_strings
      do i = 1, n
         do j = 1, r
            sort_strings(i)(j:j) = "."
            if( y(i,j) == mvcode ) sort_strings(i)(j:j) = "m"
         end do
      end do
      ! do the sort
      if( dyn_alloc(work%case_order, n, err) == RETURN_FAIL ) goto 800
      if( qsort( sort_strings, work%case_order, n, r, &
           .false., .false., err) == RETURN_FAIL ) goto 800
      ! count distinct missingness patterns
      work%npatt = 0
      tmp_string = ""
      do i = 1, n
         if( sort_strings(work%case_order(i)) /= tmp_string ) then
            work%npatt = work%npatt + 1
            tmp_string = sort_strings( work%case_order(i) )
         end if
      end do
      ! create the arrays mis, first_case_in_patt, last_case_in_patt
      if( dyn_alloc(work%mis, work%npatt, work%r, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc(work%first_case_in_patt, work%npatt, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc(work%last_case_in_patt, work%npatt, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc(work%n_in_patt, work%npatt, err) &
           == RETURN_FAIL ) goto 800
      tmp_string = ""
      patt = 0
      do i = 1, n
         if(sort_strings(work%case_order(i))/=tmp_string) then
            patt = patt + 1
            tmp_string = sort_strings(work%case_order(i))
            work%first_case_in_patt(patt) = i
            do j = 1, r
               if(tmp_string(j:j)=="m") then
                  work%mis(patt,j) = .true.
               else
                  work%mis(patt,j) = .false.
               end if
            end do
         end if
      end do
      do patt = 1, work%npatt-1
         work%last_case_in_patt(patt) = &
              work%first_case_in_patt(patt+1) - 1
      end do
      work%last_case_in_patt(work%npatt) = work%n
      work%n_in_patt(:) = work%last_case_in_patt(:) - &
           work%first_case_in_patt(:) + 1
      ! create which_patt
      if( dyn_alloc(work%which_patt, n, err) == RETURN_FAIL ) goto 800
      do patt = 1, work%npatt
         do kase = work%first_case_in_patt(patt), &
              work%last_case_in_patt(patt)
            i = work%case_order(kase)  ! actual case number
            work%which_patt(i) = patt
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function fill_sort_strings_and_sort
   !##################################################################
   integer(kind=our_int) function allocate_workspace_items( work, err) &
        result(answer)
      ! allocates all remaining workspace items
      implicit none
      ! declare args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = &
           "allocate_workspace_items"
      ! begin
      answer = RETURN_FAIL
      if( dyn_alloc(work%ysort, work%n, work%r, err) &
           == RETURN_FAIL) goto 800
      work%ysort(:,:) = 0.D0
      if( dyn_alloc(work%yimp, work%n, work%r, err) &
           == RETURN_FAIL) goto 800
      work%yimp(:,:) = 0.D0
      if( dyn_alloc(work%xsort, work%n, work%p, err) &
           == RETURN_FAIL) goto 800
      work%xsort(:,:) = 0.D0
      if( dyn_alloc(work%nobs, work%r, err) == RETURN_FAIL ) goto 800
      work%nobs(:) = 0
      if( dyn_alloc(work%ybar, work%r, err) == RETURN_FAIL ) goto 800
      work%ybar(:) = 0.D0
      if( dyn_alloc(work%ysdv, work%r, err) == RETURN_FAIL ) goto 800
      work%ysdv(:) = 0.D0
      if( dyn_alloc(work%beta, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%beta(:,:) = 0.D0
      if( dyn_alloc(work%sigma, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%sigma(:,:) = 0.D0
      if( dyn_alloc(work%oldbeta, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%oldbeta(:,:) = 0.D0
      if( dyn_alloc(work%oldsigma, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%oldsigma(:,:) = 0.D0
      if( dyn_alloc(work%oldoldbeta, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%oldoldbeta(:,:) = 0.D0
      if( dyn_alloc(work%oldoldsigma, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%oldoldsigma(:,:) = 0.D0
      if( dyn_alloc(work%ratebeta, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%ratebeta(:,:) = 0.D0
      if( dyn_alloc(work%ratesigma, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%ratesigma(:,:) = 0.D0
      if( dyn_alloc(work%sigma_ridge, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%sigma_ridge(:,:) = 0.D0
      if( dyn_alloc(work%prior_sscp, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%prior_sscp(:,:) = 0.D0
      if( dyn_alloc(work%wkppA, work%p, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkppA(:,:) = 0.D0
      if( dyn_alloc(work%wkp, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkp(:) = 0.D0
      if( dyn_alloc(work%wkprA, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%wkprA(:,:) = 0.D0
      if( dyn_alloc(work%wkrrA, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%wkrrA(:,:) = 0.D0
      if( dyn_alloc(work%wkrrB, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%wkrrB(:,:) = 0.D0
      if( dyn_alloc(work%wkrprpA, work%r*work%p, work%r*work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkrprpA(:,:) = 0.D0
      if( dyn_alloc(work%wkrpA, work%r*work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkrpA(:) = 0.D0
      if( dyn_alloc(work%wkrpB, work%r*work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkrpB(:) = 0.D0
      if( dyn_alloc(work%wknpA, work%n, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wknpA(:,:) = 0.D0
      if( dyn_alloc(work%wknpB, work%n, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wknpB(:,:) = 0.D0
      if( dyn_alloc(work%wknA, work%n, err) &
           == RETURN_FAIL ) goto 800
      work%wknA(:) = 0.D0
      if( dyn_alloc(work%wknB, work%n, err) &
           == RETURN_FAIL ) goto 800
      work%wknB(:) = 0.D0
      if( dyn_alloc(work%wknC, work%n, err) &
           == RETURN_FAIL ) goto 800
      work%wknC(:) = 0.D0
      if( dyn_alloc(work%wkpA, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkpA(:) = 0.D0
      if( dyn_alloc(work%wkpB, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkpB(:) = 0.D0
      if( dyn_alloc(work%wkpC, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%wkpC(:) = 0.D0
      if( dyn_alloc(work%iwkp, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%iwkp(:) = 0
      if( dyn_alloc(work%lwkp, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%lwkp(:) = .FALSE.
      if( dyn_alloc(work%xtxinv, work%p, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%xtxinv(:,:) = 0.D0
      if( dyn_alloc(work%xtxinvfac, work%p, work%p, err) &
           == RETURN_FAIL ) goto 800
      work%xtxinvfac(:,:) = 0.D0
      if( dyn_alloc(work%wkrA, work%r, err) == RETURN_FAIL ) goto 800
      work%wkrA(:) = 0.D0
      if( dyn_alloc(work%epsteps, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%epsteps(:,:) = 0.D0
      if( dyn_alloc(work%yhat, work%n, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%yhat(:,:) = 0.D0
      if( dyn_alloc(work%eps, work%n, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%eps(:,:) = 0.D0
      if( dyn_alloc(work%xty, work%p, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%xty(:,:) = 0.D0
      if( dyn_alloc(work%yty, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%yty(:,:) = 0.D0
      if( dyn_alloc(work%sigswp, work%r, work%r, err) &
           == RETURN_FAIL ) goto 800
      work%sigswp(:,:) = 0.D0
      if( dyn_alloc(work%swept, work%r, err) == RETURN_FAIL ) goto 800
      work%swept(:) = .false.
      work%nparam = work%p * work%r + ( work%r * (work%r + 1) ) / 2
      if( dyn_alloc(work%theta, work%nparam, err) &
          == RETURN_FAIL ) goto 800
      work%theta(:) = 0.D0
      if( dyn_alloc(work%oldtheta, work%nparam, err) &
          == RETURN_FAIL ) goto 800
      work%oldtheta(:) = 0.D0
      if( dyn_alloc(work%worst_weights, work%nparam, err) &
          == RETURN_FAIL ) goto 800
      work%worst_weights(:) = 0.D0
      work%em_store_incr = 0.D0
      work%em_store_counter = 0
      work%em_store_iter = 0
      work%next_em_store_iter = 0
      if( dyn_alloc(work%em_store_theta, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%em_store_theta(:) = 0.D0
      if( dyn_alloc(work%next_em_store_theta, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%next_em_store_theta(:) = 0.D0
      if( dyn_alloc(work%em_thetahat, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%em_thetahat(:) = 0.D0
      work%em_worst_ok = .false.
      if( dyn_alloc(work%uvec, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%uvec(:) = 0.D0
      if( dyn_alloc(work%uvec_new, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%uvec_new(:) = 0.D0
      if( dyn_alloc(work%vvec, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%vvec(:) = 0.D0
      if( dyn_alloc(work%vvec_new, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%vvec_new(:) = 0.D0
      work%worst_frac = 0.D0
      if( dyn_alloc(work%worst_linear_coef, work%nparam, err) &
           == RETURN_FAIL ) goto 800
      work%worst_linear_coef(:) = 0.D0
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function allocate_workspace_items
  !##################################################################
   integer(kind=our_int) function make_xsort_and_ysort( x, y, &
        work ) result(answer)
      implicit none
      ! declare args
      real(kind=our_dble), intent(in) :: y(:,:), x(:,:)
      type(workspace_type), intent(inout) :: work
      ! locals
      integer(kind=our_int) :: i, patt, kase
      ! and dimension n.
      character(len=*), parameter :: subname = "make_xsort_and_ysort"
      ! begin
      answer = RETURN_FAIL
      ! make xsort and ysort
      do patt = 1, work%npatt
         do kase = work%first_case_in_patt(patt), &
               work%last_case_in_patt(patt)
             i = work%case_order(kase)  
             work%ysort(kase,:) = y(i,:)
             work%xsort(kase,:) = x(i,:)
          end do
       end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
   end function make_xsort_and_ysort
   !##################################################################
   integer(kind=our_int) function find_means_and_variances( work, err) &
        result(answer)
      !  Finds the means and standard deviations of the observed responses
      !  for each variable. 
      implicit none
      ! declare args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = &
           "find_means_and_variances"
      integer(kind=our_int) :: nobs, i, j
      real(kind=our_dble) :: sum, ybar, ss
      character(len=12) :: Sint
      ! begin
      answer = RETURN_FAIL
      sInt = "???"
      do j = 1, work%r
         ! write(Sint,"(I12)") j
         ! compute mean of observed values
         nobs = 0
         sum = 0.D0
         do i = 1, work%n
            if( work%ysort(i,j) == work%mvcode ) cycle
            nobs = nobs + 1
            sum = sum + work%ysort(i,j)
         end do
         if(nobs < 2) goto 700
         work%nobs(j) = nobs
         ybar = sum / real( nobs, kind=our_dble)
         work%ybar(j) = ybar
         ! compute variance
         ss = 0.D0
         do i = 1, work%n
            if( work%ysort(i,j) == work%mvcode ) cycle
            ss = ss + ( work%ysort(i,j) - ybar )**2
         end do
         work%ysdv(j) = sqrt( ss / real(nobs-1, kind=our_dble) )
         if( work%ysdv(j) <= 0.D0 ) goto 750
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
700   Sint = adjustl(Sint)
      call err_handle(err, 1, &
           comment = "Cannot estimate variance; fewer than 2 cases" )
      call err_handle(err, 4, ivar = j )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
750   Sint = adjustl(Sint)
      call err_handle(err, 1, &
           comment = "Zero variance; observed values are identical" )
      call err_handle(err, 4, ivar = j )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function find_means_and_variances
   !##################################################################
   integer(kind=our_int) function regress_univariate( work, &
        startval_present, beta_start, sigma_start, err) result(answer)
      ! Finds the default staring values by regressing the 
      ! observed values of each Y-variable on the X's.
      ! Also computes a ridge estimate of sigma if
      ! ridge_hyperparameter is not zero.
      !
      ! New version using Householder QR decomposition
      !
      implicit none
      ! declare args
      type(workspace_type), intent(inout) :: work
      logical, intent(in) :: startval_present
      real(kind=our_dble), intent(in) :: beta_start(:,:), &
           sigma_start(:,:)
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "regress_univariate"
      integer(kind=our_int) :: nobs, i, j, jj, ii, rank
      real(kind=our_dble) :: sum, ss, mse
      character(len=12) :: Sint
      ! begin
      answer = RETURN_FAIL
      sInt = "???"
      ! get starting values for beta
      work%sigma(:,:) = 0.D0
      work%sigma_ridge(:,:) = 0.D0
      if( (.not. startval_present ) .or. &
           ( work%prior_type == "ridge"  ) ) then
         ! regress y on x's using observed cases only, store
         ! results in beta and sigma
         do jj = 1, work%r
            ! write(Sint,"(I12)") jj
            ! store x-matrix in wknpA, y-vector in wknA
            nobs = 0
            ii = 0
            do i = 1, work%n
               if( work%ysort(i,jj) == work%mvcode ) cycle
               nobs = nobs + 1
               ii = ii + 1
               do j = 1, work%p
                  work%wknpA(ii,j) = work%xsort(i,j)
               end do
               work%wknA(ii) = work%ysort(i,jj)
            end do
            if( nobs < 2 ) goto 500
            if( householder_ols( work%wknpA(1:nobs,:), &
                 work%wknA(1:nobs), rank, work%lwkp, work%wkpA, &
                 work%wknB(1:nobs), work%wknC(1:nobs), &
                 work%wkpB, work%wkpC, &
                 work%wknpB(1:nobs,:), work%wkppA, &
                 work%iwkp, err ) == RETURN_FAIL ) goto 550
            do j = 1, work%p
               work%beta(j,jj) = work%wkpA(j)
            end do
            if( nobs > work%p ) then
               ! residual variance using denominator of nobs-p
               ss = 0.D0
               do i = 1, work%n
                  if( work%ysort(i,jj) == work%mvcode ) cycle
                  sum = 0.D0
                  do j = 1, work%p
                     sum = sum + work%xsort(i,j) * work%beta(j,jj)
                  end do
                  ss = ss + ( work%ysort(i,jj) - sum )**2
               end do
               mse = ss / real( nobs - work%p, our_dble )
               if( mse <= 1.D-08 * work%ysdv(jj)**2 ) then
                  ! apparent perfect fit; set residual variance to
                  ! half the marginal variance
                  work%sigma(jj,jj) = 0.5D0 * work%ysdv(jj)**2
               else
                  work%sigma(jj,jj) = mse
               end if
            else
               ! too few observations; set residual variance to
               ! half the marginal variance
                work%sigma(jj,jj) = 0.5D0 * work%ysdv(jj)**2
            end if
            work%sigma_ridge(jj,jj) = work%sigma(jj,jj)
         end do
      end if
      ! overwrite beta and sigma if starting values were supplied
      if( startval_present ) then
         work%beta(:,:) = beta_start(:,:)
         work%sigma(:,:) = sigma_start(:,:)
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
500   Sint = adjustl(Sint)
      call err_handle(err, 1, &
           comment = "Unable to regress response on x" )
      call err_handle(err, 1, &
           comment = "Fewer than two observed cases available" )
      call err_handle(err, 4, ivar = jj )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
550   Sint = adjustl(Sint)
      call err_handle(err, 1, &
           comment = "Unable to regress response on x" )
      call err_handle(err, 4, ivar = jj )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function regress_univariate
   !##################################################################
   integer(our_int) function make_xtxinv( work, err ) result(answer)
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: j, k
      character(len=*), parameter :: subname = "make_xtxinv"
      ! begin
      answer = RETURN_FAIL
      do j = 1, work%p
         do k = 1, j
            work%wkppA(j,k) = dot_product( work%xsort(:,j), &
                 work%xsort(:,k) )
         end do
      end do
      ! invert it and put it into work%xtxinv
      if( cholesky_in_place(work%wkppA, err) == RETURN_FAIL ) goto 800
      if( invert_lower(work%wkppA, err) == RETURN_FAIL ) goto 800
      if( premult_lower_by_transpose(work%wkppA, &
           work%xtxinv, err) == RETURN_FAIL ) goto 800
      ! put cholesky factor into work%xtxinvfac
      work%xtxinvfac(:,:) = 0.D0
      do j = 1, work%p
         do k = 1, j
            work%xtxinvfac(j,k) = work%xtxinv(j,k)
         end do
      end do
      if( cholesky_in_place(work%xtxinvfac, err) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 1, &
           comment = "Predictor (X) matrix does not have full rank")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function make_xtxinv
   !##################################################################
   integer(kind=our_int) function run_estep( work, err ) &
        result(answer)
      ! Performs an E-step
      ! Computes the expectation of y
      ! and sum yi*t(yi), i=1,...,n given the current parameters
      ! stored in beta and sigma.
      ! The expectation of y is stored in ysort, and the expectation
      ! of t(eps)*eps is stored in epsteps
      ! Also computes the observed-data loglikelihood and log-prior
      ! at the current value of beta and sigma.
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "run_estep"
      integer(kind=our_int) :: patt, j, k, kase, i
      real(kind=our_dble) :: sum, logdet, pivot
      ! begin
      answer = RETURN_FAIL
      work%sigswp(:,:) = work%sigma(:,:)
      work%swept(:) = .false.
      ! compute yhat = x*beta and observed part of y - yhat
!      work%yhat = matmul( work%xsort, work%beta ) ! fixMe
      if( matmul_boundcheck( work%xsort, work%beta, &
           work%yhat, err ) == RETURN_FAIL ) goto 800
      do i = 1, work%n
         do j = 1, work%r
            if( work%ysort(i,j) /= work%mvcode ) &
                 work%eps(i,j) = work%ysort(i,j) - work%yhat(i,j)
         end do
      end do
      work%yty(:,:) = 0.D0
      !########################################
      ! compute the log-prior density
      logdet = 0.D0
      ! compute log-det of sigma and its negative inverse
      do j = 1, work%r
         pivot = work%sigswp(j,j)
         if( sweep_forward(work%sigswp, j, err) &
              == RETURN_FAIL) goto 750
         work%swept(j) = .true.
         if( pivot <= 0.D0 ) goto 700
         logdet = logdet + log(pivot)
      end do
      work%logpri = - ( work%prior_df + real( work%r+1, our_dble) ) &
           * logdet / 2.D0
      ! compute trace of -sigmainv * lambdainv
      sum = 0.D0
      do j = 1, work%r
         sum = sum + work%sigswp(j,j) * work%prior_sscp(j,j)
      end do
      work%logpri = work%logpri + sum / 2.D0
      !########################################
      ! cycle through the patterns
      work%loglik = 0.D0
      do patt = 1, work%npatt
         ! sweep and reverse-sweep, as necessary
         do j = 1, work%r
            if( (.not. work%mis(patt,j)) .and. (.not. work%swept(j)) ) then
               pivot = work%sigswp(j,j)
               if( sweep_forward(work%sigswp, j, err) &
                    == RETURN_FAIL) goto 750
               work%swept(j) = .true.
               if( pivot <= 0.D0 ) goto 700
               logdet = logdet + log(pivot)
            end if
            if( work%mis(patt,j) .and. work%swept(j) ) then
               if( sweep_reverse(work%sigswp, j, err) &
                    == RETURN_FAIL) goto 750
               work%swept(j) = .false.
               pivot = work%sigswp(j,j)
               if( pivot <= 0.D0 ) goto 700
               logdet = logdet - log(pivot)
            end if
         end do
         ! fill in the elements above the diagonal
         do j = 1, work%r
            do k = 1, j
               work%sigswp(k,j) = work%sigswp(j,k)
            end do
         end do
         ! cycle through the cases within the pattern
         do kase = work%first_case_in_patt(patt), &
              work%last_case_in_patt(patt)
            i = work%case_order(kase)  ! actual case number
            ! increment the loglikelihood
            do j = 1, work%r
               if( work%mis(patt,j) ) cycle
               sum = 0.D0
               do k = 1, work%r
                  if( work%mis(patt,k) ) cycle
                  sum = sum + work%eps(kase,k)*work%sigswp(k,j)
               end do
               work%wkrA(j) = sum
            end do
            sum = 0.D0
            do j = 1, work%r
               if( work%mis(patt,j) ) cycle
               sum = sum + work%wkrA(j)*work%eps(kase,j)
            end do
            work%loglik = work%loglik - logdet / 2.D0 + sum / 2.D0
            ! put expected value of yi_mis into ysort
            do j = 1, work%r
               if( .not.work%mis(patt,j) ) cycle ! j missing
               sum = 0.D0
               do k = 1, work%r
                  if( work%mis(patt,k) ) cycle  ! k observed
                  sum = sum + work%sigswp(j,k)*work%eps(kase,k)
               end do
               work%eps(kase,j) = sum
               work%ysort(kase,j) = sum + work%yhat(kase,j)
               work%yimp(i,j) = work%ysort(kase,j)
            end do
            ! increment lower triangle of yty
            do j = 1, work%r
               do k = 1, j
                  work%yty(j,k) = work%yty(j,k) + &
                       work%ysort(kase,j) * work%ysort(kase,k)
                  if( work%mis(patt,j) .and. work%mis(patt,k) ) then
                     work%yty(j,k) = work%yty(j,k) + work%sigswp(j,k)
                  end if
               end do
            end do
         end do
      end do
      ! fill in yty above the diagonal
      do j = 1, work%r
         do k = 1, j-1
            work%yty(k,j) = work%yty(j,k)
         end do
      end do
      ! accumulate xty
      do j = 1, work%p
         do k = 1, work%r
            work%xty(j,k) = dot_product( work%xsort(:,j), &
                 work%ysort(:,k) )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Attempted logarithm of non-positive number" )
      call err_handle(err, 1, &
           comment = "Cov. matrix became singular or negative definite")
      goto 800
750   call err_handle(err, 1, &
           comment = "Cov. matrix became singular or negative definite")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function run_estep
   !##################################################################
   integer(kind=our_int) function run_mstep( work, err ) result(answer)
      ! Performs the M-step
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "run_mstep"
      integer(kind=our_int) :: j, k
      ! begin
      answer = RETURN_FAIL
      ! update beta and sigma
!      work%beta = matmul( work%xtxinv, work%xty ) ! fixMe
      if( matmul_boundcheck( work%xtxinv, work%xty, &
           work%beta, err ) == RETURN_FAIL ) goto 800
      do j = 1, work%r
         do k = 1, j
            work%wkrrA(j,k) = dot_product( work%xty(:,j), &
                 work%beta(:,k) )
            work%wkrrA(k,j) = work%wkrrA(j,k)
         end do
      end do
      work%epsteps(:,:) = work%yty(:,:) - work%wkrrA(:,:)
      work%sigma(:,:) = ( work%epsteps(:,:) + work%prior_sscp(:,:) )  &
           / ( work%prior_df + real( work%n + work%r + 1, our_dble ) )
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function run_mstep
   !##################################################################
   integer(kind=our_int) function update_theta_and_oldtheta( work ) &
        result(answer) 
      !  Updates the vectorized parameters
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      ! locals
      character(len=*), parameter :: subname = "update_theta_and_oldtheta"
      integer(kind=our_int) :: j, k, posn
      ! begin
      answer = RETURN_FAIL
      ! vectorize theta and oldtheta
      posn = 0
      do k = 1, work%r
         do j = 1, work%p
            posn = posn + 1
            work%oldtheta( posn ) = work%oldbeta( j, k )
            work%theta( posn ) = work%beta( j, k )
         end do
      end do
      do k = 1, work%r
         do j = k, work%r
            posn = posn + 1
            work%oldtheta( posn ) = work%oldsigma( j, k )
            work%theta( posn ) = work%sigma( j, k )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
   end function update_theta_and_oldtheta
   !##################################################################
   integer(kind=our_int) function update_rates( work ) result(answer) 
      !  Updates the elementwise rates of convergence
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      ! locals
      character(len=*), parameter :: subname = "update_rates"
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: num, denom
      real(kind=our_sgle) :: dummy
      ! begin
      answer = RETURN_FAIL
      do k = 1, work%r
         do j = 1, work%p
            denom = abs( work%oldbeta(j,k) - work%oldoldbeta(j,k) )
            num = abs( work%beta(j,k) - work%oldbeta(j,k) )
            if( ( denom > tiny(dummy) ) .and. (num > tiny(dummy) ) ) then
               work%ratebeta(j,k) = num / denom
            else
               work%ratebeta(j,k) = 0.D0
            end if
         end do
      end do
      do k = 1, work%r
         do j = k, work%r
            denom = work%oldsigma(j,k) - work%oldoldsigma(j,k)
            if( denom /= 0.D0 ) then
               work%ratesigma(j,k) = &
                    ( work%sigma(j,k) - work%oldsigma(j,k) ) / denom
            else
               work%ratesigma(j,k) = 0.D0
            end if
            if( j /= k ) work%ratesigma(k,j) = work%ratesigma(j,k)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
   end function update_rates
   !##################################################################
   integer(kind=our_int) function update_theta( work ) result(answer) 
      !  Updates the vectorized parameter
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      ! locals
      character(len=*), parameter :: subname = "update_theta"
      integer(kind=our_int) :: j, k, posn
      ! begin
      answer = RETURN_FAIL
      ! vectorize theta
      posn = 0
      do k = 1, work%r
         do j = 1, work%p
            posn = posn + 1
            work%theta( posn ) = work%beta( j, k )
         end do
      end do
      do k = 1, work%r
         do j = k, work%r
            posn = posn + 1
            work%theta( posn ) = work%sigma( j, k )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
   end function update_theta
   !##################################################################
   integer(kind=our_int) function find_max_rel_diff( work ) &
        result(answer) 
      ! run this after update_theta_and_worst
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      ! locals
      character(len=*), parameter :: subname = "find_max_rel_diff"
      integer(kind=our_int) :: i
      ! begin
      answer = RETURN_FAIL
      work%reldiff = 0.D0
      do i = 1, work%nparam
         if( work%oldtheta(i) == 0.D0 ) cycle
         work%reldiff = max( work%reldiff, &
              abs( work%theta(i) - work%oldtheta(i) ) &
              / abs( work%oldtheta(i) ) )
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
   end function find_max_rel_diff
   !##################################################################
   integer(kind=our_int) function estimate_worst_frac( work, err ) &
        result(answer)
      ! Estimates the worst fraction of missing information by the
      ! power method of Fraley et al. (2007) with finite differencing
      ! instead of Richardson extrapolation.
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "estimate_worst_frac"
      real(kind=our_dble) :: sum, lambda, lambda_new, delta, &
           vscale, diff
      logical :: vvec_converged, lambda_converged
      integer(kind=our_int) :: k, i, kk, kkk, ijunk, j, posn
      real(kind=our_dble), pointer :: f1(:)=>null(), &
           f2(:)=>null(), ju(:)=>null(), ju_new(:)=>null()
      integer(kind=our_int), parameter :: maxits_power = 20, &
           maxits_richardson = 20
      real(kind=our_dble), parameter :: criterion_power = 1.D-4, &
           criterion_richardson = 1.D-4
      ! begin
      answer = RETURN_FAIL
      work%em_worst_ok = .false.
      ! put final parameters into em_thetahat
      posn = 0
      do k = 1, work%r
         do j = 1, work%p
            posn = posn + 1
            work%em_thetahat( posn ) = work%beta( j, k )
         end do
      end do
      do k = 1, work%r
         do j = k, work%r
            posn = posn + 1
            work%em_thetahat( posn ) = work%sigma( j, k )
         end do
      end do
      ! allocate local workspaces
      if( dyn_alloc( f1, work%nparam, err ) == RETURN_FAIL ) goto 800
      if( dyn_alloc( f2, work%nparam, err ) == RETURN_FAIL ) goto 800
      if( dyn_alloc( ju, work%nparam, err ) == RETURN_FAIL ) goto 800
      if( dyn_alloc( ju_new, work%nparam, err ) == RETURN_FAIL ) goto 800
      ! prepare main iteration for power method, using em trajectory
      ! as starting direction vector
      lambda_new = 0.D0
      lambda_converged = .false.
      work%vvec_new(:) = work%em_store_theta(:) - work%em_thetahat(:)
      kkk = 0
      do
         kkk = kkk + 1
         if( lambda_converged .or. ( kkk > maxits_power ) ) exit
         lambda = lambda_new
         work%vvec(:) = work%vvec_new(:)
         sum = 0.D0
         do i = 1, work%nparam
            sum = sum + work%vvec(i)**2
         end do
         if( sum <= 0.D0 ) goto 700
         sum = sqrt( sum )
         do i = 1, work%nparam
            work%uvec_new(i) = work%vvec(i) / sum
         end do
         vscale = sum
         !!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! approximate vvec_new = Jhat %*% uvec_new 
         !!!     = Jhat %*% vvec / vscale
         !!! by finite differencing, two-sided
         vvec_converged = .false.
         delta = 0.5D0
         kk = 0
         ju_new(:) = 0.2D0 ! something nonzero
         do
            kk = kk + 1
            if( vvec_converged .or. ( kk > maxits_richardson ) ) exit
            ju(:) = ju_new(:)
            delta = delta * 0.5D0
            ! take one step of em from thetahat + delta * vvec,
            ! store result in f1
            posn = 0
            do k = 1, work%r
               do j = 1, work%p
                  posn = posn + 1
                  work%beta( j, k ) = work%em_thetahat( posn ) + &
                       delta * work%vvec(posn)
               end do
            end do
            do k = 1, work%r
               do j = k, work%r
                  posn = posn + 1
                  work%sigma( j, k ) = work%em_thetahat( posn ) + &
                       delta * work%vvec(posn)
                  if( k /= j ) work%sigma( k, j ) = work%sigma( j, k )
               end do
            end do
            if( run_estep( work, err ) == RETURN_FAIL) goto 700
            if( run_mstep( work, err ) == RETURN_FAIL) goto 800
            if( update_theta( work ) == RETURN_FAIL ) goto 800
            f1(:) = work%theta(:)
            ! take one step of em from thetahat - delta * vvec,
            ! store result in f2
            posn = 0
            do k = 1, work%r
               do j = 1, work%p
                  posn = posn + 1
                  work%beta( j, k ) = work%em_thetahat( posn ) - &
                       delta * work%vvec(posn)
               end do
            end do
            do k = 1, work%r
               do j = k, work%r
                  posn = posn + 1
                  work%sigma( j, k ) = work%em_thetahat( posn ) - &
                       delta * work%vvec(posn)
                  if( k /= j ) work%sigma( k, j ) = work%sigma( j, k )
               end do
            end do
            if( run_estep( work, err ) == RETURN_FAIL) goto 750
            if( run_mstep( work, err ) == RETURN_FAIL) goto 800
            if( update_theta( work ) == RETURN_FAIL ) goto 800
            f2(:) = work%theta(:)
            !! compute ju_new
            ju_new(:) = ( f1(:) - f2(:) ) / ( 2.D0 * delta )
            !! compare ju_new and ju
            vvec_converged = .true.
            do posn = 1, work%nparam
               diff = abs( ju_new(posn) - ju(posn) )
               if( ju(posn) == 0.D0 ) cycle
               if( ( diff/abs(ju(posn)) ) >= criterion_richardson ) then
                  vvec_converged = .false.
                  exit
               end if
            end do
         end do
         !!!!!!!!!!!!!!!! end finite differencing
         work%vvec_new(:) = ju_new(:) / vscale
         if( all( work%vvec_new(:) == 0.D0 ) ) then
            lambda_new = 0.D0
            lambda_converged = .true.
            exit
         end if
         sum = 0.D0
         do i = 1, work%nparam
            sum = sum + work%uvec_new(i) * work%vvec_new(i)
         end do
         lambda_new = sum
         if( abs( lambda_new - lambda ) <= criterion_power ) &
              lambda_converged = .true.
      end do
      ! put results into workspace
      work%worst_frac = lambda_new
      sum = 0.D0
      do i = 1, work%nparam
         sum = sum + work%vvec_new(i)**2
      end do
      if( sum > 0.D0 ) then
         sum = sqrt( sum )
         do i = 1, work%nparam
            work%worst_linear_coef(i) = work%vvec_new(i) / sum
         end do
      else
         work%worst_linear_coef(:) = 0.D0
      end if
      if( lambda_converged ) then
         work%em_worst_ok = .true.
      else
         goto 760
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
700   call err_handle(err, 1, &
           comment = "Attempted division by zero; procedure aborted")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
750   call err_reset(err)
      call err_handle(err, 1, &
           comment = "Finite-differencing procedure strayed outside" )
      call err_handle(err, 1, &
           comment = "parameter space; solution at or near boundary")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
760   continue
      call err_handle(err, 1, &
           comment = "Eigen power method failed to converge")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! final cleanup
999   continue
      ijunk = dyn_dealloc( f1, err )
      ijunk = dyn_dealloc( f2, err )
      ijunk = dyn_dealloc( ju, err )
      ijunk = dyn_dealloc( ju_new, err )
      return
    end function estimate_worst_frac
   !##################################################################
   integer(kind=our_int) function run_istep( work, rand, err ) &
        result(answer)
      ! Performs an I-step.
      ! Simulates a draw of the missing parts of y given the 
      ! current parameters stored in beta and sigma.
      ! Also accumulates xty and yty.
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(random_gendata), intent(inout) :: rand
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "run_istep"
      integer(kind=our_int) :: patt, j, k, kase, i, jposn, kposn
      real(kind=our_dble) :: sum, logdet, pivot
      real(kind=our_sgle) :: z
      ! begin
      answer = RETURN_FAIL
      work%sigswp(:,:) = work%sigma(:,:)
      work%swept(:) = .false.
      ! compute yhat = x*beta and observed part of y - yhat
!      work%yhat = matmul( work%xsort, work%beta ) !fixMe
      if( matmul_boundcheck( work%xsort, work%beta, &
           work%yhat, err ) == RETURN_FAIL ) goto 800
      do i = 1, work%n
         do j = 1, work%r
            if( work%ysort(i,j) /= work%mvcode ) &
                 work%eps(i,j) = work%ysort(i,j) - work%yhat(i,j)
         end do
      end do
      work%yty(:,:) = 0.D0
      !########################################
      ! compute the log-prior density
      logdet = 0.D0
      ! compute log-det of sigma and its negative inverse
      do j = 1, work%r
         pivot = work%sigswp(j,j)
         if( sweep_forward(work%sigswp, j, err) &
              == RETURN_FAIL) goto 750
         work%swept(j) = .true.
         if( pivot <= 0.D0 ) goto 700
         logdet = logdet + log(pivot)
      end do
      work%logpri = - ( work%prior_df + real( work%r+1, our_dble) ) &
           * logdet / 2.D0
      ! compute trace of -sigmainv * lambdainv
      sum = 0.D0
      do j = 1, work%r
         sum = sum + work%sigswp(j,j) * work%prior_sscp(j,j)
      end do
      work%logpri = work%logpri + sum / 2.D0
      !########################################
      ! cycle through the patterns
      work%loglik = 0.D0
      do patt = 1, work%npatt
         ! sweep and reverse-sweep, as necessary
         do j = 1, work%r
            if( (.not. work%mis(patt,j)) .and. (.not. work%swept(j)) ) then
               pivot = work%sigswp(j,j)
               if( sweep_forward(work%sigswp, j, err) &
                    == RETURN_FAIL) goto 750
               work%swept(j) = .true.
               if( pivot <= 0.D0 ) goto 700
               logdet = logdet + log(pivot)
            end if
            if( work%mis(patt,j) .and. work%swept(j) ) then
               if( sweep_reverse(work%sigswp, j, err) &
                    == RETURN_FAIL) goto 750
               work%swept(j) = .false.
               pivot = work%sigswp(j,j)
               if( pivot <= 0.D0 ) goto 700
               logdet = logdet - log(pivot)
            end if
         end do
         ! fill in the elements above the diagonal
         do j = 1, work%r
            do k = 1, j
               work%sigswp(k,j) = work%sigswp(j,k)
            end do
         end do
         ! extract the residual covariance matrix for
         ! yimis and put its cholesky factor into lower triangle 
         ! of wkrrB
         jposn = 0
         do j = 1, work%r
            if( .not.work%mis(patt,j) ) cycle
            jposn = jposn + 1
            kposn = 0
            do k = 1, j
               if( .not.work%mis(patt,k) ) cycle
               kposn = kposn + 1
               work%wkrrA(jposn,kposn) = work%sigswp(j,k)
            end do
         end do
         if( jposn > 0 ) then
            if( cholesky_in_place( work%wkrrA(1:jposn, 1:jposn), &
                 err) == RETURN_FAIL ) goto 750
         end if
         jposn = 0
         do j = 1, work%r
            if( .not.work%mis(patt,j) ) cycle
            jposn = jposn + 1
            kposn = 0
            do k = 1, j
               if( .not.work%mis(patt,k) ) cycle
               kposn = kposn + 1
               work%wkrrB(j,k) = work%wkrrA(jposn,kposn) 
            end do
         end do
         ! cycle through the cases within the pattern
         do kase = work%first_case_in_patt(patt), &
              work%last_case_in_patt(patt)
            i = work%case_order(kase)  ! actual case number
            ! increment the loglikelihood
            do j = 1, work%r
               if( work%mis(patt,j) ) cycle
               sum = 0.D0
               do k = 1, work%r
                  if( work%mis(patt,k) ) cycle
                  sum = sum + work%eps(kase,k)*work%sigswp(k,j)
               end do
               work%wkrA(j) = sum
            end do
            sum = 0.D0
            do j = 1, work%r
               if( work%mis(patt,j) ) cycle
               sum = sum + work%wkrA(j)*work%eps(kase,j)
            end do
            work%loglik = work%loglik - logdet / 2.D0 + sum / 2.D0
            ! draw the standard normal variates
            do j = 1, work%r
               if( .not. work%mis(patt,j) ) cycle
               if( ran_snorm( rand, z, err) == RETURN_FAIL) &
                    goto 750
               work%wkrA(j) = z
            end do
            ! put expected value of yi_mis and y_imis, and add
            ! the random residual
            do j = 1, work%r
               if( .not.work%mis(patt,j) ) cycle ! j missing
               sum = 0.D0
               do k = 1, work%r
                  if( work%mis(patt,k) ) cycle  ! k observed
                  sum = sum + work%sigswp(j,k)*work%eps(kase,k)
               end do
               do k = 1, j
                  if( .not.work%mis(patt,k) ) cycle  ! k missing
                  sum = sum + work%wkrrB(j,k)*work%wkrA(k)
               end do
               work%eps(kase,j) = sum
               work%ysort(kase,j) = sum + work%yhat(kase,j)
               ! place the imputed value into yimp as well
               work%yimp(i,j) = work%ysort(kase,j)
            end do
            ! increment lower triangle of yty
            do j = 1, work%r
               do k = 1, j
                  work%yty(j,k) = work%yty(j,k) + &
                       work%ysort(kase,j) * work%ysort(kase,k)
               end do
            end do
         end do
      end do
      ! fill in yty above the diagonal
      do j = 1, work%r
         do k = 1, j-1
            work%yty(k,j) = work%yty(j,k)
         end do
      end do
      ! accumulate xty
      do j = 1, work%p
         do k = 1, work%r
            work%xty(j,k) = dot_product( work%xsort(:,j), &
                 work%ysort(:,k) )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Attempted logarithm of non-positive number" )
      call err_handle(err, 1, &
           comment = "Cov. matrix became non-positive definite")
      goto 800
750   call err_handle(err, 1, &
           comment = "Cov. matrix became non-positive definite")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function run_istep
   !##################################################################
   integer(kind=our_int) function run_pstep( work, rand, err) &
        result(answer)
      ! Performs the P-step
      implicit none
      ! args
      type(workspace_type), intent(inout) :: work
      type(random_gendata), intent(inout) :: rand
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "run_pstep"
      integer(kind=our_int) :: j, k, last, jj, kk, jposn, kposn
      real(kind=our_sgle) :: df, z
      real(kind=our_dble) :: xi
      ! begin
      answer = RETURN_FAIL
      ! compute betahat
!      work%beta = matmul( work%xtxinv, work%xty ) ! fixMe
      if( matmul_boundcheck( work%xtxinv, work%xty, &
           work%beta, err ) == RETURN_FAIL ) goto 800
      ! compute residual SSCP matrix
      do j = 1, work%r
         do k = 1, j
            work%wkrrA(j,k) = dot_product( work%xty(:,j), &
                 work%beta(:,k) )
            work%wkrrA(k,j) = work%wkrrA(j,k)
         end do
      end do
      work%epsteps(:,:) = work%yty(:,:) - work%wkrrA(:,:)
      ! add prior_sscp
      work%epsteps(:,:) = work%epsteps(:,:) + work%prior_sscp(:,:)
      ! compute cholesky factor of lambdainv
      if( cholesky_in_place( work%epsteps, err ) == RETURN_FAIL ) &
           goto 800
      ! put lower-triangular Bartlett factor into wkrrA and invert
      xi = work%prior_df + real( work%n - work%p, our_dble )
      do j = 1, work%r
         df = real( xi, kind=our_sgle) + real(-j + 1, kind=our_sgle)
         if( ran_genchi( rand, df, z, err) == RETURN_FAIL) goto 800
         work%wkrrA(j,j) = sqrt(z)
         do k = 1, j-1
            if( ran_snorm( rand, z, err) == RETURN_FAIL) goto 800
            work%wkrrA(j,k) = z
         end do
      end do
      if( invert_lower( work%wkrrA, err ) == RETURN_FAIL ) goto 800
      ! multiply wkrrA by transpose of epsteps
      do j = 1, work%r
         do k = 1, work%r
            last = min(j,k)
            work%wkrrB(j,k) = dot_product( work%wkrrA(j,1:last), &
                 work%epsteps(k,1:last) )
         end do
      end do
      ! premultiply wkrrB by its transpose, store result in sigma
      do j = 1, work%r
         do k = 1, j
            work%sigma(j,k) = dot_product( work%wkrrB(:,j), &
                 work%wkrrB(:,k) )
            if(k/=j) work%sigma(k,j) = work%sigma(j,k)
            work%wkrrA(j,k) = work%sigma(j,k)
         end do
      end do
      ! find Kronecker product of cholesky factor of sigma with
      ! cholesky factor of xtxinv
      if( cholesky_in_place( work%wkrrA, err ) == RETURN_FAIL ) &
           goto 800
      work%wkrprpA(:,:) = 0.D0
      do j = 1, work%r
         do k = 1, j
            jposn = (j-1)*work%p
            do jj = 1, work%p
               jposn = jposn + 1
               kposn = (k-1)*work%p
               do kk = 1, jj
                  kposn = kposn + 1
                  work%wkrprpA(jposn,kposn) = work%wkrrA(j,k) &
                       * work%xtxinvfac(jj,kk)
               end do
            end do
         end do
      end do
      ! add random error to betahat
      do j = 1, work%r*work%p
         if( ran_snorm( rand, z, err) == RETURN_FAIL) goto 800
         work%wkrpA(j) = z
      end do
      do j = 1, work%r*work%p
         work%wkrpB(j) = dot_product( work%wkrprpA(j,1:j), &
              work%wkrpA(1:j) )
      end do
      jposn = 0
      do k = 1, work%r
         do j = 1, work%p
            jposn = jposn + 1
            work%beta(j,k) = work%beta(j,k) + work%wkrpB(jposn)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function run_pstep
   !##################################################################
end module norm_engine
!#####################################################################
