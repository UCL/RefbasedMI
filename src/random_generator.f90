!#####################################################################
module random_generator
   ! Functions for random number generation
   ! Code adapted from RANDLIB version 1.3 by BW Brown, et al.
   use program_constants
   use error_handler
   implicit none
   private ! by default
   ! list public types
   public :: random_gendata
   ! list public functions and subroutines
   public :: ran_genchi, ran_gengam, ran_sgamma, ran_snorm, &
        ran_sexp, ran_genunf, ran_setall, ran_phrsd, ran_timeseed, &
        ran_seed, ran_seed_is_set
   ! private parameter defining size of generator arrays
   integer(kind=our_int), parameter, private :: NUM_GEN = 32
   ! private module parameters formerly stored in a common block
   integer(kind=our_int), parameter, private :: m1 = 2147483563, &
        m2 = 2147483399, a1a = 40014, a2a = 40692, a1w = 1033780774, &
        a2w = 1494757890, a1vw = 2082007225, a2vw = 784306273
   ! other private parameters
   character(len=*), parameter :: modname = "random_generator"
   !##################################################################
   type :: random_gendata
      !### Type for holding per-instance data
      sequence
      integer(kind=our_int), dimension(NUM_GEN) :: cg1 = 0
      integer(kind=our_int), dimension(NUM_GEN) :: cg2 = 0
      integer(kind=our_int), dimension(NUM_GEN) :: ig1 = 0
      integer(kind=our_int), dimension(NUM_GEN) :: ig2 = 0
      integer(kind=our_int), dimension(NUM_GEN) :: lg1 = 0
      integer(kind=our_int), dimension(NUM_GEN) :: lg2 = 0
      logical, dimension(NUM_GEN) :: qanti = .false.
      integer :: curntg = 1 ! stores current generator (1 to NUM_GEN)
      logical :: qqssd = .false. ! .true. if seed values have been set
   end type random_gendata
   !##################################################################
contains
   !##################################################################
   integer(our_int) function ran_genchi(gendata, df, val, err)
      !	Generates random chi-squared variate
      ! Requires df > 0
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(in) :: df
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables
      character(len=*), parameter :: subname = "ran_genchi"
      real(kind=our_sgle) :: val1
      ! check input arguments
      if ( df <= 0.0 ) goto 800
      if( ran_sgamma(gendata, df/2.0, val1, err) &
           == RETURN_FAIL ) goto 900
      val = 2.0*val1
      ! normal exit
      ran_genchi = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "Degrees of freedom not positive")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_genchi = RETURN_FAIL
      return
900   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_genchi = RETURN_FAIL
      return
   end function ran_genchi
   !##################################################################
   integer(our_int) function ran_gengam(gendata, b, a, val, err)
      !	Generates random gamma variate with mean a/b and 
      ! variance a/b**2
      !    b = scale parameter (must be positive)
      !    a = shape parameter (must be positive)
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(in) :: b, a
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables
      character(len=*), parameter :: subname = "ran_gengam"
      real(kind=our_sgle) :: val1
      ! check input arguments
      if ( b <= 0.0 .or. a <= 0.0 ) goto 800
      ! generate standard gamma and scale it
      if( ran_sgamma(gendata, a, val1, err) &
           == RETURN_FAIL ) goto 900
      val = val1/b
      ! normal exit
      ran_gengam = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "Shape or scale parameter not positive")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_gengam = RETURN_FAIL
      return
900   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_gengam = RETURN_FAIL
      return
   end function ran_gengam
   !##################################################################
   integer(our_int) function ran_genunf(gendata, low, high, val, err)
      !	Generates random variate uniformly distributed between 
      ! 'low' and 'high'
      ! Requires: high >= low
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(out) :: val
      real(kind=our_sgle), intent(in) :: high,low
      type(error_type), intent(inout) :: err
      ! local variables
      character(len=*), parameter :: subname = "ran_genunf"
      real(kind=our_sgle) :: val1
      ! check arguments
      if ( low > high ) goto 800
      ! generate the standard uniform and transform it
      if( ran_genreal(gendata,val1,err) == RETURN_FAIL ) goto 900
      val = low + (high-low)*val1
      ! normal exit
      ran_genunf = RETURN_SUCCESS
      return
800   continue
      call err_handle(err, 1, &
           comment = "Lower bound exceeds upper bound")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_genunf = RETURN_FAIL
      return
900   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_genunf = RETURN_FAIL
      return
   end function ran_genunf
   !##################################################################
   integer(our_int) function ran_snorm(gendata, val, err)
      !	Generates a standard normal random variate
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_snorm"
      real(kind=our_sgle) :: aa, s, tt, u, ustar, w, y
      integer(kind=our_int) :: i
      intrinsic float, int
      real(kind=our_sgle), parameter :: a(32) = (/0.0,.3917609e-1, &
           .7841241e-1,.1177699,.1573107,.1970991, &
           .2372021,.2776904,.3186394,.3601299,.4022501,.4450965, &
           .4887764,.5334097,.5791322,.6260990,.6744898,.7245144, &
           .7764218,.8305109,.8871466,.9467818,1.009990,1.077516, &
           1.150349,1.229859,1.318011,1.417797,1.534121,1.675940, &
           1.862732,2.153875/), &
           d(31) = (/0.0,0.0,0.0,0.0,0.0,.2636843,.2425085,.2255674, &
           .2116342,.1999243,.1899108,.1812252,.1736014,.1668419, &
           .1607967,.1553497,.1504094,.1459026,.1417700,.1379632, &
           .1344418,.1311722,.1281260,.1252791,.1226109,.1201036, &
           .1177417,.1155119,.1134023,.1114027,.1095039/), &
           t(31) = (/.7673828e-3,.2306870e-2,.3860618e-2,.5438454e-2, &
           .7050699e-2,.8708396e-2,.1042357e-1,.1220953e-1,.1408125e-1, &
           .1605579e-1,.1815290e-1,.2039573e-1,.2281177e-1,.2543407e-1, &
           .2830296e-1,.3146822e-1,.3499233e-1,.3895483e-1,.4345878e-1, &
           .4864035e-1,.5468334e-1,.6184222e-1,.7047983e-1,.8113195e-1, &
           .9462444e-1,.1123001,.1364980,.1716886,.2276241,.3304980, &
           .5847031/), &
           h(31) = (/.3920617e-1,.3932705e-1,.3950999e-1,.3975703e-1, &
           .4007093e-1,.4045533e-1,.4091481e-1,.4145507e-1,.4208311e-1, &
           .4280748e-1,.4363863e-1,.4458932e-1,.4567523e-1,.4691571e-1, &
           .4833487e-1,.4996298e-1,.5183859e-1,.5401138e-1,.5654656e-1, &
           .5953130e-1,.6308489e-1,.6737503e-1,.7264544e-1,.7926471e-1, &
           .8781922e-1,.9930398e-1,.1155599,.1404344,.1836142,.2790016, &
           .7010474/)
      if( ran_genreal(gendata, u, err) == RETURN_FAIL ) go to 800
      s = 0.0
      if (u.gt.0.5) s = 1.0
      u = u + u - s
      u = 32.0*u
      i = int(u)
      if (i.eq.32) i = 31
      if( i .ne. 0 ) then
         ustar = u - float(i)
         aa = a(i)
         loopa: do
            if( ustar > t(i) ) then
               w = (ustar-t(i))*h(i)
               exit loopa
            else
               if( ran_genreal(gendata, u, err) == RETURN_FAIL ) &
                    goto 800
               w = u* (a(i+1)-aa)
               tt = (0.5*w+aa)*w
               loopb: do
                  if( ustar > tt ) then
                     y = aa + w
                     val = y
                     if (s == 1.0) val = -y
                     go to 700
                  end if
                  if( ran_genreal(gendata, u, err) == RETURN_FAIL ) &
                       goto 800
                  if (ustar < u) exit loopb
                  tt = u
                  if( ran_genreal(gendata, ustar, err) == RETURN_FAIL ) &
                       goto 800
               end do loopb
               if( ran_genreal(gendata, ustar, err) == RETURN_FAIL ) &
                    goto 800
            end if
         end do loopa
      else
         i = 6
         aa = a(32)
         loopc: do
            u = u + u
            if( u >= 1.0 ) exit loopc
            aa = aa + d(i)
            i = i + 1
         end do loopc
         u = u - 1.0
         loopd: do
            w = u*d(i)
            tt = (0.5*w+aa)*w
            loope: do
               if( ran_genreal(gendata, ustar, err) == RETURN_FAIL ) &
                    goto 800
               if (ustar > tt) then
                  y = aa + w
                  val = y
                  if (s == 1.0) val = -y
                  go to 700
               end if
               if( ran_genreal(gendata, u, err) == RETURN_FAIL ) &
                    goto 800
               if (ustar < u) exit loope
               tt = u
            end do loope
            if( ran_genreal(gendata, u, err) == RETURN_FAIL ) &
                 goto 800
         end do loopd
      end if
      y = aa + w
      val = y
      if (s == 1.0) val = -y
      ! normal exit
700   continue
      ran_snorm = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_snorm = RETURN_FAIL
      return
   end function ran_snorm
   !##################################################################
   integer(our_int) function ran_sgamma(gendata, a, val, err)
      !	Generates a standard gamma random variate
      implicit none
      ! Required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(in) :: a
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_sgamma"
      real(kind=our_sgle) :: b, b0, c, d, e, &
           p, q, q0, r, s, s2, si, t, u, v, w, x
      intrinsic abs, alog, exp, sign, sqrt
      real(kind=our_sgle), parameter :: q1=.04166669, q2=.02083148, &
           q3=.00801191, q4=.00144121, q5=-.00007388, q6=.00024511, &
           q7=.00024240, a1=.3333333, a2=-.2500030, a3=.2000062, &
           a4=-.1662921, a5=.1423657, a6=-.1367177, a7=.1233795, &
           e1=1.0, e2=.4999897, e3=.1668290, e4=.0407753, e5=.0102930
      real(kind=our_sgle), parameter :: sqrt32=5.656854249
      real(kind=our_sgle) :: val1
      ! do it
      if( a < 1.0 ) then
         b0 = 1.0 + .3678794*a
         loopa: do
            if( ran_genreal(gendata,val1,err) == RETURN_FAIL ) goto 800
            p = b0*val1
            if( p < 1.0 ) then
               val = exp(alog(real(p))/a)
               if( ran_sexp(gendata,val1,err) == RETURN_FAIL ) goto 800
               if( val <= val1 ) go to 700
            else
               val = -alog(real((b0-p)/a))
               if( ran_sexp(gendata,val1,err) == RETURN_FAIL ) goto 800
               if( (1.0-a)*alog(real(val)) <= val1 ) goto 700
            end if
         end do loopa
      end if
      !	1) Recalculations of s2,s,d if a has changed
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
      
      !	2) t=standard normal deviate, x=(s,1/2)-normal deviate.
      if( ran_snorm(gendata, t, err) == RETURN_FAIL ) go to 800
      x = s + 0.5*t
      val = x*x
      if( t >= 0.0 ) go to 700
      
      !	3) u= 0,1 -uniform sample. squeeze acceptance (s)
      if( ran_genreal(gendata,u,err) == RETURN_FAIL ) go to 800
      if( d*u <= t*t*t) go to 700
      
      !	4) recalculations of q0,b,si,c
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
      !	Approximation depending on size of parameter a.
      !	The constants in the expressions for b, si and
      !	c were established by numerical experiments.
      if( a > 3.686 ) then
         if( a > 13.022 ) then
            b = 1.77
            si = .75
            c = .1515/s
         else
            b = 1.654 + .0076*s2
            si = 1.68/s + .275
            c = .062/s + .024
         end if
      else
         b = .463 + s + .178*s2
         si = 1.235
         c = .195/s - .079 + .16*s
      end if
      
      !	5) no quotient test if x not positive
      if( x > 0 ) then
         !	6) calculation of v and quotient q
         v = t/ (s+s)
         if( abs(v) > 0.25 ) then
            q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(real(1.0+v))
         else
            q = q0 + 0.5*t*t* &
                 ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
         end if
         !	7) quotient acceptance (q)
         if( alog(real(1.0-u)) <= q ) go to 700
      end if
      
      !	8) e=standard exponential deviate
      !		u= 0,1 -uniform deviate
      !		t=(b,si)-double exponential (laplace) sample
      loopb: do
         if( ran_sexp(gendata, e, err) == RETURN_FAIL ) go to 800
         !e = ran_sexp(gendata,err)
         if( ran_genreal(gendata,u,err) == RETURN_FAIL ) go to 800
         u = u + u - 1.0
         t = b + sign(si*e,u)
         !	9) rejection if t .lt. tau(1) = -.71874483771719
         if( t >= (-.7187449) ) then
            !	10) calculation of v and quotient q
            v = t/ (s+s)
            if( abs(v) > 0.25 ) then
               q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(real(1.0+v))
            else
               q = q0 + 0.5*t*t* &
                    ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
            end if
            !	11) Hat acceptance (h) (if q not positive go to step 8)
            if( q > 0.0 ) then
               if( q > 0.5 ) then
                  if( q >= 15.0 ) then
                     if( (q+e-0.5*t*t) > 87.49823 ) exit loopb
                     if( c*abs(u) <= exp(q+e-0.5*t*t) ) exit loopb
                  else
                     w = exp(q) - 1.0
                     if( c*abs(u) <= w*exp(e-0.5*t*t) ) exit loopb
                  end if
               else
                  w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
                  if( c*abs(u) <= w*exp(e-0.5*t*t) ) exit loopb
               end if
            end if
         end if
      end do loopb
      x = s + 0.5*t
      val = x*x
      ! normal exit
700   continue
      ran_sgamma = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_sgamma = RETURN_FAIL
      return
   end function ran_sgamma
   !####################################################################
   integer(our_int) function ran_genreal(gendata, val, err)
      !	Returns a random floating point number from a uniform 
      !	distribution over 0 - 1 (endpoints of this interval are 
      !	not returned) using the current generator.
      !	4.656613057E-10 is 1/M1.  M1 is set in inrgcm to 2147483563.
      !	If M1 changes, change this also.
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      integer(kind=our_int) :: ival
      ! call integer generator and translate result to (0,1)
      if( ran_genint(gendata, ival, err) == RETURN_SUCCESS ) then
         val = ival*4.656613057e-10 
         ran_genreal = RETURN_SUCCESS
      else
         val = 0.0
         ran_genreal = RETURN_FAIL
      end if
   end function ran_genreal
   !####################################################################
   integer(our_int) function ran_genint(gendata, ival, err)
      !	Returns a random integer following a uniform distribution over
      !	(1, 2147483562) using the current generator.  Uses the 
      !	following module variables:  a1a, a1vw, a1w, a2a, a2vw, a2w, m1,
      !	m2, cg1, cg2.
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      integer(kind=our_int), intent(out) :: ival
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_genint"
      integer(kind=our_int) :: currentg, k, s1, s2, z
      ! do it
      if (.not. (gendata%qqssd) ) goto 700
      currentg = gendata%curntg
      s1 = gendata%cg1(currentg)
      s2 = gendata%cg2(currentg)
      k = s1/53668
      s1 = a1a* (s1-k*53668) - k*12211
      if ( s1 < 0 ) s1 = s1 + m1
      k = s2/52774
      s2 = a2a* (s2-k*52774) - k*3791
      if ( s2 < 0 ) s2 = s2 + m2
      gendata%cg1(currentg) = s1
      gendata%cg2(currentg) = s2
      z = s1 - s2
      if ( z < 1 ) z = z + m1 - 1
      if (gendata%qanti(currentg)) z = m1 - z
      ! normal exit
      ival = z
      ran_genint = RETURN_SUCCESS
      return
      ! error traps
700   continue
      call err_handle(err, 1, &
           comment = "Random generator seeds have not been set")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ival = 0
      ran_genint = RETURN_FAIL
      return
   end function ran_genint
   !####################################################################
   integer(our_int) function ran_setall(gendata, iseed1, iseed2, err)
      !	Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
      !	initial seeds of the other generators are set accordingly, and
      !	all generators states are set to these seeds.
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      integer(kind=our_int), intent(in) :: iseed1,iseed2
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_setall"
      integer(kind=our_int) :: g, ocgn
      ! do it
      gendata%qqssd = .true.
      ocgn = gendata%curntg   ! save current generator
      gendata%ig1(1) = iseed1
      gendata%ig2(1) = iseed2
      if( initgn(gendata, -1,err) == RETURN_FAIL) goto 800
      do g = 2, NUM_GEN
         if( mltmod(a1vw,gendata%ig1(g-1),m1,gendata%ig1(g),err) &
              == RETURN_FAIL ) go to 800
         if( mltmod(a2vw,gendata%ig2(g-1),m2,gendata%ig2(g),err) &
              == RETURN_FAIL ) go to 800
         gendata%curntg = g   ! Set current generator number
         if( initgn(gendata, -1, err) == RETURN_FAIL) goto 800
      end do
      gendata%curntg = ocgn   ! reset current generator
      ! normal exit
      ran_setall = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ran_setall = RETURN_FAIL
      return
   end function ran_setall
   !####################################################################
   integer(our_int) function mltmod(a, s, m, ival, err)
      !   Requires: 
      !       0 < a <= m
      !       0 < s <= m
      !	Returns (a*s) mod m
      !
      !	  H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      !	   machine. On a different machine recompute H
      implicit none
      ! required arguments
      integer(kind=our_int), intent(in) :: a, m, s
      integer(kind=our_int), intent(out) :: ival
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "mltmod"
      integer(kind=our_int), parameter :: h=32768
      integer(kind=our_int) :: a0, a1, k, p, q, qh, rh
      ! do it
      if ( a<=0 .or. a>=m .or. s<=0 .or. s>=m ) goto 800
      if ( a < h ) then
         a0 = a
         p = 0
      else
         a1 = a/h
         a0 = a - h*a1
         qh = m/h
         rh = m - h*qh
         if ( a1 >= h ) then
            a1 = a1 - h
            k = s/qh
            p = h* (s-k*qh) - k*rh
            loopa: do while( p < 0 )
               p = p + m
            end do loopa
         else
            p = 0
         end if
         if ( a1 /= 0 ) then
            q = m/a1
            k = s/q
            p = p - k* (m-a1*q)
            if (p.gt.0) p = p - m
            p = p + a1* (s-k*q)
            loopb: do while( p < 0 )
               p = p + m
            end do loopb
         end if
         k = p/qh
         p = h* (p-k*qh) - k*rh
         loopc: do while( p < 0 )
            p = p + m
         end do loopc
      end if
      if( a0 /= 0 ) then
         q = m/a0
         k = s/q
         p = p - k* (m-a0*q)
         if (p.gt.0) p = p - m
         p = p + a0* (s-k*q)
         loopd: do while( p < 0 )
            p = p + m
         end do loopd
      end if
      ival = p
      ! normal exit
      mltmod = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      call err_handle(err, 1, &
           comment="Parameters out of order")
      ival = 0
      mltmod = RETURN_FAIL
      return
   end function mltmod
   !####################################################################
   integer(our_int) function initgn(gendata, isdtyp, err)
      !	  Reinitializes the state of the current generator
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      integer(kind=our_int), intent(in) :: isdtyp
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "initgn"
      if ( isdtyp == -1 ) then
         gendata%lg1(gendata%curntg) = gendata%ig1(gendata%curntg)
         gendata%lg2(gendata%curntg) = gendata%ig2(gendata%curntg)
      else if ( isdtyp == 1 ) then
         if( mltmod(a1w,gendata%lg1(gendata%curntg),m1,&
              gendata%lg1(gendata%curntg),err) &
              == RETURN_FAIL ) go to 800
         if( mltmod(a2w,gendata%lg2(gendata%curntg),m2,&
              gendata%lg2(gendata%curntg),err) &
              == RETURN_FAIL ) go to 800
      end if
      gendata%cg1(gendata%curntg) = gendata%lg1(gendata%curntg)
      gendata%cg2(gendata%curntg) = gendata%lg2(gendata%curntg)
      ! normal exit
      initgn = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment="Could not initialize random seeds")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      initgn = RETURN_FAIL
      return
   end function initgn
   !####################################################################
   integer(our_int) function ran_sexp(gendata, val, err)
      !	Generates a value from standard exponential distribution
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      real(kind=our_sgle), intent(out) :: val
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_sexp"
      real(kind=our_sgle) :: a, u, umin, ustar
      integer(kind=our_int) :: i
      real(kind=our_sgle), dimension(8), parameter :: q = (/.6931472, &
           .9333737,.9888778,.9984959,.9998293,.9999833,.9999986,&
           .9999999/)

      a = 0.0
      if( ran_genreal(gendata,u,err) == RETURN_FAIL ) goto 800
      
      loopa: do
         u = u + u
         if( u >= 1.0 ) exit loopa
         a = a + q(1)
      end do loopa
      
      u = u - 1.0
      if( u <= q(1) ) then
         val = a + u
         goto 700
      end if
      i = 1
      if( ran_genreal(gendata,ustar,err) == RETURN_FAIL ) goto 800
      umin = ustar
      
      loopb: do
         if( ran_genreal(gendata,ustar,err) == RETURN_FAIL ) goto 800
         if (ustar < umin) umin = ustar
         i = i + 1
         if ( u <= q(i) ) exit loopb
      end do loopb
      val = a + umin*q(1)
      
      ! normal exit
700   continue
      ran_sexp = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment="Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      val = 0.0
      ran_sexp = RETURN_FAIL
      return
   end function ran_sexp
   !####################################################################
   integer(our_int) function ran_timeseed(gendata, err)
      !	Yields two seed values based on the intrinsic function 
      ! "date_and_time".  The first is the number of 
      ! milliseconds since 0:00:00.000 on 1/1/2000, and the 
      ! second is a permutation of the first, obtained by 
      ! swapping the top 15 bits with the bottom 16.
      !	Range: 1 to (2**31-1).  The generated seeds will repeat
      !	every 24 days, 20 hours, 31 minutes, 23.647 seconds.
      !	Created by David R. Lemmon, Penn State University, 
      ! Sept 2000.
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_timeseed"
      integer(kind=our_int) :: isd1, isd2
      integer(kind=our_int), dimension(8) :: ivalues
      integer(kind=our_int) :: i
      real(kind=our_sgle) :: ms, monthdays(1:12)
      real(kind=our_sgle), parameter :: themax = (2.0**31 - 1.0)
      intrinsic date_and_time, int, mod
      monthdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      !	Get date/time info
      call date_and_time( values = ivalues )
      ms = 0
      ! accumulate ms since 0:00:00.000 on 1/1/2000
      ! Years:
      do i=2001,ivalues(1)
         if(is_leap_year(i)) then
            ms = ms + 366.0*24.0*3600000.0
         else
            ms = ms + 365.0*24.0*3600000.0
         end if
      end do
      ! Months:
      if(is_leap_year(ivalues(1))) then
         monthdays(2) = 29
      end if
      do i=1,(ivalues(2)-1)
         ms = ms + monthdays(i)*24.0*3600000.0
      end do
      ! Days:
      ms = ms + (ivalues(3)-1)*24.0*3600000.0
      ! Hours, minutes, seconds, milliseconds:
      ms = ms + ((ivalues(5)*60 + ivalues(6))*60 + &
           ivalues(7))*1000 + ivalues(8)
      ms = mod(ms, themax)
      if(ms == 0) ms = themax   ! highly unlikely, but possible
      isd1 = int(ms, our_int)
      ! For the second seed, take the top 15 bits and swap
      ! with the bottom 16 bits.
      isd2 = mod(isd1,2**16)*(2**15) + int(isd1/(2**15))
      ! set the seeds
      if( ran_setall(gendata, isd1, isd2, err) &
           == RETURN_FAIL ) goto 900
      ! normal exit
      ran_timeseed = RETURN_SUCCESS
      return
      ! error traps
900   continue
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      call err_handle(err, 1, &
           comment = "Operation failed")
      ran_timeseed = RETURN_FAIL
      return
   end function ran_timeseed
   !####################################################################
   integer(our_int) function ran_phrsd(gendata, phrase, err)
      !	Uses a phrase (character string) to generate two integer seeds
      !	which are then used to set the random number generator.
      ! Calling program should use "trim" to remove trailing spaces.
	  ! WARNING:  This algorithm produces the same seed values for strings
	  !  that are simply re-arranged version of each other.  For 
	  !  example, "test123" and "test321" produce the same seed values.
	  !
      implicit none
      ! required arguments
      type(random_gendata), intent(inout) :: gendata
      character(len=*), intent(in) :: phrase
      type(error_type), intent(inout) :: err
      ! local variables and parameters
      character(len=*), parameter :: subname = "ran_phrasd"
      integer(kind=our_int) :: seed1 ,seed2
      character(len=*), parameter :: table='abcdefghijklmnopqrstuvwxyz'//  &
           'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//'0123456789'// &
           '!@#$%^&*()_+[];:''"<>?,./'
      integer(kind=our_int), parameter :: twop30=1073741824
      integer(kind=our_int) :: i, ichr, j, lphr
      integer(kind=our_int) :: values(5)
      integer(kind=our_int), parameter :: &
           shift(0:4)=(/1,64,4096,262144,16777216/)
      intrinsic index,mod
      ! make sure that the string is not blank
      if( len(phrase) == 0 ) goto 800
      ! define the integer seeds
      seed1 = 1234567890
      seed2 = 123456789
      lphr = len(phrase)
      do i = 1,lphr
         ichr = mod(index(table,phrase(i:i)),64)
         if (ichr == 0) ichr = 63
         do j = 1,5
            values(j) = ichr - j
            if (values(j).lt.1) values(j) = values(j) + 63
         end do
         do j = 1,5
            seed1 = mod(seed1+shift(j-1)*values(j),twop30)
            seed2 = mod(seed2+shift(j-1)*values(6-j),twop30)
         end do
      end do
      ! set the seeds
      if( ran_setall(gendata, seed1, seed2, err) &
           == RETURN_FAIL ) goto 900
      ! normal exit
      ran_phrsd = RETURN_SUCCESS
      return
      ! error traps
800   continue
      call err_handle(err, 1, &
           comment = "A zero-length string is invalid")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ran_phrsd = RETURN_FAIL
      return
900   continue
      call err_handle(err, 1, &
           comment = "Operation failed")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ran_phrsd = RETURN_FAIL
      return
   end function ran_phrsd
   !####################################################################
   integer(kind=our_int) function ran_seed()
      !	Returns a positive-valued, randomly generated, 31-bit seed value
      !   by calling the intrinsic function random_seed
      implicit none
      integer, dimension(:), allocatable :: theseed
      integer :: idim
      call random_seed()
      call random_seed( size = idim )
      allocate (theseed(idim))
      call random_seed( get = theseed(1:idim) )
      ran_seed = theseed(1)
      deallocate (theseed)
   end function ran_seed
   !####################################################################
   logical function is_leap_year(iyear)
      !	Returns .true. if iyear is a leap year.  Should work for the
      !	next several centuries.
      implicit none
      integer(kind=our_int), intent(in) :: iyear
      is_leap_year = ((mod(iyear, 4)== 0) &
           .and. (mod(iyear,100) /= 0) &
           .or. (mod(iyear,400) == 0))
   end function is_leap_year
   !####################################################################
   logical function ran_seed_is_set(gendata)
      ! Returns .true. if random generator seed values have been set
      implicit none
      type(random_gendata) :: gendata
      ran_seed_is_set = gendata%qqssd
   end function ran_seed_is_set
   !####################################################################
end module random_generator
!#######################################################################
