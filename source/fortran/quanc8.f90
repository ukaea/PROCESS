      subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
!
!      double precision fun, a, b, abserr, relerr, result, errest, flag
      double precision a, b, abserr, relerr, result, errest, flag
      integer nofun
      external fun
!
!   estimate the integral of fun(x) from a to b
!   to a user provided tolerance.
!   an automatic adaptive routine based on
!   the 8-panel newton-cotes rule.
!
!   input ..
!
!   fun     the name of the integrand function subprogram fun(x).
!   a       the lower limit of integration.
!   b       the upper limit of integration.(b may be less than a.)
!   relerr  a relative error tolerance. (should be non-negative)
!   abserr  an absolute error tolerance. (should be non-negative)
!
!   output ..
!
!   result  an approximation to the integral hopefully satisfying the
!           least stringent of the two error tolerances.
!   errest  an estimate of the magnitude of the actual error.
!   nofun   the number of function values used in calculation of result.
!   flag    a reliability indicator.  if flag is zero, then result
!           probably satisfies the error tolerance.  if flag is
!           xxx.yyy , then  xxx = the number of intervals which have
!           not converged and  0.yyy = the fraction of the interval
!           left to do when the limit on  nofun  was approached.
!
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),f(16),x(16),fsave(8,30),xsave(8,30)
      double precision dabs,dmax1
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
!
!   ***   stage 1 ***   general initialization
!   set constants.
!
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
!
!   trouble when nofun reaches nofin
!
      w0 =   3956.0d0 / 14175.0d0
      w1 =  23552.0d0 / 14175.0d0
      w2 =  -3712.0d0 / 14175.0d0
      w3 =  41984.0d0 / 14175.0d0
      w4 = -18160.0d0 / 14175.0d0
!
!   initialize running sums to zero.
!
      flag = 0.0d0
      result = 0.0d0
      cor11  = 0.0d0
      errest = 0.0d0
      area   = 0.0d0
      nofun = 0
      if (a .eq. b) return
!
!   ***   stage 2 ***   initialization for first interval
!
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev  = 0.0d0
      f0 = fun(x0)
      stone = (b - a) / 16.0d0
      x(8)  =  (x0  + x(16)) / 2.0d0
      x(4)  =  (x0  + x(8))  / 2.0d0
      x(12) =  (x(8)  + x(16)) / 2.0d0
      x(2)  =  (x0  + x(4))  / 2.0d0
      x(6)  =  (x(4)  + x(8))  / 2.0d0
      x(10) =  (x(8)  + x(12)) / 2.0d0
      x(14) =  (x(12) + x(16)) / 2.0d0
      do 25 j = 2, 16, 2
         f(j) = fun(x(j))
   25 continue
      nofun = 9
!
!   ***   stage 3 ***   central calculation
!   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
!   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
!
   30 x(1) = (x0 + x(2)) / 2.0d0
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.0d0
         f(j) = fun(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0d0
      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6)) &
       + w3*(f(3)+f(5))  +  w4*f(4)) * step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14)) &
       + w3*(f(11)+f(13)) + w4*f(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
!
!   ***   stage 4 *** interval convergence test
!
      esterr = dabs(qdiff) / 1023.0d0
      tolerr = dmax1(abserr,relerr*dabs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
!
!   ***   stage 5   ***   no convergence
!   locate next interval.
!
   50 nim = 2*nim
      lev = lev+1
!
!   store right hand elements for future use.
!
      do 52 i = 1, 8
         fsave(i,lev) = f(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
!
!   assemble left hand elements for immediate use.
!
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         f(2*j+18) = f(j+9)
         x(2*j+18) = x(j+9)
   55 continue
      go to 30
!
!   ***   stage 6   ***   trouble section
!   number of function values is about to exceed limit.
!
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
!
!   current level is levmax.
!
   62 flag = flag + 1.0d0
!
!   ***   stage 7   ***   interval converged
!   add contributions into running sums.
!
   70 result = result + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0d0
!
!   locate next interval.
!
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
!
!   assemble elements required for the next interval.
!
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8
         f(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      go to 30
!
!   ***   stage 8   ***   finalize and return
!
   80 result = result + cor11
!
!   make sure errest not less than roundoff level.
!
      if (errest .eq. 0.0d0) return
   82 temp = dabs(result) + errest
      if (temp .ne. dabs(result)) return
      errest = 2.0d0*errest
      go to 82
      end
