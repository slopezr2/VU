!###############################################################################
!
! Chemistry solver: TwoStep
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Chem_Solver_TwoStep

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  TwoStep


  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'LE_Chem_Solver_TwoStep'


contains



  !======================================================================
  !
  !     AUTHOR: Jan Verwer, Centre for Mathematics and Computer Science
  !     (CWI), Kruislaan 413, 1098 SJ Amsterdam, The Netherlands
  !
  !     EMAIL ADDRESS: janv@cwi.nl
  !
  !----------------------------------------------------------------------
  !
  !     PURPOSE: TWOSTEP is designed for the numerical solution of
  !     stiff ODE systems
  !
  !       dy(t)/dt = yp(t,y) - yl(y,t)*y(t)
  !
  !     originating from atmospheric chemistry. Its underlying integration
  !     method is the implicit, 2-nd order, 2-step BDF formula. A simple
  !     explicit Gauss-Seidel technique is used for approximating the
  !     implicitly defined BDF solutions, rather than the usual modified
  !     Newton method [1]. By this approach TWOSTEP is explicit, whereas
  !     the excellent stability of the BDF method is maintained. Also the
  !     Jacobi technique may be appropriate. In general, however, its use
  !     generally decreases the stability of TWOSTEP and hence also the
  !     efficiency.
  !
  !                     y(t+dt) - y(t)
  !     Euler implicit: -------------- = yp - yl*y(t+dt)
  !                          dt
  !
  !      y(t+dt) = y(t) + dt*yp - dt*yl*y(t+dt)
  !                                                        y(t) + dt*yp
  !      y(t+dt)*[1 + dt*yl] = y(t) + dt*yp <=> y(t+dt) = -------------- .
  !                                                           1 + dt*yl
  !
  !     The general solution of the differential equation dy/dt = f((t,y), using BDF2,
  !     is: y(t+dt) = {a1*y(t) + a2*y(t-dt)} + g*dt*f(t,y).
  !
  !     Define ysum = a1*y(t) + a2*y(t-dt), then the BDF2-solution is
  !     (analoguous to the Euler implicit solution):
  !
  !               ysum + g*dt*yp
  !     y(t+dt) = --------------- .
  !                 1 + g*dt*yl
  !
  !     Taylor expansion in dt and setting all terms of order < 3 to zero
  !     leads to a set of 3 equations for the coefficients a1, a2 and g,
  !     with solution a1 = (c+1)**2/(c**2 + 2c), a2 = -1/(c**2 + 2c), g = (c+1)/c+2),
  !     and c = dtold/dt, the ratio of the last two time steps.
  !
  !     While in the prototype solver discussed in [1] the number of
  !     iterations is determined by a convergence criterion, TWOSTEP
  !     works with a fixed, a priori given number of iterations. This
  !     leads to a somewhat simpler code and in various experiments
  !     this more simple strategy has turned out to be as efficient as
  !     the iteration stategy [2]. In fact, TWOSTEP may work with
  !     only a single iteration per integration step. In this situation
  !     TWOSTEP uses one ODE evaluation per time step, which means that
  !     its workload then becomes even less than popular solvers based
  !     on the quasi-steady-state-approach (QSSA) [3]. Note that when
  !     only one or a few Gauss-Seidel iterations are used, the order
  !     of the components influences the accuracy. To our experience,
  !     this influence is minor.
  !
  !     [1] J.G. Verwer, Gauss-Seidel iteration for stiff ODEs from
  !     chemical kinetics, Report NM-R9315, CWI, Amsterdam, 1993
  !     (to appear in SIAM Journal on Scientific Computing)
  !
  !     [2] J.G. Verwer & D. Simpson, Explicit
  !     methods for stiff ODEs from atmospheric chemistry. Report
  !     NM-R9409, CWI, Amsterdam, 1994.
  !
  !     [3] J.G. Verwer & M. van Loon,  An evaluation of explicit pseudo-
  !     steady-state approximation schemes for stiff ODE systems from
  !     chemical kinetics, Report NM-R9312, CWI, Amsterdam, 1993 (to
  !     appear in J. Comput. Phys.).
  !
  !----------------------------------------------------------------------
  !
  !     MEANING OF PARAMETERS:
  !
  !     n       - Integer. Number of components.
  !     t       - Real. The independent variable time.
  !     te      - Real. The endpoint of time.
  !     dt      - Real. The stepsize.
  !     dtmin   - Real. A minimum for dt.
  !     dtmax   - Real. A maximum for dt.
  !     yold    - Real array (n). Solution at previous time point.
  !     y       - Real array (n). Solution at current time point.
  !     ynew    - Real array (n). Solution at forward time point.
  !     yp      - Real array (n). Storage for the production term.
  !     yl      - Real array (n). Storage for the loss term.
  !     Q       - Real array (n). emission array.
  !     ysum    - Real array (n). Workarray for the BDF2 method.
  !     atol    - Real array (n). Absolute tolerances.
  !     rtol    - Real array (n). Relative tolerances.
  !     method  - Integer. Determines whether Gauss-Seidel or Jacobi
  !               iteration is used; method only occurs as parameter
  !               of subroutine ITER.
  !     numit   - Integer. Number of Gauss-Seidel or Jacobi iterations.
  !     nfcn    - Integer. The total number of (function) calls of ITER.
  !     naccpt  - Integer. The number of accepted integration steps.
  !     nrejec  - Integer. The number of rejected integration steps.
  !     nstart  - Integer. The number of restarts + 1.
  !     startdt - Real. The last computed starting stepsize.
  !
  !----------------------------------------------------------------------
  !
  !     STORAGE: 9 real arrays of dimension n.
  !
  !----------------------------------------------------------------------
  !
  !     INPUT:
  !     n   - Number of components.
  !     t       - Initial time; t is changed.
  !     te      - The endpoint of time.
  !     dtmin   - The minimal stepsize that TWOSTEP is allowed to use.
  !     dtmax   - The maximal stepsize that TWOSTEP is allowed to use.
  !
  !               If on input dtmin = dtmax, then dt:=dtmin and dt is kept
  !               fixed throughout the integration. In this case, the
  !               length  of the integration interval te-t must be an
  !               integer multiple >=2 of dt.
  !
  !               If dtmin = dtmax, the stepsize control is switched off.
  !               The user thus should ascertain that the integration
  !               process remains stable for the selected stepsize.
  !
  !               If dtmax > dtmin, then stepsize control is carried out
  !               and TWOSTEP determines an initial stepsize itself.
  !
  !     y       - Initial solution vector; y is changed.
  !     atol    - Absolute tolerances.
  !     rtol    - Relative tolerances.
  !     method  - method = 0 for the Gauss-Seidel technique.
  !               method = 1 for the Jacobi technique.
  !
  !               The Gauss-Seidel technique is recommended. The user
  !               must program the selected technique in his own
  !               subroutine ITER which defines also the ODE system
  !               (see below for an example).
  !
  !     numit   - The number of the Gauss-Seidel or Jacobi iterations used
  !               per time step. This number thus is a priori described
  !               for the whole integration. A low number is recommended,
  !               numit = 1,2 or 3, say. In the rural case of [2], numit = 1
  !               has been chosen. In the urban case of [2], numit = 2
  !               turned out to be more efficient.
  !
  !               Note that for the Gauss-Seidel technique the order of
  !               the components within ITER can play an important role.
  !               It is recommended to order the components in decreasing
  !               of the loss rates. Hence first the component with the
  !               largest loss rate, etc. If loss rates are close, then
  !               the ordering need not be strict.
  !
  !----------------------------------------------------------------------
  !
  !     OUTPUT:
  !
  !     t       - t = te.
  !     dtold   - The last stepsize value used.
  !     y       - The computed solution at t = te.
  !     nfcn, naccpt, nrejec, nstart, startdt - See MEANING of parameters.
  !
  !----------------------------------------------------------------------
  !
  !     SUBROUTINES:
  !
  !     TWOSTEP calls three subroutines, viz. NEWDT, FIT and ITER. Only
  !     subroutine ITER is to be defined by the user.
  !
  !     NEWDT  - Computes the new stepsize. NEWDT itself calls FIT.
  !     FIT    - May slightly adjust dt to guarantee that the remainder
  !              of the integration interval is an integer multiple of
  !              the current stepsize. The adjustment is carried as soon
  !              as (te-t)/dt <= 10.0. Hence the adjustment may lead to
  !              a stepsize slightly smaller than dtmin for approximately
  !              ten integration steps.
  !     ITER   - A user defined routine for the ODE system. Within ITER
  !              also the Gauss-Seidel or Jacobi technique is to be
  !              implemented by the user, as exemplified below:
  !
  !              subroutine ITER(t,method,n,y,ysum,gdt,yp,yl,Q)
  !              integer    method,n,i,
  !              real*8     t,y(n),ysum(n),gdt,yp(n),yl(n),Q(n)
  !
  !              if (method.eq.0) then
  !              --------------------------------------------------
  !              The following statements show how to program the
  !              Gauss-Seidel technique. Note that array y is
  !              changed. The ispec integers define the order
  !              in which the components are processed.
  !              -------------------------------------------------
  !              i=ispec1
  !              yp(i)= USER DEFINED(t,y)
  !              yl(i)= USER DEFINED(t,y)
  !              y(i)=(ysum(i)+gdt*yp(i))/(1.0+gdt*yl(i))
  !
  !              i=ispec2
  !              yp(i)= USER DEFINED(t,y)
  !              yl(i)= USER DEFINED(t,y)
  !              y(i)=(ysum(i)+gdt*yp(i))/(1.0+gdt*yl(i))
  !
  !              etcetera
  !
  !              else if (method.eq.1) then
  !              ------------------------------------------------
  !              This shows how to implement the Jacobi technique.
  !              Now the order of the components is irrelevant.
  !              Note that again y is changed, but only in the
  !              second do-loop.
  !              ------------------------------------------------
  !              for i=1,n
  !              yp(i)= USER DEFINED(t,y)
  !              yl(i)= USER DEFINED(t,y)
  !              enddo
  !
  !              for i=1,n
  !              y(i)=(ysum(i)+gdt*yp(i))/(1.0+gdt*yl(i))
  !              enddo
  !              endif
  !              end of ITER
  !
  !----------------------------------------------------------------------

  subroutine TWOSTEP( n, t, te, dtmin, dtmax, &
                      y, Q, &
#ifdef with_vbs
                      drog, &
#endif
                      nreac, rk, &
                      naux, aux, &
                      atol, rtol, method, numit, &
                      nfcn, naccpt, nrejec, nstart, startdt, &
                      Iter, status )
    ! --- in/out -----------------------------

    integer, intent(in)       ::  n
    real, intent(inout)       ::  t
    real, intent(in)          ::  te
    real, intent(in)          ::  dtmin, dtmax
    real, intent(inout)       ::  y(n)
    real, intent(in)          ::  Q(n)
#ifdef with_vbs
    real, intent(out)         ::  drog(n)
#endif
    integer, intent(in)       ::  nreac
    real(8), intent(in)       ::  rk(nreac)
    integer, intent(in)       ::  naux
    real(8), intent(in)       ::  aux(naux)
    real, intent(in)          ::  atol(n), rtol(n)
    integer, intent(in)       ::  method
    integer, intent(in)       ::  numit
    integer, intent(out)      ::  nfcn, naccpt, nrejec, nstart
    real, intent(out)         ::  startdt
    integer, intent(out)      ::  status

    ! --- interfaces -----------------------------

    interface
      pure  subroutine Iter( n, y, ysum, dt, yp, yl, Q, &
#ifdef with_vbs
                              drog_loss_temp, drog_prod_neg_temp, &
#endif
                              nreac, rk, &
                              naux, aux )
        integer, intent(in)    ::  n
        real, intent(inout)    ::  y(n)
        real, intent(in)       ::  ysum(n)
        real, intent(in)       ::  dt
        real, intent(out)      ::  yp(n), yl(n)
        real, intent(in)       ::  Q(n)
#ifdef with_vbs
        real, intent(out)      ::  drog_loss_temp(n)
        real, intent(out)      ::  drog_prod_neg_temp(n)
#endif        
        integer, intent(in)    ::  nreac
        real(8), intent(in)    ::  rk(nreac)
        integer, intent(in)    ::  naux
        real(8), intent(in)    ::  aux(naux)
      end subroutine Iter
    end interface

    ! --- const ----------------------------------

    character(len=*), parameter  :: rname = mname//'/TWOSTEP'

    ! --- local ----------------------------------

    integer   ::  i, j
    real      ::  dt
    real      ::  yold(n), ynew(n)
    real      ::  ysum(n)
    real      ::  yp(n), yl(n)
    real      ::  ytol
    real      ::  ratio
    real      ::  dtold
    real      ::  a1, a2
    real      ::  c, cp1
    real      ::  gdt
    real      ::  errlte
    real(8)   ::  dy
    logical   ::  accept,failer,restart
    logical   ::  goto10
#ifdef with_vbs
    real      ::  drog_loss_temp(n)
    real      ::  drog_prod_neg_temp(n)
    real      ::  drog_old(n)
#endif

    ! --- begin ----------------------------------

    ! Initialization of counters, etc.

    naccpt=0
    nrejec=0
    nfcn=0
    nstart=0
    startdt = 0.0
    status = 0

    failer=.false.
    restart=.false.
    accept=.true.

    ! init temporary concentration arrays:
    yold = y
    ynew = y
    ysum = y
    yp = 0.0
    yl = 0.0
    
#ifdef with_vbs
    drog_loss_temp     = 0.0    
    drog_prod_neg_temp = 0.0
    drog_old           = 0.0
    drog               = 0.0
#endif    
    ! loop '10' ...
    do

      ! Initial stepsize computation.
      if ( dtmin == dtmax ) then
        nstart=1
        dt=dtmin
      else
        nstart=nstart+1
        dt=0.0
        call iter( n, y, ysum, dt, yp, yl, Q, &
#ifdef with_vbs
                     drog_loss_temp, drog_prod_neg_temp, &
#endif                    
                     nreac, rk, naux, aux )
        nfcn = nfcn+1
        dt=te-t
        do i=1,n
          ytol=atol(i)+rtol(i)*abs(y(i))
          dy=yp(i)-y(i)*yl(i)
          if (dy.ne.0.0) dt=min(dt,ytol/abs(dy))
        end do
        if (restart) dt=dt/10.0
        restart=.true.
        dt=max(dtmin,min(dt,dtmax))
        call FIT(t,te,dt)
        if (dt.eq.te-t) dt=dt/2.0
        startdt=dt
      endif

      ! The starting step is carried out, using the implicit Euler method.
      ! and with an extra 2 iterations
      do i=1,n
        ynew(i)=y(i)
        yold(i)=y(i)
        ysum(i)=y(i)
      end do

      ! iteration loop:
      do i=1,numit+2
        call Iter( n, ynew, ysum, dt, yp, yl, Q, &
#ifdef with_vbs
                     drog_loss_temp, drog_prod_neg_temp, &
#endif                    
                     nreac, rk, naux, aux )
        nfcn=nfcn+1
      end do

#ifdef with_vbs
      ! y_reacted = yold - ynew
      !
      !                      yold + dt*yp    
      !           = yold -  ------------- = 
      !                      1 + dt*yl      
      !
      !                      yl*yold + yp
      !           =   dt *  --------------
      !                      1 + dt*yl 
      !
      drog = drog + dt * ( drog_loss_temp * yold + drog_prod_neg_temp) / ( 1 + dt * drog_loss_temp)
      drog_old = drog
#endif
      naccpt=naccpt+1
      t=t+dt
      do j=1,n
        y(j)=ynew(j)
      end do

      ! Subsequent steps are carried out with the two-step BDF method.
      dtold=dt
      ratio=1.0

      ! Start of time loop
      goto10 = .false.
      do
        c=1.0/ratio
        cp1=c+1.0
        a1=((c+1.0)**2)/(c*c+2.0*c)
        a2=-1.0/(c*c+2.0*c)
        gdt=dt*(1.0+c)/(2.0+c)
        do j=1,n
          ysum(j)=a1*y(j)+a2*yold(j)
          ynew(j)=max(0.0,y(j)+ratio*(y(j)-yold(j)))
        end do

        ! iteration loop:
        do i=1,numit
          call iter( n, ynew, ysum, gdt, yp, yl, Q, &
#ifdef with_vbs
                       drog_loss_temp, drog_prod_neg_temp, &
#endif                    
                       nreac, rk, naux, aux )
          nfcn=nfcn+1
        end do

        ! If stepsizes should remain equal, stepsize control is omitted.
        if ( dtmin /= dtmax ) then
          ! Otherwise stepsize control is carried out.
          errlte=0.0
          do i=1,n
            ytol=atol(i)+rtol(i)*abs(y(i))
            errlte=max(errlte,abs(c*ynew(i)-cp1*y(i)+yold(i))/ytol)
          end do
          errlte=2.0*errlte/(c+c*c)
          call NEWDT(t,te,dt,dtold,ratio,errlte,accept,dtmin,dtmax)
        end if

        ! Here the step has been accepted. The endpoint check
        ! is carried out.
        if (accept) then
          failer=.false.
          restart=.false.
          t=t+dtold
          naccpt=naccpt+1
          do j=1,n
            yold(j)=y(j)
            y(j)=ynew(j)
          end do
#ifdef with_vbs
          drog_old = drog
          drog     = drog + dt * ( drog_loss_temp * yold + drog_prod_neg_temp) / (1 + dt * drog_loss_temp )
#endif
          ! if time >= end time  -> jump out of time loop
          if ( t >= te ) exit

          ! if accepted, return to begin of time loop
          cycle
        endif

        ! if not accepted, restart from label 10
        if (failer) then
          nrejec=nrejec+1
          failer=.false.
          naccpt=naccpt-1
          t=t-dtold
#ifdef with_vbs
          drog = drog_old
#endif          
          do j=1,n
            y(j)=yold(j)
          end do
          goto10 = .true.
          exit
        endif

        ! Here the step has been rejected;
        ! return to begin of time loop (label 60)
        nrejec=nrejec+1
        failer=.true.
      end do

      if ( .not. goto10 ) exit
    end do

    ! End of TWOSTEP.

  end subroutine TWOSTEP


  ! ***


  pure subroutine NEWDT(t,te,dt,dtold,ratio,errlte, &
                   accept,dtmin,dtmax)

    ! --- in/out -----------------------------

    real, intent(in)        ::  t, te, dtmin, dtmax
    real, intent(inout)     ::  dt, dtold, errlte
    real, intent(out)       ::  ratio
    logical, intent(out)    ::  accept

    ! --- local ------------------------------

    real                    ::  ts

    ! --- begin ------------------------------

    if (errlte.gt.1.0.and.dt.gt.dtmin) then
      accept=.false.
      ts=t
    else
      accept=.true.
      dtold=dt
      ts=t+dtold
    endif
    errlte=max(errlte,1e-10)
    dt=max(0.5,min(2.0,0.8/sqrt(errlte)))*dt
    dt=max(dtmin,min(dt,dtmax))
    call FIT(ts,te,dt)
    ratio=dt/dtold

  end subroutine NEWDT


  ! ***


  pure subroutine FIT(t,te,dt)

    ! --- in/out -----------------------------

    real, intent(in)      ::  t,te
    real, intent(inout)   ::  dt

    ! --- local ------------------------------

    real      ::  rns
    integer   ::  ns

    ! --- begin ------------------------------

    rns=(te-t)/dt
    if (rns.gt.10.0) return
    ns=int(rns)+1
    dt=(te-t)/ns
    dt=(dt+t)-t

  end subroutine FIT


end module LE_Chem_Solver_TwoStep
