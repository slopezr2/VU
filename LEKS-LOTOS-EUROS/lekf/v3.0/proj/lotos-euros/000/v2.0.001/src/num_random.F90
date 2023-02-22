!###############################################################################
!
! numerical tools : random numbers
!
!
! USAGE
!
!    ! 
!    ! Trick to define multiple number generators.
!    ! Every time a call to 'Random_Get_Normal' the seed is reset actually,
!    ! since the state of the Fortran random genererator cannot be saved.
!    !
!
!    ! random number generator:
!    type(T_Random)   ::  rnd
!
!    ! initialize new random generator with seed number:
!    call rnd%Init( 1234, status )
!
!    ! draw normal distributed numbers:
!    call rnd%Get_Normal( x, status )
!    call rnd%Get_Normal( x, status )
!    call rnd%Get_Normal( x, status )
!     :
!
!    ! done:
!    call rnd%Done( status )
!  
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################



module Num_Random

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  gasdev
  
  public  ::  T_Random


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'Num_Random'
  
  
  ! --- types ------------------------------------
  
  type T_Random
    ! number of calls:
    integer                   ::  n
    ! user defined seed:
    integer                   ::  seed0
    ! seed size:
    integer                   ::  nseed
    ! actual seed:
    integer, pointer          ::  seed(:)  ! (nseed)
  contains
    procedure ::  Init          =>  Random_Init
    procedure ::  Done          =>  Random_Done
    procedure ::  Get_Normal    =>  Random_Get_Normal
  end type T_Random
  
  
contains


  ! -------------------------------------------
  ! Returns a normally distributed deviate with zero mean and unit variance
  ! From Numerical Recipes in Fortran, 2nd edition, p.280
  !
  ! Copied from Var4D code.
  ! -------------------------------------------

  real function gasdev()

    integer, save  ::  iset = 0
    real, save     ::  gset
    real           ::  fac,rsq,v1,v2,rnd1,rnd2

    ! generate new values ?
    if ( iset == 0 ) then
      do 
        ! draw two uniform distributed random numbers:
        call random_number(rnd1)  ! [0,1]
        call random_number(rnd2)  ! [0,1]
        ! convert range:
        v1 = 2.0 * rnd1 - 1.0   ! [-1,1]
        v2 = 2.0 * rnd2 - 1.0   ! [-1,1]
        ! squared radius in 2D plane:
        rsq = v1**2 + v2**2   ! [0,2]
        ! ok if in (0,1)
        if ( (rsq > 0.0) .and. (rsq < 1.0) ) exit
      end do
      ! magic ...
      fac = sqrt( -2.0*log(rsq)/rsq )
      ! two Gaussian distributed random numbers
      gset   = v1*fac   ! store for later use
      gasdev = v2*fac   ! return value
      ! reset flag:
      iset = 1
    else
      ! now return the other one:
      gasdev = gset
      ! reset flag:
      iset = 0
    end if
    
  end function gasdev
  
  
  ! ====================================================================
  ! ===
  ! === reproducible random number generator
  ! ===
  ! ====================================================================
  
  !
  ! Inititalize new T_Random type with integer seed.
  !
  
  subroutine Random_Init( self, seed0, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Random), intent(out)    ::  self
    integer, intent(in)             ::  seed0
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Random_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get size:
    call Random_Seed( size=self%nseed )
    ! storage:
    allocate( self%seed(self%nseed), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! no call set:
    self%n = 0
    
    ! user defined seed:
    self%seed0 = seed0
    
    ! ok
    status = 0
    
  end subroutine Random_Init
  
  !
  ! done
  !
  
  subroutine Random_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Random), intent(inout)  ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Random_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear
    deallocate( self%seed, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Random_Done
  
  !
  ! Reset seed given current number of calls
  !
  
  subroutine Random_ResetSeed( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Random), intent(inout)  ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Random_ResetSeed'
    
    ! --- local ----------------------------------
    
    integer        ::  i
    real           ::  x
    
    ! --- begin ----------------------------------
    
    ! increase number of calls:
    self%n = self%n + 1
    
    !! info ...
    !write (gol,'("RANDOM: reset seed to state ",2i6)') self%seed0, self%n; call goPr
    
    ! fill seed values, fill with user-provided number plus value index:
    do i = 1, self%nseed
      self%seed(i) = i + self%seed0
    end do
    ! reset:
    call Random_Seed( put=self%seed )
    
    ! call the random generator for the past calls:
    do i = 1, self%n
      call random_number( x )
    end do
    
    ! increase number of calls:
    self%n = self%n + 1
    
    ! ok
    status = 0
    
  end subroutine Random_ResetSeed
  
  
  !
  ! Draw number from Normal distribution ;
  ! copied from 'gasdev', uses algorithm that draws two uniform numbers,
  ! but in this case the second is not saved but simply ignored
  !
  ! Eventually provide mean and std.dev. for N(mu,sigma)
  !
  
  subroutine Random_Get_Normal( self, x, status, &
                                  mu, sigma )
  
    ! --- in/out ---------------------------------
    
    class(T_Random), intent(inout)  ::  self
    real, intent(out)               ::  x
    integer, intent(out)            ::  status
    
    real, intent(in), optional      ::  mu
    real, intent(in), optional      ::  sigma

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Random_Get_Normal'
    
    ! --- local ----------------------------------

    real           ::  self1, self2
    real           ::  v1, v2
    real           ::  rsq
    real           ::  fac
    
    ! --- begin ----------------------------------
    
    ! reset seed:
    call Random_ResetSeed( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! draw pairs of random numbers in [0,1]x[0,1]
    ! until distance is < 1.0 :
    do 
      ! draw two uniform distributed random numbers:
      call random_number(self1)  ! [0,1]
      call random_number(self2)  ! [0,1]
      ! convert range:
      v1 = 2.0 * self1 - 1.0   ! [-1,1]
      v2 = 2.0 * self2 - 1.0   ! [-1,1]
      ! squared radius in 2D plane:
      rsq = v1**2 + v2**2   ! [0,2]
      ! ok if in (0,1)
      if ( (rsq > 0.0) .and. (rsq < 1.0) ) exit
    end do
    ! magic ...
    fac = sqrt( -2.0*log(rsq)/rsq )
    ! Gaussian distributed random number:
    x = v1*fac
    
    ! std.dev.?
    if ( present(sigma) ) x = x * sigma
    ! mean?
    if ( present(mu) ) x = mu + x
    
    !! info ...
    !write (gol,'("RANDOM: draw normal : ",3f12.4)') self1, self1, x; call goPr
    
    ! ok
    status = 0
    
  end subroutine Random_Get_Normal


end module Num_Random

