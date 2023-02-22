!#######################################################################
!
! LE_Emis_Point 
!   routines to add emissions at single point ;
!   location, component and release rate defined in rcfile
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Point

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_Point


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Point'


  ! --- local --------------------------------

  ! emission data base
  type T_Emis_Point
    ! label assigned to this emission:
    character(len=32)                 ::  label
    ! grid cell:
    integer                           ::  i, j
    ! tracer:
    integer                           ::  ispec
    ! release rate:
    real                              ::  rate
    real                              ::  stackheight
    character(len=32)                 ::  units
    character(len=32)                 ::  heightunits
    ! time range:
    type(TDate)                       ::  t1, t2
  contains
    procedure     ::  Init    =>  LE_Emis_Point_Init
    procedure     ::  Done    =>  LE_Emis_Point_Done
    procedure     ::  Setup   =>  LE_Emis_Point_Setup
  end type T_Emis_Point

contains


  ! ===============================================================


  subroutine LE_Emis_Point_Init( self, label, rcF, rckey, status )

    use GO     , only : TrcFile, ReadRc
    use GO     , only : goMatchValue, goReadFromLine
    use GO     , only : operator(<), wrtgol
    use Grid   , only : In_Domain, GetLocation
    use LE_Grid, only : ugg
    use indices, only : specname
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    class(T_Emis_Point), intent(out)    ::  self
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Point_Init'

    ! --- local ------------------------------
    
    character(len=32)           ::  name
    real                        ::  lon, lat
    logical                     ::  indomain
    character(len=64)           ::  line

    ! --- begin -------------------------------

    ! store label:
    self%label = trim(label)
             
    ! read location:
    call ReadRc( rcF, trim(rckey)//'.lon', lon, status )
    IF_NOTOK_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.lat', lat, status )
    IF_NOTOK_RETURN(status=1)

    ! check if it is the local domain:
    call ugg%InDomain( lon, lat, indomain, status )
    IF_NOTOK_RETURN(status=1)
    ! in domain?
    if ( indomain ) then

      ! grid cell indices:
      call ugg%GetLocation( lon, lat, self%i, self%j, status )
      IF_NOTOK_RETURN(status=1)

      ! read name of tracer:
      call ReadRc( rcF, trim(rckey)//'.tracer', name, status )
      IF_NOTOK_RETURN(status=1)    
      ! match tracer name to global tracer index
      call goMatchValue( trim(name), specname, self%ispec, status )
      IF_NOTOK_RETURN(status=1)

      ! emission rate and units:
      call ReadRc( rcF, trim(rckey)//'.rate', line, status )
      IF_NOTOK_RETURN(status=1)
      ! value:
      call goReadFromLine( line, self%rate, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! remainder is units:
      self%units = trim(line)

      ! emission height:
      call ReadRc( rcF, trim(rckey)//'.height', line, status )
      IF_NOTOK_RETURN(status=1)
      ! value:
      call goReadFromLine( line, self%stackheight, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! remainder is units:
      self%heightunits = trim(line)

      ! time range:
      !   1900-01-01 00:00:00 ; 2099-01-01 00:00:00
      call ReadRc( rcF, trim(rckey)//'.trange', line, status )
      IF_NOTOK_RETURN(status=1)
      ! value:
      call goReadFromLine( line, self%t1, status, sep=';' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( line, self%t2, status, sep=';' )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( self%t2 < self%t1 ) then
        call wrtgol( 'invalid time range : ', self%t1, ' - ', self%t2 ); call goErr
        TRACEBACK; status=1; return
      end if

      ! enable data:    
      call LE_Data_Enable( 'area', status )
      IF_NOTOK_RETURN(status=1)
      
    else
    
      ! dummy ...
      self%i = -999
      self%j = -999
      
    end if
    
    ! ok
    status = 0

  end subroutine LE_Emis_Point_Init


  ! ***


  subroutine LE_Emis_Point_Done( self, status )
    
    ! --- in/out ------------------------------

    class(T_Emis_Point), intent(inout)    ::  self
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Point_Done'

    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! ok
    status = 0

  end subroutine LE_Emis_Point_Done


  ! ***


  subroutine LE_Emis_Point_Setup( self, emis_a, t1, t2, status )

    use GO     , only : TDate, rTotal, operator(-), operator(<), operator(<=), wrtgol
    use Dims   , only : nx, ny, nz, nspec
    use Indices, only : specunit, specname, specmolm
    use LE_Data, only : LE_Data_GetPointer
 
    ! --- in/out ---------------------------

    class(T_Emis_Point), intent(inout)    ::  self
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Point_Setup'
    
    real, parameter   ::  hour_per_min = 1.0/60.0  ! hour/min
    real, parameter   ::  sec_per_min  = 60.0      ! sec/min
    real, parameter   ::  ug_per_kg    = 1.0e9     ! ug/kg
    
    ! --- local ----------------------------
    
    real           ::  dt_sec
    real           ::  demis
    integer        ::  iz
    real, pointer  ::  area(:,:,:)   ! (lon,lat,1)    
    real, pointer  ::  h_m(:,:,:)   ! (lon,lat,lev)

    ! --- begin ----------------------------
    
    ! in domain?
    if ( (self%i > 0) .and. (self%j > 0) ) then

      ! pointers to meteo fields:
      call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
      IF_NOTOK_RETURN(status=1)    
      call LE_Data_GetPointer( 'h', h_m, status, check_units=self%heightunits )
      IF_NOTOK_RETURN(status=1)

      ! check for time range:
      if ( (t2 <= self%t1) .or. (self%t2 <= t1) ) then

        ! outside range, no emission

      else if ( (self%t1 <= t1) .and. (t2 <= self%t2) ) then

        ! time step in seconds:
        dt_sec = rTotal( t2 - t1, 'sec' )

        ! switch on target units:
        select case ( specunit(self%ispec) )
          !
          !~ volume mixing ratio's:
          case ( 'ppb' )
            ! unit emis_a for gases   : mol/min
            select case ( trim(self%units) )
              !
              case ( 'mol/min' )
                ! mol/min =  mol/min 
                demis     =  self%rate
              !
              case ( 'mol/s' )
                ! mol/min =   mol/s        s/min
                demis     = self%rate * sec_per_min
              !
              case ( 'kg/hour' )
                ! mol/min =  kg/hour /     (kg/mol)            hour/min
                demis     = self%rate / specmolm(self%ispec) * hour_per_min
              !
              case default
                write (gol,'("could not convert point emis units `",a,"` to `",a,"`")') &
                         trim(self%units), 'mol/min'; call goErr
                TRACEBACK; status=1; return
            end select
          !
          !~ mass concentrations:
          case ( 'ug/m3' )
            !
            ! unit emis_a for gases   : ug/min
            select case ( trim(self%units) )
              !
              case ( 'ug/min' )
                ! ug/min =  ug/min
                demis    = self%rate
              !
              case ( 'kg/hour' )
                ! ug/min =  kg/hour    ug/kg        hour/min
                demis    = self%rate * ug_per_kg * hour_per_min
              !
              case ( 'ug/m2/s' )
                ! mol/min =  ug/m2/s          m2               s/min
                demis     = self%rate * area(self%i,self%j,1) * sec_per_min
              !
              case default
                write (gol,'("could not convert point emis units `",a,"` to `",a,"`")') &
                         trim(self%units), 'ug/min'; call goErr
                TRACEBACK; status=1; return
            end select
          !
          !~ unknown ..
          case default
            write (gol,'("unsupported unit `",a,"` for tracer ",i4," (",a,")")') &
                    trim(specunit(self%ispec)), self%ispec, trim(specname(self%ispec)); call goErr
            TRACEBACK; status=1; return
        end select

        ! determine in which layer mass is injected
        do iz=2,nz
          if (self%stackheight >= h_m(self%i,self%j,iz-1) .and. self%stackheight < h_m(self%i,self%j,iz) ) then
            ! add contribution:
            emis_a(self%i,self%j,iz,self%ispec) = emis_a(self%i,self%j,iz,self%ispec) + demis 
           end if
        enddo
        if (self%stackheight > h_m(self%i,self%j,nz)) then
          write( gol, '("Test point source emission above model top")' ) ; call goErr
          write( gol, '("Model top: ", f6.3)' ) h_m(self%i,self%j,nz) ; call goErr
          write( gol, '("Stackheight: ", f6.3)' ) self%stackheight ; call goErr
          TRACEBACK;status=1;return             
        else if  (self%stackheight < h_m(self%i,self%j,1) ) then
          ! add contribution to lowest layer:
          emis_a(self%i,self%j,1,self%ispec) = emis_a(self%i,self%j,1,self%ispec) + demis
        end if

      else

        ! not yet ...
        write (gol,'("partly overlapping time ranges not supported yet:")'); call goErr
        call wrtgol( '  emission time range : ', self%t1, ' - ', self%t2 ); call goErr
        call wrtgol( '  model time range    : ', t1, ' - ', t2 ); call goErr
        TRACEBACK; status=1; return

      end if  ! check time range
      
    end if  ! in local domain

    ! ok
    status = 0

  end subroutine LE_Emis_Point_Setup
  
  
  ! ***
 
end module LE_Emis_Point


 
