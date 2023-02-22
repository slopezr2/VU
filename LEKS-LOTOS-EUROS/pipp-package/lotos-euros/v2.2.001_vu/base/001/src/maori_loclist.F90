!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_LocList

  use GO, only : gol, goPr, goErr

  use MAORI_Param, only : MAORI_LEN_NAME, MAORI_LEN_LONGNAME, MAORI_LEN_LINE

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  T_MAORI_Loc, T_MAORI_LocList
  public  ::  MAORI_LocList_Init, MAORI_LocList_Done
  public  ::  MAORI_LocList_Inq
  public  ::  MAORI_LocList_Find
  public  ::  MAORI_Loc_Inq


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_LocList'


  ! --- types ------------------------------

  ! single location:
  type T_MAORI_Loc
    ! unique location id (record number in list?)
    integer                             ::  loc_id
    ! coordinate:
    real                                ::  lon, lat   ! degree
    ! height:
    real                                ::  alt        ! m (above sea level)
    ! other (name, code, etc)
    character(len=MAORI_LEN_LONGNAME)   ::  name
    character(len=MAORI_LEN_NAME)       ::  code
    character(len=MAORI_LEN_LINE)       ::  meta
  end type T_MAORI_Loc

  ! list
  type T_MAORI_LocList
    ! number of elements:
    integer                         ::  n
    ! list:
    type(T_MAORI_Loc), pointer      ::  loc(:)
  end type T_MAORI_LocList


contains


  ! ==================================================================


  subroutine MAORI_LocList_Init( loclist, key, status, bounds )

    use GO, only : goVarValue
    use GO, only : goGetFU

    use MAORI_LocList_LML     , only : T_MAORI_LocList_LML, Init, Done, ReadRecord
    use MAORI_LocList_EMEP    , only : T_MAORI_LocList_EMEP, Init, Done, ReadRecord
    use MAORI_LocList_EARLINET, only : T_MAORI_LocList_EARLINET, Init, Done, ReadRecord
    use MAORI_LocList_AERONET , only : T_MAORI_LocList_AERONET, Init, Done, ReadRecord
    use MAORI_LocList_AirBase , only : T_MAORI_LocList_AirBase, Init, Done, ReadRecord
    use MAORI_LocList_AQORD   , only : T_MAORI_LocList_AQORD, Init, Done, ReadRecord
    use MAORI_LocList_MACC_EVA, only : T_MAORI_LocList_MACC_EVA, Init, Done, ReadRecord
    use MAORI_LocList_CSV     , only : T_MAORI_LocList_CSV

    use MAORI_Param, only : MAORI_LEN_FILE, MAORI_LEN_NAME

    ! --- in/out ---------------------------------

    type(T_MAORI_LocList), intent(out)    ::  loclist
    character(len=*), intent(in)          ::  key
    integer, intent(out)                  ::  status
    real, intent(in), optional            ::  bounds(4)  ! (/west,east,south,north/) in degrees

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_LocList_Init'

    ! --- local -----------------------------------

    character(len=MAORI_LEN_FILE)   ::  fname
    character(len=MAORI_LEN_NAME)   ::  ftype
    type(T_MAORI_LocList_LML)       ::  lml
    type(T_MAORI_LocList_EMEP)      ::  emep
    type(T_MAORI_LocList_EARLINET)  ::  earlinet
    type(T_MAORI_LocList_AERONET)   ::  aeronet
    type(T_MAORI_LocList_AirBase)   ::  airbase
    type(T_MAORI_LocList_AQORD)     ::  aqord
    type(T_MAORI_LocList_MACC_EVA)  ::  macc_eva
    type(T_MAORI_LocList_CSV)       ::  csv
    integer                         ::  i
    integer                         ::  nrec
    integer                         ::  iloc
    real                            ::  lon, lat
    real                            ::  alt
    character(len=MAORI_LEN_LINE)   ::  name
    character(len=MAORI_LEN_LINE)   ::  code
    character(len=MAORI_LEN_LINE)   ::  meta

    ! --- begin -----------------------------------

    ! key format:
    !    file=<fname>;type=<ftype>
    fname = 'no-loc-list-file.txt'
      call goVarValue( key, ';', 'file', '=', fname, status )
      IF_ERROR_RETURN(status=1)
    ftype = 'none'
      call goVarValue( key, ';', 'type', '=', ftype, status )
      IF_ERROR_RETURN(status=1)

    ! use file type depended routines:
    select case ( ftype )
      ! LML station list:
      case ( 'lml' )
        ! open file, count records:
        call Init( lml, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! EMEP station list:
      case ( 'emep' )
        ! open file, count records:
        call Init( emep, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! EARLINET station list:
      case ( 'earlinet' )
        ! open file, count records:
        call Init( earlinet, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! AERONET station list:
      case ( 'aeronet' )
        ! open file, count records:
        call Init( aeronet, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! AirBase station list:
      case ( 'airbase' )
        ! open file, count records:
        call Init( airbase, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! aqord station list:
      case ( 'aqord' )
        ! open file, count records:
        call Init( aqord, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! MACC EVA station list:
      case ( 'macc_eva' )
        ! open file, count records:
        call Init( macc_eva, trim(fname), nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! generic csv files, specification of header names in keys:
      case ( 'csv' )
        ! open file, count records:
        call csv%Init( trim(fname), key, nrec, status )
        IF_NOTOK_RETURN(status=1)
      ! unknown ...
      case default
        write (gol,'("unsupported loclist file type : ",a)') trim(ftype); call goErr
        TRACEBACK; status=1; return
    end select

    ! check ...
    if ( nrec <= 0 ) then
      write (gol,'("no records found in loclist described by : ",a)') trim(key); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( loclist%loc(nrec), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over records:
    iloc = 0
    do i = 1, nrec

      ! use file type depended routines:
      select case ( ftype )
        ! LML station list:
        case ( 'lml' )
          ! fill record:
          call ReadRecord( lml, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! EMEP station list:
        case ( 'emep' )
          ! fill record:
          call ReadRecord( emep, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! EARLINET station list:
        case ( 'earlinet' )
          ! fill record:
          call ReadRecord( earlinet, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! AERONET station list:
        case ( 'aeronet' )
          ! fill record:
          call ReadRecord( aeronet, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! AirBase station list:
        case ( 'airbase' )
          ! fill record:
          call ReadRecord( airbase, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! AQORD station list:
        case ( 'aqord' )
          ! fill record:
          call ReadRecord( aqord, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! MACC/EVA station list:
        case ( 'macc_eva' )
          ! fill record:
          call ReadRecord( macc_eva, name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! CSV station list:
        case ( 'csv' )
          ! fill record:
          call csv%ReadRecord( name, code, lon, lat, alt, meta, status )
          IF_NOTOK_RETURN(status=1)
        ! unknown ...
        case default
          write (gol,'("unsupported loclist file type : ",a)') trim(ftype); call goErr
          TRACEBACK; status=1; return
      end select

      ! bounding box provided ?
      if ( present(bounds) ) then
        ! outside box ?
        if ( (lon < bounds(1)) .or. (lon > bounds(2)) .or. &
             (lat < bounds(3)) .or. (lat > bounds(4)) ) then
          ! next:
          cycle
        end if
      end if

      ! next record:
      iloc = iloc + 1
      ! store:
      loclist%loc(iloc)%loc_id = i
      loclist%loc(iloc)%name   = trim(name)
      loclist%loc(iloc)%code   = trim(code)
      loclist%loc(iloc)%lon    = lon
      loclist%loc(iloc)%lat    = lat
      loclist%loc(iloc)%alt    = alt
      loclist%loc(iloc)%meta   = trim(meta)
      
    end do ! locations
    
    ! accepted number of locations:
    loclist%n = iloc

    ! use file type depended routines:
    select case ( ftype )
      ! LML station list:
      case ( 'lml' )
        ! close file:
        call Done( lml, status )
        IF_NOTOK_RETURN(status=1)
      ! EMEP station list:
      case ( 'emep' )
        ! close file:
        call Done( emep, status )
        IF_NOTOK_RETURN(status=1)
      ! EARLINET station list:
      case ( 'earlinet' )
        ! close file:
        call Done( earlinet, status )
        IF_NOTOK_RETURN(status=1)
      ! AERONET station list:
      case ( 'aeronet' )
        ! close file:
        call Done( aeronet, status )
        IF_NOTOK_RETURN(status=1)
      ! AirBase station list:
      case ( 'airbase' )
        ! close file:
        call Done( airbase, status )
        IF_NOTOK_RETURN(status=1)
      ! aqord station list:
      case ( 'aqord' )
        ! close file:
        call Done( aqord, status )
        IF_NOTOK_RETURN(status=1)
      ! MACC EVA station list:
      case ( 'macc_eva' )
        ! close file:
        call Done( macc_eva, status )
        IF_NOTOK_RETURN(status=1)
      ! CSV station list:
      case ( 'csv' )
        ! close file:
        call csv%Done( status )
        IF_NOTOK_RETURN(status=1)
      ! unknown ...
      case default
        write (gol,'("unsupported loclist file type : ",a)') trim(ftype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! info ...
    write (gol,'("  MAORI selected stations:")'); call goPr
    write (gol,'("    ---------------------------------------------------------------------")'); call goPr
    write (gol,'("      nr      location        altitude [code] name                       ")'); call goPr
    write (gol,'("    ---- -------------------- -------- ----------------------------------")'); call goPr
    do i = 1, loclist%n
      write (gol,'("    ",i4," (",f8.2,",",f8.2,") ",f8.2," [",a,"] ",a)') &
              i, loclist%loc(i)%lon , loclist%loc(i)%lat, loclist%loc(i)%alt, &
              trim(loclist%loc(i)%code), trim(loclist%loc(i)%name); call goPr
    end do
    write (gol,'("    ---------------------------------------------------------------------")'); call goPr

    ! ok
    status = 0

  end subroutine MAORI_LocList_Init


  ! ***


  subroutine MAORI_LocList_Done( loclist, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_LocList), intent(inout)    ::  loclist
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_LocList_Done'

    ! --- begin -----------------------------------

    ! clear:
    if ( loclist%n > 0 ) then
      ! clear:
      deallocate( loclist%loc, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! reset:
      loclist%n = 0
    end if

    ! ok
    status = 0

  end subroutine MAORI_LocList_Done



  ! ***


  subroutine MAORI_LocList_Inq( loclist, status, nloc, loc_id_range )

    ! --- in/out ---------------------------------

    type(T_MAORI_LocList), intent(in)       ::  loclist
    integer, intent(out)                    ::  status

    integer, intent(out), optional          ::  nloc
    integer, intent(out), optional          ::  loc_id_range(2)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_LocList_Inq'

    ! --- local ----------------------------------
    
    integer     ::  i

    ! --- begin ----------------------------------

    ! return values:
    if ( present(nloc) ) nloc = loclist%n
    
    ! min/max range of location ids?
    if ( present(loc_id_range) ) then
      ! any obs?
      if ( loclist%n > 0 ) then
        ! loop:
        do i = 1, loclist%n
          ! first?
          if ( i == 1 ) then
            ! fill from first value:
            loc_id_range = (/ loclist%loc(i)%loc_id, loclist%loc(i)%loc_id /)
          else
            ! update:
            loc_id_range(1) = min( loc_id_range(1), loclist%loc(i)%loc_id )
            loc_id_range(2) = max( loc_id_range(2), loclist%loc(i)%loc_id )
          end if
        end do ! i
      else
        ! dummy:
        loc_id_range = (/9999,-9999/)
      end if
    end if

    ! ok
    status = 0

  end subroutine MAORI_LocList_Inq



  ! ***
  
  
  !
  ! Search location matching for optional argument:
  !   code   : station code
  ! Return status:
  !  -1 : not found
  !   0 : found
  !  >0 : error
  !

  subroutine MAORI_LocList_Find( loclist, iloc, status, &
                                    code, quiet )

    ! --- in/out ---------------------------------

    type(T_MAORI_LocList), intent(in)         ::  loclist
    integer, intent(out)                      ::  iloc
    integer, intent(out)                      ::  status

    character(len=*), intent(in), optional    ::  code
    logical, intent(in), optional             ::  quiet

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_LocList_Find'

    ! --- local ----------------------------------
    
    logical         ::  found
    logical         ::  shutup

    ! --- begin ----------------------------------
    
    ! shut up?
    shutup = .false.
    if ( present(quiet) ) shutup = quiet
    
    ! loop:
    do iloc = 1, loclist%n
      ! init flag:
      found = .true.
      ! match?
      if ( present(code) ) found = found .and. (trim(loclist%loc(iloc)%code) == trim(code))
      ! still correct?
      if ( found ) exit
    end do
    ! check ...
    if ( .not. found ) then
      if ( .not. shutup ) then
        write (gol,'("could not find any location matching:")'); call goErr
        if ( present(code) ) then
          write (gol,'("  code == ",a)') trim(code); call goErr
        end if
        TRACEBACK
      end if
      ! return with warning:
      status = -1; return
    end if

    ! ok
    status = 0

  end subroutine MAORI_LocList_Find


  ! ***


  subroutine MAORI_Loc_Inq( loclist, iloc, status, &
                                 loc_id, name, code, lon, lat, alt )

    ! --- in/out ---------------------------------

    type(T_MAORI_LocList), intent(in)           ::  loclist
    integer, intent(in)                         ::  iloc
    integer, intent(out)                        ::  status

    integer, intent(out), optional              ::  loc_id
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  code
    real, intent(out), optional                 ::  lon
    real, intent(out), optional                 ::  lat
    real, intent(out), optional                 ::  alt

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Loc_Inq'

    ! --- begin ----------------------------------

    ! check:
    if ( (iloc < 1) .or. (iloc > loclist%n) ) then
      write (gol,'("loc index ",i6," outside valid range 1, ..",i6)') iloc, loclist%n; call goErr
      TRACEBACK; status=1; return
    end if

    ! return values:
    if ( present(loc_id) ) loc_id = loclist%loc(iloc)%loc_id
    if ( present(name  ) ) name   = loclist%loc(iloc)%name
    if ( present(code  ) ) code   = loclist%loc(iloc)%code
    if ( present(lon   ) ) lon    = loclist%loc(iloc)%lon
    if ( present(lat   ) ) lat    = loclist%loc(iloc)%lat
    if ( present(alt   ) ) alt    = loclist%loc(iloc)%alt

    ! ok
    status = 0

  end subroutine MAORI_Loc_Inq


end module MAORI_LocList
