!#################################################################
!
! NAME
!
!   LE_Output_Emis_Label  -  LOTOS-EUROS output of labelled emissions fields
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
!
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!
!#################################################################

module LE_Output_Emis_Label

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Emis_Label

  public  ::  LE_Output_Emis_Label_Init, LE_Output_Emis_Label_Done
  public  ::  LE_Output_Emis_Label_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Emis_Label'
  
  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10

  ! --- types ------------------------------

  type T_LE_Output_Emis_label
    ! name for this file:
    character(len=16)                     ::  typ
    character(len=16)                     ::  name
    ! common stuff:
    type(T_LE_Output_Common)              ::  com
    ! adhoc convention, might change per file ...
    character(len=32)                     ::  convention
    ! replace existing files ?
    logical                               ::  replace
    ! file opened ?
    logical                               ::  opened
    ! current time range:
    type(TDate)                           ::  tr(2)
    ! time resolution:
    real                                  ::  dhour
    ! collect: daily, instant
    character(len=32)                     ::  collect
    ! time record counter
    integer                               ::  itrec
    ! state name:
    character(len=16)                     ::  state
    ! file name:
    character(len=512)                    ::  fname
    ! file handle:
    integer                               ::  ncid
    ! dimension handles:
    integer                               ::  dimid_lon      
    integer                               ::  dimid_lat      
    integer                               ::  dimid_lev      
    integer                               ::  dimid_time     
    integer                               ::  dimid_label    
    integer                               ::  dimid_labelname_len 
    ! dimension variables:
    integer                               ::  varid_lon      
    integer                               ::  varid_lat      
    integer                               ::  varid_lev      
    integer                               ::  varid_time     
    integer                               ::  varid_time_dtg     
    integer                               ::  varid_labelnames    
    ! tracer variables:
    integer                               ::  ntr         
    integer, pointer                      ::  itr_lab(:)      
    integer, pointer                      ::  itr_glob(:)      
    character(len=32), pointer            ::  name_tr(:)  
    character(len=32), pointer            ::  unit_tr(:)
    real, pointer                         ::  unitconv(:)
    integer, pointer                      ::  varid_tr(:) 
    ! level selection:
    character(len=16)                     ::  levtype
    integer                               ::  nlev        
    integer, pointer                      ::  ilev(:)    
    real                                  ::  heights(maxlev)
    ! grads ctl file ?
    logical                               ::  grads_ctl
    character(len=256)                    ::  grads_ctl_file
    character(len=256)                    ::  grads_ctl_base
    integer                               ::  grads_ctl_nt
    type(TDate)                           ::  grads_ctl_t1
    type(TIncrDate)                       ::  grads_ctl_dt
    ! bounding box
    integer                               ::  i1, i2, ni
    integer                               ::  j1, j2, nj 
    real                                  ::  westb, southb
    
  end type T_LE_Output_Emis_label


contains


  ! ====================================================


  subroutine LE_Output_Emis_Label_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues, goSplitString
    use GO     , only : AnyDate
    use Dims   , only : nx, ny, nz
    use Dims   , only : nspec_all
    use Indices, only : specname, specunit
    use LE_Grid, only : ugg
    use SA_Labeling, only : nspec_labelled, labelled_specs_names

    ! --- in/out --------------------------------

    type(T_LE_Output_Emis_Label), intent(out)   ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Emis_Label_Init'

    ! --- local ---------------------------------
    
    character(len=64)               ::  basekey
    character(len=512)              ::  output_tracer_names
    character(len=32)               ::  output_level_numbers
    
    character(len=512)              ::  key
    real                            ::  west, east, south, north
    real, pointer                   ::  ff(:,:)

    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name
    
    ! files not open yet:
    leo%opened = .false.
    
    ! init common stuff:
    call leo%com%Init( rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! no time range set yet:
    leo%tr(1) = AnyDate()
    leo%tr(2) = AnyDate()    
    
    ! Read in output properties from rc file
    call rcF%Get( trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)
    
    ! write GrADS ctl file ?
    call rcF%Get( trim(rckey)//'.ctl', leo%grads_ctl, status )
    IF_NOTOK_RETURN(status=1)
    
    ! base key:
    write (basekey,'(a,".",a,".",a)') trim(rckey), trim(typ), trim(name)
    
    ! output convention:
    call rcF%Get( trim(basekey)//'.convention', leo%convention, status )
    IF_NOTOK_RETURN(status=1)

    ! collect daily or instant
    call rcF%Get( trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call rcF%Get( trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! level type:
    call rcF%Get( trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer names :
    call rcF%Get( trim(basekey)//'.specs', output_tracer_names, status )
    IF_NOTOK_RETURN(status=1)

    ! leveldescriptions
    call rcF%Get( trim(basekey)//'.levels', output_level_numbers, status )
    IF_NOTOK_RETURN(status=1)       

    ! state name:
    leo%state = trim(state)
        
    ! Allocate Type variables for output
    allocate( leo%itr_lab (nspec_all), stat = status ) 
    IF_NOTOK_RETURN(status=1)    
    allocate( leo%itr_glob(nspec_all), stat = status ) 
    IF_NOTOK_RETURN(status=1)    
    allocate( leo%name_tr (nspec_all), stat = status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%unit_tr (nspec_all), stat = status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%unitconv (nspec_all), stat = status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%varid_tr(nspec_all), stat = status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%ilev(1+nz+1), stat = status )    !(surface + nz-levels + upperboundary)
    IF_NOTOK_RETURN(status=1)
    
    ! Match labeled emitted tracer names with selected output
    call goMatchValues( output_tracer_names, labelled_specs_names(1:nspec_labelled), &
                        leo%ntr, leo%name_tr, leo%itr_lab, &
                        status )
    IF_NOTOK_RETURN(status=1)
        
    ! Match labeled emitted tracer names with global tracer names
    call goMatchValues( output_tracer_names, specname, &
                        leo%ntr, leo%name_tr, leo%itr_glob, &
                        status )
    IF_NOTOK_RETURN(status=1)
        
    ! Match levelnumbers with selected output
    call goMatchValues( output_level_numbers, 0, nz+1, &
                        leo%nlev, leo%ilev, &
                        status )
    IF_NOTOK_RETURN(status=1) 
    
    ! init GrADS stuff:
    if ( leo%grads_ctl ) then
      ! no times written yet:
      leo%grads_ctl_nt = 0
    end if
    
    ! bounding box
    call rcF%Get( trim(basekey)//'.bounding_box', key, status )
    IF_NOTOK_RETURN(status=1)   
    ! empty?     
    if (len_trim(key) == 0) then
      ! full domain
      leo%i1 = 1
      leo%i2 = nx
      leo%ni = nx
      leo%j1 = 1
      leo%j2 = ny
      leo%nj = ny
      leo%westb = ugg%longitude_bnds(1,1,1)
      leo%southb = ugg%latitude_bnds(1,1,1)
    else
      ! not yet ...
      write (gol,'("no output subset supported for domain decomposition yet")'); call goErr
      TRACEBACK; status=1; return
      
      select case ( trim(ugg%type) ) 
        
        case ( 'cartesian-regular') 
          ! read domain from key
          read(key,*,iostat=status) west, east, south, north
          if(status/=0) then
            write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
          endif

          ! Check if bounding box is in run domain
          if ( west < ugg%longitude_bnds_1d(1,1) .or. east > ugg%longitude_bnds_1d(2,ugg%nlon) .or. &
               south < ugg%latitude_bnds_1d(1,1) .or. north > ugg%latitude_bnds_1d(2,ugg%nlat) ) then
            write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
            write( gol, '("Run domain: ", 4f8.2)' ) ugg%longitude_bnds_1d(1,1),ugg%longitude_bnds_1d(2,ugg%nlon),ugg%latitude_bnds_1d(1,1),ugg%latitude_bnds_1d(2,ugg%nlat); call goErr
            write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
            TRACEBACK;status=1;return
          endif

          ! for safety
          nullify(ff)
          ! get cell range covered by box
          call ugg%GetDistribution(west,east,south,north,leo%i1,leo%i2,leo%j1,leo%j2,ff,status)
          IF_NOTOK_RETURN(status=1)
          !clear, fractions not used
          if ( associated(ff) ) deallocate(ff)
          ! set shape
          leo%ni = leo%i2-leo%i1+1
          leo%nj = leo%j2-leo%j1+1
          ! set west/south bounds
          leo%westb  = ugg%longitude_bnds_1d(1,leo%i1)
          leo%southb = ugg%latitude_bnds_1d(1,leo%j1)
        case default 
          write( gol, '("Definition of bounding box not clear for grid-type: ", a)' ) trim(ugg%type) ; call goErr
          TRACEBACK;status=1;return
      end select
    end if

    ! ok
    status = 0

  end subroutine LE_Output_Emis_Label_Init


  ! ***


  subroutine LE_Output_Emis_Label_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Emis_Label), intent(inout)   ::  leo
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Label_Done'

    character(len=256) :: commandline
    ! --- begin ---------------------------------

    ! file opened ?
    if ( leo%opened ) then
      ! close:
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NF90_NOTOK_RETURN(status=1)
    end if 
    
    ! deallocate type output variables
    deallocate( leo%itr_lab  )
    deallocate( leo%itr_glob )
    deallocate( leo%name_tr  )
    deallocate( leo%unit_tr  )
    deallocate( leo%unitconv )
    deallocate( leo%varid_tr )
    deallocate( leo%ilev     )

    ! done with common stuff ...
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Output_Emis_Label_Done


  ! ***


  subroutine LE_Output_Emis_Label_PutOut( leo, t, c, cg, status )
    
    use NetCDF , only : NF90_Open
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att 

    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF , only : NF90_NETCDF4
    use NetCDF , only : NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR
    use NetCDF , only : NF90_Write, NF90_NoWrite
    
    use LE_CF_Conventions   , only : LE_CF_names

    use LE_Output_Tools , only : LE_Output_Define_Dims_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO, only : TIncrDate, NewDate, IncrDate, Get, MidNight
    use LE_Comm    , only : goc_dom
    
    use Dims       , only : nx, ny, nz
    
    use LE_Grid    , only : ugg, glb_ugg
    use C3PO       , only : T_Grid_NcDef
    use Indices    , only : nspec
    use Indices    , only : specname, specunit, specmolm
    
    use SA_Labeling, only : SA_nlabel
    use SA_Labeling, only : SA_Label_Names, labelname_len
    use SA_Labeling, only : SA_Emis

    ! --- in/out --------------------------------

    type(T_LE_Output_Emis_Label), intent(inout)   ::  leo
    type(TDate), intent(in)                 ::  t
    real, intent(in)                        ::  c (nx,ny,nz,nspec)
    real, intent(in)                        ::  cg(nx,ny,nspec)
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Emis_Label_PutOut'

    ! --- local ---------------------------------
    
    integer                                 ::  time6(6)
    integer                                 ::  time
    integer                                 ::  cmode
    real, allocatable                       ::  pat(:,:,:)
    real, allocatable                       ::  field(:,:,:,:)  ! nx,ny,0:nz+1, SA_nlabel
    
    integer                                 ::  ix,iy,iz,ispec
    integer                                 ::  i, j
    integer                                 ::  k, l
    integer                                 ::  itr, ilev
    integer                                 ::  varid
    type(T_Grid_NcDef)                      ::  gncd
    type(TDate)                             ::  t0
    real                                    ::  inv_denum
    character(len=32)                       ::  varname
    character(len=256)                      ::  cf_standard_name, cf_long_name, cf_units, cf_tracer_name
    character(len=512)                      ::  comment
    integer                                 ::  SA_ilabel
    integer                                 ::  labellen

    character(len=256) :: commandline
    ! --- begin ---------------------------------
        
    ! This time in current NC-file?
    if( (t < leo%tr(1)) .or. t > (leo%tr(2) ) ) then
      
      ! file opened ?
      if (leo%opened ) then
        ! close
        status = NF90_Close( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! reset flag:
        leo%opened = .false.
      endif
      
      ! extract time fields:
      call Get( t, time6=time6 )

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily' )
          ! set time range [00,24) for this day:
          leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
          leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:3)
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! collect daily for (00,24]
        case ( 'daily24' )
          ! set time range (00,24] for this day:
          if ( Midnight(t) ) then
            leo%tr(1) = t - IncrDate(day=1)
            leo%tr(2) = t
          else
            leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
            leo%tr(2) = leo%tr(1) + IncrDate( day=1 )
            leo%tr(1) = leo%tr(1) + IncrDate(mili=1)
          end if
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), leo%tr(1)%year, leo%tr(1)%month, leo%tr(1)%day
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! files with instant fields:
        case ( 'instant' )
          ! set time range for current instant time:
          leo%tr(1) = t
          leo%tr(2) = t
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:5)
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

      ! root only:
      if ( goc_dom%root ) then
        
        ! set creation mode flag:
        if (leo%replace) then
          cmode = NF90_CLOBBER    ! overwrite existing file
        else
          cmode = NF90_NOCLOBBER  ! do not overwrite existing file
        endif
        
        ! enable large file support:
        cmode = or( cmode, NF90_NETCDF4 )

        ! create file
        status = NF90_Create(trim(leo%fname), cmode, leo%ncid)
        if ( status /= NF90_NOERR ) then
          write (gol,'("creating file : ")'); call goErr
          write (gol,'("  file name  : ",a)') trim(leo%fname); call goErr
          write (gol,'("  nf90 error : ",a)') trim(nf90_strerror(status)); call goErr
          TRACEBACK; status=1; return
        end if 

        ! reset flag:
        leo%opened = .true.

        ! write global attributes:
        call leo%com%PutOut_GlobalAttributes( leo%ncid, status )
        IF_NOTOK_RETURN(status=1)
      
        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                 dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
        IF_NOTOK_RETURN(status=1)
      
        ! level dimension
        call LE_Output_Define_Dims_Lev(leo%ncid, leo%dimid_lev, leo%nlev, trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time dimensions
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        ! label dimension
        status = NF90_Def_Dim( leo%ncid, 'label', SA_nlabel, leo%dimid_label )
        IF_NF90_NOTOK_RETURN(status=1)

        ! len_labelname dimension
        status = NF90_Def_Dim( leo%ncid, 'labelname_len', labelname_len, leo%dimid_labelname_len )
        IF_NF90_NOTOK_RETURN(status=1)

        !
        ! define variables:
        !

        ! level variables
        call LE_Output_Define_Vars_Lev(leo%ncid, leo%varid_lev, leo%dimid_lev, &
                                         trim(leo%levtype), trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)
      
        ! labelname variable
        status = NF90_Def_Var( leo%ncid, 'labelnames', NF90_CHAR, (/leo%dimid_labelname_len, leo%dimid_label/), varid )
        IF_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'Names of the label defnitions')
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_labelnames = varid
      end if ! root

      ! loop over tracers to be written:
      do l = 1, leo%ntr

        ! global tracer index
        itr = leo%itr_glob(l)

        ! no comment yet
        comment = ''

        ! set names etc following some convention:
        select case ( trim(leo%convention) )
          
          case ( 'cf' )

            ! standard variable name:
            varname = trim(leo%name_tr(l))

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status, cf_tracer_name=cf_tracer_name )
            IF_NOTOK_RETURN(status=1)
            
            ! reset to new unit ...
            cf_standard_name = 'emission_of_'//trim(cf_tracer_name)//'_in_air'
            cf_long_name     = trim(cf_standard_name)
            cf_units         = 'kg m-2 s-1'
            leo%unitconv(l)  = -999.9    ! this enables conversion further on ...
          
          case default
          
            write (gol,'("unsupported output convention : ",a)') trim(leo%convention); call goErr
            TRACEBACK; status=1; return
            
        end select
        
        ! store units for later usage (GrADS ctl file, conversions):
        leo%unit_tr(l) = trim(cf_units)

        if (goc_dom%root) then
          ! define variable
          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                              (/leo%dimid_lon, leo%dimid_lat, leo%dimid_lev &
                              , leo%dimid_label, leo%dimid_time/), varid, &
                                 deflate_level=leo%com%deflate_level )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'long_name', trim(cf_long_name)  )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'units', trim(cf_units)  )
          ! molemass:
          call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'molemass', specmolm(itr) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'molemass_unit', 'kg mole-1' )
          IF_NF90_NOTOK_RETURN(status=1)
          ! optional comment:
          if ( len_trim(comment) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
          leo%varid_tr(l) = varid

        end if  ! root

      end do  ! written tracers
      
      ! root?
      if (goc_dom%root ) then
        ! end definition mode:
        status = NF90_EndDef(leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if ! root

      ! no records written yet:
      leo%itrec = 0

    end if

    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root ?
    if ( goc_dom%root ) then

      ! GrADS time counter:
      if ( leo%grads_ctl ) then
        ! increase counter:
        leo%grads_ctl_nt = leo%grads_ctl_nt + 1
        ! set times if necessary:
        if ( leo%grads_ctl_nt == 1 ) then
          leo%grads_ctl_t1 = t
          leo%grads_ctl_dt = IncrDate(day=1)   ! dummy ...
        end if
        if ( leo%grads_ctl_nt == 2 ) then
          leo%grads_ctl_dt = t - leo%grads_ctl_t1
        end if
      end if

      ! write dimension data only once ...
      if ( leo%itrec == 1 ) then

        ! write grid to netCDF file
        call glb_ugg%PutGrid_NetCDF( gncd, status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
        ! write level indices
        call LE_Output_Put_Var_Lev(leo%ncid, leo%varid_lev, leo%nlev, &
                                     trim(leo%levtype), leo%ilev, leo%heights, status)
        IF_NOTOK_RETURN(status=1)
#endif
        do SA_ilabel = 1, SA_nlabel
          labellen = len(trim( SA_Label_Names(SA_ilabel)) )
          status = NF90_Put_Var(leo%ncid, leo%varid_labelnames, trim(SA_Label_Names(SA_ilabel)), start=(/1,SA_ilabel/), count=(/labellen,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do

      end if  ! first record

      ! date up to seconds:
      call Get( t, time6=time6 )

      ! time since t0
      t0 = leo%com%t0      
      time = iTotal( t - t0, 'sec' )

      ! write time record:
#ifdef with_netcdf
      call LE_Output_Put_Var_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                   time, time6, trim(leo%com%CF_convention), leo%itrec, status )
      IF_NOTOK_RETURN(status=1)
#endif

    end if  ! root
    ! 4d field from surface to aloft with labels
    allocate( field(1:nx,1:ny,0:nz+1,1:SA_nlabel) )
    ! 3d field with labels
    allocate( pat(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel) )

    ! loop over all tracer to be written:
    do l = 1, leo%ntr

      ! global tracer index:
      itr = leo%itr_glob(l)

      do k = 1, leo%nlev
        ilev = leo%ilev(k)

        ! init 2d field to zero:
        pat = 0.0
        
        if ( ilev == 0 ) then ! total column output
          pat = sum( SA_Emis(leo%i1:leo%i2,leo%j1:leo%j2,:,leo%itr_lab(l),:), dim=3)
        else                  ! output per model layer
          pat = SA_Emis(leo%i1:leo%i2,leo%j1:leo%j2,ilev,leo%itr_lab(l),:)         
        end if

        ! unit conversion:
        if ( leo%unitconv(l) > 0.0 ) then
          ! single factor ...
          pat = pat * leo%unitconv(l)
        else
          ! specials ...
          select case ( 'tracer in '//trim(specunit(itr))//' -> emis in '//trim(leo%unit_tr(l)) )
            !~ other notations ...
            case ( 'tracer in ug/m3 -> emis in kg m-2 s-1' )
              ! aerosols; emis_a_label is in ug/min
              !     ug/min * kg/ug / s/min
              pat =   pat  * 1e-9  / 60.0  ! kg/s/cell
              ! per area:
              do SA_ilabel = 1, SA_nlabel
                call ugg%AreaOper( pat(1:nx,1:ny, SA_ilabel), '/', 'm2', status )  ! kg/m2/s
                IF_NOTOK_RETURN(status=1)
              end do
            !~ volume mixing ratio to mass concentration:
            case ( 'tracer in ppb -> emis in kg m-2 s-1' )
             ! gasses: emis_a_label is in mol/min
              !                 mole/min      * (kg tr/mole tr) / s/min
              pat(1:nx,1:ny,:) = pat(1:nx,1:ny,:) *  specmolm(itr)  / 60.0  ! kg/s/cell
              ! per area:
              do SA_ilabel = 1, SA_nlabel
                call ugg%AreaOper( pat(1:nx,1:ny, SA_ilabel), '/', 'm2', status )  ! kg/m2/s
                IF_NOTOK_RETURN(status=1)
              end do
            !~ unkown ...
            case default
              write (gol,'("do not know how to convert ",a," emissions from `tracer in ",a,"` to `",a,"`")') &
                      trim(specname(itr)), trim(specunit(itr)), trim(leo%unit_tr(l)); call goErr
              TRACEBACK; status=1; return
          end select
        end if

        ! store:
        field(leo%i1:leo%i2,leo%j1:leo%j2,ilev,:) = pat

      end do   ! model layers
      
      ! * vertical
      
      ! which output levels ?
      select case ( trim(leo%levtype) )

        !~ selected model levels:
        case ( 'levels' )

          ! loop over selected layer:
          do k = 1, leo%nlev
            ! global level index:
            ilev = leo%ilev(k)
            
            do SA_ilabel = 1, SA_nlabel

              ! write 2d field with level and time index
              call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_tr(l), k, leo%itrec, &
                                              field(leo%i1:leo%i2,leo%j1:leo%j2,ilev,SA_ilabel), &
                                              status, ilab=SA_ilabel )
              IF_NOTOK_RETURN(status=1)
            end do
          end do
      
        !~ unknown ...
        case default
          write (gol,'("unsupported level type : ",a)') trim(leo%levtype); call goErr
          TRACEBACK; status=1; return

      end select

    end do   ! tracers
    
    deallocate(pat)
    deallocate(field)

    ! root?
    if ( goc_dom%root ) then
      ! next time exceeds interval ?
      if ( t+IncrDate(hour=int(leo%dhour)) > leo%tr(2) ) then
        ! close:
#ifdef with_netcdf
        status = NF90_Close( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
        ! reset flag:
        leo%opened = .false.
        ! write GrADS ctl file if necessary:
        call Write_GrADS_Ctl( leo, status )
        IF_NOTOK_RETURN(status=1)
      end if
    end if ! root

    ! ok
    status = 0

  end subroutine LE_Output_Emis_Label_PutOut

  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use SA_Labeling, only : SA_nlabel
    use SA_Labeling, only : SA_Label_Names, SA_Short_Label_Names
    use LE_Grid, only : glb_ugg
    

    ! --- in/out ---------------------------------

    type(T_LE_Output_Emis_Label), intent(inout)   ::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l
    character(len=512)    ::  line
    integer               ::  SA_ilabel
    integer               ::  Nvars
    character(len=512)    ::  dim_descr
    character(len=512)    ::  comment_line

    ! --- begin ----------------------------------

    ! write ctl file ?
    if ( leo%grads_ctl ) then

      ! ctl file name:
      write (leo%grads_ctl_file,'(a,"_",a,"_",a)') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_file,'(a,"_",a)') trim(leo%grads_ctl_file), trim(leo%state)
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24) or [00,24] :
        case ( 'daily', 'daily24' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:
        case ( 'instant' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2_%h2%n2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

      ! open ctl file:
      call GrADS_Ctl_Init( ctl, trim(leo%com%outdir), trim(leo%grads_ctl_file), status )
      IF_NOTOK_RETURN(status=1)
      ! comment ...
      call GrADS_Ctl_Comment( ctl, '', status )
      call GrADS_Ctl_Comment( ctl, 'GrADS Data Descriptor File', status )
      call GrADS_Ctl_Comment( ctl, '', status )
      call GrADS_Ctl_Comment( ctl, 'Label definitions', status )
      call GrADS_Ctl_Comment( ctl, '', status )
      do SA_ilabel = 1, SA_nlabel
        write( comment_line , '(a," : ", a)' ) trim(SA_Short_Label_names(SA_ilabel)), trim(SA_Label_Names(SA_ilabel))
        call GrADS_Ctl_Comment( ctl, trim(comment_line) , status )
      end do
      call GrADS_Ctl_Comment( ctl, '', status )      
      ! data set:
      call GrADS_Ctl_DSet( ctl, trim(leo%grads_ctl_base), status )
      IF_NOTOK_RETURN(status=1)
      ! title:
      write (line,'("model: ",a,"; expid: ",a)') trim(leo%com%model), trim(leo%com%expid)
      call GrADS_Ctl_Title( ctl, trim(line), status )
      IF_NOTOK_RETURN(status=1)
      ! write xdef/ydef from grid/projection definition
      call glb_ugg%WriteCtlProjection( ctl, status)
      IF_NOTOK_RETURN(status=1)
      ! zdef:
      call GrADS_Ctl_ZDef( ctl, leo%ilev(1:leo%nlev), status )
      IF_NOTOK_RETURN(status=1)
      ! tdef:
      call GrADS_Ctl_TDef( ctl, leo%grads_ctl_nt, &
                               (/leo%grads_ctl_t1%year,leo%grads_ctl_t1%month,leo%grads_ctl_t1%day,leo%grads_ctl_t1%hour,leo%grads_ctl_t1%min/), &
                               (/                    0,                     0,leo%grads_ctl_dt%day,leo%grads_ctl_dt%hour,leo%grads_ctl_dt%min/), &
                               status )
      IF_NOTOK_RETURN(status=1)
      ! number of variables in ctl file, number of variables x number of labels
      Nvars = leo%ntr * SA_nlabel
      call GrADS_Ctl_Vars( ctl, Nvars, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over tracers to be written:
      do l = 1, leo%ntr
        ! loop over labels to be written
        do SA_ilabel = 1, SA_nlabel
          ! set variable description:
          write (line,'(a," [",a,"]")') trim(leo%name_tr(l)), 'fraction'
          ! set dimension description t, labelnr, z,y,x
          write (dim_descr, '("t, ",i2,", z,y,x")') SA_ilabel-1  
          ! add variable description:
          call GrADS_Ctl_Var( ctl, trim(leo%name_tr(l)), leo%nlev, trim(dim_descr), trim(line), status, label_name = SA_Short_Label_Names(SA_ilabel) )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      ! end of variables section:
      call GrADS_Ctl_EndVars( ctl, status )
      IF_NOTOK_RETURN(status=1)
      ! close ctl file:
      call GrADS_Ctl_Done( ctl, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! ok
    status = 0

  end subroutine Write_GrADS_Ctl

end module LE_Output_Emis_Label
