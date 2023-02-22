!#################################################################
!
! NAME
!
!   LE_Output_Label  -  LOTOS-EUROS output of concentration fields
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

module LE_Output_Labelled_Conc

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Labelled_Conc

  public  ::  LE_Output_Labelled_Conc_Init, LE_Output_Labelled_Conc_Done
  public  ::  LE_Output_Labelled_Conc_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Labelled_Conc'
  
  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10

  ! --- types ------------------------------

  type T_LE_Output_labelled_Conc
    ! name for this file:
    character(len=16)                     ::  typ
    character(len=32)                     ::  name
    ! adhoc convention, might change per file ...
    character(len=32)                     ::  convention
    ! common stuff:
    type(T_LE_Output_Common)              ::  com
    !replace existing file?
    logical                               ::  replace
    ! file opened ?
    logical                               ::  opened
    ! current time range :   
    type(TDate)                           ::  tr(2)
    ! time resolution:
    real                                  ::  dhour
    ! collect: daily, instant
    character(len=32)                     ::  collect
    ! time record counter
    integer                               ::  itrec
    ! state name:
    character(len=16)                     ::  state
    ! file name;
    character(len=512)                    ::  fname
    ! file handle:
    integer                               ::  ncid
    ! Dimension handles:
    integer                               ::  dimid_lon      
    integer                               ::  dimid_lat      
    integer                               ::  dimid_lev      
    integer                               ::  dimid_time     
    integer                               ::  dimid_spec     
    integer                               ::  dimid_label    
    integer                               ::  dimid_labelname_len 
    ! Dimension variables:
    integer                               ::  varid_lon      
    integer                               ::  varid_lat      
    integer                               ::  varid_lev      
    integer                               ::  varid_time     
    integer                               ::  varid_time_dtg     
    integer                               ::  varid_labelnames    
    ! tracer variables
    integer                               ::  ntr         
    integer, pointer                      ::  itr(:)      
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
    ! indices
    integer                               ::  i1, i2, ni
    integer                               ::  j1, j2, nj 
    ! bounding box
    integer                               ::  bb_i1, bb_i2, bb_ni
    integer                               ::  bb_j1, bb_j2, bb_nj
    real                                  ::  westb, southb
    
  end type T_LE_Output_labelled_Conc


contains


  ! ====================================================


  subroutine LE_Output_Labelled_Conc_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues, goSplitString
    use GO     , only : AnyDate
    use Dims   , only : nx, ny, nz
    use Dims   , only : nspec_all
    use LE_Data, only : LE_Data_Enable
    use Indices, only : specname, specunit
    use LE_Grid, only : glb_ugg
    use SA_Labeling     , only : nspec_labelled
    use SA_Labeling     , only : labelled_specs_names

    ! --- in/out --------------------------------

    type(T_LE_Output_Labelled_Conc), intent(out)   ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Labelled_Conc_Init'

    ! --- local ---------------------------------
    
    character(len=64)               ::  basekey
    character(len=512)              ::  output_tracer_names
    character(len=32)               ::  output_level_numbers
    
    character(len=512)              ::  key
    real                            ::  west, east, south, north
    real, pointer                   ::  ff(:,:)
    integer                         ::  i, itr

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
        
    ! Allocate
    allocate( leo%itr     (nspec_all) )
    allocate( leo%name_tr (nspec_all) )
    allocate( leo%varid_tr(nspec_all) )
    allocate( leo%unit_tr (nspec_all) )
    allocate( leo%unitconv(nspec_all) )
    allocate( leo%ilev(1+nz+1) )    !(surface + nz-levels + upperboundary)
    
    ! Match tracer indices
    call goMatchValues( output_tracer_names, specname, &
                        leo%ntr, leo%name_tr, leo%itr, &
                        status )
    IF_NOTOK_RETURN(status=1) 
   
    ! info ...
    write (gol,'("selected tracers for labelled conc output:")'); call goPr
    do i = 1, leo%ntr
      itr = leo%itr(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_tr(i), &
                  itr, trim(specname(itr)), trim(specunit(itr)); call goPr
    end do
                                
    ! Match levelnumbers with selected output
    if ( trim(output_level_numbers) == 'all' ) then
      ! all
      leo%nlev = nz
      do i = 1, leo%nlev
        leo%ilev(i) = i
      end do
    else
      call goMatchValues( output_level_numbers, 0, nz+1, &
                          leo%nlev, leo%ilev, &
                          status )
      IF_NOTOK_RETURN(status=1) 
    end if
    
    ! init GrADS stuff:
    if ( leo%grads_ctl ) then
      ! no times written yet:
      leo%grads_ctl_nt = 0
    end if
    
    ! full domain
    leo%i1 = 1
    leo%i2 = nx
    leo%ni = nx
    leo%j1 = 1
    leo%j2 = ny
    leo%nj = ny

    ! bounding box
    call rcF%Get( trim(basekey)//'.bounding_box', key, status )
    IF_NOTOK_RETURN(status=1)   
    ! empty?     
    if (len_trim(key) == 0) then
      ! full domain
      leo%bb_i1 = 1
      leo%bb_i2 = glb_ugg%nlon
      leo%bb_ni = glb_ugg%nlon
      leo%bb_j1 = 1
      leo%bb_j2 = glb_ugg%nlat
      leo%bb_nj = glb_ugg%nlat      
      leo%westb = glb_ugg%longitude_bnds(1,1,1)
      leo%southb = glb_ugg%latitude_bnds(1,1,1)      
    else
      
      select case ( trim(glb_ugg%type) ) 
        
        case ( 'cartesian-regular') 
          ! read domain from key
          read(key,*,iostat=status) west, east, south, north
          if(status/=0) then
            write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
          endif

          ! Check if bounding box is in run domain
          if ( west < glb_ugg%longitude_bnds_1d(1,1) .or. east > glb_ugg%longitude_bnds_1d(2,glb_ugg%nlon) .or. &
               south < glb_ugg%latitude_bnds_1d(1,1) .or. north > glb_ugg%latitude_bnds_1d(2,glb_ugg%nlat) ) then
            write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
            write( gol, '("Run domain: ", 4f8.2)' ) glb_ugg%longitude_bnds_1d(1,1),glb_ugg%longitude_bnds_1d(2,glb_ugg%nlon),glb_ugg%latitude_bnds_1d(1,1),glb_ugg%latitude_bnds_1d(2,glb_ugg%nlat); call goErr
            write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
            TRACEBACK;status=1;return
          endif

          ! for safety
          nullify(ff)
          ! get cell range covered by box
          call glb_ugg%GetDistribution(west,east,south,north,leo%bb_i1,leo%bb_i2,leo%bb_j1,leo%bb_j2,ff,status)
          IF_NOTOK_RETURN(status=1)
          !clear, fractions not used
          if ( associated(ff) ) deallocate(ff)
          ! set shape
          leo%bb_ni = leo%bb_i2-leo%bb_i1+1
          leo%bb_nj = leo%bb_j2-leo%bb_j1+1
          ! set west/south bounds
          leo%westb  = glb_ugg%longitude_bnds_1d(1,leo%bb_i1)
          leo%southb = glb_ugg%latitude_bnds_1d(1,leo%bb_j1)
        case default 
          write( gol, '("Definition of bounding box not clear for grid-type: ", a)' ) trim(glb_ugg%type) ; call goErr
          TRACEBACK;status=1;return
      end select
    end if

    ! enable data:
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Labelled_Conc_Init


  ! ***


  subroutine LE_Output_Labelled_Conc_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Labelled_Conc), intent(inout)   ::  leo
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Labelled_Conc_Done'
    
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
      IF_NOTOK_RETURN(status=1)
    end if 
    
    ! deallocate type output variables
    deallocate( leo%itr )
    deallocate( leo%name_tr  )
    deallocate( leo%varid_tr )
    deallocate( leo%unit_tr  )
    deallocate( leo%unitconv )

    ! done with common stuff ...
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Output_Labelled_Conc_Done


  ! ***


  subroutine LE_Output_Labelled_Conc_PutOut( leo, t, c, cg, status )
    
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

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    use LE_Data      , only : LE_Data_GetPointer

    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO, only : TIncrDate, NewDate, IncrDate, Get
    use GO, only : GoMatchValue
    use LE_Comm    , only : goc_dom
    
    use Binas      , only : xm_air, xm_nh4, xm_no3, xm_so4
    use Dims       , only : nx, ny, nz
    use LE_Grid    , only : glb_ugg
    use C3PO       , only : T_Grid_NcDef
    
    use Indices    , only : nspec, specname, specunit, specmolm
    use Indices    , only : accum_n, accum_ww, accum_ii, accum_ppb_to_ugm3
    use Indices    , only : ispec_no3a_f, ispec_so4a_f, ispec_nh4a_f
    use SA_Labeling, only : nspec_labelled
    use SA_Labeling, only : SA_frac, SA_nlabel
    use SA_Labeling, only : SA_Label_Names, labelname_len
    use SA_Labeling, only : labelled_specs_names

    ! --- in/out --------------------------------

    type(T_LE_Output_Labelled_Conc), intent(inout)   ::  leo
    type(TDate), intent(in)                 ::  t
    real, intent(in)                        ::  c (nx,ny,nz,nspec)
    real, intent(in)                        ::  cg(nx,ny,nspec)
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Labelled_Conc_PutOut'

    ! --- local ---------------------------------
    
    integer                                 ::  time6(6)
    integer                                 ::  time
    integer                                 ::  cmode
    
    integer                                 ::  ix,iy,iz,ispec, itr
    integer                                 ::  i, j
    integer                                 ::  k, l, la
    integer                                 ::  itr_conc, itr_frac, ilev
    integer                                 ::  varid
    type(T_Grid_NcDef)                      ::  gncd
    type(TDate)                             ::  t0
    real                                    ::  inv_denum
    integer                                 ::  SA_ilabel
    integer                                 ::  labellen
    real, allocatable                       ::  pat(:,:,:) ! nx,ny,SA_nlabel
    real, allocatable                       ::  convfact(:,:) ! nx,ny
    
    
    ! specials for combination of NH4NO3 and NH4_on_SO4in labelled output
    real, allocatable                       ::  mol_nh4no3(:,:)
    real, allocatable                       ::  mass_nh4no3(:,:)
    real, allocatable                       ::  frac_nh4(:,:,:)
    real, allocatable                       ::  frac_no3(:,:,:)
    real, allocatable                       ::  mass_nh4no3_lab(:,:,:)
    
    real, allocatable                       ::  mol_nh4(:,:)
    real, allocatable                       ::  mol_nh4_on_so4(:,:)
    real, allocatable                       ::  mass_nh4_on_so4(:,:)
    real, allocatable                       ::  mass_nh4_on_so4_lab(:,:,:)
    
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units, cf_tracer_name
    character(len=512)    ::  comment
    character(len=32)     ::  varname
    character(len=32)     ::  afac

    character(len=256) :: commandline
    
    ! meteo data:
    real, pointer        ::  dens(:,:,:)   ! (lon,lat,lev)    
    ! --- begin ---------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)
    
    !$OMP parallel &
    !$OMP default( none ) &
    !$OMP shared( SA_frac ) &
    !$OMP shared( nx, ny, nz, nspec_labelled, SA_nlabel ) &
    !$OMP private( ix,iy,iz,SA_ilabel,ispec ) &
    !$OMP private( inv_denum )
    !$OMP do SCHEDULE(DYNAMIC)
    ! Fractions should add up to one, otherwise scale the fractions
    do ix = 1, nx
    do iy = 1, ny
    do iz = 1, nz
      do ispec = 1, nspec_labelled
        if (abs( sum(SA_frac(ix,iy,iz,ispec,:)) - 1.0) > 1.e-6 .and. sum(SA_frac(ix,iy,iz,ispec,:)) > 0.0 ) then
          inv_denum = 1.0 / sum(SA_frac(ix,iy,iz,ispec,:))
          do SA_ilabel = 1, SA_nlabel
            SA_frac(ix,iy,iz,ispec,SA_ilabel) = SA_frac(ix,iy,iz,ispec,SA_ilabel) * inv_denum
          enddo
        endif
      enddo
    enddo
    enddo
    enddo
    !$OMP end do
    !$OMP end parallel
    
    ! 2d field + label dimension
    allocate( pat(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! conversion field
    allocate( convfact(nx,ny) )    
    
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
      
      call Get(t, time6 = time6)
      
      ! set time range [00,24) for this day:
      leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)
      
      write (leo%fname, '(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                trim(leo%name), time6(1:3)
      write (leo%fname,'(a,".nc")') trim(leo%fname)            
      
      ! root only :
      if ( goc_dom%root ) then
        if (leo%replace) then
          cmode = NF90_CLOBBER    ! overwrite existing file
        else
          cmode = NF90_NOCLOBBER  ! do not overwrite existing file
        endif

        ! enable large file support:
        cmode = or( cmode, NF90_NETCDF4 )

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
                                    dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat, &
                                    subset=(/leo%bb_i1,leo%bb_i2,leo%bb_j1,leo%bb_j2/) )

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
      
      end if  ! root
                  
      ! loop over tracers to be written:
      do l = 1, leo%ntr

        ! global tracer index
        itr = leo%itr(l)

        ! CF standard name for concentration/mixing ratio/column:

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
                         status )
            IF_NOTOK_RETURN(status=1)

          case ( 'megapoli' )

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
            cf_standard_name = 'mass_concentration_of_'//trim(cf_tracer_name)//'_in_air'
            cf_long_name     = trim(cf_standard_name)
            cf_units = 'ug m-3'
            leo%unitconv(l) = -999.9    ! this enables conversion from ppb to ug/m3 below ...

          case default

            write (gol,'("unsupported output convention : ",a)') trim(leo%convention); call goErr
            TRACEBACK; status=1; return

        end select

        ! store units for later usage (GrADS ctl file, conversions):
        leo%unit_tr(l) = trim(cf_units)

        ! accumulated specs ?
        if ( accum_n(itr) > 1 ) then
          if ( len_trim(comment) > 0 ) comment = trim(comment)//' ; '
          do la = 1, accum_n(itr)
            if ( la > 1 ) comment = trim(comment)//' + '
            if ( accum_ppb_to_ugm3(itr,la) ) then
              write (afac,'(f5.2)') accum_ww(itr,la)/xm_air
              comment = trim(comment)//trim(specname(accum_ii(itr,la)))//'*'//trim(afac)//'*dens'
            else
              write (afac,'(f5.2)') accum_ww(itr,la)
              comment = trim(comment)//trim(specname(accum_ii(itr,la)))//'*'//trim(afac)
            end if
          end do
        end if
        
        ! root?
        if ( goc_dom%root ) then
          ! define variable
          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                              (/leo%dimid_lon, leo%dimid_lat, leo%dimid_lev &
                              , leo%dimid_label, leo%dimid_time/), varid, &
                                 deflate_level=leo%com%deflate_level )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write attributes:
          if ( len_trim(cf_standard_name) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
          status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
          IF_NF90_NOTOK_RETURN(status=1)
          ! molemass:
          status = nf90_put_att( leo%ncid, varid, 'molemass', specmolm(itr) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'molemass_unit', 'kg mole-1' )
          IF_NF90_NOTOK_RETURN(status=1)
          ! optional comment:
          if ( len_trim(comment) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if

          ! store variable id:
          leo%varid_tr(l) = varid
          
        end if  ! root
      enddo

      ! root?
      if ( goc_dom%root ) then
        ! end defintion mode:
        status = NF90_EndDef(leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
        
      ! no records written yet
      leo%itrec = 0
    
    endif
      
    ! record the time step
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
    
    ! write tracer variables
    do l = 1, leo%ntr
      itr = leo%itr(l)
      do k = 1, leo%nlev
        ilev = leo%ilev(k)
        
        ! nullify
        pat = 0.0
        
        if ( ilev == 0 ) then
          do la = 1, accum_n(itr) 
            
            ! special cases for nh4a_f_no3a_f and nh4a_f_on_so4a_f
            if ( trim( specname(accum_ii(itr,la)) ) == 'nh4a_f_no3a_f' ) then              
              
              ! Combine labels for NH4 and NO3, with combined aerosol 50% of NH4 origin ann 50% NO3 origin
              ! Note original distribution is based on mass ( so NO3 origin has to large contribution)                            
              allocate( mol_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_nh4(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_no3(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4no3_lab(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              
              ! # mol NH4NO3 (equal to # mol NO3)
              mol_nh4no3  = cg(leo%i1:leo%i2,leo%j1:leo%j2,ispec_no3a_f) / xm_no3
              mass_nh4no3 = mol_nh4no3 * (xm_no3+xm_nh4)
              
              ! Labelled fractions of NH4
              call GoMatchValue( 'nh4a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_nh4 = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,1,itr_frac,:)

              ! Labelled fractions of NO3
              call GoMatchValue( 'no3a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_no3 = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,1,itr_frac,:)
              
              ! Divide equally to total mass of NH4NO3              
              do SA_ilabel = 1, SA_nlabel
               mass_nh4no3_lab(:,:,SA_ilabel) = &
                 ( 0.5*frac_nh4(:,:,SA_ilabel) + 0.5*frac_no3(:,:,SA_ilabel) ) * mass_nh4no3
              end do
              
              ! add to total PM
              pat = pat + mass_nh4no3_lab
              
              ! clear
              deallocate( mol_nh4no3 )
              deallocate( mass_nh4no3 )
              deallocate( frac_nh4 )
              deallocate( frac_no3 )
              deallocate( mass_nh4no3_lab )
                            
            else if ( trim(specname(accum_ii(itr,la)) ) == 'nh4a_f_on_so4a_f' ) then     
              
              allocate( mol_nh4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mol_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mol_nh4_on_so4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4_on_so4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_nh4(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4_on_so4_lab(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
                            
              ! total moles of NH4
              mol_nh4    = cg(leo%i1:leo%i2,leo%j1:leo%j2,ispec_nh4a_f) / xm_nh4
              ! amount of moles attached to NO3
              mol_nh4no3 = cg(leo%i1:leo%i2,leo%j1:leo%j2,ispec_no3a_f) / xm_no3
              
              ! Rest of NH4 attached on SO4
              mol_nh4_on_so4  = mol_nh4 - mol_nh4no3
              mass_nh4_on_so4 = mol_nh4_on_so4 * xm_nh4
              
              ! Labelled fractions of NH4
              call GoMatchValue( 'nh4a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_nh4 = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,1,itr_frac,:)
              
              ! mass of NH4 on SO4 per label
              do SA_ilabel = 1, SA_nlabel
                mass_nh4_on_so4_lab(:,:,SA_ilabel) = mass_nh4_on_so4 * frac_nh4(:,:,SA_ilabel)
              end do
              
              ! add to total PM
              pat = pat + mass_nh4_on_so4_lab
              
              deallocate( mol_nh4 )
              deallocate( mol_nh4no3 )
              deallocate( mol_nh4_on_so4 )
              deallocate( mass_nh4_on_so4 )
              deallocate( frac_nh4 )
              deallocate( mass_nh4_on_so4_lab )
              
            else          
              ! match with labelled tracer number to get fractions
              call GoMatchValue( specname(accum_ii(itr,la)), labelled_specs_names(1:nspec_labelled), itr_frac, status )
              if ( status /= 0 ) then
                write( gol, '("****")' ); call goErr
                write( gol, '("No matching labelfractions, you probably try to put out total PM, while (at least) one of the subtracers is not labelled")' ) ; call goErr
                write( gol, '("****")' ); call goErr
                TRACEBACK;status=1;return
              end if

              ! conversion needed?
              convfact = 1.0
              if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(1:nx,1:ny,1)/xm_air
              ! Surface level concentrations, take fractions from first layer
              do SA_ilabel = 1, SA_nlabel
                pat(:,:,SA_ilabel) = pat(:,:,SA_ilabel) + &
                  SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,1,itr_frac,SA_ilabel) * &
                    cg(leo%i1:leo%i2,leo%j1:leo%j2,accum_ii(itr,la) ) * &
                    convfact(leo%i1:leo%i2,leo%j1:leo%j2) * accum_ww(itr,la)
              end do ! labels
            end if  ! check special subtracers
          end do ! subspecs
        else
          ! add contributions
          do la = 1, accum_n(itr) 

            ! special cases for nh4a_f_no3a_f and nh4a_f_on_so4a_f
            if ( trim( specname(accum_ii(itr,la)) ) == 'nh4a_f_no3a_f' ) then              
              ! Combine labels for NH4 and NO3, with combined aerosol 50% of NH4 origin ann 50% NO3 origin
              ! Note original distribution is based on mass ( so NO3 origin has to large contribution)                            
              allocate( mol_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_nh4(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_no3(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4no3_lab(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              
              ! # mol NH4NO3 (equal to # mol NO3)
              mol_nh4no3 = c(leo%i1:leo%i2,leo%j1:leo%j2,ilev,ispec_no3a_f) / xm_no3
              mass_nh4no3 = mol_nh4no3 * (xm_no3+xm_nh4)
              
              ! Labelled fractions of NH4
              call GoMatchValue( 'nh4a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_nh4 = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,ilev,itr_frac,:)

              ! Labelled fractions of NO3
              call GoMatchValue( 'no3a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_no3 = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,ilev,itr_frac,:)
              
              ! Divide equally to total mass of NH4NO3              
              do SA_ilabel = 1, SA_nlabel
               mass_nh4no3_lab(:,:,SA_ilabel) = &
                 ( 0.5*frac_nh4(:,:,SA_ilabel) + 0.5*frac_no3(:,:,SA_ilabel) ) * mass_nh4no3
              end do
              
              ! add to total PM
              pat = pat + mass_nh4no3_lab
              
              ! clear
              deallocate( mol_nh4no3 )
              deallocate( mass_nh4no3 )
              deallocate( frac_nh4 )
              deallocate( frac_no3 )
              deallocate( mass_nh4no3_lab )
                            
            else if ( trim(specname(accum_ii(itr,la)) ) == 'nh4a_f_on_so4a_f' ) then     
              
              allocate( mol_nh4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mol_nh4no3(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mol_nh4_on_so4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4_on_so4(leo%i1:leo%i2,leo%j1:leo%j2), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( frac_nh4(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
              allocate( mass_nh4_on_so4_lab(leo%i1:leo%i2,leo%j1:leo%j2,1:SA_nlabel), stat=status )
              IF_NOTOK_RETURN(status=1)
                            
              ! total moles of NH4
              mol_nh4    = c(leo%i1:leo%i2,leo%j1:leo%j2,ilev,ispec_nh4a_f) / xm_nh4
              ! amount of moles attached to NO3
              mol_nh4no3 = c(leo%i1:leo%i2,leo%j1:leo%j2,ilev,ispec_no3a_f) / xm_no3
              
              ! Rest of NH4 attached on SO4
              mol_nh4_on_so4  = mol_nh4 - mol_nh4no3
              mass_nh4_on_so4 = mol_nh4_on_so4 * xm_nh4
              
              ! Labelled fractions of NH4
              call GoMatchValue( 'nh4a_f', labelled_specs_names(1:nspec_labelled), itr_frac, status )
              IF_NOTOK_RETURN(status=1)
              frac_nh4(:,:,:) = SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,ilev,itr_frac,:)
              
              ! mass of NH4 on SO4 per label
              do SA_ilabel = 1, SA_nlabel
                mass_nh4_on_so4_lab(:,:,SA_ilabel) = mass_nh4_on_so4 * frac_nh4(:,:,SA_ilabel)
              end do
              
              ! add to total PM
              pat = pat + mass_nh4_on_so4_lab
              
              deallocate( mol_nh4 )
              deallocate( mol_nh4no3 )
              deallocate( mol_nh4_on_so4 )
              deallocate( mass_nh4_on_so4 )
              deallocate( frac_nh4 )
              deallocate( mass_nh4_on_so4_lab )
              
            else          
              ! match with labelled tracer number to get fractions
              call GoMatchValue( specname(accum_ii(itr,la)), labelled_specs_names(1:nspec_labelled), itr_frac, status )
              if ( status /= 0 ) then
                write( gol, '("****")' ); call goErr
                write( gol, '("No matching labelfractions, you probably try to put out total PM, while (at least) one of the subtracers is not labelled")' ) ; call goErr
                write( gol, '("****")' ); call goErr
                TRACEBACK;status=1;return
              end if

              ! conversion needed ?
              convfact = 1.0
              if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(1:nx,1:ny,ilev)/xm_air
              ! higher level concentrations
              do SA_ilabel = 1, SA_nlabel
                pat(:,:,SA_ilabel) = pat(:,:,SA_ilabel) + &
                  SA_frac(leo%i1:leo%i2,leo%j1:leo%j2,ilev,itr_frac,SA_ilabel) * &
                    c(leo%i1:leo%i2,leo%j1:leo%j2,ilev,accum_ii(itr,la) ) *&
                  convfact(leo%i1:leo%i2,leo%j1:leo%j2) * accum_ww(itr,la)
              end do ! labels
              
            end if  ! special tracer check
          end do ! subspecs
        end if ! level check
        
        ! unit conversion:
        if ( leo%unitconv(l) > 0.0 ) then
          ! single factor ...
          pat = pat * leo%unitconv(l)
        else
          ! specials ...
          select case ( trim(specunit(itr))//' -> '//trim(leo%unit_tr(l)) )
            !~ other notations ...
            case ( 'ug/m3 -> ug m-3' )
              ! no further conversion
            !~ volume mixing ratio to mass concentration:
            case ( 'ppb -> ug m-3' )
              ! 1-nz, index in density array
              iz = min( max( 1, ilev ), nz )
              do SA_ilabel = 1, SA_nlabel
                ! ug/m3            =  ppb                * ((mole tr/mole air)/ppb) *  (kg tr/mole tr) *             (kg air/m3 air)         /  (kg air/mole air) * ug/kg
                pat(:,:,SA_ilabel) = pat(:,:,SA_ilabel)  *       1e-9               *   specmolm(itr)  * dens(leo%i1:leo%i2,leo%j1:leo%j2,iz) /     xm_air        * 1e9
              end do
            !~ unkown ...
            case default
              write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                      trim(specname(itr)), trim(specunit(itr)), trim(leo%unit_tr(l)); call goErr
              TRACEBACK; status=1; return
          end select
        end if ! unit conversion
        
        ! Write out per label
        do SA_ilabel = 1, SA_nlabel                          
          ! write 2D field with level and time index:
          call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_tr(l), k, leo%itrec, &
                                          pat(:,:,SA_ilabel), status, bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/), ilab=SA_ilabel )
          IF_NOTOK_RETURN(status=1)
        end do
        
      enddo ! levels
    enddo ! tracers
    
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
    end if  ! root?

    ! Done
    deallocate(pat)
    deallocate(convfact)
    
    ! ok
    status = 0

  end subroutine LE_Output_Labelled_Conc_PutOut

  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use SA_Labeling, only : SA_nlabel
    use SA_Labeling, only : SA_Label_Names, SA_Short_Label_Names
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Labelled_Conc), intent(inout)   ::  leo
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
        ! collect daily for [00,24)
        case ( 'daily' )
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
      select case ( trim(leo%levtype) )
        case ( 'levels' )
          call GrADS_Ctl_ZDef( ctl, leo%ilev(1:leo%nlev), status )
          IF_NOTOK_RETURN(status=1)          
        case ( 'heights', 'elevations' )
          call GrADS_Ctl_ZDef( ctl, nint(leo%heights(1:leo%nlev)), status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("do not know how to write zdef for levtyp : ",a)') trim(leo%levtype); call goPr
          TRACEBACK; status=1; return
      end select
     ! tdef:
      call GrADS_Ctl_TDef( ctl, leo%grads_ctl_nt, &
                               (/leo%grads_ctl_t1%year,leo%grads_ctl_t1%month,leo%grads_ctl_t1%day,leo%grads_ctl_t1%hour,leo%grads_ctl_t1%min/), &
                               (/                    0,                     0,leo%grads_ctl_dt%day,leo%grads_ctl_dt%hour,leo%grads_ctl_dt%min/), &
                               status )
      IF_NOTOK_RETURN(status=1)
      
      ! number of variables in ctl file, number of variables time number of labels
      Nvars = leo%ntr * SA_nlabel      
      call GrADS_Ctl_Vars( ctl, Nvars, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over tracers to be written:
      do l = 1, leo%ntr
        ! loop over labels to be written
        do SA_ilabel = 1, SA_nlabel
          ! set variable description:
          write (line,'(a," [",a,"]")') trim(leo%name_tr(l)), trim(leo%unit_tr(l))
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

end module LE_Output_Labelled_Conc
