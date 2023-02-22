!###############################################################################
!
! landuse  -  Read landuse database from smiatek
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

module LE_LandUse_Smiatek

  use GO, only : gol, goPr, goErr
  
  use LE_Landuse_Data

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  LE_LandUse_Smiatek_Init
  
  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_Landuse_Smiatek'

  ! number of landuse classes:
  integer, parameter  ::  nlu_smtk = 13
  !
  ! the  lotos landuse database consists of the following "Smiatek" clases:
  integer, parameter  ::  ilu_smtk_urban_areas        =  1
  integer, parameter  ::  ilu_smtk_agriculture        =  2
  integer, parameter  ::  ilu_smtk_grassland          =  3
  integer, parameter  ::  ilu_smtk_deciduous_forest   =  4
  integer, parameter  ::  ilu_smtk_coniferous_forest  =  5
  integer, parameter  ::  ilu_smtk_mixed_forest       =  6
  integer, parameter  ::  ilu_smtk_water              =  7
  integer, parameter  ::  ilu_smtk_marsh_or_wetland   =  8
  integer, parameter  ::  ilu_smtk_sand_or_bare_rocks =  9
  integer, parameter  ::  ilu_smtk_tundra             = 10
  integer, parameter  ::  ilu_smtk_permanent_ice      = 11
  integer, parameter  ::  ilu_smtk_tropical_forest    = 12
  integer, parameter  ::  ilu_smtk_woodland_scrub     = 13
  
  ! --- var --------------------------------------

contains


  ! ========================================================================


  subroutine LE_Landuse_Smiatek_Init( query, status )

    use GO, only : goVarValue
    use dims, only : nx, ny
    
    ! --- in/out ---------------------------------

    character(len=*), intent(in)  ::  query
    integer, intent(out)          ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Smiatek_Init'
    ! --- local ----------------------------------

    character(len=512)    ::  fname
    real, allocatable     ::  lu_smtk(:,:,:)              ! (nx,ny,nlu_smtk)
    real, allocatable     ::  veg_smtk(:,:,:)             ! (nx,ny,nveg)
    real, allocatable     ::  veg_smtk_to_depac(:,:,:)    ! (nx,ny,nlu)

    ! --- begin ----------------------------------
    
    ! storage:
    allocate( lu_smtk(nx,ny,nlu_smtk), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( veg_smtk(nx,ny,nveg), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( veg_smtk_to_depac(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! * landuse map, vegetation map:    
    lu_smtk  = 0.0
    veg_smtk = 0.0
    veg_smtk_to_depac = 0.0
    
    ! name of landuse file:
    call goVarValue( trim(query), ';', 'file', '=', fname, status )
    IF_NOTOK_RETURN(status=1)

    ! Read land use database:
    call get_lu( trim(fname), lu_smtk, veg_smtk, status )
    IF_NOTOK_RETURN(status=1)
   
    ! setup DEPAC landuse classes for deposition:
    call Translate_LU_to_Depac( lu_smtk, veg_smtk, veg_smtk_to_depac, status)
    IF_NOTOK_RETURN(status=1)
    
    ! fill in vegetation data in common array
    ! For Smiatek this is identical
    veg_fracs(:,:,1:ntree_type)      = veg_smtk(:,:,1:ntree_type)
    veg_fracs(:,:,ntree_type+1:nveg) = veg_smtk_to_depac(:,:,:)
    
    ! clear:
    deallocate( lu_smtk, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( veg_smtk, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( veg_smtk_to_depac, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Landuse_Smiatek_Init


  ! ***


  subroutine get_lu( fname, lu_smtk, veg_smtk, status )

    use GO          , only : pathsep
    use GO          , only : goGetFU
    use Grid        , only : TllGridInfo, Init, Done, IndexFractions

    use dims        , only : runF
    use dims        , only : nx, ny

    use LE_Grid     , only : lli
    use LE_Grid     , only : indomain
    use LE_LandUse_File_nc, only : T_Landuse_File_nc, Init, Done, Get

    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  fname
    real, intent(inout)             ::  lu_smtk(:,:,:)
    real, intent(inout)             ::  veg_smtk(:,:,:)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/get_lu'

    ! --- local ------------------------------------------

    integer                 ::  l, k
    character(len=3)        ::  ext

    ! landuse and grid:
    real                    ::  west_bound, east_bound, dlon
    real                    ::  south_bound, north_bound, dlat
    integer                 ::  nlon, nlat
    integer, allocatable    ::  landuseX(:,:)
    integer, allocatable    ::  vegetationX(:,:)
    type(TllGridInfo)       ::  lliX

    ! read from nc file:
    type(T_Landuse_File_nc) ::  luf

    ! read from .txt file:
    integer                 ::  fu
    character(len=256)      ::  line
    integer                 ::  iline, ilon, ilat
    real                    ::  cell_lon, cell_lat
    integer                 ::  cell_landuse, cell_forest

    ! --- begin ------------------------------------------

    ! info ...
    write (gol,'("read landuse data base ...")'); call goPr

    ! extract extension:
    l = len_trim(fname)
    k = index( fname(1:l), '.', back=.true. )
    if ( (k == 0) .or. (l-k < 1) .or. (l-k > 3) ) then
      write (gol,'("could not extract 3-character extension of : ",a)') trim(fname); call goErr
      IF_NOTOK_RETURN(status=1)
    end if
    ext = fname(k+1:l)

    ! open using different routines given extension:
    select case ( ext )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'txt' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("    WARNING - reading landuse from a txt file will take a while ...")'); call goPr
        write (gol,'("    WARNING - to avoid getting borred, try reading from a nc file in future!")'); call goPr

        !
        ! Example lines:
        ! ------------------------
        !  lon      lat     landuse    tree
        !  -35.992  29.008    7    0
        !   :
        ! ------------------------
        !

        ! free file unit:
        call goGetFU( fu, status )
        IF_NOTOK_RETURN(status=1)

        ! open file:
        open( fu, file=trim(fname), status='old', form='formatted', iostat=status )
        if (status/=0) then
          write (gol,'("opening land use file :")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          IF_NOTOK_RETURN(status=1)
        endif

        ! read header:
        iline = 1
        read (fu,'(a)',iostat=status) line
        if (status/=0) then
          write (gol,'("reading header line from :")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          IF_NOTOK_RETURN(status=1)
        endif

        ! read data:
        ilon = 1
        ilat = 0
        do

          ! next line:
          iline = iline + 1

          ! read line:
          read (fu,'(2f8.3,2i5)',iostat=status) cell_lon, cell_lat, cell_landuse, cell_forest
          if (status < 0) exit  ! eof
          if (status/=0) then
            write (gol,'("reading line from :")'); call goErr
            write (gol,'("  file : ",a)') trim(fname); call goErr
            write (gol,'("  line : ",i6)') iline; call goErr
            IF_NOTOK_RETURN(status=1)
          endif

          ! guess grid ...
          if ( iline == 2 ) then
            if ( (cell_lon==-35.992) .and. (cell_lat==29.008) ) then
              west_bound = -36.0 ; dlon = 1.0/60.0 ; nlon = 5520
              south_bound = 29.0 ; dlat = 1.0/60.0 ; nlat = 2760
            else
              write (gol,'("could not guess the grid from the first point:")'); call goErr
              write (gol,'("  lon, lat : ",2f8.3)') cell_lon, cell_lat; call goErr
              IF_NOTOK_RETURN(status=1)
            end if
            ! setup storage:
            allocate( landuseX(nlon,nlat) )
            allocate( vegetationX(nlon,nlat) )
          end if

          ! next cell index:
          ilat = ilat + 1
          if ( ilat > nlat ) then
            ilon = ilon + 1
            ilat = 1
          end if

          ! store:
          landuseX(ilon,ilat) = cell_landuse
          vegetationX (ilon,ilat) = cell_forest

        end do

        ! close:
        close( fu, iostat=status )
        if (status/=0) then
          write (gol,'("closing land use file")'); call goErr
          IF_NOTOK_RETURN(status=1)
        endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'nc' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! open database:
        call Init( luf, trim(fname), status )
        IF_NOTOK_RETURN(status=1)

        ! extract grid parameters:
        call Get( luf, status, &
                    west_bound=west_bound, east_bound=east_bound, nlon=nlon, &
                    south_bound=south_bound, north_bound=north_bound, nlat=nlat )
        IF_NOTOK_RETURN(status=1)

        ! spacing:
        dlon = (  east_bound -  west_bound ) / nlon
        dlat = ( north_bound - south_bound ) / nlat

        ! setup storage:
        allocate( landuseX(nlon,nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( vegetationX(nlon,nlat), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read:
        call Get( luf, status, landuse=landuseX )
        IF_NOTOK_RETURN(status=1)
        call Get( luf, status, vegetation=vegetationX )
        IF_NOTOK_RETURN(status=1)

        ! close:
        call Done( luf, status )
        IF_NOTOK_RETURN(status=1)

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported extension for landuse file : ",a)') trim(ext); call goErr
        IF_NOTOK_RETURN(status=1)

    end select

         !! dump
         !call nc_dump( 'lux-'//trim(ext)//'.nc', landuseX, 'landuseX', (/'x','y'/), status )
         !IF_NOTOK_RETURN(status=1)

    write (gol,'("  convert to LE grid ...")'); call goPr
    
    ! setup input grid:
    call Init( lliX, west_bound+0.5*dlon, dlon, nlon, &
                    south_bound+0.5*dlat, dlat, nlat, status )
    IF_NOTOK_RETURN(status=1)

    ! convert from 2D field with indices to 3D field with fraction per cell for each index:
    call IndexFractions( lliX, landuseX, lli, 1, nlu_smtk, lu_smtk, status )
    IF_NOTOK_RETURN(status=1)
    call IndexFractions( lliX, vegetationX, lli, 1, nveg, veg_smtk, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( lliX, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    if ( allocated( landuseX   ) ) then
      deallocate( landuseX, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    if ( allocated( vegetationX) ) then
      deallocate( vegetationX, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine get_lu


  ! ***


  subroutine Translate_LU_to_Depac(lu_smtk, veg_smtk, veg_smtk_to_depac, status )
  
    use dims            , only :  nx, ny

    ! --- in/out ---------------------------------
    
    real, intent(in)              ::  lu_smtk(nx,ny,nlu_smtk)
    real, intent(in)              ::  veg_smtk(nx,ny,nveg)
    real, intent(inout)           ::  veg_smtk_to_depac(nx,ny,nlu)
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   :: rname = mname//'/Translate_LU_to_Depac'
    
    ! --- local ----------------------------------
    
    integer                       ::  i, j
    
    ! --- begin ------------
    
    ! storage:        
    do i =1,nx
      do j = 1, ny
      
        ! assign depac classes from smiatek classes
        lu_fracs(i,j,ilu_grass)              =       lu_smtk(i,j,ilu_smtk_grassland)
        lu_fracs(i,j,ilu_arable)             = 0.9 * lu_smtk(i,j,ilu_smtk_agriculture)
        lu_fracs(i,j,ilu_permanent_crops)    = 0.1 * lu_smtk(i,j,ilu_smtk_agriculture)
        lu_fracs(i,j,ilu_coniferous_forest)  =       lu_smtk(i,j,ilu_smtk_coniferous_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_mixed_forest) 
        lu_fracs(i,j,ilu_deciduous_forest)   =       lu_smtk(i,j,ilu_smtk_deciduous_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_mixed_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_marsh_or_wetland) + &
                                                     lu_smtk(i,j,ilu_smtk_woodland_scrub) + &
                                                     lu_smtk(i,j,ilu_smtk_tropical_forest)                                                                                               
        lu_fracs(i,j,ilu_water_sea)          =       lu_smtk(i,j,ilu_smtk_water) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_marsh_or_wetland) + &
                                                     lu_smtk(i,j,ilu_smtk_permanent_ice)   
        lu_fracs(i,j,ilu_urban)              =       lu_smtk(i,j,ilu_smtk_urban_areas)
        lu_fracs(i,j,ilu_other)              =       lu_smtk(i,j,ilu_smtk_tundra)
        lu_fracs(i,j,ilu_desert)             =       lu_smtk(i,j,ilu_smtk_sand_or_bare_rocks)
        lu_fracs(i,j,ilu_ice)                =       0.0  
!        lu_fracs(i,j,ilu_ice)                =       lu_smtk(i,j,ilu_smtk_permanent_ice)   
        lu_fracs(i,j,ilu_savanna)            =       0.0
        lu_fracs(i,j,ilu_tropical_forest)    =       0.0
!        lu_fracs(i,j,ilu_tropical_forest)    =       lu_smtk(i,j,ilu_smtk_tropical_forest)
        lu_fracs(i,j,ilu_water_inland)       =       0.0
        lu_fracs(i,j,ilu_mediterrean_scrub)  =       0.0
        lu_fracs(i,j,ilu_semi_natural_veg)   =       0.0
        if ( with_ozone_specials ) then
          lu_fracs(i,j,ilu_wheat)            =       0.0
          lu_fracs(i,j,ilu_beech)            =       0.0
          lu_fracs(i,j,ilu_spruce)           =       0.0
        end if
        
        
        ! translation of common vegations classes
        veg_smtk_to_depac(i,j,ilu_grass)              =       veg_smtk(i,j,ntree_type+ilu_smtk_grassland)
        veg_smtk_to_depac(i,j,ilu_arable)             = 0.9 * veg_smtk(i,j,ntree_type+ilu_smtk_agriculture)
        veg_smtk_to_depac(i,j,ilu_permanent_crops)    = 0.1 * veg_smtk(i,j,ntree_type+ilu_smtk_agriculture)
        veg_smtk_to_depac(i,j,ilu_coniferous_forest)  =       veg_smtk(i,j,ntree_type+ilu_smtk_coniferous_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_mixed_forest) 
        veg_smtk_to_depac(i,j,ilu_deciduous_forest)   =       veg_smtk(i,j,ntree_type+ilu_smtk_deciduous_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_mixed_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_marsh_or_wetland) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_woodland_scrub) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_tropical_forest)                                                                                                   
        veg_smtk_to_depac(i,j,ilu_water_sea)          =       veg_smtk(i,j,ntree_type+ilu_smtk_water) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_marsh_or_wetland) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_permanent_ice)
        veg_smtk_to_depac(i,j,ilu_urban)              =       veg_smtk(i,j,ntree_type+ilu_smtk_urban_areas)
        veg_smtk_to_depac(i,j,ilu_other)              =       veg_smtk(i,j,ntree_type+ilu_smtk_tundra)
        veg_smtk_to_depac(i,j,ilu_desert)             =       veg_smtk(i,j,ntree_type+ilu_smtk_sand_or_bare_rocks)
        veg_smtk_to_depac(i,j,ilu_ice)                =       0.0
!        veg_smtk_to_depac(i,j,ilu_ice)                =       veg_smtk(i,j,200+ilu_smtk_permanent_ice)   
        veg_smtk_to_depac(i,j,ilu_savanna)            =       0.0
        veg_smtk_to_depac(i,j,ilu_tropical_forest)    =       0.0
!        veg_smtk_to_depac(i,j,ilu_tropical_forest)    =       veg_smtk(i,j,200+ilu_smtk_tropical_forest)
        veg_smtk_to_depac(i,j,ilu_water_inland)       =       0.0
        veg_smtk_to_depac(i,j,ilu_mediterrean_scrub)  =       0.0
        veg_smtk_to_depac(i,j,ilu_semi_natural_veg)   =       0.0
        if ( with_ozone_specials ) then
          veg_smtk_to_depac(i,j,ilu_wheat)              =       0.0
          veg_smtk_to_depac(i,j,ilu_beech)              =       0.0
          veg_smtk_to_depac(i,j,ilu_spruce)             =       0.0
        end if

        ! check 
      	if (abs(sum(lu_fracs(i,j,:))-1.0) > 0.001) then
          print *, 'lu_depac = ', lu_fracs(i, j, :)
          print *, 'abs(sum(lu_depac(i,j,:))-1.0) = ', abs(sum(lu_fracs(i,j,:))-1.0)
          print *, 'landuse fractions do not count up to 1 for cell (', i, j, ')'
          stop
        end if
      
      end do
    end do
    
    ! ok
    status = 0
           
  end subroutine Translate_LU_to_Depac 


end module LE_Landuse_Smiatek
