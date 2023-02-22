!###############################################################################
!
! LE_Emis_SeaSalt - LOTOS-EUROS seasalt emission routines
!
! History
!   2010-10, Astrid Manders, TNO
!     Extension and update of old sea salt routine.
!     Implemented Monahan and Martensson source functions.
!     and reduced source strength in Baltic sea.
!     Originally: 4 sea salt classes.
!     Adapted for standard fine and coarse sea salt class.
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

module LE_Emis_SeaSalt

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  LE_Emis_SeaSalt_Init, LE_Emis_SeaSalt_Done
  public  ::  LE_Emis_SeaSalt_Add


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_SeaSalt'

  !dry particle density of NaCl in ug/um3
  real, parameter       ::  dens_NaCl = 2.17e-6


  ! --- local ------------------------------------

  ! bound_ss: boundaries of sea salt particle classes. Number of classes nclass
  ! must correspond with sea salt class definitition in dims.f90
  ! rp=particle radius, Dp=particle diameter (80%rh), in um
  ! rdr=dry particle radius, Ddr=dry particle diameter, in um
  ! wind velocity at 10m 
  ! boundaries of particle classes: diameter in um at 80%rh
  ! nclass particle size classes
  integer, parameter          ::  nclass=5
  real, dimension(nclass,2)   ::  bound_ss
  real, dimension(nclass)     ::  Smon_Dp, Smar_Dp
  integer                     ::  iclass

  ! Parameters specific for Martensson
  integer, parameter                ::  Martstep=50
  real, dimension (nclass,Martstep) ::  Ak, Bk, Ddr
  real, dimension (nclass)          ::  Ddstep
  real, parameter                   ::  Martbound=1.00
  

contains


  ! ====================================================================
  
  
  !Initialisation, sets boundaries of particle size classes 
  !calculates Monahan emission per particle size class
  !generates Martensson lookup tables for Ak and Bk

  subroutine LE_Emis_SeaSalt_Init( status )

    use LE_Data, only : LE_Data_Enable

    ! --- in/out ---------------------------------
    
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Emis_SeaSalt_Init'
    
    ! --- local ------------------------------------------
    
    real    :: Ddup, Ddlow,  rhconv
    real    :: emismon
    integer :: istep

   

    ! --- begin ------------------------------------------

    ! enable data:
    call LE_Data_Enable( 'sstk', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'wspd_surf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lon', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)

    rhconv=2.0 !conversion for the particle  diameter at 80% rh to dry

    bound_ss(1,1) = 0.14
    bound_ss(1,2) = 1.00
    bound_ss(2,1) = 1.00
    bound_ss(2,2) = 2.50 
    bound_ss(3,1) = 2.50
    bound_ss(3,2) = 4.00 
    bound_ss(4,1) = 4.00
    bound_ss(4,2) = 7.0
    bound_ss(5,1) = 7.00
    bound_ss(5,2) = 10.0

    !initialisation of emission arrays and construction of Martensson look-up table
    !and Monahan function
    do iclass = 1,nclass
      Smon_Dp(iclass)=0.
      Smar_Dp(iclass)=0.

      if (bound_ss(iclass,1).lt.0.14) bound_ss(iclass,1) = 0.14
      if (bound_ss(iclass,2).lt.0.14) write(*,*) 'WARNING: NACL emission not valid for lowest size class'


      if (bound_ss(iclass,2).le.martbound) then 
        !Initialisation of Martensson lookup table for Ak, Bk, Ddr
        !Computation of Martensson flux through calc_seasalt, depends on Tsea
        Ddlow=bound_ss(iclass,1)/rhconv
        Ddup=bound_ss(iclass,2)/rhconv
        Ddstep(iclass)=(log10(Ddup)-log10(Ddlow))/Martstep
        do istep=1,Martstep
          Ddr(iclass,istep)=10**(log10(Ddlow)+(istep-1)*Ddstep(iclass))
        enddo
        call akbkmartensson
      else if (bound_ss(iclass,2).gt.martbound) then
        !Monahan
        call emismonahan(emismon)
        smon_Dp(iclass) =  emismon
        !write(*,*) smon_Dp(iclass)
        !write(*,*) 'per hour:', 3600*smon_Dp(iclass)
        !write(*,*) 'per hour per percent sea:', 3600/100*smon_Dp(iclass)
      endif
    enddo

    ! ok
    status = 0
    
  end subroutine LE_Emis_SeaSalt_Init
  
  
  ! ***
  
  
  subroutine LE_Emis_SeaSalt_Done( status )
  
    use GO, only : LUT_Done
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_SeaSalt_Done'
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_SeaSalt_Done


  ! ***
  

  ! apply the sea salt emissions dependency on windspeed, e.g. multiply with U10**3.41

  subroutine LE_Emis_SeaSalt_Add( emis, status )

    use Binas        , only : pi
    use dims         , only : nx, ny, nz
    use dims         , only : runF, outF 
    use Indices      , only : nspec 
    use Indices      , only : n_seasalt, ispecs_seasalt
    use Indices      , only : ispec_Na_f,  ispec_Na_ff, ispec_Na_ccc, ispec_Na_cc, ispec_Na_c
    use Indices      , only : n_basecation, ispecs_basecation, n_seasalt
    use Indices      , only : ispec_Ca_f, ispec_Ca_c
    use Indices      , only : ispec_Mg_f, ispec_Mg_c
    use Indices      , only : ispec_K_f , ispec_K_c
    use LE_Data, only : LE_Data_GetPointer
    use LE_Landuse_Data, only : lu_fracs, ilu_water_sea
  
    ! --- in/out ---------------------------------
    
    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_SeaSalt_Add'
    
    ! Baltic sea, lower salinity
    real   , parameter :: baltwest=14.0
    real   , parameter :: balteast=32.0
    real   , parameter :: baltnorth=66.0
    real   , parameter :: baltsouth=54.0

    ! IJsselmeer, Wadden Sea, no sea salt emissions
    ! ijsselmeer/waddenzee lon-lat domein
    real   , parameter :: Ywest1=4.8 
    real   , parameter :: Yeast1=6.8  
    real   , parameter :: Ynorth1=53.15
    real   , parameter :: Ysouth1=52.0 

    real   , parameter :: Ywest2=5.0 
    real   , parameter :: Yeast2=7.0  
    real   , parameter :: Ynorth2=53.35
    real   , parameter :: Ysouth2=53.15 

    real   , parameter :: Ywest3=5.3 
    real   , parameter :: Yeast3=7.0  
    real   , parameter :: Ynorth3=53.48
    real   , parameter :: Ysouth3=53.35 

    !mass percentages of cations for sea salt:
    real   , parameter :: nafrac=0.306
    real   , parameter :: mgfrac=0.037
    real   , parameter :: cafrac=0.012
    real   , parameter :: kfrac =0.011
    
    ! conversion from ug/m3 to ug/cm3 and ug/cm3 to molec/cm3
    real, parameter:: convcm=1.e-6
    real, parameter:: kgtoug=1.e9
    
    ! --- local ----------------------------------

    real      ::  emismart
    integer   ::  ix,iy               ! grid cell index ix, iy
    real      ::  u_fact, conv, Tss   ! fraction sea and unit conversion sec to min
    real      ::  xcoor, ycoor
    integer   ::  itr, ispec
    real      ::  delta_emis

    real      ::  radavmass
    real      ::  volavmass

    real, pointer   ::  sst  (:,:,:)   ! (nlon,nlat,1)
    real, pointer   ::  wsurf(:,:,:)   ! (nlon,nlat,1)
    real, pointer   ::  area(:,:,:)   ! (lon,lat,1)
    real, pointer   ::  lons(:,:,:)   ! (lon,lat,1)
    real, pointer   ::  lats(:,:,:)   ! (lon,lat,1)


    ! --- begin ----------------------------------

    call LE_Data_GetPointer( 'sstk' , sst  , status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'wspd_surf', wsurf, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lon', lons, status, check_units ='degrees_east' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units ='degrees_north' )
    IF_NOTOK_RETURN(status=1)
    

    ! info ...
    write (gol,'("<constructing sea salt emissions>")'); call goPr

    do ix = 1, nx
      do iy = 1, ny
        if (lu_fracs(ix,iy,ilu_water_sea) > 0.0) then

          if (wsurf(ix,iy,1)<12.5) then
            u_fact = wsurf(ix,iy,1)**3.41  !wind factor
          else
            u_fact=12.5**3.41  !no further increase of whitecap beyond 12.5 m/s
          end if

          ! loop over classes:
          do iclass=1,nclass
	          ! Calculate Tsea-dependent Martensson emission. 
            ! Monahan flux already calculated in init-seasalt.
            if (bound_ss(iclass,2).le.Martbound) then
	            ! Martensson	  
              Tss = sst(ix,iy,1) !sea surface temperature
              ! dry particle volume flux   
              call emismartensson( emismart, Tss )
              ! store:
              Smar_Dp(iclass) = emismart
            end if  
          end do

          !lotos grid size in Km, translate to m, multiply by 60 to get per minute
          conv =lu_fracs(ix,iy,ilu_water_sea)*60 *area(ix,iy,1)

          !Check for points with lower sea salt emissions:                            
          !Baltic Sea, Lake IJssel, Wadden Sea
          !check if point in Baltic Sea, if yes then emission factor 10 lower
          !crude approximation for less emission due to low salinity
          !Martensson et al 2003 find factor 10 lower number of aerosol when
          !salinity lower than 10 psu
          xcoor=lons(ix,iy,1)
          ycoor=lats(ix,iy,1)
          if(baltwest.lt.xcoor.and.xcoor.lt.balteast.and. &
             baltsouth.lt.ycoor.and.ycoor.lt.baltnorth) then
             conv=conv/10.
             !write(*,*) 'ix,iy,xcoor,ycoor', ix,iy,xcoor,ycoor
          elseif(Ywest1.lt.xcoor.and.xcoor.lt.Yeast1.and. &
             Ysouth1.lt.ycoor.and.ycoor.lt.Ynorth1) then
             conv=0.0     
             !print *,' ijsselmeer within domain ,ix,iy,xcoor,ycoor :',ix,iy,xcoor,ycoor,seafraction(ix,iy)
          elseif(Ywest2.lt.xcoor.and.xcoor.lt.Yeast2.and. &
             Ysouth2.lt.ycoor.and.ycoor.lt.Ynorth2) then
             conv=0.0     
             !print *,' ijsselmeer within domain ,ix,iy,xcoor,ycoor :',ix,iy,xcoor,ycoor,seafraction(ix,iy)
          elseif(Ywest3.lt.xcoor.and.xcoor.lt.Yeast3.and. &
             Ysouth3.lt.ycoor.and.ycoor.lt.Ynorth3) then
             conv=0.0     
            ! print *,' ijsselmeer within domain ,ix,iy,xcoor,ycoor :',ix,iy,xcoor,ycoor,seafraction(ix,iy)
          end if

          !~ fine mode increment:
          delta_emis = u_fact*(smar_Dp(1)+smon_Dp(1))*conv 
          ! tracer enabled ?
          
          ispec = ispec_Na_ff
          if ( ispec > 0 ) then
            emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + delta_emis*Nafrac
            !print *, 'emis_ispec', ispec, delta_emis
          end if
          !~ fine mode increment:
          delta_emis = u_fact*(smar_Dp(2)+smon_Dp(2))*conv 
          ! tracer enabled ?
          ispec = ispec_Na_f
          if ( ispec > 0 ) then
            emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + delta_emis*Nafrac
          end if

          !~ coarse mode increment:
          delta_emis = u_fact*(smar_Dp(3)+smon_Dp(3))*conv
          ! tracer enabled ?
          ispec = ispec_Na_ccc
          if ( ispec > 0 ) then
            ! add:
            emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + delta_emis*Nafrac
            !print *, 'emis_ispec', ispec, delta_emis
          end if
          
          delta_emis = u_fact*(smar_Dp(4)+smon_Dp(4))*conv
          ! tracer enabled ?
          ispec = ispec_Na_cc
          if ( ispec > 0 ) then
            ! add:
            emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + delta_emis*Nafrac
            !print *, 'emis_ispec', ispec, delta_emis
          end if
          
          
          
          delta_emis = u_fact*(smar_Dp(5)+smon_Dp(5))*conv
          
          ! tracer enabled ?
          ispec = ispec_Na_c
          if ( ispec > 0 ) then
            ! add:
            emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + delta_emis*Nafrac
          end if

          ! add sea salt contribution, conversion to base cations;
          ! loop over basecation tracers:
          do itr = 1, n_basecation
            ! tracer index:
            ispec = ispecs_basecation(itr)
            ! switch:
            select case ( ispec )
              case(ispec_Na_ff)  ! already done as part of seasalt
              case(ispec_Na_f)  ! already done as part of seasalt
              case(ispec_Na_ccc)  ! already done as part of seasalt
              case(ispec_Na_cc)  ! already done as part of seasalt
              case(ispec_Na_c)  ! already done as part of seasalt
              
              case(ispec_Ca_f); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(1)+smon_Dp(1))*conv*Cafrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(2)+smon_Dp(2))*conv*Cafrac
              case(ispec_Ca_c); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(3)+smon_Dp(3))*conv*Cafrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(4)+smon_Dp(4))*conv*Cafrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(5)+smon_Dp(5))*conv*Cafrac
              case(ispec_Mg_f); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(1)+smon_Dp(1))*conv*Mgfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(2)+smon_Dp(2))*conv*Mgfrac
              case(ispec_Mg_c); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(3)+smon_Dp(3))*conv*Mgfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(4)+smon_Dp(4))*conv*Mgfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(5)+smon_Dp(5))*conv*Mgfrac
                                
              case(ispec_K_f ); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(1)+smon_Dp(1))*conv*Kfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(2)+smon_Dp(2))*conv*Kfrac
              case(ispec_K_c ); emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(3)+smon_Dp(3))*conv*Kfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(4)+smon_Dp(4))*conv*Kfrac
                                emis(ix,iy,1,ispec) = emis(ix,iy,1,ispec) + u_fact*(smar_Dp(5)+smon_Dp(5))*conv*Kfrac
              case default
                write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                TRACEBACK; status=1; return
            end select
          end do

        endif
        
      end do   ! i
    end do    ! j
    
    ! ok
    status = 0

  end subroutine LE_Emis_SeaSalt_Add

  
  ! ***
  

  ! emission function holds for 80% RH
  ! at 80 % RH  r_dry = 0.5*r_wet, emission is assumed to be wet, 
  ! but for mass cqalc o NaCl dry mass  is taken: thus masscor = 1/6 *Pi*rho*  r**2

  ! do the integration of the Monahan sea salt flux without wind factor
  ! per unit area per second
  ! 1.373 * r**-3 * (1+0.057r**1.05)*10**(1.19exp(-B**2)) *massfactor
  !	B = (0.380 - log10(r))/0.650

  subroutine emismonahan( emismon )
  
    use Binas, only : pi

    real, intent(out)      ::  emismon

    integer       ::  istep, nsteps
    real          ::  stepsize, emism, rp, B

    !determine the stepsize and number of steps:
    stepsize=0.01
    nsteps = nint((bound_ss(iclass,2) - bound_ss(iclass,1))/stepsize)

    emismon = 0.0
    do istep = 1,nsteps
      rp = 0.5* (bound_ss(iclass,1) + (istep-0.5)*stepsize)
      B = (0.380 - log10(rp))/0.65
      B = 1.19*exp(-B**2)
      emisM= 1.373*(1.0+0.057*rp**1.05)*10**B*pi*(1./6.)*dens_nacl 
      emismon = emismon+ emisM * stepsize/2 ! stepsize in diameter, thus for radius divide by 2
    enddo

  end subroutine 


  ! ***
  
  
  ! Calculate seasalt emission for smallest particles, according
  ! to Martensson et al. 2003 eqn 6
  ! variant met andere integratie

  ! Dry particle diameter Ddr in micrometer
  ! Tsea: sea water temperature in K, from ECMWF
  ! Ak, Bk polynomials (martenssonpolynomial)

  subroutine emisMartensson( emismart, Tsea )
  
    use Binas, only : pi

    real, intent(out)   ::  emismart
    real, intent(in)    ::  Tsea     

    integer   ::  i

    ! Martensson emission depends on sea water temperature
    ! Min in feb, max in aug, valid for North Sea
    ! Could be improved by implementing region-dependent Tsea 
    ! or even more accurately by ECMWF-fields

    emismart=0.0

    ! integration scheme: mean of forward and backward integration      
    emismart = emismart + (Ak(iclass,       1)*Tsea+Bk(iclass,       1))*(Ddr(iclass,       1)**3)     
    emismart = emismart + (Ak(iclass,Martstep)*Tsea+Bk(iclass,Martstep))*(Ddr(iclass,Martstep)**3)
    emismart=emismart*0.5
    do i=2,Martstep-1                 
      emismart=emismart+(Ak(iclass,i)*Tsea+Bk(iclass,i))*(Ddr(iclass,i)**3)
    enddo

    emismart=emismart*Ddstep(iclass)*dens_nacl*(1./6.)*pi*3.84e-6

  end subroutine



  ! ***


  ! Calculates Martensson polynomials Ak and Bk
  ! based on coefficients cn and dn specified for different particle size intervals
  ! according to Martensson et al 2003.
  ! These size intervals do not correspond to the size 
  ! intervals defined in lotos-euros.

  subroutine akbkmartensson

    integer :: interval, j,k
    real :: ck(3,5)
    real :: dk(3,5)
    real :: sizebound(4)

    sizebound=(/0.02, 0.145, 0.419, 2.8/)
    ck(1,:)=(/-2.881e6,-3.003e13,-2.867e21,5.932e28,-2.576e35/)
    ck(2,:)=(/-6.743e6,1.183e14,-8.148e20,2.404e27,-2.452e33/)
    ck(3,:)=(/2.181e6,-4.165e12,3.132e18,-9.841e23,1.085e29/)
    dk(1,:)=(/7.609e8 ,1.829e16,6.791e23,-1.616e31,7.188e37/)
    dk(2,:)=(/2.279e9,-3.787e16,2.528e23,-7.310e29,7.368e35/)
    dk(3,:)=(/-5.800e8,1.105e15,-8.297e20,2.601e26,-2.859e31/)  



    do j=1,Martstep
      !selection of right interval in ck, dk     
      if(Ddr(iclass,j).lt.sizebound(2)) then
            interval=1
        elseif(sizebound(2).lt.Ddr(iclass,j).and.&
               Ddr(iclass,j).lt.sizebound(3)) then
            interval=2
        elseif(sizebound(3).lt.Ddr(iclass,j).and.&
               Ddr(iclass,j).lt.sizebound(4))then
             interval=3
        else
             write(*,*) 'invalid particle size for Martensson!!'
      endif

      !calculation of Ak and Bk.
      !coefficients for Dd in m, but Ddr here in um->factor 10e6
      Ak(iclass,j)=0.
      Bk(iclass,j)=0.

      do k=1,5
        Ak(iclass,j)=Ak(iclass,j)&
                    +ck(interval,k)*(Ddr(iclass,j)/1e6)**(k-1)
        Bk(iclass,j)=Bk(iclass,j)&
                    +dk(interval,k)*(Ddr(iclass,j)/1e6)**(k-1)
      enddo       
    enddo

  end subroutine

end module LE_Emis_SeaSalt
