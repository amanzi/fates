module ATSFatesInterfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library/API with the CLM/ALM/ATS/etc model driver.  
   ! All connections between the two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see FATES memory 
   ! structures.
   ! The routines here, that call FATES library routines, will not pass any types defined
   ! by the driving land model (HLM).
   ! 
   ! either native type arrays (int,real,log, etc) or packed into ED boundary condition
   ! structures.
   !
   ! Note that CLM/ALM does use Shared Memory Parallelism (SMP), where processes such as 
   ! the update of state variables are forked.  However, IO is not assumed to be 
   ! threadsafe and therefore memory spaces reserved for IO must be continuous vectors,
   ! and moreover they must be pushed/pulled from history IO for each individual 
   ! bounds_proc memory space as a unit.
   !
   ! Therefore, the state variables in the alm_fates communicator is vectorized by
   ! threadcount, and the IO communication arrays are not.
   !
   !
   ! Conventions:
   ! keep line widths within 90 spaces
   ! HLM acronym = Host Land Model
   !
   ! -------------------------------------------------------------------------------------


   use, intrinsic :: ISO_C_BINDING 
  
   !  use ed_driver_interface, only: 
   
   ! Used CLM Modules
   ! use VegetationType    , only : veg_pp
   use shr_kind_mod      , only : r8 => shr_kind_r8
   ! use decompMod         , only : bounds_type
   ! use WaterStateType    , only : waterstate_type
   ! use WaterFluxType     , only : waterflux_type
   ! use CanopyStateType   , only : canopystate_type
   ! use TemperatureType   , only : temperature_type
   ! use EnergyFluxType    , only : energyflux_type

   ! use SoilStateType     , only : soilstate_type 
   ! use clm_varctl        , only : iulog
   ! use clm_varctl        , only : use_vertsoilc 
   ! use clm_varctl        , only : use_fates_spitfire
   ! use clm_varctl        , only : use_fates_planthydro
   ! use clm_varctl        , only : use_fates_ed_st3
   ! use clm_varctl        , only : use_fates_ed_prescribed_phys
   ! use clm_varctl        , only : use_fates_logging
   ! use clm_varctl        , only : use_fates_inventory_init
   ! use clm_varctl        , only : fates_inventory_ctrl_filename
   ! use clm_varcon        , only : tfrz
   ! use clm_varcon        , only : spval 
   ! use clm_varcon        , only : denice
   ! use clm_varcon        , only : ispval

   ! use clm_varpar        , only : natpft_size
   ! use clm_varpar        , only : numrad
   ! use clm_varpar        , only : ivis
   ! use clm_varpar        , only : inir
   ! use clm_varpar        , only : nlevgrnd
   ! use clm_varpar        , only : nlevdecomp
   ! use clm_varpar        , only : nlevdecomp_full
   ! use clm_varpar        , only : i_met_lit, i_cel_lit, i_lig_lit
   ! use PhotosynthesisType , only : photosyns_type
   ! use atm2lndType       , only : atm2lnd_type
   ! use SurfaceAlbedoType , only : surfalb_type
   ! use SolarAbsorbedType , only : solarabs_type
   ! use CNCarbonFluxType  , only : carbonflux_type
   ! use CNCarbonStateType , only : carbonstate_type
   ! use FrictionVelocityType , only : frictionvel_type
   ! use clm_time_manager  , only : is_restart
   ! use ncdio_pio         , only : file_desc_t, ncd_int, ncd_double
   ! use restUtilMod,        only : restartvar
   ! use clm_time_manager  , only : get_days_per_year, &
   !                                get_curr_date,     &
   !                                get_ref_date,      &
   !                                timemgr_datediff,  &
   !                                is_beg_curr_day,   &
   !                                get_step_size,     &
   !                                get_nstep
   ! use spmdMod           , only : masterproc
   ! use decompMod         , only : get_proc_bounds,   &
   !                                get_proc_clumps,   &
   !                                get_clump_bounds
   ! use GridCellType      , only : grc_pp
   ! use ColumnType        , only : col_pp
   ! use LandunitType      , only : lun_pp
   ! use landunit_varcon   , only : istsoil
   ! use abortutils        , only : endrun
   use FatesGlobals        , only : endrun => fates_endrun
   use shr_log_mod       , only : errMsg => shr_log_errMsg    
   ! use clm_varcon        , only : dzsoi_decomp
   ! use FuncPedotransferMod, only: get_ipedof
   

   ! Used FATES Modules
   use FatesInterfaceMod     , only : fates_interface_type
   use FatesInterfaceMod     , only : allocate_bcin
   use FatesInterfaceMod     , only : allocate_bcout
   use FatesInterfaceMod     , only : SetFatesTime
   use FatesInterfaceMod     , only : set_fates_ctrlparms
   use FatesInterfaceMod     , only : InitPARTEHGlobals
   
   use FatesHistoryInterfaceMod, only : fates_history_interface_type
   use FatesRestartInterfaceMod, only : fates_restart_interface_type

!   use ChecksBalancesMod     , only : SummarizeNetFluxes, FATES_BGC_Carbon_BalanceCheck
   use EDTypesMod            , only : ed_patch_type
   use FatesInterfaceMod     , only : hlm_numlevgrnd
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site
   use EDInitMod             , only : zero_site
   use EDInitMod             , only : init_site_vars
   use EDInitMod             , only : init_patches
   use EDInitMod             , only : set_site_properties
   use EDPftVarcon           , only : EDpftvarcon_inst
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs, ED_Norman_Radiation
   use EDBtranMod            , only : btran_ed, &
                                      get_active_suction_layers
   use EDCanopyStructureMod  , only : canopy_summarization, update_hlm_dynamics
   use FatesPlantRespPhotosynthMod, only : FatesPlantRespPhotosynthDrive
   use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
   use EDPhysiologyMod       , only : FluxIntoLitterPools
   use FatesPlantHydraulicsMod, only : hydraulics_drive
   use FatesPlantHydraulicsMod, only : HydrSiteColdStart
   use FatesPlantHydraulicsMod, only : InitHydrSites
   use FatesPlantHydraulicsMod, only : UpdateH2OVeg
   use FatesInterfaceMod      , only : bc_in_type, bc_out_type

   implicit none
   
   type, public :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:) 

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

   end type f2hmap_type

   type, bind (C) :: SiteInfo
      integer (C_INT) :: nlevbed
      integer (C_INT) :: nlevdecomp
      integer (C_INT) :: patchno
      integer (C_INT) :: altmax_lastyear_indx_col
      REAL(C_DOUBLE) :: temp_veg24_patch
      REAL(C_DOUBLE) :: latdeg, londeg      
   end type SiteInfo


   type, bind (C) :: PhotoSynthesisInput
      REAL(C_DOUBLE) :: dayl_factor      ! scalar (0-1) for daylength
      REAL(C_DOUBLE) :: esat_tv     ! saturation vapor pressure at t_veg (Pa)
      REAL(C_DOUBLE) :: eair        ! vapor pressure of canopy air (Pa)
      REAL(C_DOUBLE) :: oair        ! Atmospheric O2 partial pressure (Pa)
      REAL(C_DOUBLE) :: cair        ! Atmospheric CO2 partial pressure (Pa)
      REAL(C_DOUBLE) :: rb          ! boundary layer resistance (s/m)
      REAL(C_DOUBLE) :: t_veg       ! vegetation temperature (Kelvin)
      REAL(C_DOUBLE) :: tgcm        ! air temperature at agcm reference height (Kelvin)
   end type PhotoSynthesisInput
   
   type, bind (C) :: TimeInput
   ! -------------------------------------------------------------------------------------
   integer (C_INT)  :: current_year    ! Current year
   integer (C_INT)  :: current_month   ! month of year
   integer (C_INT)  :: current_day     ! day of month
   integer (C_INT)  :: current_tod     ! time of day (seconds past 0Z)
   integer (C_INT)  :: current_date    ! YYYYMMDD
   integer (C_INT)  :: reference_date  ! YYYYMMDD
   REAL(C_DOUBLE)   :: model_day       ! elapsed days between current date and ref
   integer (C_INT)  :: day_of_year     ! The integer day of the year
   integer (C_INT)  :: days_per_year   ! The HLM controls time, some HLMs may 
                                           ! include a leap
   end type TimeInput

   type, public,  bind (C) :: bounds_type
      integer (C_INT) :: begg, endg       ! beginning and ending gridcell index
      integer (C_INT):: begl, endl       ! beginning and ending landunit index
      integer (C_INT):: begc, endc       ! beginning and ending column index
      integer (C_INT):: begp, endp       ! beginning and ending pft index
      integer (C_INT):: lbj, ubj
      integer (C_INT):: begCohort, endCohort ! beginning and ending cohort indices

      integer (C_INT):: level            ! whether defined on the proc or clump level
      integer (C_INT):: clump_index      ! if defined on the clump level, this gives the clump index
   end type bounds_type

!   type, bind (C) :: CanopyState_type

     ! integer(C_INT) :: frac_veg_nosno_patch     (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 
     ! integer(C_INT) :: frac_veg_nosno_alb_patch (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 

     ! real(r8) , pointer :: tlai_patch               (:)   ! patch canopy one-sided leaf area index, no burying by snow
     ! real(r8) , pointer :: tsai_patch               (:)   ! patch canopy one-sided stem area index, no burying by snow
     ! real(r8) , pointer :: elai_patch               (:)   ! patch canopy one-sided leaf area index with burying by snow
     ! real(r8) , pointer :: esai_patch               (:)   ! patch canopy one-sided stem area index with burying by snow
     ! real(r8) , pointer :: elai_p_patch             (:)   ! patch canopy one-sided leaf area index with burying by snow average over timestep 
     ! real(r8) , pointer :: laisun_patch             (:)   ! patch patch sunlit projected leaf area index  
     ! real(r8) , pointer :: laisha_patch             (:)   ! patch patch shaded projected leaf area index  
     ! real(r8) , pointer :: laisun_z_patch           (:,:) ! patch patch sunlit leaf area for canopy layer 
     ! real(r8) , pointer :: laisha_z_patch           (:,:) ! patch patch shaded leaf area for canopy layer 
     ! real(r8) , pointer :: mlaidiff_patch           (:)   ! patch difference between lai month one and month two (for dry deposition of chemical tracers)
     ! real(r8) , pointer :: annlai_patch             (:,:) ! patch 12 months of monthly lai from input data set (for dry deposition of chemical tracers) 
     ! real(r8) , pointer :: htop_patch               (:)   ! patch canopy top (m)
     ! real(r8) , pointer :: hbot_patch               (:)   ! patch canopy bottom (m)
     ! real(r8) , pointer :: displa_patch             (:)   ! patch displacement height (m)
     ! real(r8) , pointer :: fsun_patch               (:)   ! patch sunlit fraction of canopy         
     ! real(r8) , pointer :: fsun24_patch             (:)   ! patch 24hr average of sunlit fraction of canopy 
     ! real(r8) , pointer :: fsun240_patch            (:)   ! patch 240hr average of sunlit fraction of canopy

     ! real(r8) , pointer :: alt_col                  (:)   ! col current depth of thaw 
     ! integer  , pointer :: alt_indx_col             (:)   ! col current depth of thaw 
     ! real(r8) , pointer :: altmax_col               (:)   ! col maximum annual depth of thaw 
     ! real(r8) , pointer :: altmax_lastyear_col      (:)   ! col prior year maximum annual depth of thaw 
     ! integer  , pointer :: altmax_indx_col          (:)   ! col maximum annual depth of thaw 
     ! integer  , pointer :: altmax_lastyear_indx_col (:)   ! col prior year maximum annual depth of thaw 
     
     ! real(r8) , pointer :: dewmx_patch              (:)   ! patch maximum allowed dew [mm] 

     ! real(r8) , pointer :: dleaf_patch              (:)   ! patch characteristic leaf width (diameter) [m]
     !                                                      ! for non-ED/FATES this is the same as pftcon%dleaf()
     ! real(r8),  pointer :: lbl_rsc_h2o_patch        (:)   ! laminar boundary layer resistance for water over dry leaf (s/m)

!  end type CanopyState_type
  
   ! hlm_bounds_to_fates_bounds is not currently called outside the interface.
   ! Although there may be good reasons to, I privatized it so that the next
   ! developer will at least question its usage (RGK)
   !private :: hlm_bounds_to_fates_bounds

   logical :: DEBUG  = .false.
   character(len=*), parameter, private :: sourcefile = &
        __FILE__

   
   type(fates_interface_type), allocatable :: fates (:)
   type(f2hmap_type), allocatable  :: f2hmap(:)
   ! fates_hist is the interface class for the history output
   type(fates_history_interface_type) :: fates_hist  
   ! fates_restart is the inteface calss for restarting the model
   type(fates_restart_interface_type) :: fates_restart   
   integer                            :: iulog   
  
 contains

  ! ====================================================================================

   subroutine init_ats_fates(num_sites, site_ats) BIND(C)
     
      ! ---------------------------------------------------------------------------------
      !
      ! sites is the root of the ED state hierarchy (instantaneous info on 
      ! the state of the ecosystem).  As such, it governs the connection points between
      ! the host (which also dictates its allocation) and its patch structures.
      !
      ! sites may associate with different scales in different models. In
      ! ATS, it is being designed to relate to column scale.
      !
      ! This global may become relegated to this module. 
      !
      ! Note: ATS/ALM currently wants sites to be allocated even if ed
      ! is not turned on
      ! ---------------------------------------------------------------------------------
     
      use FatesInterfaceMod, only : FatesInterfaceInit 
      use FatesInterfaceMod, only : FatesReportParameters
      use FatesParameterDerivedMod, only : param_derived
      use FatesInterfaceMod, only : numpft_fates => numpft



      implicit none
      
      !Input Arguments
      ! class(hlm_fates_interface_type), intent(inout) :: this
      ! type(bounds_type),intent(in)                   :: bounds_proc
      integer (C_INT)                            :: num_sites
      type(SiteInfo),intent(in)                   :: site_ats(num_sites)

      !temprorary local variables
      real(r8) , parameter                           :: spval = 1.e36_r8  ! special value for real data
!      real(r8)                                       :: latdeg, londeg
      integer                                        :: nlevgrnd
      integer                                        :: max_patch_per_site
      integer                                        :: soilwater_ipedof
      character(len=256)                             :: fates_inventory_ctrl_filename = ''
      integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
      integer, parameter :: ngases      =   3     ! CH4, O2, & CO2
      integer, parameter :: nlevcan     =   1     ! number of leaf layers in canopy layer
      integer, parameter :: numwat      =   5     ! number of water types (soil, ice, 2 lakes, wetland)
      integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
      integer, parameter :: ivis        =   1     ! index for visible band
      integer, parameter :: inir        =   2     ! index for near-infrared band
      integer, parameter :: numsolar    =   2     ! number of solar type bands: direct, diffuse
      integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
      integer, parameter :: dst_src_nbr =   3     ! number of size distns in src soil (BGC only)
      integer, parameter :: sz_nbr      = 200     ! number of sub-grid bins in large bin of dust size distribution (BGC only)
      integer, parameter :: mxpft       =  24     ! maximum number of PFT's for any mode;
      ! FIX(RF,032414) might we set some of these automatically from reading pft-physiology?
      integer, parameter :: numveg      =  16     ! number of veg types (without specific crop)
      integer, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang


      

      real(r8), allocatable                          :: dzsoi_decomp(:)
      real(r8), allocatable                          :: zi(:,:), dz(:,:), z(:,:)
      integer, allocatable                           :: nlevbed(:)

      ! local variables
      integer                                        :: nclumps   ! Number of threads
      logical                                        :: verbose_output
      integer                                        :: pass_masterproc
      integer                                        :: pass_vertsoilc
      integer                                        :: pass_spitfire     
      integer                                        :: pass_ed_st3
      integer                                        :: pass_logging
      integer                                        :: pass_ed_prescribed_phys
      integer                                        :: pass_planthydro
      integer                                        :: pass_inventory_init
      integer                                        :: pass_is_restart
      integer                                        :: fates_parteh_mode
      integer                                        :: nc        ! thread index
      integer                                        :: s         ! FATES site index
      integer                                        :: c         ! HLM column index
      integer                                        :: l         ! HLM LU index
      integer                                        :: g         ! HLM grid index
      integer                                        :: pi,pf
      integer, allocatable                           :: collist (:)
      ! type(bounds_type)                              :: bounds_clump
      integer                                        :: nmaxcol
      integer                                        :: ndecomp


     
      !write(*, *) "Enter Init fates"     
      DEBUG = .true.
            
      ! Initialize the FATES communicators with the HLM
      ! This involves to stages
      ! 1) allocate the vectors
      ! 2) add the history variables defined in clm_inst to the history machinery
      call param_derived%Init( numpft_fates )

      iulog = 10
      verbose_output = .false.
      pass_is_restart = 0
      pass_vertsoilc = 0
      pass_spitfire = 0
      pass_ed_st3 = 0
      pass_logging = 0
      pass_ed_prescribed_phys = 0
      fates_parteh_mode = 1
      
      pass_planthydro = 0
      pass_inventory_init = 0
      pass_masterproc = 0
      fates_inventory_ctrl_filename=""
     
      call FatesInterfaceInit(iulog, verbose_output)

      nclumps = 1 ! get_proc_clumps()
      allocate(fates(nclumps))
      allocate(f2hmap(nclumps))

      ! ! ---------------------------------------------------------------------------------
      ! ! Send dimensions and other model controling parameters to FATES.  These
      ! ! are obviously only those parameters that are dictated by the host
      ! ! ---------------------------------------------------------------------------------
      
      ! Force FATES parameters that are recieve type, to the unset value
      call set_fates_ctrlparms('flush_to_unset')
      
      ! Send parameters individually

      nlevgrnd = 1
      soilwater_ipedof = 0
      max_patch_per_site =site_ats(1)%patchno
      
      call set_fates_ctrlparms('num_sw_bbands',ival=numrad)
      call set_fates_ctrlparms('vis_sw_index',ival=ivis)
      call set_fates_ctrlparms('nir_sw_index',ival=inir)
      
      call set_fates_ctrlparms('num_lev_ground',ival=nlevgrnd)
      call set_fates_ctrlparms('hlm_name',cval='ATS')
      call set_fates_ctrlparms('hio_ignore_val',rval=spval)
      call set_fates_ctrlparms('soilwater_ipedof',ival=soilwater_ipedof)
      call set_fates_ctrlparms('max_patch_per_site',ival=max_patch_per_site) ! FATES IGNORES
                                                                          ! AND DOESNT TOUCH
                                                                          ! THE BARE SOIL PATCH      
      call set_fates_ctrlparms('parteh_mode', ival=1)

      call set_fates_ctrlparms('is_restart',ival=pass_is_restart)
      call set_fates_ctrlparms('use_vertsoilc',ival=pass_vertsoilc)     
      call set_fates_ctrlparms('use_spitfire',ival=pass_spitfire)
      call set_fates_ctrlparms('use_ed_st3',ival=pass_ed_st3)
      call set_fates_ctrlparms('use_logging',ival=pass_logging)
      call set_fates_ctrlparms('use_ed_prescribed_phys',ival=pass_ed_prescribed_phys)
      call set_fates_ctrlparms('use_planthydro',ival=pass_planthydro)
      call set_fates_ctrlparms('use_inventory_init',ival=pass_inventory_init)
      call set_fates_ctrlparms('inventory_ctrl_file',cval=fates_inventory_ctrl_filename)
      call set_fates_ctrlparms('masterproc',ival=pass_masterproc)
      ! ! Check through FATES parameters to see if all have been set
      call set_fates_ctrlparms('check_allset')

      if(DEBUG)then
         write(iulog,*) 'ats_fates%init():  allocating for ',nclumps,' threads'
      end if

      
      ! nclumps = get_proc_clumps()

      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,nmaxcol,s,c,l,g,collist,pi,pf)
      do nc = 1, nclumps
         
         if(DEBUG)then
            write(iulog,*) 'alm_fates%init(): thread',nc,': allocated ', num_sites,' sites'
         end if

         ! Allocate vectors that match FATES sites with HLM columns
         ! RGK: Sites and fcolumns are forced as args during clm_driv() as of 6/4/2016
         ! We may have to give these a dummy allocation of 1, which should
         ! not be a problem since we always iterate on nsites.
        
         ! Set the number of FATES sites
         fates(nc)%nsites = num_sites

         ! Allocate the FATES sites
         allocate ( fates(nc)%sites( fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (in)
         allocate( fates(nc)%bc_in( fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (out)
         allocate(fates(nc)%bc_out(fates(nc)%nsites))         
       
         ! Allocate and Initialize the Boundary Condition Arrays
         ! These are staticaly allocated at maximums, so
         ! No information about the patch or cohort structure is needed at this step
         
         do s = 1,  fates(nc)%nsites

            if (pass_vertsoilc) then
               ndecomp = site_ats(s)%nlevbed
            else
               ndecomp = 1
            end if


            call allocate_bcin( fates(nc)%bc_in(s), site_ats(s)%nlevbed, ndecomp)
            call allocate_bcout( fates(nc)%bc_out(s), site_ats(s)%nlevbed, ndecomp)
            
            call  fates(nc)%zero_bcs(s)

      !       ! Pass any grid-cell derived attributes to the site
      !       ! ---------------------------------------------------------------------------

            fates(nc)%sites(s)%lat = site_ats(s)%latdeg
            fates(nc)%sites(s)%lon = site_ats(s)%londeg

         end do


         ! Initialize site-level static quantities dictated by the HLM
         ! currently ground layering depth

      !   call  init_soil_depths(nc, zi, dz, z, dzsoi_decomp)
         
      !    if (use_fates_planthydro) then
      !       call InitHydrSites( fates(nc)%sites, fates(nc)%bc_in)
      !    end if

         if(  fates(nc)%nsites == 0 ) then
            write(iulog,*) 'Clump ',nc,' had no valid FATES sites'
            write(iulog,*) 'This will likely cause problems until code is improved'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


      end do
      !$OMP END PARALLEL DO

      !write(*, *) "Exit Init fates"

    end subroutine init_ats_fates


    subroutine get_nlevsclass(nlevel_class) BIND(C)

      use FatesInterfaceMod, only : nlevsclass

      integer(C_INT),intent(out)            :: nlevel_class

      nlevel_class = nlevsclass

    end subroutine get_nlevsclass


    subroutine init_soil_depths(nc, site_id, site_info, zi, dz, z, dzsoi_decomp) BIND(C)

      use FatesInterfaceMod, only : FatesReportParameters
      
      ! Input Arguments
      integer(C_INT),intent(in)            :: nc   ! Clump
      integer(C_INT),intent(in)            :: site_id   ! Site_id
      type(SiteInfo),intent(in)            :: site_info
      real(C_DOUBLE), intent(in)           :: dzsoi_decomp(site_info%nlevdecomp)
      real(C_DOUBLE), intent(in)           :: zi(site_info%nlevbed + 1), dz(site_info%nlevbed), z(site_info%nlevbed)

      ! Locals
      integer :: s  ! site index
      integer :: c  ! column index
      integer :: j  ! Depth index
      integer :: nlevsoil
      integer :: nlevdecomp
      logical :: masterproc

      !write(*, *) "Enter Init_soil_depth"

      if (site_id.le.0.or.site_id.gt.fates(nc)%nsites) then
         write(iulog,*) "Incorrect site_id is provides", site_id
         write(iulog,*) "Number of site is ", fates(nc)%nsites
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      nlevsoil =  fates(nc)%bc_in(site_id)%nlevsoil
      if (nlevsoil.gt.site_info%nlevbed) then
         write(iulog,*) 'Incorrect number of soil levels is provided ', site_info%nlevbed
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      nlevdecomp =  fates(nc)%bc_in(site_id)%nlevdecomp
      if (nlevdecomp.gt.site_info%nlevdecomp) then
         write(iulog,*) 'Incorrect number of decomposed soil levels is provided ',  site_info%nlevdecomp
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      
      fates(nc)%bc_in(site_id)%zi_sisl(0:nlevsoil)    = zi(0:nlevsoil)
      fates(nc)%bc_in(site_id)%dz_sisl(1:nlevsoil)    = dz(1:nlevsoil)
      fates(nc)%bc_in(site_id)%z_sisl(1:nlevsoil)     = z(1:nlevsoil)
      fates(nc)%bc_in(site_id)%dz_decomp_sisl(1:nlevdecomp) = dzsoi_decomp(1:nlevdecomp)
      
      do j=1,nlevsoil
             fates(nc)%bc_in(site_id)%decomp_id(j) = 1
      end do     

      !call init_history_io(bounds_proc)      
      ! Report Fates Parameters (debug flag in lower level routines)
      ! masterproc = .true.  !cx: need to be passed from the function
      ! call FatesReportParameters(masterproc)
      
      !write(*, *) "Exit Init_soil_depth fates"
      return
    end subroutine init_soil_depths


    subroutine dynamics_driv_per_site( &
         nc, site_id, site_info, time_input, dtime,&
         h2osoi_vol_col, temp_veg24_patch, prec24_patch, rh24_patch, wind24_patch) BIND(C)
    
!       ! This wrapper is called daily from ATS_driver for each site
!       ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
!       ! ed_driver is not a hlm_fates_inst_type procedure because we need an extra step 
!       ! to process array bounding information 
      
     implicit none
      ! Input Arguments
      integer(C_INT),intent(in)                 :: nc   ! Clump
      integer(C_INT),intent(in)                 :: site_id   ! site index
      type(SiteInfo),intent(in)                 :: site_info
      type(TimeInput),intent(in)                :: time_input
      real(C_DOUBLE), intent(in)                :: dtime
      real(C_DOUBLE), intent(in)                :: h2osoi_vol_col(site_info%nlevbed)
      real(C_DOUBLE), intent(in)                :: temp_veg24_patch(1)
      real(C_DOUBLE), intent(in)                :: prec24_patch(1)
      real(C_DOUBLE), intent(in)                :: rh24_patch(1)
      real(C_DOUBLE), intent(in)                :: wind24_patch(1) !site_info%patchno
      ! real(C_DOUBLE), intent(inout)             :: decomp_cpools_sourcesink_met(site_info%nlevdecomp)
      ! real(C_DOUBLE), intent(inout)             :: decomp_cpools_sourcesink_cel(site_info%nlevdecomp)
      ! real(C_DOUBLE), intent(inout)             :: decomp_cpools_sourcesink_lig(site_info%nlevdecomp)
      
      
!       ! !LOCAL VARIABLES:
       integer  :: ifp                      ! patch index


       integer  :: yr                       ! year (0, ...)
       integer  :: mon                      ! month (1, ..., 12)
       integer  :: day                      ! day of month (1, ..., 31)
       integer  :: sec                      ! seconds of the day
       integer  :: nlevsoil                 ! number of soil layers at the site
       integer  :: current_year             
       integer  :: current_month
       integer  :: current_day
       integer  :: current_tod
       integer  :: current_date
       integer  :: jan01_curr_year
       integer  :: reference_date
       integer  :: days_per_year
       real(r8) :: model_day
       real(r8) :: day_of_year
       !-----------------------------------------------------------------------

       ! ---------------------------------------------------------------------------------
       ! Part I.
       ! Prepare input boundary conditions for FATES dynamics
       ! Note that timing information is the same across all sites, this may
       ! seem redundant, but it is possible that we may have asynchronous site simulations
       ! one day.  The cost of holding site level boundary conditions is minimal
       ! and it keeps all the boundaries in one location
       ! ---------------------------------------------------------------------------------

!       days_per_year = get_days_per_year()
!       call get_curr_date(current_year,current_month,current_day,current_tod)
!       current_date = current_year*10000 + current_month*100 + current_day
!       jan01_curr_year = current_year*10000 + 100 + 1

!       call get_ref_date(yr, mon, day, sec)
        reference_date = yr*10000 + mon*100 + day

!       call timemgr_datediff(reference_date, sec, current_date, current_tod, model_day)

!       call timemgr_datediff(jan01_curr_year,0,current_date,sec,day_of_year)
       
        day_of_year = time_input%day_of_year
        call SetFatesTime(time_input%current_year, time_input%current_month, &
                         time_input%current_day, time_input%current_tod, &
                         time_input%current_date, time_input%reference_date, &
                         time_input%model_day, floor(day_of_year), &
                         time_input%days_per_year, 1.0_r8/dble(time_input%days_per_year))

      if (site_id.le.0.or.site_id.gt.fates(nc)%nsites) then
         write(iulog,*) "Incorrect site_id is provided", site_id
         write(iulog,*) "Number of site is ", fates(nc)%nsites
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if


      nlevsoil = fates(nc)%bc_in(site_id)%nlevsoil

      fates(nc)%bc_in(site_id)%h2o_liqvol_sl(1:nlevsoil)  = h2osoi_vol_col(1:nlevsoil) 
      fates(nc)%bc_in(site_id)%t_veg24_si = site_info%temp_veg24_patch

      fates(nc)%bc_in(site_id)%max_rooting_depth_index_col = &
               min(nlevsoil, site_info%altmax_lastyear_indx_col)     
      
      if ( fates(nc)%sites(site_id)%youngest_patch%patchno.gt.site_info%patchno ) then
         write(iulog,*) "Incorrect patchno is provided", site_info%patchno, "<", fates(nc)%sites(site_id)%youngest_patch%patchno
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      
      do ifp = 1, fates(nc)%sites(site_id)%youngest_patch%patchno
         
         fates(nc)%bc_in(site_id)%t_veg24_pa(ifp) = &
              temp_veg24_patch(1)

         fates(nc)%bc_in(site_id)%precip24_pa(ifp) = &
              prec24_patch(1)

         fates(nc)%bc_in(site_id)%relhumid24_pa(ifp) = &
              rh24_patch(1)

         fates(nc)%bc_in(site_id)%wind24_pa(ifp) = &
              wind24_patch(1)

      end do

         
!          if(use_fates_planthydro)then
!             this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil)  = soilstate_inst%hksat_col(c,1:nlevsoil)
!             this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = soilstate_inst%watsat_col(c,1:nlevsoil)
!             this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = spval !soilstate_inst%watres_col(c,1:nlevsoil)
!             this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = soilstate_inst%sucsat_col(c,1:nlevsoil)
!             this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)    = soilstate_inst%bsw_col(c,1:nlevsoil)
!             this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) =  waterstate_inst%h2osoi_liq_col(c,1:nlevsoil)
!          end if
         



!       ! ---------------------------------------------------------------------------------
!       ! Part II: Call the FATES model now that input boundary conditions have been
!       ! provided.
!       ! ---------------------------------------------------------------------------------
      


      call ed_ecosystem_dynamics(fates(nc)%sites(site_id),    &
           fates(nc)%bc_in(site_id))
            
      call ed_update_site(fates(nc)%sites(site_id), &
           fates(nc)%bc_in(site_id))
            

      
       ! call subroutine to aggregate ED litter output fluxes and 
       ! package them for handing across interface
        call FluxIntoLitterPools(fates(nc)%nsites, &
              fates(nc)%sites,  &
              fates(nc)%bc_in,  &
              fates(nc)%bc_out)

        ! ---------------------------------------------------------------------------------       
	! Part III: Process FATES output into the dimensions and structures that are part
        ! of the HLMs API.  (column, depth, and litter fractions)
        ! ---------------------------------------------------------------------------------
       ! call UpdateLitterFluxes_per_site(nc, site_id, dtime, &
       !      decomp_cpools_sourcesink_met, &
       !      decomp_cpools_sourcesink_cel, &
       !      decomp_cpools_sourcesink_lig)         


!       ! ---------------------------------------------------------------------------------
!       ! Part III.2 (continued).
!       ! Update diagnostics of the FATES ecosystem structure that are used in the HLM.
!       ! ---------------------------------------------------------------------------------
!       call this%wrap_update_hlmfates_dyn(nc,               &
!                                          bounds_clump,     &
!                                          waterstate_inst,  &
!                                          canopystate_inst, &
!                                          frictionvel_inst)
       
        call wrap_update_atsfates_dyn(nc)
       ! ---------------------------------------------------------------------------------
       ! Part IV: 
       ! Update history IO fields that depend on ecosystem dynamics
       ! ---------------------------------------------------------------------------------
       ! call fates_hist%update_history_dyn( nc,                    &
       !                                     fates(nc)%nsites, &
       !                                     fates(nc)%sites) 
     
      return
    end subroutine dynamics_driv_per_site
    

!    ! ------------------------------------------------------------------------------------
    subroutine UpdateLitterFluxes_per_site(nc, site_id, dtime, &
         decomp_cpools_sourcesink_met, &
         decomp_cpools_sourcesink_cel, &
         decomp_cpools_sourcesink_lig)         

       implicit none
       integer  :: site_id                        ! site index
       integer  , intent(in) :: nc                       ! clump index
       real(r8) , intent(in) :: dtime
       real(r8) , intent(inout) :: decomp_cpools_sourcesink_met(:)
       real(r8) , intent(inout) :: decomp_cpools_sourcesink_cel(:)
       real(r8) , intent(inout) :: decomp_cpools_sourcesink_lig(:)

!       ! !LOCAL VARIABLES:
       integer  :: nld_si

       nld_si = fates(nc)%bc_in(site_id)%nlevdecomp
       
       ! decomp_cpools_sourcesink_met(1:nld_si) = &
       !      fates(nc)%bc_out(site_id)%FATES_c_to_litr_lab_c_col(1:nld_si) * dtime
       ! decomp_cpools_sourcesink_cel(1:nld_si) = &
       !      fates(nc)%bc_out(site_id)%FATES_c_to_litr_cel_c_col(1:nld_si) * dtime
       ! decomp_cpools_sourcesink_lig(1:nld_si) = &
       !      fates(nc)%bc_out(site_id)%FATES_c_to_litr_lig_c_col(1:nld_si) * dtime


    end subroutine UpdateLitterFluxes_per_site

!    !--------------------------------------------------------------------------------------

    subroutine wrap_update_atsfates_dyn(nc)
!          waterstate_inst, canopystate_inst, frictionvel_inst )

       ! ---------------------------------------------------------------------------------
       ! This routine handles the updating of vegetation canopy diagnostics, (such as lai)
       ! that either requires HLM boundary conditions (like snow accumulation) or
       ! provides boundary conditions (such as vegetation fractional coverage)
       ! ---------------------------------------------------------------------------------

      implicit none

      integer                 , intent(in)           :: nc
 !     type(waterstate_type)   , intent(inout)        :: waterstate_inst
 !     type(canopystate_type)  , intent(inout)        :: canopystate_inst
 !     type(frictionvel_type)  , intent(inout)        :: frictionvel_inst
     
      integer :: npatch  ! number of patches in each site
      integer :: ifp     ! index FATES patch 
      integer :: p       ! HLM patch index
      integer :: s       ! site index
      integer :: c       ! column index

!      associate(                                &
!          tlai => canopystate_inst%tlai_patch , &
!          elai => canopystate_inst%elai_patch , &
!          tsai => canopystate_inst%tsai_patch , &
!          esai => canopystate_inst%esai_patch , &
!          htop => canopystate_inst%htop_patch , &
!          hbot => canopystate_inst%hbot_patch , & 
!          z0m  => frictionvel_inst%z0m_patch  , & ! Output: [real(r8) (:)   ] momentum roughness length (m)      
!          displa => canopystate_inst%displa_patch, &
!          dleaf_patch => canopystate_inst%dleaf_patch, &
!          snow_depth => waterstate_inst%snow_depth_col, &
!          frac_sno_eff => waterstate_inst%frac_sno_eff_col, &
!          frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch)


!        ! Process input boundary conditions to FATES
!        ! --------------------------------------------------------------------------------
       do s=1,fates(nc)%nsites
          fates(nc)%bc_in(s)%snow_depth_si   = 0.0_r8  !snow_depth(c)
          fates(nc)%bc_in(s)%frac_sno_eff_si = 0.0_r8  !frac_sno_eff(c)
       end do
       
       ! Canopy diagnostics for FATES
       call canopy_summarization(fates(nc)%nsites, &
            fates(nc)%sites,  &
            fates(nc)%bc_in)

        ! Canopy diagnostic outputs for HLM
        ! call update_hlm_dynamics(fates(nc)%nsites, &
        !      fates(nc)%sites,  &
        !      fates(nc)%bc_out )
        

     end subroutine wrap_update_atsfates_dyn



     subroutine calculate_biomass(nc, ats_biomass_array, nsites, num_scls) BIND(C)

       use FatesInterfaceMod, only : nlevsclass
       use EDtypesMod          , only : nfsc
       use FatesLitterMod      , only : ncwd
       use EDtypesMod          , only : ican_upper
       use EDtypesMod          , only : ican_ustory
       use FatesSizeAgeTypeIndicesMod, only : get_sizeage_class_index
       use FatesSizeAgeTypeIndicesMod, only : get_sizeagepft_class_index
       use FatesSizeAgeTypeIndicesMod, only : get_agepft_class_index
       use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
       use FatesSizeAgeTypeIndicesMod, only : get_height_index
       use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
       use EDTypesMod        , only : nlevleaf
       use EDParamsMod,           only : ED_val_history_height_bin_edges
       use EDtypesMod               , only : ed_cohort_type
       use EDtypesMod               , only : ed_patch_type
       use EDtypesMod               , only : AREA
       use EDtypesMod               , only : AREA_INV
       use FatesConstantsMod        , only : g_per_kg
       use PRTGenericMod            , only : leaf_organ, fnrt_organ, sapw_organ
       use PRTGenericMod            , only : struct_organ, store_organ, repro_organ
       use PRTGenericMod            , only : all_carbon_elements
       
       implicit none

       integer(C_INT),intent(in)                   :: nc   ! Clump
       real (C_DOUBLE),dimension(*), intent(inout) :: ats_biomass_array
       integer (C_INT), value :: nsites
       integer (C_INT), value :: num_scls


       integer :: s, scpf, scls
       real(r8) :: sapw_c, struct_c, leaf_c
       real(r8) :: fnrt_c, store_c
       real(r8) :: total_c, alive_c
       integer  :: ft               ! functional type index
       type(ed_patch_type),pointer  :: cpatch
       type(ed_cohort_type),pointer :: ccohort
       real(r8) :: n_density   ! individual of cohort per m2.
       real(r8) :: n_perm2     ! individuals per m2 for the whole column
       integer :: io_id, i_pa

       if (nsites.ne.fates(1)%nsites) then
          write(iulog,*) 'Number of sites provided by ATS does not match FATES'
          write(iulog,*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       if (num_scls.ne.nlevsclass) then
          write(iulog,*) 'Number of size classes provided by ATS does not match FATES'
          write(iulog,*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
          
       do s=1,nsites
          cpatch => fates(nc)%sites(s)%oldest_patch
          i_pa = 0
          do while(associated(cpatch))
             i_pa = i_pa + 1
             ccohort => cpatch%shortest
             do while(associated(ccohort))
                if(isnan(ccohort%n))then
                   write(iulog,*) 'cohort n become nan'
                   write(iulog,*) 'Aborting'	     
                   write(iulog,*)  ccohort%n, ccohort%gpp_acc
                   call endrun(msg=errMsg(sourcefile, __LINE__)) 
                endif
                if(isnan(ccohort%dbh))then
                   write(iulog,*) 'cohort dbh become nan'
                   write(iulog,*) 'Aborting'	     
                   write(iulog,*)  ccohort%n, ccohort%gpp_acc
                   call endrun(msg=errMsg(sourcefile, __LINE__)) 
                endif

                ft = ccohort%pft
                call sizetype_class_index(ccohort%dbh, ccohort%pft, ccohort%size_class, ccohort%size_by_pft_class)


                if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then

                   ! for quantities that are at the CLM patch level, because of the way 
                   ! that CLM patches are weighted for radiative purposes this # density needs 
                   ! to be over either ED patch canopy area or ED patch total area, whichever is less
                   !n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                   n_density = ccohort%n/cpatch%area 

                   ! for quantities that are natively at column level, calculate plant 
                   ! density using whole area
                   n_perm2   = ccohort%n * AREA_INV
                   ! write(100+i_pa,*) cpatch%area,cpatch%total_canopy_area,AREA_INV
                   ! write(200+i_pa,*) n_density

                else
                   n_density = 0.0_r8
                   n_perm2   = 0.0_r8
                endif



                associate(scpf => ccohort%size_by_pft_class, &
                     scls => ccohort%size_class)

                  ! Mass pools [kgC]
                sapw_c   = ccohort%prt%GetState(sapw_organ, all_carbon_elements)
                struct_c = ccohort%prt%GetState(struct_organ, all_carbon_elements)
                leaf_c   = ccohort%prt%GetState(leaf_organ, all_carbon_elements)
                fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)
                store_c  = ccohort%prt%GetState(store_organ, all_carbon_elements)

                alive_c  = leaf_c + fnrt_c + sapw_c
                total_c  = alive_c + store_c + struct_c

                io_id = (scls-1)*nsites + s
                ats_biomass_array(io_id) = ats_biomass_array(io_id) +  &
                      total_c * ccohort%n * AREA_INV * g_per_kg
                !ats_biomass_array(io_id) = ats_biomass_array(io_id) +  &
                !     total_c * n_density * g_per_kg

              end associate
              ccohort => ccohort%taller
           end do !cohort loop
           cpatch => cpatch%younger
        end do !patch loop


       end do
       
     end subroutine calculate_biomass

     

!    ! ====================================================================================

!    subroutine restart( this, bounds_proc, ncid, flag, waterstate_inst, &
!                              canopystate_inst, frictionvel_inst )

!       ! ---------------------------------------------------------------------------------
!       ! The ability to restart the model is handled through three different types of calls
!       ! "Define" the variables in the restart file, we "read" those variables into memory
!       ! or "write" data into the file from memory.  This subroutine accomodates all three
!       ! of those modes through the "flag" argument.  FATES as an external model also
!       ! requires an initialization step, where we set-up the dimensions, allocate and
!       ! flush the memory space that is used to transfer data in and out of the file.  This
!       ! Only occurs once, where as the define step occurs every time a file is opened.
!       !
!       ! Note: waterstate_inst and canopystate_inst are arguments only because following
!       ! the reading of variables, it is necessary to update diagnostics of the canopy
!       ! throug the interface call alm_fates%wrap_update_hlmfates_dyn() which requires
!       ! this information from the HLM.
!       ! ---------------------------------------------------------------------------------


!      use FatesConstantsMod, only : fates_long_string_length
!      use FatesIODimensionsMod, only: fates_bounds_type
!      use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int
!      use EDMainMod, only :        ed_update_site
!      use FatesInterfaceMod, only:  fates_maxElementsPerSite

!       implicit none

!       ! Arguments

!       class(hlm_fates_interface_type), intent(inout) :: this
!       type(bounds_type)              , intent(in)    :: bounds_proc
!       type(file_desc_t)              , intent(inout) :: ncid    ! netcdf id
!       character(len=*)               , intent(in)    :: flag
!       type(waterstate_type)          , intent(inout) :: waterstate_inst
!       type(canopystate_type)         , intent(inout) :: canopystate_inst
!       type(frictionvel_type)         , intent(inout) :: frictionvel_inst
      
!       ! Locals
!       type(bounds_type) :: bounds_clump
!       integer           :: nc
!       integer           :: nclumps
!       type(fates_bounds_type) :: fates_bounds
!       type(fates_bounds_type) :: fates_clump
!       integer                 :: c   ! HLM column index
!       integer                 :: s   ! Fates site index
!       integer                 :: g   ! HLM grid index
!       integer                 :: dk_index
!       character(len=fates_long_string_length) :: ioname
!       integer                 :: nvar
!       integer                 :: ivar
!       logical                 :: readvar

!       logical, save           :: initialized = .false.


!       nclumps = get_proc_clumps()

!       ! ---------------------------------------------------------------------------------
!       ! note (rgk: 11-2016) The history and restart intialization process assumes
!       ! that the number of site/columns active is a static entity.  Thus
!       ! we only allocate the mapping tables for the column/sites we start with.
!       ! If/when we start having dynamic column/sites (for reasons uknown as of yet)
!       ! we will need to re-evaluate the allocation of the mapping tables so they
!       ! can be unallocated,reallocated and set every time a new column/site is spawned
!       ! ---------------------------------------------------------------------------------

!       ! ---------------------------------------------------------------------------------
!       ! Only initialize the FATES restart structures the first time it is called
!       ! Note that the allocations involved with initialization are static.
!       ! This is because the array spaces for IO span the entire column, patch and cohort
!       ! range on the proc.
!       ! With DYNAMIC LANDUNITS or SPAWNING NEW OR CULLING OLD SITES:
!       ! we will in that case have to de-allocate, reallocate and then re-set the mapping
!       ! tables:  this%fates_restart%restart_map(nc)
!       ! I think that is it...
!       ! ---------------------------------------------------------------------------------

!       if(.not.initialized) then

!          initialized=.true.
      
!          ! ------------------------------------------------------------------------------
!          ! PART I: Set FATES DIMENSIONING INFORMATION
!          ! ------------------------------------------------------------------------------
         
!          call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)
         
!          call this%fates_restart%Init(nclumps, fates_bounds)
         
!          ! Define the bounds on the first dimension for each thread
!          !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
!          do nc = 1,nclumps
!             call get_clump_bounds(nc, bounds_clump)
            
!             ! thread bounds for patch
!             call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
!             call this%fates_restart%SetThreadBoundsEach(nc, fates_clump)
!          end do
!          !$OMP END PARALLEL DO
         
!          !$OMP PARALLEL DO PRIVATE (nc,s,c,g)
!          do nc = 1,nclumps
            
!             allocate(this%fates_restart%restart_map(nc)%site_index(this%fates(nc)%nsites))
!             allocate(this%fates_restart%restart_map(nc)%cohort1_index(this%fates(nc)%nsites))            
!             do s=1,this%fates(nc)%nsites
!                c = this%f2hmap(nc)%fcolumn(s)
!                this%fates_restart%restart_map(nc)%site_index(s)   = c
!                g = col_pp%gridcell(c)
!                this%fates_restart%restart_map(nc)%cohort1_index(s) = (g-1)*fates_maxElementsPerSite + 1
!             end do
            
!          end do
!          !$OMP END PARALLEL DO
         
!          ! ------------------------------------------------------------------------------------
!          ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
!          ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
!          ! ------------------------------------------------------------------------------------
!          call this%fates_restart%assemble_restart_output_types()
         
         
!          ! ------------------------------------------------------------------------------------
!          ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
!          ! HLM ACCORDING TO THEIR TYPES
!          ! ------------------------------------------------------------------------------------
!          call this%fates_restart%initialize_restart_vars()
         
!       end if

!       ! ---------------------------------------------------------------------------------
!       ! If we are writing, we must loop through our linked list structures and transfer the
!       ! information in the linked lists (FATES state memory) to the output vectors.
!       ! ---------------------------------------------------------------------------------

!       if(flag=='write')then
!          !$OMP PARALLEL DO PRIVATE (nc)
!          do nc = 1, nclumps
!             if (this%fates(nc)%nsites>0) then
!                call this%fates_restart%set_restart_vectors(nc,this%fates(nc)%nsites, &
!                                                            this%fates(nc)%sites)
!             end if
!          end do
!          !$OMP END PARALLEL DO
!       end if

!       ! ---------------------------------------------------------------------------------
!       ! In all cases, iterate through the list of variable objects
!       ! and either define, write or read to the NC buffer
!       ! This seems strange, but keep in mind that the call to restartvar()
!       ! has a different function in all three cases.
!       ! ---------------------------------------------------------------------------------

!       nvar = this%fates_restart%num_restart_vars()
!       do ivar = 1, nvar
            
!          associate( vname => this%fates_restart%rvars(ivar)%vname, &
!               vunits      => this%fates_restart%rvars(ivar)%units,   &
!               vlong       => this%fates_restart%rvars(ivar)%long )

!            dk_index = this%fates_restart%rvars(ivar)%dim_kinds_index
!            ioname = trim(this%fates_restart%dim_kinds(dk_index)%name)
        
!            select case(trim(ioname))
!            case(cohort_r8)

!               call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
!                     xtype=ncd_double,dim1name=trim('cohort'),long_name=trim(vlong), &
!                     units=trim(vunits),interpinic_flag='interp', &
!                     data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)
              
!            case(site_r8)
              
!               call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
!                     xtype=ncd_double,dim1name=trim('column'),long_name=trim(vlong), &
!                     units=trim(vunits),interpinic_flag='interp', &
!                     data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)
              
!            case(cohort_int)
              
!               call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
!                     xtype=ncd_int,dim1name=trim('cohort'),long_name=trim(vlong), &
!                     units=trim(vunits),interpinic_flag='interp', &
!                     data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)
              
!            case(site_int)
           
!               call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
!                     xtype=ncd_int,dim1name=trim('column'),long_name=trim(vlong), &
!                     units=trim(vunits),interpinic_flag='interp', &
!                     data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)
              
!            case default
!               write(iulog,*) 'A FATES iotype was created that was not registerred'
!               write(iulog,*) 'in ATS.:',trim(ioname)
!               call endrun(msg=errMsg(sourcefile, __LINE__))
!            end select
           
!          end associate
!       end do
      
!       ! ---------------------------------------------------------------------------------
!       ! If we are in a read mode, then we have just populated the sparse vectors
!       ! in the IO object list. The data in these vectors needs to be transferred
!       ! to the linked lists to populate the state memory.
!       ! ---------------------------------------------------------------------------------

!       if(flag=='read')then
         
!          !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s)
!          do nc = 1, nclumps
!             if (this%fates(nc)%nsites>0) then

!                call get_clump_bounds(nc, bounds_clump)

!                ! ------------------------------------------------------------------------
!                ! Convert newly read-in vectors into the FATES namelist state variables
!                ! ------------------------------------------------------------------------
!                call this%fates_restart%create_patchcohort_structure(nc, &
!                     this%fates(nc)%nsites, this%fates(nc)%sites, this%fates(nc)%bc_in)
               
!                call this%fates_restart%get_restart_vectors(nc, this%fates(nc)%nsites, &
!                     this%fates(nc)%sites )

!                ! I think ed_update_site and update_hlmfates_dyn are doing some similar
!                ! update type stuff, should consolidate (rgk 11-2016)
!                do s = 1,this%fates(nc)%nsites
!                   call ed_update_site( this%fates(nc)%sites(s), &
!                         this%fates(nc)%bc_in(s) )
!                end do

!                ! ------------------------------------------------------------------------
!                ! Update diagnostics of FATES ecosystem structure used in HLM.
!                ! ------------------------------------------------------------------------
!                call this%wrap_update_hlmfates_dyn(nc,bounds_clump, &
!                      waterstate_inst,canopystate_inst,frictionvel_inst)
               
!                ! ------------------------------------------------------------------------
!                ! Update history IO fields that depend on ecosystem dynamics
!                ! ------------------------------------------------------------------------
!                call this%fates_hist%update_history_dyn( nc, &
!                      this%fates(nc)%nsites,                 &
!                      this%fates(nc)%sites) 

               
!             end if
!          end do
!          !$OMP END PARALLEL DO
         
!       end if
      
!       return
!    end subroutine restart

!    !=====================================================================================

    
    subroutine init_coldstart(nc) BIND(C)

      !      ! Arguments
      implicit none
      ! Input Arguments
      integer(C_INT),intent(in)                 :: nc   ! Clump
      
!      ! locals
      integer :: s


!      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s,c,j,vol_ice,eff_porosity)

      call InitPARTEHGlobals()      
        
      if ( fates(nc)%nsites>0 ) then

         do s = 1, fates(nc)%nsites
            call init_site_vars(fates(nc)%sites(s), fates(nc)%bc_in(s))
            call zero_site(fates(nc)%sites(s))
         end do
           
         call set_site_properties(fates(nc)%nsites, fates(nc)%sites)


         call init_patches(fates(nc)%nsites, fates(nc)%sites, fates(nc)%bc_in)

         do s = 1, fates(nc)%nsites
            call ed_update_site(fates(nc)%sites(s), fates(nc)%bc_in(s))
         end do

         ! ------------------------------------------------------------------------
         ! Update diagnostics of FATES ecosystem structure used in HLM.
         ! ------------------------------------------------------------------------
         call wrap_update_atsfates_dyn(nc)

         ! ------------------------------------------------------------------------
         ! Update history IO fields that depend on ecosystem dynamics
         ! ------------------------------------------------------------------------
         ! call fates_hist%update_history_dyn(nc, fates(nc)%nsites, fates(nc)%sites) 
           
      end if

!      !$OMP END PARALLEL DO

    end subroutine init_coldstart
    
    subroutine DegCtoF(degC, degF) bind(c,name='DegCtoF')
      use, intrinsic :: iso_c_binding, only : c_double
      implicit none
      real(c_double), intent(in), dimension(:) :: degC
      real(c_double), intent(out), dimension(*) :: degF

      degF(1:SIZE(degC)) = degC*1.8+32
    end subroutine DegCtoF

    

    
    ! ======================================================================================
   
    subroutine wrap_sunfrac(nc, array_size, forc_solad, forc_solai) BIND(C)
         
       ! ---------------------------------------------------------------------------------
       ! This interface function is a wrapper call on ED_SunShadeFracs. The only
       ! returned variable is a patch vector, fsun_patch, which describes the fraction
       ! of the canopy that is exposed to sun.
       ! ---------------------------------------------------------------------------------
      
       implicit none
      
!       ! Input Arguments
!       class(hlm_fates_interface_type), intent(inout) :: this
!       type(bounds_type)              , intent(in)    :: bounds_clump
      
!       ! direct and diffuse downwelling radiation (W/m2)
!       type(atm2lnd_type),intent(in)        :: atm2lnd_inst
        integer(C_INT),intent(in)             :: nc                              !clump index
	integer(C_INT),intent(in)             :: array_size
        real(C_DOUBLE), intent(in)            :: forc_solad(array_size)          ! direct radiation (W/m**2); 1=visible lights; 2=near infrared radition
        real(C_DOUBLE), intent(in)            :: forc_solai(array_size)          ! diffuse radiation (W/m**2);  1=visible lights; 2=near infrared radition 

!       ! Input/Output Arguments to ATS
!       type(canopystate_type),intent(inout) :: canopystate_inst
      
!       ! Local Variables
!       integer  :: p                           ! global index of the host patch
!       integer  :: g                           ! global index of the host gridcell
!       integer  :: c                           ! global index of the host column

        integer  :: s                           ! FATES site index
!        integer  :: nc

!        nc = 1
        integer  :: ifp                         ! FATEs patch index
!                                               ! this is the order increment of patch
!                                               ! on the site
!       integer  :: nc                          ! clump index
      
!       type(ed_patch_type), pointer :: cpatch  ! c"urrent" patch  INTERF-TODO: SHOULD
!                                               ! BE HIDDEN AS A FATES PRIVATE

!       associate( forc_solad => atm2lnd_inst%forc_solad_grc, &
!                  forc_solai => atm2lnd_inst%forc_solai_grc, &
!                  fsun       => canopystate_inst%fsun_patch, &
!                  laisun     => canopystate_inst%laisun_patch, &               
!                  laisha     => canopystate_inst%laisha_patch )

!         nc = bounds_clump%clump_index
!         ! -------------------------------------------------------------------------------
!         ! Convert input BC's
!         ! The sun-shade calculations are performed only on FATES patches
!         ! -------------------------------------------------------------------------------
          do s = 1, fates(nc)%nsites
!            c = this%f2hmap(nc)%fcolumn(s)
!            g = col_pp%gridcell(c)

            do ifp = 1, fates(nc)%sites(s)%youngest_patch%patchno
!            !do ifp = 1, this%fates(nc)%bc_in(s)%npatches

!               p = ifp+col_pp%pfti(c)

               fates(nc)%bc_in(s)%solad_parb(ifp,:) = forc_solad(:)
               fates(nc)%bc_in(s)%solai_parb(ifp,:) = forc_solai(:)

            end do
         end do

         ! -------------------------------------------------------------------------------
         ! Call FATES public function to calculate internal sun/shade structures
         ! as well as total patch sun/shade fraction output boundary condition
         ! -------------------------------------------------------------------------------

         call ED_SunShadeFracs(fates(nc)%nsites, &
              fates(nc)%sites,  &
              fates(nc)%bc_in,  &
              fates(nc)%bc_out)

!         ! -------------------------------------------------------------------------------
!         ! Transfer the FATES output boundary condition for canopy sun/shade fraction
!         ! to the HLM
!         ! -------------------------------------------------------------------------------

!         do s = 1, this%fates(nc)%nsites
!            c = this%f2hmap(nc)%fcolumn(s)
!            do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
!               p = ifp+col_pp%pfti(c)
!               fsun(p)   = this%fates(nc)%bc_out(s)%fsun_pa(ifp)
!               laisun(p) = this%fates(nc)%bc_out(s)%laisun_pa(ifp)
!               laisha(p) = this%fates(nc)%bc_out(s)%laisha_pa(ifp)
!            end do
!         end do

!       end associate

    end subroutine wrap_sunfrac
   
!    ! ===================================================================================

    subroutine prep_canopyfluxes(nc) BIND(C)

!      ! ----------------------------------------------------------------------
!      ! the main function for calculating photosynthesis is called within a
!      ! loop based on convergence.  Some intitializations, including 
!      ! canopy resistance must be intitialized before the loop
!      ! The photosyns_ structure is currently unused, leaving it for now
!      ! in case we want to do any value initializing in future.
!      ! ----------------------------------------------------------------------

       implicit none  
       
!      ! Arguments
       integer(C_INT),intent(in)             :: nc                              !clump index
!      class(hlm_fates_interface_type), intent(inout) :: this
!      type(bounds_type)              , intent(in)    :: bounds_clump

      ! locals
       integer                                        :: s
!      integer                                        :: nc

!      nc = bounds_clump%clump_index
       do s = 1, fates(nc)%nsites
         ! filter flag == 1 means that this patch has not been called for photosynthesis
         fates(nc)%bc_in(s)%filter_photo_pa(:) = 1
        ! set transpiration input boundary condition to zero. The exposed
        ! vegetation filter may not even call every patch.
!        if (use_fates_planthydro) then
!           this%fates(nc)%bc_in(s)%qflx_transp_pa(:) = 0._r8
!        end if
      end do
    end subroutine prep_canopyfluxes

!    ! ====================================================================================
   
    subroutine wrap_btran(nc, array_size, t_soil, h2osoi_liqvol, eff_porosity, watsat, soil_suc, salinity) BIND(C)
      ! ,soilstate_inst, waterstate_inst, &
      !                    temperature_inst, energyflux_inst,  &
      !                    soil_water_retention_curve)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates btran for FATES, this will be an input boundary
      ! condition for FATES photosynthesis/transpiration.
      !
      ! This subroutine also calculates rootr
      ! 
      ! ---------------------------------------------------------------------------------

      ! use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

      implicit none
      
      ! Arguments
      integer(C_INT),intent(in)             :: nc                              !clump index
      integer(C_INT),intent(in)             :: array_size
      real(C_DOUBLE), intent(in)            :: t_soil(array_size)          !soil temperature (Kelvin)
      real(C_DOUBLE), intent(in)            :: h2osoi_liqvol(array_size)   !liquid volumetric moisture, will be used for BeTR
      real(C_DOUBLE), intent(in)            :: eff_porosity(array_size)    !effective porosity = porosity - vol_ice 
      real(C_DOUBLE), intent(in)            :: watsat(array_size)          !volumetric soil water at saturation (porosity)
      real(C_DOUBLE), intent(in)            :: soil_suc(array_size)        !soil suction,  negative, [mm]
      real(C_DOUBLE), intent(in)            :: salinity(array_size)        !salt concentration [ppt]
      
                                                                        ! columns with exposed veg
      ! type(soilstate_type)   , intent(inout)         :: soilstate_inst
      ! type(waterstate_type)  , intent(in)            :: waterstate_inst
      ! type(temperature_type) , intent(in)            :: temperature_inst
      ! type(energyflux_type)  , intent(inout)         :: energyflux_inst
      ! class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

      ! local variables
      real(r8) :: smp_node ! Soil suction potential, negative, [mm]
      real(r8) :: s_node
      integer  :: s
      integer  :: c
      integer  :: j
      integer  :: k
      integer  :: ifp
      integer  :: p
      integer  :: nlevsoil
      !integer  :: nc

      !nc = 1

      ! -------------------------------------------------------------------------------
      ! Convert input BC's
      ! Critical step: a filter is being passed in that dictates which columns have
      ! exposed vegetation (above snow).  This is necessary, because various hydrologic
      ! variables like h2osoi_liqvol are not calculated and will have uninitialized
      ! values outside this list.
      !
      ! bc_in(s)%filter_btran      (this is in, but is also used in this subroutine)
      !
      ! We also filter a second time within this list by determining which soil layers
      ! have conditions for active uptake based on soil moisture and temperature. This
      ! must be determined by FATES (science stuff).  But the list of layers and patches
      ! needs to be passed back to the interface, because it then needs to request
      ! suction on these layers via ATS/ALM functions.  We cannot wide-swath calculate
      ! this on all layers, because values with no moisture or low temps will generate
      ! unstable values and cause sigtraps.
      ! -------------------------------------------------------------------------------
      k = 0
      do s = 1, fates(nc)%nsites
         nlevsoil = fates(nc)%bc_in(s)%nlevsoil

         fates(nc)%bc_in(s)%filter_btran = .true.

         do j = 1,nlevsoil
            k=k+1
            fates(nc)%bc_in(s)%tempk_sl(j)         = t_soil(k)
            fates(nc)%bc_in(s)%h2o_liqvol_sl(j)    = h2osoi_liqvol(k)
            fates(nc)%bc_in(s)%eff_porosity_sl(j)  = eff_porosity(k)
            fates(nc)%bc_in(s)%watsat_sl(j)        = watsat(k)
            fates(nc)%bc_in(s)%smp_sl(j)           = soil_suc(k)
            fates(nc)%bc_in(s)%salinity_sl(j)      = salinity(k)
         end do

      end do


      ! -------------------------------------------------------------------------------
      ! Suction and active uptake layers calculated, lets calculate uptake (btran)
      ! This will calculate internals, as well as output boundary conditions: 
      ! btran, rootr
      ! -------------------------------------------------------------------------------

      call btran_ed(fates(nc)%nsites, &
           fates(nc)%sites,  &
           fates(nc)%bc_in,  &
           fates(nc)%bc_out)


   end subroutine wrap_btran

!    ! ====================================================================================

   subroutine wrap_photosynthesis(nc, dtime, p_atm, array_size, t_soil, photosys_in) BIND(C)


   
     !ARGUMENTS:
     integer(C_INT),intent(in)             :: nc                              !clump index
     real(C_DOUBLE), intent(in)            :: p_atm, dtime
     integer(C_INT),intent(in)             :: array_size
     real(C_DOUBLE), intent(in)            :: t_soil(array_size)
     type(PhotoSynthesisInput),intent(in)  :: photosys_in


     !LOCAL:
     integer                               :: nlevsoil
     integer                               :: s,c,p,ifp,j,icp,k

!     call t_startf('edpsn')     

!     nc = 1
     k = 0
     do s = 1, fates(nc)%nsites

        nlevsoil = fates(nc)%bc_in(s)%nlevsoil

        do j = 1,nlevsoil
           k=k+1
           fates(nc)%bc_in(s)%t_soisno_sl(j)   = t_soil(k)  ! soil temperature (Kelvin)
        end do
        fates(nc)%bc_in(s)%forc_pbot           = p_atm   ! atmospheric pressure (Pa)
     end do

     do s = 1, fates(nc)%nsites

        ! This filter is flushed to 1 before the canopyflux stability iterator
        ! It is set to status 2 if it is an active patch within the iterative loop
        ! After photosynthesis is called, it is upgraded to 3 if it was called.
        ! After all iterations we can evaluate which patches have a final flag
        ! of 3 to check if we missed any.
        
        fates(nc)%bc_in(s)%filter_photo_pa(:) = 2
        fates(nc)%bc_in(s)%dayl_factor_pa(:) = photosys_in%dayl_factor ! scalar (0-1) for daylength
        fates(nc)%bc_in(s)%esat_tv_pa(:)     = photosys_in%esat_tv     ! saturation vapor pressure at t_veg (Pa)
        fates(nc)%bc_in(s)%eair_pa(:)        = photosys_in%eair        ! vapor pressure of canopy air (Pa)
        fates(nc)%bc_in(s)%oair_pa(:)        = photosys_in%oair        ! Atmospheric O2 partial pressure (Pa)
        fates(nc)%bc_in(s)%cair_pa(:)        = photosys_in%cair        ! Atmospheric CO2 partial pressure (Pa)
        fates(nc)%bc_in(s)%rb_pa(:)          = photosys_in%rb          ! boundary layer resistance (s/m)
        fates(nc)%bc_in(s)%t_veg_pa(:)       = photosys_in%t_veg       ! vegetation temperature (Kelvin)     
        fates(nc)%bc_in(s)%tgcm_pa(:)        = photosys_in%tgcm        ! air temperature at agcm 
                                                                       ! reference height (kelvin)
     end do


      
    ! Call photosynthesis
    
    call FatesPlantRespPhotosynthDrive (fates(nc)%nsites, &
         fates(nc)%sites,  &
         fates(nc)%bc_in,  &
         fates(nc)%bc_out, &
         dtime)

    ! Perform a double check to see if all patches on naturally vegetated columns
    ! were activated for photosynthesis
    ! ---------------------------------------------------------------------------------
    do s = 1, fates(nc)%nsites	
       do p=1, fates(nc)%sites(s)%youngest_patch%patchno
         if(fates(nc)%bc_in(s)%filter_photo_pa(p) /= 2)then
                write(iulog,*) 'Not all patches on the natveg column in the photosynthesis'
                write(iulog,*) 'filter ran photosynthesis'
                call endrun(msg=errMsg(sourcefile, __LINE__))
         else  
	    fates(nc)%bc_in(s)%filter_photo_pa(p) = 3
	 end if
       end do	  
    end do
    !      do icp = 1,fn
    !         p = filterp(icp)
    !         c = veg_pp%column(p)
    !         s = this%f2hmap(nc)%hsites(c)
    !         ! do if structure here and only pass natveg columns
    !         ifp = p-col_pp%pfti(c)
    !         if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 2)then
    !            write(iulog,*) 'Not all patches on the natveg column in the photosynthesis'
    !            write(iulog,*) 'filter ran photosynthesis'
    !            call endrun(msg=errMsg(sourcefile, __LINE__))
    !         else
    !            this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 3
    !            rssun(p) = this%fates(nc)%bc_out(s)%rssun_pa(ifp)
    !            rssha(p) = this%fates(nc)%bc_out(s)%rssha_pa(ifp)

    !            ! These fields are marked with a bad-value flag
    !            photosyns_inst%psnsun_patch(p)   = spval
    !            photosyns_inst%psnsha_patch(p)   = spval
    !         end if
    !      end do

    !    end associate
!    call t_stopf('edpsn')

    end subroutine wrap_photosynthesis

!  ! ======================================================================================

  subroutine wrap_accumulatefluxes(nc, dtime) BIND(C)

!    ! !ARGUMENTS:
!    class(hlm_fates_interface_type), intent(inout) :: this
!    type(bounds_type)              , intent(in)    :: bounds_clump
!    integer                        , intent(in)    :: fn                   ! size of pft filter
!    integer                        , intent(in)    :: filterp(fn)          ! pft filter
     integer(C_INT),intent(in)             :: nc    ! Clump index
     real(C_DOUBLE), intent(in)            :: dtime !time step, seconds
   
!    ! Locals
!    integer                                        :: s,c,p,ifp,icp
!    real(r8)                                       :: dtime
!    integer                                        :: nc

!    nc = bounds_clump%clump_index
!     ! Run a check on the filter
!     do icp = 1,fn
!        p = filterp(icp)
!        c = veg_pp%column(p)
!        s = this%f2hmap(nc)%hsites(c)
!        ifp = p-col_pp%pfti(c)
!        if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 3)then
!           call endrun(msg=errMsg(sourcefile, __LINE__))
!        end if
!     end do


!     dtime = get_step_size()

     call  AccumulateFluxes_ED(fates(nc)%nsites,  &
                                fates(nc)%sites, &
                                fates(nc)%bc_in,  &
                                fates(nc)%bc_out, &
                                dtime)

    
!     call this%fates_hist%update_history_prod(nc, &
!                                this%fates(nc)%nsites,  &
!                                this%fates(nc)%sites, &
!                                dtime)

  end subroutine wrap_accumulatefluxes
  
 
! ======================================================================================

  subroutine wrap_canopy_radiation(nc, jday,array_size, albgrd,albgri) BIND(C)


!     ! Arguments
!     class(hlm_fates_interface_type), intent(inout) :: this
!     type(bounds_type),  intent(in)             :: bounds_clump
!     ! filter for vegetated pfts with coszen>0
!     integer            , intent(in)            :: num_vegsol                 
!     integer            , intent(in)            :: filter_vegsol(num_vegsol)    
!     ! cosine solar zenith angle for next time step
!     real(r8)           , intent(in)            :: coszen( bounds_clump%begp: )        
!     type(surfalb_type) , intent(inout)         :: surfalb_inst 
    
!     ! locals
!     integer                                    :: s,c,p,ifp,icp,nc
      integer(C_INT),intent(in)             :: nc                              !clump index
      real   (C_DOUBLE), intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
      integer(C_INT),intent(in)  :: array_size   !radiation index size = 2
      real(C_DOUBLE), intent(in) :: albgrd(array_size)  !ground albedo (direct) 1=visiable; 2=near infrared (nir)
      real(C_DOUBLE), intent(in) :: albgri(array_size)  !ground albedo (diffuse) 1=visiable; 2=near infrared (nir)
      real (r8) :: lat    ! Centered latitude degrees(-90,...90)
      real (r8) :: lon    ! Centered longitude degrees(0...360)
      real (r8) :: coszen   ! cosine solar zenith angle for next time step
      integer  :: s         ! the site index
!     integer  :: nc         ! the thread index

!      nc = 1
!     associate(&
!          albgrd_col   =>    surfalb_inst%albgrd_col         , & !in
!          albgri_col   =>    surfalb_inst%albgri_col         , & !in
!          albd         =>    surfalb_inst%albd_patch         , & !out
!          albi         =>    surfalb_inst%albi_patch         , & !out
!          fabd         =>    surfalb_inst%fabd_patch         , & !out
!          fabi         =>    surfalb_inst%fabi_patch         , & !out
!          ftdd         =>    surfalb_inst%ftdd_patch         , & !out
!          ftid         =>    surfalb_inst%ftid_patch         , & !out
!          ftii         =>    surfalb_inst%ftii_patch)            !out

!     nc = bounds_clump%clump_index

     do s = 1, fates(nc)%nsites

!        c = this%f2hmap(nc)%fcolumn(s)
!        do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
          
!           p = ifp+col_pp%pfti(c)
          
!           if( any(filter_vegsol==p) )then !cxu: need to updated for non-vegetation land;
              lat= fates(nc)%sites(s)%lat 
              lon= fates(nc)%sites(s)%lon
              coszen = shr_orb_cosz(jday,lat,lon)
              fates(nc)%bc_in(s)%filter_vegzen_pa(:) = .true.
              fates(nc)%bc_in(s)%coszen_pa(:)  = coszen
              fates(nc)%bc_in(s)%albgr_dir_rb(:) = albgrd(:)
              fates(nc)%bc_in(s)%albgr_dif_rb(:) = albgri(:)

!           else
             
!              this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) = .false.

!           end if

!        end do
     end do

     call ED_Norman_Radiation(fates(nc)%nsites,  &
          fates(nc)%sites, &
          fates(nc)%bc_in,  &
          fates(nc)%bc_out)
    
!     ! Pass FATES BC's back to HLM
!     ! -----------------------------------------------------------------------------------
!     do icp = 1,num_vegsol
!        p = filter_vegsol(icp)
!        c = veg_pp%column(p)
!        s = this%f2hmap(nc)%hsites(c)
!        ! do if structure here and only pass natveg columns
!        ifp = p-col_pp%pfti(c)

!        if(.not.this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) )then
!           write(iulog,*) 'Not all patches on the natveg column were passed to canrad'
!           call endrun(msg=errMsg(sourcefile, __LINE__))
!        else
!           albd(p,:) = this%fates(nc)%bc_out(s)%albd_parb(ifp,:)
!           albi(p,:) = this%fates(nc)%bc_out(s)%albi_parb(ifp,:)
!           fabd(p,:) = this%fates(nc)%bc_out(s)%fabd_parb(ifp,:)
!           fabi(p,:) = this%fates(nc)%bc_out(s)%fabi_parb(ifp,:)
!           ftdd(p,:) = this%fates(nc)%bc_out(s)%ftdd_parb(ifp,:)
!           ftid(p,:) = this%fates(nc)%bc_out(s)%ftid_parb(ifp,:)
!           ftii(p,:) = this%fates(nc)%bc_out(s)%ftii_parb(ifp,:)
!        end if
!     end do
    
!   end associate

  end subroutine wrap_canopy_radiation

!  ! ======================================================================================

!  subroutine wrap_bgc_summary(this, bounds_clump, carbonflux_inst, carbonstate_inst)

   

!     ! Arguments
!     class(hlm_fates_interface_type), intent(inout) :: this
!     type(bounds_type),  intent(in)                 :: bounds_clump
!     type(carbonflux_type), intent(in)              :: carbonflux_inst
!     type(carbonstate_type), intent(in)             :: carbonstate_inst

!     ! locals
!     real(r8) :: dtime
!     integer  :: nstep
!     logical  :: is_beg_day
!     integer  :: s,c,nc

!     associate(& 
!         hr            => carbonflux_inst%hr_col,      & ! (gC/m2/s) total heterotrophic respiration
!         totsomc       => carbonstate_inst%totsomc_col, & ! (gC/m2) total soil organic matter carbon
!         totlitc       => carbonstate_inst%totlitc_col)   ! (gC/m2) total litter carbon in BGC pools
      
!       nc = bounds_clump%clump_index
      
!       ! Summarize Net Fluxes
!       do s = 1, this%fates(nc)%nsites
!          c = this%f2hmap(nc)%fcolumn(s)
!          this%fates(nc)%bc_in(s)%tot_het_resp = hr(c)
!          this%fates(nc)%bc_in(s)%tot_somc     = totsomc(c)
!          this%fates(nc)%bc_in(s)%tot_litc     = totlitc(c)
!       end do
      
!       is_beg_day = is_beg_curr_day()
!       dtime = get_step_size()
!       nstep = get_nstep()
      
!       call SummarizeNetFluxes(this%fates(nc)%nsites,  &
!                              this%fates(nc)%sites,    &
!                              this%fates(nc)%bc_in,    &
!                              is_beg_day)
      

!       call FATES_BGC_Carbon_Balancecheck(this%fates(nc)%nsites,  &
!                                          this%fates(nc)%sites, &
!                                          this%fates(nc)%bc_in,  &
!                                          is_beg_day,            &
!                                          dtime, nstep)
      

!       ! Update history variables that track these variables
!       call this%fates_hist%update_history_cbal(nc, &
!                                this%fates(nc)%nsites,  &
!                                this%fates(nc)%sites)

      
!     end associate
!  end subroutine wrap_bgc_summary

!  ! ======================================================================================


!  subroutine TransferZ0mDisp(this,bounds_clump,frictionvel_inst,canopystate_inst)

!     ! Arguments
!     class(hlm_fates_interface_type), intent(inout) :: this
!     type(bounds_type),intent(in)                   :: bounds_clump
!     type(canopystate_type)  , intent(inout)        :: canopystate_inst
!     type(frictionvel_type)  , intent(inout)        :: frictionvel_inst

!     ! Locals
!     integer :: ci   ! Current clump index
!     integer :: s    ! Site index
!     integer :: c    ! Column index
!     integer :: ifp  ! Fates patch index
!     integer :: p    ! ATS patch index

!     ci = bounds_clump%clump_index

!     do s = 1, this%fates(ci)%nsites
!        c = this%f2hmap(ci)%fcolumn(s)

!        frictionvel_inst%z0m_patch(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8
!        canopystate_inst%displa_patch(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8

!        do ifp = 1, this%fates(ci)%sites(s)%youngest_patch%patchno
!           p = ifp+col_pp%pfti(c)
!           frictionvel_inst%z0m_patch(p) = this%fates(ci)%bc_out(s)%z0m_pa(ifp)
!           canopystate_inst%displa_patch(p) = this%fates(ci)%bc_out(s)%displa_pa(ifp)
!        end do
!     end do

!     return
!  end subroutine TransferZ0mDisp

!  ! ======================================================================================

  subroutine init_history_io()

   !use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 

   use FatesConstantsMod, only : fates_short_string_length, fates_long_string_length
   use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
   use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
   use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
   use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
   use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
   use FatesIOVariableKindMod, only : site_height_r8
   use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
   use FatesIODimensionsMod, only : fates_bounds_type


 !   ! Arguments
 !   class(hlm_fates_interface_type), intent(inout) :: this
 !  type(bounds_type),intent(in)                   :: bounds_proc  ! Currently "proc"
   
   
   ! Locals
 !  type(bounds_type)                              :: bounds_clump
   integer :: nvar  ! number of IO variables found
   integer :: ivar  ! variable index 1:nvar
   integer :: nc    ! thread counter 1:nclumps
   integer :: nclumps ! number of threads on this proc
   integer :: s     ! FATES site index
   integer :: c     ! ALM/ATS column index
   character(len=fates_short_string_length) :: dim2name
   character(len=fates_long_string_length) :: ioname
   integer :: d_index, dk_index
   
   type(fates_bounds_type) :: fates_bounds
   type(fates_bounds_type) :: fates_clump

 !   ! This routine initializes the types of output variables
 !   ! not the variables themselves, just the types
 !   ! ---------------------------------------------------------------------------------

    nclumps = 1   !get_proc_clumps()

 !   ! ------------------------------------------------------------------------------------
 !   ! PART I: Set FATES DIMENSIONING INFORMATION
 !   !       
 !   ! -------------------------------------------------------------------------------
 !   ! Those who wish add variables that require new dimensions, please
 !   ! see FATES: FatesHistoryInterfaceMod.F90.  Dimension types are defined at the top of the
 !   ! module, and a new explicitly named instance of that type should be created.
 !   ! With this new dimension, a new output type/kind can contain that dimension.
 !   ! A new type/kind can be added to the dim_kinds structure, which defines its members
 !   ! in created in init_dim_kinds_maps().  Make sure to increase the size of fates_num_dim_kinds.
 !   ! A type/kind of output is defined by the data type (ie r8,int,..)
 !   ! and the dimensions.  Keep in mind that 3D variables (or 4D if you include time)
 !   ! are not really supported in ATS/ALM right now.  There are ways around this
 !   ! limitations by creating combined dimensions, for instance the size+pft dimension
 !   ! "scpf"
 !   ! ------------------------------------------------------------------------------------
   
 !   call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)

    call fates_hist%Init(nclumps, fates_bounds)

 !   ! Define the bounds on the first dimension for each thread
 !   !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
 !   do nc = 1,nclumps
      
 !      call get_clump_bounds(nc, bounds_clump)
      
 !      ! thread bounds for patch
 !      call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
 !      call this%fates_hist%SetThreadBoundsEach(nc, fates_clump)
 !   end do
 !   !$OMP END PARALLEL DO

 !   ! ------------------------------------------------------------------------------------
 !   ! PART I.5: SET SOME INDEX MAPPINGS SPECIFICALLY FOR SITE<->COLUMN AND PATCH 
 !   ! ------------------------------------------------------------------------------------
   
 !   !$OMP PARALLEL DO PRIVATE (nc,s,c)
 !   do nc = 1,nclumps
      
       ! allocate(fates_hist%iovar_map(nc)%site_index(this%fates(nc)%nsites))
       ! allocate(fates_hist%iovar_map(nc)%patch1_index(this%fates(nc)%nsites))
      
       ! do s=1,fates(nc)%nsites
       !    c = f2hmap(nc)%fcolumn(s)
       !    this%fates_hist%iovar_map(nc)%site_index(s)   = c
       !    this%fates_hist%iovar_map(nc)%patch1_index(s) = col_pp%pfti(c)+1
       ! end do
       
 !   end do
 !   !$OMP END PARALLEL DO
   
   ! ------------------------------------------------------------------------------------
   ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
   ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
   ! ------------------------------------------------------------------------------------
   call fates_hist%assemble_history_output_types()
   
   ! ------------------------------------------------------------------------------------
   ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
   ! HLM ACCORDING TO THEIR TYPES
   ! ------------------------------------------------------------------------------------
   call fates_hist%initialize_history_vars()
   nvar = fates_hist%num_history_vars()
   
   do ivar = 1, nvar
      
      associate( vname    => fates_hist%hvars(ivar)%vname, &
                 vunits   => fates_hist%hvars(ivar)%units,   &
                 vlong    => fates_hist%hvars(ivar)%long, &
                 vdefault => fates_hist%hvars(ivar)%use_default, &
                 vavgflag => fates_hist%hvars(ivar)%avgflag)

      dk_index = fates_hist%hvars(ivar)%dim_kinds_index
      ioname = trim(fates_hist%dim_kinds(dk_index)%name)
        
      ! select case(trim(ioname))
      ! case(patch_r8)
      !    call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
      !         avgflag=trim(vavgflag),long_name=trim(vlong), &
      !         ptr_patch=this%fates_hist%hvars(ivar)%r81d,    &
      !         default=trim(vdefault),                       &
      !         set_lake=0._r8,set_urb=0._r8)
           
 !        case(site_r8)
 !           call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r81d,      & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)

 !        case(patch_ground_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         & ! <--- addfld2d
 !                              type2d=trim(dim2name),                        & ! <--- type2d
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_patch=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
           
 !        case(patch_size_pft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_patch=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_ground_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,      & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_size_pft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,      & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_size_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_pft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_age_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)

 !        case(site_height_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,     & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)

 !        case(site_scagpft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,     & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)

 !        case(site_agepft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,     & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)

 !        case(site_fuel_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_cwdsc_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_can_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_cnlf_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_cnlfpft_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
 !        case(site_scag_r8)
 !           d_index = this%fates_hist%dim_kinds(dk_index)%dim2_index
 !           dim2name = this%fates_hist%dim_bounds(d_index)%name
 !           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
 !                              type2d=trim(dim2name),                        &
 !                              avgflag=trim(vavgflag),long_name=trim(vlong), &
 !                              ptr_col=this%fates_hist%hvars(ivar)%r82d,    & 
 !                              default=trim(vdefault),                       &
 !                              set_lake=0._r8,set_urb=0._r8)
           

 !        case default
 !           write(iulog,*) 'A FATES iotype was created that was not registerred'
 !           write(iulog,*) 'in ATS.:',trim(ioname)
 !           call endrun(msg=errMsg(sourcefile, __LINE__))
 !        end select
          
    end associate
   end do
  end subroutine init_history_io

 ! ======================================================================================
 


!  ! ======================================================================================

!  subroutine ComputeRootSoilFlux(this, bounds_clump, num_filterc, filterc, &
!        soilstate_inst, waterflux_inst)

!     class(hlm_fates_interface_type), intent(inout) :: this
!     type(bounds_type),intent(in)                   :: bounds_clump
!     integer,intent(in)                             :: num_filterc
!     integer,intent(in)                             :: filterc(num_filterc)
!     type(soilstate_type), intent(inout)            :: soilstate_inst
!     type(waterflux_type), intent(inout)            :: waterflux_inst
    
!     ! locals
!     integer :: s
!     integer :: c 
!     integer :: l
!     integer :: nc
!     integer :: num_filter_fates
!     integer :: nlevsoil


!     if( .not. use_fates_planthydro ) return
       
!     nc = bounds_clump%clump_index
    
!     ! Perform a check that the number of columns submitted to fates for 
!     ! root water sink is the same that was expected in the hydrology filter
!     num_filter_fates = 0
!     do s = 1,num_filterc
!        l = col_pp%landunit(filterc(s))
!        if (lun_pp%itype(l) == istsoil ) then
!           num_filter_fates = num_filter_fates + 1
!        end if
!     end do
    
!     if(num_filter_fates .ne. this%fates(nc)%nsites )then
!        write(iulog,*) 'The HLM list of natural veg columns during root water transfer'
!        write(iulog,*) 'is not the same size as the fates site list?'
!        call endrun(msg=errMsg(sourcefile, __LINE__))
!     end if
    
!     do s = 1, this%fates(nc)%nsites
!        c = this%f2hmap(nc)%fcolumn(s)
!        nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
!        waterflux_inst%qflx_rootsoi_col(c,1:nlevsoil) = &
!             this%fates(nc)%bc_out(s)%qflx_soil2root_sisl(1:nlevsoil)
!     end do
    
!  end subroutine ComputeRootSoilFlux

!  ! ======================================================================================
! !
! ! THIS WAS MOVED TO WRAP_HYDRAULICS_DRIVE()
! !
! ! subroutine TransferPlantWaterStorage(this, bounds_clump, nc, waterstate_inst)
! !   
! !   implicit none
! !   class(hlm_fates_interface_type), intent(inout) :: this
! !   type(bounds_type),intent(in)                   :: bounds_clump
! !   integer,intent(in)                             :: nc
! !   type(waterstate_type)   , intent(inout)        :: waterstate_inst
! !
! !   ! locals
! !   integer :: s
! !   integer :: c 
! !   
! !   if (.not. (use_fates .and. use_fates_planthydro) ) return
! !   
! !   do s = 1, this%fates(nc)%nsites
! !      c = this%f2hmap(nc)%fcolumn(s)
! !      waterstate_inst%total_plant_stored_h2o_col(c) = &
! !            this%fates(nc)%bc_out(s)%plant_stored_h2o_si
! !   end do
! !   return
! !end subroutine TransferPlantWaterStorage




!  ! ======================================================================================

!  subroutine wrap_hydraulics_drive(this, bounds_clump, &
!                                  soilstate_inst, waterstate_inst, waterflux_inst, &
!                                  solarabs_inst, energyflux_inst)


!    implicit none
!    class(hlm_fates_interface_type), intent(inout) :: this
!    type(bounds_type),intent(in)                   :: bounds_clump
!    type(soilstate_type)    , intent(inout)        :: soilstate_inst
!    type(waterstate_type)   , intent(inout)        :: waterstate_inst
!    type(waterflux_type)    , intent(inout)        :: waterflux_inst
!    type(solarabs_type)     , intent(in)           :: solarabs_inst
!    type(energyflux_type)   , intent(inout)        :: energyflux_inst

!     ! locals
!    integer :: s
!    integer :: c 
!    integer :: j
!    integer :: ifp
!    integer :: p
!    integer :: nc
!    real(r8) :: dtime
!    integer  :: nlevsoil


!    if ( .not.use_fates_planthydro ) return

!    nc = bounds_clump%clump_index
!    dtime = get_step_size()

!    ! Prepare Input Boundary Conditions
!    ! ------------------------------------------------------------------------------------

!    do s = 1, this%fates(nc)%nsites
!       c = this%f2hmap(nc)%fcolumn(s)
!       nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

!       this%fates(nc)%bc_in(s)%smpmin_si                 = &
!             soilstate_inst%smpmin_col(c)
!       this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil)    = &
!             soilstate_inst%watsat_col(c,1:nlevsoil) 
!       this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil)    = spval !&
! !            soilstate_inst%watres_col(c,1:nlevsoil)
!       this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil)     = &
!             soilstate_inst%sucsat_col(c,1:nlevsoil)
!       this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)        = &
!             soilstate_inst%bsw_col(c,1:nlevsoil)
!       this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil)    = &
!             waterstate_inst%h2osoi_liq_col(c,1:nlevsoil)
!       this%fates(nc)%bc_in(s)%eff_porosity_sl(1:nlevsoil) = &
!             soilstate_inst%eff_porosity_col(c,1:nlevsoil)

!       do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno 
!          p = ifp+col_pp%pfti(c)
!          this%fates(nc)%bc_in(s)%swrad_net_pa(ifp) = solarabs_inst%fsa_patch(p)
!          this%fates(nc)%bc_in(s)%lwrad_net_pa(ifp) = energyflux_inst%eflx_lwrad_net_patch(p)
!          this%fates(nc)%bc_in(s)%qflx_transp_pa(ifp) = waterflux_inst%qflx_tran_veg_patch(p)
!       end do
!    end do

!    ! Call Fates Hydraulics
!    ! ------------------------------------------------------------------------------------


!    call hydraulics_drive(this%fates(nc)%nsites, &
!             this%fates(nc)%sites,  &
!             this%fates(nc)%bc_in,  &
!             this%fates(nc)%bc_out, &
!             dtime)

!    ! Prepare Output Boundary Conditions
!    ! ------------------------------------------------------------------------------------

!    do s = 1, this%fates(nc)%nsites
!       c = this%f2hmap(nc)%fcolumn(s)
!       waterstate_inst%total_plant_stored_h2o_col(c) = &
!             this%fates(nc)%bc_out(s)%plant_stored_h2o_si
               
!    end do
   
   


!    ! Update History Buffers that need to be updated after hydraulics calls

!    call this%fates_hist%update_history_hydraulics(nc, &
!          this%fates(nc)%nsites, &
!          this%fates(nc)%sites, &
!          dtime)


!    return
!  end subroutine wrap_hydraulics_drive

!  ! ======================================================================================

!  subroutine hlm_bounds_to_fates_bounds(hlm, fates)

!    use FatesIODimensionsMod, only : fates_bounds_type
!    use FatesInterfaceMod, only : nlevsclass_fates => nlevsclass
!    use FatesInterfaceMod, only : nlevage_fates    => nlevage
!    use FatesInterfaceMod, only : nlevheight_fates => nlevheight
!    use EDtypesMod,        only : nfsc_fates       => nfsc
!    use EDtypesMod,        only : ncwd_fates       => ncwd
!    use EDtypesMod,        only : nlevleaf_fates   => nlevleaf
!    use EDtypesMod,        only : nclmax_fates     => nclmax
!    use clm_varpar,        only : nlevgrnd
!    use FatesInterfaceMod, only : numpft_fates     => numpft

!    implicit none

!    type(bounds_type), intent(in)        :: hlm
!    type(fates_bounds_type), intent(out) :: fates

!    fates%cohort_begin = hlm%begcohort
!    fates%cohort_end = hlm%endcohort
   
!    fates%patch_begin = hlm%begp
!    fates%patch_end = hlm%endp
   
!    fates%column_begin = hlm%begc
!    fates%column_end = hlm%endc
   
!    fates%ground_begin = 1
!    fates%ground_end = nlevgrnd
   
!    fates%sizepft_class_begin = 1
!    fates%sizepft_class_end = nlevsclass_fates * numpft_fates
   
!    fates%size_class_begin = 1
!    fates%size_class_end = nlevsclass_fates

!    fates%pft_class_begin = 1
!    fates%pft_class_end = numpft_fates

!    fates%age_class_begin = 1
!    fates%age_class_end = nlevage_fates

!    fates%sizeage_class_begin = 1
!    fates%sizeage_class_end   = nlevsclass_fates * nlevage_fates
   
!    fates%fuel_begin = 1
!    fates%fuel_end = nfsc_fates
   
!    fates%cwdsc_begin = 1
!    fates%cwdsc_end = ncwd_fates
   
!    fates%can_begin = 1
!    fates%can_end = nclmax_fates
   
!    fates%cnlf_begin = 1
!    fates%cnlf_end = nlevleaf_fates * nclmax_fates
   
!    fates%cnlfpft_begin = 1
!    fates%cnlfpft_end = nlevleaf_fates * nclmax_fates * numpft_fates

!    fates%height_begin = 1
!    fates%height_end = nlevheight_fates

!    fates%agepft_class_begin = 1
!    fates%agepft_class_end   = nlevage_fates * numpft_fates
   
!    fates%sizeagepft_class_begin = 1
!    fates%sizeagepft_class_end   = nlevsclass_fates * nlevage_fates * numpft_fates

   
!  end subroutine hlm_bounds_to_fates_bounds

  real(r8) FUNCTION shr_orb_cosz(jday,lat,lon)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the cosine of the solar zenith angle.
   ! Assumes 365.0 days/year.
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: Brian Kauffman
   ! Date:            Jan/98
   ! History:         adapted from statement FUNCTION in share/orb_cosz.h
   !
   !----------------------------------------------------------------------------

   real   (r8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
   real   (r8),intent(in) :: lat    ! Centered latitude (degrees, -90,...,90)
   real   (r8),intent(in) :: lon    ! Centered longitude (degrees, 0, ...,360)
   real   (r8) :: latr    ! Centered latitude (radians)
   real   (r8) :: lonr    ! Centered longitude (radians)
   real   (r8) :: declin ! Solar declination (radians)
   real   (r8), parameter :: dg2radian =0.0174533_r8 !degree to radian converion factor
   !----------------------------------------------------------------------------
   latr = lat * dg2radian
   lonr = lon * dg2radian
   declin = -23.44_r8*dg2radian*cos(360/365*(jday+10)*dg2radian)  
   shr_orb_cosz = sin(latr)*sin(declin) - &
   &              cos(latr)*cos(declin)*cos(jday*2.0_r8*3.1415_r8 + lonr)

  END FUNCTION shr_orb_cosz

end module ATSFatesInterfaceMod
