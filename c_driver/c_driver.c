/* *******************************************

This is a model problem which simulates a 1D column of soils.
The column has nlevel=5 soil layears and nlevdecomp=1 decomp layers.

************************************************************ */

#include "ISO_Fortran_binding.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mathimf.h>




typedef struct {  
  int nlevbed;
  int nlevdecomp;
  int patchno;
  int altmax_lastyear_indx_col;
  double temp_veg24_patch;
  double latdeg, londeg;
} site_info;


typedef struct {  
  double dayl_factor;  // scalar (0-1) for daylength
  double esat_tv;      // saturation vapor pressure at t_veg (Pa)
  double eair;         // vapor pressure of canopy air (Pa)
  double oair;         // Atmospheric O2 partial pressure (Pa)
  double cair;         // Atmospheric CO2 partial pressure (Pa)
  double rb;           // boundary layer resistance (s/m)
  double t_veg;        // vegetation temperature (Kelvin)
  double tgcm;         // air temperature at agcm reference height (Kelvin)
  double solad[2]; //direct radiation (W/m**2); 1=visible lights; 2=near infrared radition
  double solai[2]; //diffuse radiation (W/m**2); 1=visible lights; 2=near infrared radition
  double albgrd[2]; //!ground albedo (direct) 1=visiable; 2=near infrared (nir)
  double albgri[2]; //ground albedo (diffuse) 1=visiable; 2=near infrared (nir)
} PhotoSynthesisInput;

void init_ats_fates(int*, site_info*);
void init_soil_depths(int*, int*, site_info*, double*, double*, double*, double*);
void init_coldstart(int* );
extern void fatessetmasterproc(int*);
extern void fatessetinputfiles(CFI_cdesc_t * clm, CFI_cdesc_t * fates);
extern void fatesreadparameters();
extern void fatesreadpfts();
extern void set_fates_global_elements();
extern void get_nlevsclass(int*);
extern void wrap_btran(int*, double*, double*, double*, double*, double*);
extern void wrap_photosynthesis(double*, double*, int*, double*, PhotoSynthesisInput*);
  

extern void dynamics_driv_per_site(int*, int*, site_info*, double*,
                            double*, double*, double*, double*, double*);




main()
{

  int iulog=10;
  int masterproc = 1;

  CFI_cdesc_t fatesdesc, clmdesc;
  int retval;
 
  
  char *fates_file = "/turquoise/usr/projects/veg/cxu/fates-ats2/c_driver/fates_params_default_c20191007.nc";
  char *clm_file = "/turquoise/usr/projects/veg/cxu/fates-ats2/c_driver/clm_params_c180301.nc";

  retval = CFI_establish(&fatesdesc, fates_file, CFI_attribute_other, CFI_type_char, strlen(fates_file), 0, NULL);
  retval = CFI_establish(&clmdesc, clm_file, CFI_attribute_other, CFI_type_char, strlen(clm_file), 0, NULL); 
  fatessetmasterproc(&masterproc);
  fatessetinputfiles(&clmdesc, &fatesdesc);    
  fatesreadparameters();
    
 
  /* Initialization */
  /*  ! Read in FATES parameter values early in the call sequence as well
    ! The PFT file, specifically, will dictate how many pfts are used
    ! in fates, and this will influence the amount of memory we
    ! request from the model, which is relevant in set_fates_global_elements()*/
  
  fatesreadpfts();

  /*   ------------------------------------------------------------------------
     Ask Fates to evaluate its own dimensioning needs.
     This determines the total amount of space it requires in its largest
     dimension.  We are currently calling that the "cohort" dimension, but
     it is really a utility dimension that captures the models largest
     size need.
     Sets:
     fates_maxElementsPerPatch
     fates_maxElementsPerSite (where a site is roughly equivalent to a column)
     
     (Note: fates_maxELementsPerSite is the critical variable used by CLM
     to allocate space)
     ------------------------------------------------------------------------*/

  set_fates_global_elements();


  int nlevelclass;
  get_nlevsclass(&nlevelclass);
  
  site_info site[10];
  int num_sites = 1;
  int clump = 1;
  int nlevel = 5;
  
  /* Define site information */
  for (int i=0; i<num_sites; i++){ 
    site[i].nlevbed = nlevel;
    site[i].nlevdecomp = 1;
    site[i].patchno = 1;
    site[i].temp_veg24_patch = 283;
    site[i].altmax_lastyear_indx_col = 1;
    site[i].latdeg = 30.;
    site[i].londeg = 30.;
  }

  /* Preliminary initialization of FATES */
  init_ats_fates(&num_sites, site);

  double zi[6], z[5], dz[5], dzsoil_decomp[5];
  /* Define soil layers */
  zi[0] = 0;
  for (int i=0; i<nlevel; i++){
    zi[i+1] = zi[i] + 1;
    z[i] = 0.5*(zi[i+1] + zi[i]);
    dz[i] = 1;
  }  
  dzsoil_decomp[0] = 1;
  
  /* Initialize soil layers in FATES*/
  for (int i=0; i<num_sites; i++){
    int s = i+1;    
    init_soil_depths(&clump, &s,  &(site[i]), zi, dz, z, dzsoil_decomp);
  }


  /* Init cold start of FATES */
  init_coldstart(&clump);
   

  /* One time step */
  

  double h2osoi_vol_col[nlevel];

  for (int i=0;i<nlevel;i++) h2osoi_vol_col[i]=1;
  
  double temp_veg24_patch[site[0].patchno],
         prec24_patch[site[0].patchno],
         rh24_patch[site[0].patchno],
         wind24_patch[site[0].patchno];

  double dtime = 1800.;
  int array_size = nlevel * num_sites;
  double t_soil[array_size], poro[array_size], eff_poro[array_size];
  double sat[array_size], soil_suc[array_size];
  double bsw = 0.5;
  double soil_suc_min = 10;
  PhotoSynthesisInput photosys_in;


  double p_atm = 101325;
  double Time = 0;
  int   radnum = 2; //number of radiation bands
  double jday; //julian days (1-365)

  for (int n=0; n<48; n++){
    for (int i=0;i<array_size;i++) t_soil[i] = 290;
    for (int i=0;i<array_size;i++) poro[i] = 0.7;
    for (int i=0;i<array_size;i++) eff_poro[i] = 0.7;
    for (int i=0;i<array_size;i++) sat[i] = 1.;
    
    for (int i=0;i<array_size;i++){
      soil_suc[i] = -soil_suc_min * pow(sat[i], bsw);
    }

    wrap_btran(&array_size, t_soil, sat, eff_poro, poro, soil_suc);

    
    photosys_in.dayl_factor = 0.7;
    photosys_in.esat_tv = 2300.;     // Saturated vapor pressure in leaves (Pa)
    photosys_in.eair = 2000.;        // Air water vapor pressure (Pa)
    photosys_in.oair = 21280;        // Oxygen partial pressure
    photosys_in.cair = 5985;       // CO2 partial pressure
    photosys_in.rb = 3.;            // Boundary layer resistance (s/m)
    photosys_in.t_veg = 305;        // Leaf temperature (K)
    photosys_in.tgcm = 305;         // Air temperature (K)
    photosys_in.albgrd[0] = 0.15;
    photosys_in.albgrd[1] = 0.15;
    photosys_in.albgri[0] = 0.1;
    photosys_in.albgri[1] = 0.1;
    photosys_in.solad[0] = 200.0;
    photosys_in.solad[1] = 40.0;
    photosys_in.solai[0] = 50.0;
    photosys_in.solai[1] = 10.0;    
    jday = 1.0+n/48.0;

    wrap_sunfrac(&radnum, photosys_in.solad, photosys_in.solai);

    wrap_canopy_radiation(&jday,&radnum, photosys_in.albgrd, photosys_in.albgri);

    wrap_photosynthesis(&dtime, &p_atm, &array_size, t_soil, &photosys_in);

    Time += dtime;
                        
  }

  
  double day_sec = 86400.;  
  for (int i=0; i<num_sites; i++){
    int s = i+1;

    temp_veg24_patch[0] = 305;
    prec24_patch[0] = 3.472222e-05;
    rh24_patch[0] = 80;
    wind24_patch[0] = 2.0;
    
    dynamics_driv_per_site( &clump,  &s,  &(site[i]), &dtime,
                            h2osoi_vol_col, temp_veg24_patch,
                            prec24_patch, rh24_patch, wind24_patch);
  }

  printf("Finished successfully.\n");
  
  return 0;
}
