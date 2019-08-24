# Synthetic turbulence interface

This utility provides interfaces to different synthetic turbulence models and to
export them to NetCDF format.

## NetCDF file format

```
netcdf turbulence {
types:
  byte enum wind_profile_t {Constant = 0, Power_Law = 1, Log_Law = 2} ;
dimensions:
	ndim = 3 ;
	nx = 256 ;
	ny = 256 ;
	nz = 256 ;
variables:
	wind_profile_t wind_profile ;
	double uref ;
		uref:units = "m/s" ;
	double href ;
		href:units = "m" ;
	double hoffset ;
		hoffset:units = "m" ;
	double roughness_height ;
		roughness_height:units = "m" ;
	int mean_wind_added ;
		mean_wind_added:comment = "1 = Mean wind added to uvel; 0 = No mean component added" ;
	double box_lengths(ndim) ;
		box_lengths:units = "m" ;
	double scale_factors(ndim) ;
	double uvel(nx, ny, nz) ;
		uvel:units = "m/s" ;
	double vvel(nx, ny, nz) ;
		vvel:units = "m/s" ;
	double wvel(nx, ny, nz) ;
		wvel:units = "m/s" ;

// global attributes:
		:title = "WindSim: Mann turbulence file" ;
}
```
