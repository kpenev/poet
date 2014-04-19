Usage {#usage}
=====
SubPixPhot [options]

Supported Options
-----------------
<dl><dt>-K, --low-mass-K=<\<double\>></dt>
 <dd>The wind strength for low mass stars in units of Msun*Rsun^2*day^2/(rad^2*Gyr). To identify this quantity in --input-columns use '__K__' (identical to the --high-mass-K since only one K can be applicable to a run). Default: 0.17.</dd>

<dt>--high-mass-K=<\<double\>></dt>
 <dd>The wind strength for high mass stars in units of Msun*Rsun^2*day^2/(rad^2*Gyr). To identify this quantity in --input-columns use '__K__' (identical to the --low-mass-K since only one K can be applicable to a run). Default: 0.17.</dd>

<dt>--low-mass-wind-sat-w=<\<double\>></dt>
 <dd>The frequency at which the wind saturates for low mass stars in units of rad/day. To identify this quantity in --input-columns use '__wsat__' (identical to the --high-mass-wind-sat-w since only one saturation frequency can be applicable to a run). Default: 2.45.</dd>

<dt>--high-mass-wind-sat-w=<\<double\>></dt>
 <dd>The frequency at which the wind saturates for high mass stars in units of rad/day. To identify this quantity in --input-columns use '__wsat__' (identical to the --high-mass-wind-sat-w since only one saturation frequency can be applicable to a run). Default: 0.</dd>

<dt>--low-mass-wind-sat-p=<\<double\>></dt>
 <dd>The period at which the wind saturates for low mass stars in days. This is an alternative to --low-mass-wind-sat-w. If both options are specified, --low-mass-wind-sat-w is used. In --input-columns identified by '__psat__'.</dd>

<dt>--high-mass-wind-sat-p=<\<double\>></dt>
 <dd>The period at which the wind saturates for high mass stars in days. This is an alternative to --high-mass-wind-sat-w. If both options are specified, --high-mass-wind-sat-w is used. In --input-columns identified by '__psat__'.</dd>

<dt>--core-env-coupling-timescale=<\<double\>></dt>
 <dd>The timescale on which the core end envelope are coupled in Myr. In --input-columns identified by '__tcoup__'. Default: 55.</dd>

<dt>--lgQ=<\<double\>></dt>
 <dd>Log base 10 of the tidal quality factor of the star. In --input-columns identified by __lgQ__'. Default: 6.</dd>

<dt>-M, --Mstar=<\<double\>></dt>
 <dd>Mass of the star in solar masses. In --input-columns identified by '__M__'. Default: 1</dd>

<dt>-m, --Mplanet=<\<double\>></dt>
 <dd>Mass of the planet in jupiter masses. In --input-columns identified by '__m__'. Default: 1.</dd>

<dt>-r, --Rplanet=<\<double\>></dt>
 <dd>Radius of the planet in jupiter radii. In --input-columns identified by '__r__'. Default: 1.</dd>

<dt>--planet-formation-age=<\<double\>></dt>
 <dd>The age at which the planet forms. If it is smaller than the disk dissipation age, the planet forms at the disk dissipation age. In --input-columns identifid by '__tform__'. Default 0.</dd>

<dt>--w-disk=<\<double\>></dt>
 <dd>The spin frequency at which the star is locked while the disk is present in rad/day. In --input-columns identifid by '__wdisk__'. Default: 0.898.</dd>

<dt>--p-disk=<\<double\>></dt>
 <dd>The spin period at which the star is locked while the disk is present in days. This is an alternative to --w-disk. If both options are specified, --w-disk is used. In --input-columns identified by '__pdisk__'.</dd>

<dt>--t-disk=<\<double\>></dt>
 <dd>The age at which the disk dissipates in Myr. In --input-columns identified by '__tdisk__'. Default: 5.</dd>

<dt>-a, --init-semimajor=<\<double\>></dt>
 <dd>The semimajor axis at which the planet first appears in AU. In --input-columns identified by '__aform__'. Default: 0.05.</dd>

<dt>-p, --init-orbit-period=<\<double\>></dt>
 <dd>The orbital period at which the planet first appears in days. This is an alternative to --init-semimajor. In--input-columns identified by '__pform__'.</dd>

<dt>--t0=<\<double\>></dt>
 <dd>The starting age for the evolution in Gyr. If this argument is used, --w-disk, --w-disk-col, --t-disk and --t-disk-col are ignored and --init-semimajor or --init-semimajor-col is the semimajor axis that the planet has at the age specified by this argument. Further, if this argument is  used, --init-wsurf or --init-wsurf-col must also be specified and if the evolution of low mass stars is to be computed --init-wrad or --init-wrad-col is also required. In --input-columns identified by '__t0__'.</dd>

<dt>--tmax=<\<double\>></dt>
 <dd>The maximum end age for the evolution in Gyr. If a negative value is passed, the evolution stops at absolute value  of this parameter of if the star moves off the main sequence In --input-columns identified by '__tmax__'. Default: -inf.</dd>

<dt>--init-wrad=<\<double\>></dt>
 <dd>Initial spin of the stellar core in rad/day for low mass stars. This argument is ignored, unless --t0 or --t0-col is also specified. In --input-columns identified by '__wrad0__'.</dd>

<dt>--init-wsurf=<\<double\>></dt>
 <dd>Initial spin of the stellar surface in rad/day. For low mass stars this is the rotation of the convective zone, and for high mass stars it is the unique rotation of the star. This argument is ignored, unless --t0 or --t0-col is also specified. In --input-columns identified by '__wsurf0__'.</dd>

<dt>--max-step=<\<double\>></dt>
 <dd>An upper limit to impose on the sover timestep. In --input-columns identified by '__maxdt__'. Default: inf.</dd>

<dt>--precision=<\<double\>></dt>
 <dd>The number of digits of precision to require of the solution. Need not be an integer. In --input-columns identified by '__prec__'. Default: 6.</dd>

<dt>-o, --output=<\<file\>|index></dt>
 <dd>If nothing is read from an input \<file\>, then this argument should be the name of the file to use to output the evolution. If at least one quantity is read from a list file, that file should also contain a column (named '__outf__') specifying the input file name for each parameter set, and this column should be specified via the --input-columns option. In the latter case, this argument is ignored. Default: poet.evol</dd>

<dt>--start-locked</dt>
 <dd>Whether the planet should start locked to the star. If true, causes the value of --init-wsurf-col to be ignored.</dd>

<dt>--require-ages=<\<comma separated list\>></dt>
 <dd>A list of ages to include in the calculated evolution. This argument can be overwritten or appended to by the '__trequired__' column in the input \<file\>. If the entry in the file starts with a comma, the ages listed there are combined with the value specified by this argument, otherwise only the ages listed in the file are used.</dd>

<dt>--input-columns=<\<comma separated list\>></dt>
 <dd>Allows multiple evolutions to be calculated with a single run. This argument should be a list of the columns in an input \<file\> that includes all quantities that change, and also a column specifying the file to output the corresponding evolution to. Any argument not included in this list of columns is held constant at the value specified at the corresponding option value or its default. Columns that should be skipped should have no identifier listed (i.e. there should be two consecutive commas). Not all columns in the input file need to be listed. Any content past the last column identified in this argument on each line is ignored.</dd>

<dt>--output-columns=<\<comma separated list\>></dt>
 <dd>Specifies what to output. This argument should be a \<comma separated list\> containing some of the following columns:
	* __t__: Age in Gyr
	* __Iconv__: Convective zone moment of inertia Msun Rsun^2 (low mass stars only else NaN).
	* __Irad__: Radiative zone moment of inertia Msun Rsun^2 (low mass stars only else NaN).
	* __I__: Total moment of inertia of the star Msun Rsun^2.
	* __R__: Stellar radius in solar radii
	* __Lum__: Stellar luminosity in solar luminosities
	* __Rrad__: Radius of the stellar radiative zone in solar radii (low mass stars only, else NaN).
	* __Mrad__: Mass of the stellar radiative zone in solar masses (low mass stars only, else NaN).
	* __DIconv__: First order derivative with respect to age of the moment of inertia of the stallar surface convective zone in Msun Rsun^2/Gyr (low mass stars only, else NaN).
	* __DIrad__: First order derivative with respect to age of the moment of inertia of the stallar radiative core in Msun Rsun^2/Gyr (low mass stars only, else NaN).
	* __DI__: First order derivative with respect to age of the moment of inertia of the entire star in Msun Rsun^2/Gyr.
	* __DR__: First order derivative with respect to age of the stellar radius in solar radii per Gyr.
	* __DRrad__: First order derivative with respect to age of the radius of the stellar radiative core in solar radii per Gyr (low mass stars only, else NaN).
	* __DMrad__: 
	* __DDIconv__: Second order derivative with respect to age of the moment of inertia of the stallar surface convective zone in Msun Rsun^2/Gyr^2 (low mass stars only, else NaN).
	* __DDIrad__: Second order derivative with respect to age of the moment of inertia of the stallar radiative core in Msun Rsun^2/Gyr^2 (low mass stars only, else NaN).
	* __DDI__: Second order derivative with respect to age of the moment of inertia of the entire star in Msun Rsun^2/Gyr^2
	* __DDRrad__: First order derivative with respect to age of the radius of the stellar radiative core in solar radii per Gyr^2.
	* __a__: Semimajor axis in AU
	* __Worb__: Orbital frequency in rad/day
	* __Porb__: Orbital period in days
	* __Lconv__: Convective zone angular momentum in Msun Rsun^2 rad/day (low mass stars only else NaN).
	* __Lrad__: Radiative zone angular momentum in Msun Rsun^2 rad/day (low mass stars only else NaN).
	* __L__: Total angular momentum of the star in Msun Rsun^2 rad/day.
	* __Wsurf__: Stellar surface angular velocity of the star in rad/day.
	* __Wrad__: Angular velocity of the radiative zone in rad/day (low mass stars only, else NaN)
	* __Psurf__: Stellar surface spin period in days.
	* __Prad__: Spin period of the radiative zone in days (low mass stars only, else NaN)
	* __mode__: Evolution mode of the step starting at this age.
	* __wind__: The wind state (saturated/not saturated) of the step starting at this age.
Default: t,a,Lconv,Lrad,L,Iconv,Irad,I,mode</dd>

<dt>-i, --input=<\<file\>></dt>
 <dd>The \<file\> to read the parameters of the planet-star systems to calculate evolutions for. If omitted, standard input is used instead. Any lines beginning with '#' are ignored. The other lines should start with at least the number of columns specified in the --input-columns option with white space only between columns.</dd>

<dt>--serialized-stellar-evol=<\<file\>></dt>
 <dd>The \<file\> to read previously serialized stellar evolution from or write one if the file does not exist. Default: './serialized_evolution'.</dd>

<dt>--custom-stellar-evolution=<\<file\>></dt>
 <dd>A single stellar evolution track from which to construct the stellar evolution to use (assumed to apply to all input systems regardless of the stellar mass). It should contain the quantities specified by --custom-stellar-evolution-format as columns in the precise order specified there. If this option is used, all input stars are treated as low mass stars, since that way the high mass star behavior can be reproduced  by simply assigning the entirestar to a surface  convective zone.</dd>

<dt>--custom-stellar-evolution-format=<col1,col2,...></dt>
 <dd>A \<comma separated list\> of the columns in the custom stellar evolution track specified by the --custom-stellar-evolution option. The recognized column names are:
	* __Iconv__: Convective zone moment of inertia in Msun Rsun^2
	* __Irad__: Radiative zone moment of inertia in Msun Rsun^2 
	* __R__: Stellar radius in solar radii
	* __Rrad__: Radius of the radiative zone in solar radii.
	* __Mrad__: Mass of the radiative zone in solar masses.
	* __Lum__: Stellar luminosity in solar luminosities
	* __t__: Stellar age in Gyr
</dd>

<dt>--Iconv-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Convective zone moment of inertia when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--Iconv-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Convective zone moment of inertia of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --Iconv-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>--Irad-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Radiative zone moment of inertia when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--Irad-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Radiative zone moment of inertia of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --Irad-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>--R-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Stellar radius when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--R-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Stellar radius of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --R-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>--Rrad-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Radius of the radiative zone when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--Rrad-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Radius of the radiative zone of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --Rrad-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>--Mrad-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Mass of the radiative zone when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--Mrad-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Mass of the radiative zone of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --Mrad-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>--Lum-smoothing=<\<double\>></dt>
 <dd>Smoothing to apply to the Stellar luminosity when interpolating the custom stellar evolution track. If no smoothing should be applied use 'NaN'. Default: nanThis option is ignored, unless --custom-stellar-evolution is specified.</dd>

<dt>--Lum-nodes=<integer></dt>
 <dd>The number of nodes to use when smoothing the Stellar luminosity of the custom stellar evolution track. This value is ignored unless both --custom-stellar-evolution is used and --Lum-smoothing is not NaN. Negative values result in using a number of nodes which is the smaller of the absolute of the given value and three times the number of age entries in the track. Default: -1000.</dd>

<dt>-h, --help</dt>
 <dd>Print this help and exit.</dd>

<dt>--doxygen-help</dt>
 <dd>Print this help formatted suitably for Doxygen and exit.</dd>

</dl>