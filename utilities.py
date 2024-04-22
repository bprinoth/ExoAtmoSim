__all__ = [
    'paramget',
    'BB', 
    'generate_config'
    
    
]
######### HELPER FUNCTIONS - DEFINED OUTSIDE
# # paramget allows you to get parameters from a config file instead of typing all the stuff into a list

def paramget(keyword,dp,full_path=False,force_string = False):
    """This code queries a planet system parameter from a config file located in the folder
    specified by the path dp; or run configuration parameters from a file speciefied by the full
    path dp, if full_path is set to True. It is taken from tayph @Hoeijmakers

    Parameters
    ----------
    keyword : str
        A keyword present in the cofig file.

    dp : str, Path
        Output filename/path.

    full_path: bool
        If set, dp refers to the actual file, not the location of a folder with a config.dat;
        but the actual file itself.

    Returns
    -------
    value : int, float, bool, str
        The value corresponding to the requested keyword.

    """
    import pathlib
    import distutils.util
    import pdb


    #This is probably the only case where I don't need obs_times and config to exist together...
    #dp=check_path(dp)
    #typetest(keyword,str,'keyword in paramget()')

    if isinstance(dp,str) == True:
        dp=pathlib.Path(dp)
    try:
        if full_path == False:
            dp = dp/'config'
        f = open(dp, 'r')
    except FileNotFoundError:
        raise FileNotFoundError('Parameter file does not exist at %s' % str(dp)) from None
    
    x = f.read().splitlines()
    f.close()
    n_lines=len(x)
    keywords={}
    
    for i in range(0,n_lines):
        line=x[i].split()
        if len(line) > 1:
            if force_string:
                value=(line[1])
            else:
                try:
                    value=float(line[1])
                except ValueError:
                    try:
                        value=bool(distutils.util.strtobool(line[1]))
                    except ValueError:
                        value=(line[1])
            keywords[line[0]] = value
    try:
        return(keywords[keyword])
    
    except KeyError:
        # print(keywords)
        raise Exception('Keyword %s is not present in parameter file at %s' % (keyword,dp)) from None

def BB(temperature): # BB function from frei @ Brett Morris.
    import numpy as np
    from astropy import constants as const
    
    """
    Compute the blackbody flux

    Parameters
    ----------
    temperature : ~astropy.units.Quantity
        Temperature of the blackbody

    Returns
    -------
    bb : function
        Returns a function which takes the wavelength as an
        `~astropy.units.Quantity` and returns the Planck flux
    """
    # h = 6.62607015e-34  # J s
    # c = 299792458.0  # m/s
    # k_B = 1.380649e-23  # J/K
    return lambda wavelength: (
            2 * const.h * const.c ** 2 / np.power(wavelength, 5) /
            np.expm1(const.h * const.c / (wavelength * const.k_B * temperature))
    )
    
    
def generate_config(planetname, resolution=140_000, air=True, print_params=True, savepath=''):
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    planet_data = NasaExoplanetArchive.query_criteria(table='pscomppars', where=f"pl_name='{planetname}'")
    
    if len(planet_data) != 1:
        print(f'[ERROR] Uh, something is wrong. I retrieved {len(planet_data)} planets. Check your spelling please.')
    
    import astropy.units as u
    import astropy.constants as const
    import numpy as np
    
    # P, a, aRstar, Rp, Mp, K, RpRstar, vsys, RA, DEC, Tc, vsini
    P = planet_data['pl_orbper']
    a = planet_data['pl_orbsmax']
    aRstar = planet_data['pl_ratdor']    
    Rp = planet_data['pl_radj']
    Mp = planet_data['pl_bmassj']
    K = (planet_data['pl_rvamp']).to("km/s")
    vsini = planet_data['st_vsin']
    RpRstar = planet_data['pl_ratror']
    vsys = planet_data['st_radv']
    Tc = planet_data['pl_tranmid']
    T14 = planet_data['pl_trandur']
    inclination = planet_data['pl_orbincl']
    Rs = planet_data['st_rad']
    Ms = planet_data['st_mass']

    loggstar = planet_data['st_logg']
    metallicity = planet_data['st_met']
    
    #print(Rs)
    #print(planet_data.columns)
    omega = planet_data['pl_orblper']
    lampoo = planet_data['pl_projobliq']
    
    # position in the sky
    RA = planet_data['ra']
    DEC = planet_data['dec']
    
    teq = planet_data['pl_eqt']
    teff = planet_data['st_teff']
    
    ecc = planet_data['pl_orbeccen']
    impact = planet_data['pl_imppar']
    
    loggplanet = np.log10((const.G * Mp / Rp**2).to('cm/s2').value)
    
    
    if print_params:
        print(f"Parameters for config of {planetname}")

        print('\tP', P[0].value, P[0].unit)
        print('\ta', a[0].value, a[0].unit)
        print('\taRs', aRstar[0])
        print('\tRp', Rp[0].value, Rp[0].unit)
        print('\tMp', Mp[0].value, Mp[0].unit)
        print('\tK', K[0].value, K[0].unit)
        print('\tvsini', vsini[0].value, vsini[0].unit)
        print('\tRpRstar', RpRstar[0])
        print('\tvsys', vsys[0].value, vsys[0].unit)
        print('\tinclination', inclination[0].value, inclination[0].unit)
        print('\tTc', Tc[0].value, Tc[0].unit)
        print('\tRA', RA[0].value)
        print('\tDEC', DEC[0].value)
        print('\tteq', teq[0].value, teq[0].unit)
        print('\tteff', teff[0].value, teff[0].unit)
        print('\tresolution', resolution)
        print('\tair', air)
        print('\tecc', ecc[0])
        print('\timpact', impact[0])
        print('\tomega', omega[0])
        print('\tloggstar', loggstar[0])
        print('\tmet', metallicity[0])
        print('\tloggp', loggplanet[0])
        print('\tlampoo', lampoo[0].value, lampoo[0].unit)
        #print('\tecc', ecc[0])
    
    print(f"[INFO] Generating config file for {planetname}")
    
    lines = [
        f"P\t{P[0].value}",
        f"a\t{a[0].value}",
        f"aRstar\t{aRstar[0]}",
        f"Rp\t{Rp[0].value}",
        f"Mp\t{Mp[0].value}",
        f"K\t{K[0].value}",
        f"vsini\t{vsini[0].value}",
        f"RpRstar\t{RpRstar[0]}",
        f"vsys\t{vsys[0].value}",
        f"inclination\t{inclination[0].value}",
        f"Tc\t{Tc[0].value}",
        f"RA\t{RA[0].value}",
        f"DEC\t{DEC[0].value}",
        f"Teq\t{teq[0].value}",
        f"Teff\t{teff[0].value}",
        f"resolution\t{resolution}",
        f"air\t{air}",
        f"ecc\t{ecc[0]}",
        f"impact\t{impact[0]}",
        f"omega\t{omega[0].value}",
        f"Rs\t{Rs[0].value}",
        f"Ms\t{Ms[0].value}",
        f"logg\t{loggstar[0]}",
        f"metallicity\t{metallicity[0].value}",
        f"loggp\t{loggplanet[0]}",
        f"lampoo\t{lampoo[0].value}"
        
    ]
    
    with open(f"{savepath}/config", 'w') as f:
        f.write('\n'.join(lines))
    
    return planet_data