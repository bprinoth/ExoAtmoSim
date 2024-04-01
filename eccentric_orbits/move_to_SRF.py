import os
import numpy as np
import tayph.util as ut
import tayph.operations as ops
import astropy.constants as const
import astropy.units as u
from scipy.interpolate import interp1d

def move_to_SRF(dp, planet):
    dp = Path(dp)
    
    if not os.path.exists(f'{planet}_SRF'):
        os.mkdir(f'{planet}_SRF')
    
    # copy the orders into the SRF folder
    os.system(f"scp -pr {planet}/*.fits {planet}_SRF/")
    
    # lists for the data
    list_of_wls=[]
    list_of_orders=[]
    list_of_sigmas=[]

    # figuring out how many orders we have
    filelist_orders= [str(i) for i in Path(dp).glob('order_*.fits')]
    order_numbers = [int(i.split('order_')[1].split('.')[0]) for i in filelist_orders]
    order_numbers.sort()
    n_orders = len(order_numbers)

    # some corrections (E.g. negative values and NaNs)
    trigger2 = 0#These triggers are used to limit the generation of output in the forloop.
    trigger3 = 0
    n_negative_total = 0#This will hold the total number of pixels that were set to NaN because
    #they were zero when reading in the data.
    
    for i in order_numbers:
        wavepath = dp/f'wave_{i}.fits'
        orderpath= dp/f'order_{i}.fits'
        sigmapath= dp/f'sigma_{i}.fits'
           
        wave_order = ut.readfits(wavepath)#2D or 1D?
        order_i = ut.readfits(orderpath)

        list_of_wls.append(ops.airtovac(wave_order))
            

        #Now test for negatives, set them to NaN and track them.
        n_negative = len(order_i[order_i <= 0])
        if trigger3 == 0 and n_negative > 0:
            print("------Setting negative values to NaN.")
            trigger3 = -1
        n_negative_total+=n_negative
        order_i[order_i <= 0] = np.nan #This is very important for later when we are computing
        #average spectra and the like, to avoid divide-by-zero cases.
        postest(order_i,f'order {i} in run_instance().')#make sure whatever comes out here is
        #strictly positive.
        list_of_orders.append(order_i)

        #Try to get a sigma file. If it doesn't exist, we raise a warning. If it does, we test
        #its dimensions and append it.
        try:
            sigma_i = ut.readfits(sigmapath)
            dimtest(sigma_i,[n_exp,n_px],f'order {i} in run_instance().')
            list_of_sigmas.append(sigma_i)
        except FileNotFoundError:
            if trigger2 == 0:
                ut.tprint('------WARNING: Sigma (flux error) files not provided. '
                'Assuming sigma = sqrt(flux). This is standard practise for HARPS data, but '
                'e.g. ESPRESSO has a pipeline that computes standard errors on each pixel more '
                'accurately. The sqrt approximation is known to be inappropriate for data '
                'of CARMENES or GIANO, because these pipeline renormalise the spectra. Make '
                'sure that the errors are treated properly, otherwise the error-bars on the '
                'output CCFs will be invalid.')
                trigger2=-1
            list_of_sigmas.append(np.sqrt(order_i))
            
    
    list_of_orders_t, list_of_sigmas_t = list_of_orders, list_of_sigmas    

    # here comes the actual velocity correction!

    rv_corr = -PlanetSpectrum.RV_star_ecc(f'{dp}/') - sp.paramget('vsys', f'{dp}/')
    gamma = 1.0 + (rv_corr * u.km / u.s / const.c).decompose()
    
    list_of_shifted_orders = []

    # velocity shift
    for i in order_numbers:
        shifted_order = np.zeros_like(list_of_orders[i])
        
        for k in range(len(shifted_order)):
            
            shifted_order[k] = interp1d(list_of_wls[i][k]*gamma[k], list_of_orders_t[i][k], bounds_error=False)(list_of_wls[i][k])
    
        list_of_shifted_orders.append(shifted_order)
    
    # saving
    for i in order_numbers:
        ut.writefits(dp/f'order_{i}.fits', list_of_shifted_orders[i])
        ut.writefits(dp/f'wave_{i}.fits', ops.vactoair(list_of_wls[i]))
    
    return list_of_shifted_orders, list_of_orders#,list_of_sigmas_t,list_of_Ts_t
