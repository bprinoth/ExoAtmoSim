import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# WAVELENGTH REGIONS OF THE CRIRES+ BANDS

V = [[ 356.081 , 376.5509999999814 ],
[ 376.5609999999814 , 397.0309999999628 ],
[ 397.04099999996276 , 417.51099999994415 ],
[ 417.52099999994414 , 437.9909999999255 ],
[ 438.0009999999255 , 458.4709999999069 ],
[ 458.4809999999069 , 478.95099999988827 ],
[ 478.96099999988826 , 499.43099999986964 ],
[ 499.44099999986963 , 519.910999999851 ],
[ 519.920999999851 , 540.3909999998324 ],
[ 540.4009999998324 , 560.8709999998138 ],
[ 560.8809999998138 , 581.3509999997951 ],
[ 581.3609999997951 , 601.8309999997765 ],
[ 601.8409999997765 , 622.3109999997579 ],
[ 622.3209999997579 , 642.7909999997393 ],
[ 642.8009999997392 , 663.2709999997206 ],
[ 663.2809999997206 , 683.750999999702 ],
[ 683.760999999702 , 704.2309999996834 ],
[ 704.2409999996834 , 724.7109999996648 ],
[ 724.7209999996647 , 745.1909999996461 ],
[ 745.2009999996461 , 765.6709999996275 ],
[ 765.6809999996275 , 786.1509999996089 ],
[ 786.1609999996089 , 806.6309999995902 ],
[ 806.6409999995902 , 827.1109999995716 ],
[ 827.1209999995716 , 847.590999999553 ],
[ 847.600999999553 , 868.0709999995344 ],
[ 868.0809999995344 , 888.5509999995157 ],
[ 888.5609999995157 , 909.0309999994971 ],
[ 909.0409999994971 , 929.5109999994785 ],
[ 929.5209999994785 , 949.9909999994599 ]]

# this is just an estimate to close up to the Y band!

# Y band

Y = [[949.898 , 968.753],
[965.933 , 985.463],
[982.886 , 1002.754],
[1000.440 ,  1020.655],
[1018.747 ,  1039.200],
[1037.493 ,  1058.426],
[1057.032 ,  1078.362],
[1077.326 ,  1099.053],
[1098.400 ,  1120.709]]

# J band

J = [[1122.121 , 1144.645],
[1145.040 , 1168.019],
[1168.906 , 1192.356],
[1193.779 , 1217.720],
[1219.734 , 1244.182],
[1246.813 , 1271.600],
[1275.136 , 1300.661],
[1304.738 , 1330.921]]

# H band

H = [[1445.816 , 1474.143],
[1483.886 , 1512.948],
[1523.999 , 1554.147],
[1566.338 , 1596.995],
[1611.072 , 1642.598],
[1658.417 , 1690.860],
[1708.663 , 1742.235],
[1761.918 , 1796.339],
[1818.644 , 1854.250]]

# K band 
K = [[1968.499 , 2004.852],
[2038.829 , 2076.268],
[2114.306 , 2153.390],
[2195.651 , 2235.941],
[2283.444 , 2325.318],
[2378.509 , 2421.842],
[2481.614 , 2527.165]]


dict_order_length = {
    'B+V': len(V),
    'Y': len(Y),
    'J': len(J),
    'H': len(H),
    'K': len(K)
}

centers_of_bands = np.array([
            445,  # B
            551,  # V
            1220, # J
            1630, # H
            2190  # K
        ]) # in nm


dict_get_band = {
    'B+V': V,
    'Y': Y,
    'J': J, 
    'H', H,
    'K', K
    
}


magnitudes_AB_conv = np.array([0.09, 0.55, 0.91, 1.39, 1.85]) # conversion factor between magnitude systems


class ANDESSimulator:
    
    def __init__(self, bands=['BV', 'Y', 'J', 'H', 'K'], sampling=0.01):
        self.bands = bands 

        # Let's check here if the bands are well defined:
        for band in self.bands:
            if band is not in ['B+V', 'Y', 'J', 'H', 'K']:
                print(f"[ERROR] Sorry, but I don't know the {band} band.")
        
        self.sampling = 0.01 # default for CRIRES+ or ANDES!
        self.wavelengths = self.generate_wavelengths()
        
    def generate_order_def(self):

        self.total_number_of_orders = 0
        for band in self.bands:
            total_number_of_orders += dict_order_length[band]

        self.wavelengths = np.zeros((self.total_number_of_orders, 2048))

        len_prv = 0
        
        for band in self.bands:

            wavelengths_band = np.zeros((dict_order_length[band], 2048)) # all CRIRES+ orders have the same shape + I forced the B+V order to be the same!
            
            # Now we generate the wavelength according to the band and add it to the wavelengths of the band
            
            for i in range(len(wavelengths_band)):
                band_ranges = dict_get_band[band]
                wavelengths_band[i] = np.linspace(band_ranges[i][0], band_ranges[i][0], 2048)            
            
            self.wavelengths[len_prv:dict_order_length[band]+len_prv] = wavelengths_band

            # overwrite the length of the previous band(s)
            len_prv = dict_order_length[band]
            
        np.save('order_def_ANDES.npy', self.wavelengths)
            

    def generate_snr_file(self, magnitudes_vega, SNR_ANDES):
        magnitudes_AB = magnitudes_vega + self.magnitudes_AB_conv
        interpolated_SNR = interp1d(self.centers_of_bands, SNR_ANDES, bounds_error=False,
                                    fill_value=(SNR_ANDES[0], SNR_ANDES[-1]), kind='linear')(self.wavelengths)
        
        np.save('SNR_ANDES.npy', interpolated_SNR)


# Example usage:
if __name__ == "__main__":
    magnitudes_vega = np.array([10.35, 9.82, 8.722, 8.470, 8.429])
    SNR_ANDES = np.array([379.7, 487.6, 843.8, 806.1, 689.3])

    simulator = ANDESSimulator()
    simulator.generate_order_definition_file()
    simulator.generate_snr_file(magnitudes_vega, SNR_ANDES)
