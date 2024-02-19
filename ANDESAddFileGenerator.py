import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class ANDESSimulator:
    # Static attributes
    centers_of_bands = np.array([445, 551, 1220, 1630, 2190])
    magnitudes_AB_conv = np.array([0.09, 0.55, 0.91, 1.39, 1.85])

    def __init__(self, resolving_power=100_000, wavelength_start=350, number_of_orders=127):
        self.resolving_power = resolving_power
        self.wavelength_start = wavelength_start
        self.number_of_orders = number_of_orders
        self.wavelengths = self.generate_wavelengths()

    def generate_wavelengths(self):
        wavelengths = [self.wavelength_start]
        while wavelengths[-1] < 2400:
            wavelengths.append(wavelengths[-1] + wavelengths[-1] / self.resolving_power)
        return np.asarray(wavelengths)

    def generate_order_definition_file(self):
        len_order = int(len(self.wavelengths) / self.number_of_orders)
        wave_orders = np.zeros((self.number_of_orders, len_order))
        for i in range(self.number_of_orders):
            wave_orders[i] = self.wavelengths[i * len_order:(i + 1) * (len_order)]
        np.save('order_def.npy', wave_orders)

    def generate_snr_file(self, magnitudes_vega, SNR_ANDES):
        magnitudes_AB = magnitudes_vega + self.magnitudes_AB_conv
        interpolated_SNR = interp1d(self.centers_of_bands, SNR_ANDES, bounds_error=False,
                                    fill_value=(SNR_ANDES[0], SNR_ANDES[-1]), kind='linear')(self.wavelengths)
        np.save('SNR.npy', interpolated_SNR)

        plt.plot(self.centers_of_bands, SNR_ANDES)
        plt.plot(self.wavelengths, interpolated_SNR)
        plt.show()


# Example usage:
if __name__ == "__main__":
    magnitudes_vega = np.array([10.35, 9.82, 8.722, 8.470, 8.429])
    SNR_ANDES = np.array([379.7, 487.6, 843.8, 806.1, 689.3])

    simulator = ANDESSimulator()
    simulator.generate_order_definition_file()
    simulator.generate_snr_file(magnitudes_vega, SNR_ANDES)
