# Pressure
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

# Wild type pressure file
wild_type_data = np.loadtxt('/content/gdrive/MyDrive/1aki/npt_press_dens.xvg')

# Mutant pressure file
mutant_data = np.loadtxt('/content/gdrive/MyDrive/1aki/npt_press_dens.xvg')

plt.title('Pressure during 1000 ps Equilibration (NPT)')
plt.xlabel('Time (ps)')
plt.ylabel('Pressure [bar]')
plt.ylim(-500,500)

# Smoothing using Savitzky-Golay
from scipy.signal import savgol_filter
wild_type_smoothed = savgol_filter(wild_type_data[:,1], 21, 5)
mutant_smoothed = savgol_filter(mutant_data[:,1], 21, 5)

# Plotting raw data and smoothed data for wild type
plt.plot(wild_type_data[:,0], wild_type_data[:,1], linestyle='solid', linewidth='2', color='red', label='Raw - Wild Type')
plt.plot(wild_type_data[:,0], wild_type_smoothed, linestyle='solid', linewidth='2', color='blue', label='Smoothed - Wild Type')

# Plotting raw data and smoothed data for mutant
plt.plot(mutant_data[:,0], mutant_data[:,1], linestyle='solid', linewidth='2', color='green', label='Raw - Mutant')
plt.plot(mutant_data[:,0], mutant_smoothed, linestyle='solid', linewidth='2', color='purple', label='Smoothed - Mutant')

plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig('pressure_plot.png')

# Show the plot
plt.show()
