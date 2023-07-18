# Temperature
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

# Wild type temperature file
wild_type_data = np.loadtxt('/content/gdrive/MyDrive/1aki/nvt_temp.xvg')

# Mutant temperature file
mutant_data = np.loadtxt('/content/gdrive/MyDrive/1aki/nvt_temp.xvg')

plt.title('Temperature during 1000 ps Equilibration (NVT)')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature [K]')

# Plotting wild type temperature
plt.plot(wild_type_data[:,0], wild_type_data[:,1], linestyle='solid', linewidth='2', color='red', label='Wild Type')

# Plotting mutant temperature
plt.plot(mutant_data[:,0], mutant_data[:,1], linestyle='solid', linewidth='2', color='blue', label='Mutant')

plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig('temperature_plot.png')

# Show the plot
plt.show()
