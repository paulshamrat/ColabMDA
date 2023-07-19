#==================================
# 01 Potential Energy
#=================================

# Potential Energy
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

# Wild type potential energy file
wild_type_data = np.loadtxt('/content/gdrive/MyDrive/1aki/em_potential.xvg')

# Mutant potential energy file
mutant_data = np.loadtxt('/content/gdrive/MyDrive/1aki/em_potential.xvg')

plt.title('Potential Energy during Minimization')
plt.xlabel('Energy Minimization Step')
plt.ylabel(r'E$_P$ [kJ•mol$^{-1}]$')

# Plotting wild type potential energy
plt.plot(wild_type_data[:,0], wild_type_data[:,1], linestyle='solid', linewidth='2', color='red', label='Wild Type')

# Plotting mutant potential energy
plt.plot(mutant_data[:,0], mutant_data[:,1], linestyle='solid', linewidth='2', color='blue', label='Mutant')

plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig('potential_plot.png')

# Show the plot
plt.show()


#=================================
# 02 Temperature
#=================================

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


#=================================
# 03 Pressure
#=================================

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



#=================================
# 02 Density
#=================================

# Density
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

# Wild type density file
wild_type_data = np.loadtxt('/content/gdrive/MyDrive/1aki/npt_press_dens.xvg')

# Mutant density file
mutant_data = np.loadtxt('/content/gdrive/MyDrive/1aki/npt_press_dens.xvg')

plt.title('Density during 1000 ps Equilibration (NPT)')
plt.xlabel('Time (ps)')
plt.ylabel('Density [kg•m$^{-3}$]')
plt.ylim(1000,1025)

# Plotting wild type density
plt.plot(wild_type_data[:,0], wild_type_data[:,2], linestyle='solid', linewidth='2', color='red', label='Wild Type')

# Plotting mutant density
plt.plot(mutant_data[:,0], mutant_data[:,2], linestyle='solid', linewidth='2', color='blue', label='Mutant')

plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig('density_plot.png')

# Show the plot
plt.show()
