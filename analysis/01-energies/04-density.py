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
