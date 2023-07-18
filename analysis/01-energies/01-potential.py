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
