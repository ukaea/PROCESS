import matplotlib.pyplot as plt
import numpy as np

plt.style.use("ggplot")

# Define the ion_temperature range
ion_temperature = np.linspace(0, 500, 5000)

# Calculate the values
values = (
    1.02934
    - 8.3264e-3 * ion_temperature
    + 1.7631e-4 * ion_temperature**2
    - 1.8201e-6 * ion_temperature**3
    + 6.9855e-9 * ion_temperature**4
)

# Plot the values
plt.plot(ion_temperature, values)
plt.xlabel("Ion Temperature [keV]")
plt.ylabel(r"$\langle \sigma v \rangle_{DDP} \ / \ \langle \sigma v \rangle_{DDN}$")
plt.title("Deuterium - Deuterium Branching Ratio")
plt.ylim(0.8, 1.1)
plt.xscale("log")
plt.grid(True)
plt.minorticks_on()
plt.grid(which="minor", axis="x", linestyle=":", linewidth="0.5")
plt.savefig("deuterium_branching_plot.png")
