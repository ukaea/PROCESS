import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

import process.superconductors as superconductors

temp_c0max = 16.34  # Critical temperature (K) at zero field and strain
b_c20max = (
    33.24  # Upper critical field (T) for superconductor at zero temperature and strain
)
epsilon = 0.00  # Strain in superconductor

# Create a grid of temperature and field values
temp_conductor = np.linspace(1, 10, 50)  # Temperature range (K)
b_conductor = np.linspace(1, 14, 50)  # Magnetic field range (T)
temp_grid, b_grid = np.meshgrid(temp_conductor, b_conductor)


j_scaling = np.zeros_like(temp_grid)
for i in range(temp_grid.shape[0]):
    for j in range(temp_grid.shape[1]):
        j_scaling[i, j], _, _ = superconductors.western_superconducting_nb3sn(
            temp_grid[i, j], b_grid[i, j], epsilon, b_c20max, temp_c0max
        )
        # Convert from A/m² to kA/mm² (1 A/m² = 1e-6 A/mm²)
        j_scaling[i, j] *= 1e-9
        print(f"j_scaling[{i}, {j}] = {j_scaling[i, j]} A/mm²")

# Plot the critical current density as a function of field and temperature
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(temp_grid, b_grid, j_scaling, cmap="viridis", edgecolor="none")


ax.set_xlabel("Temperature of superconductor (K)")
ax.set_ylabel("Magnetic field at superconductor (T)")
ax.set_zlabel("$J_{\\text{crit}}$ (kA/mm²)")
ax.set_title(
    "Western Superconducting Ltd. Nb\u2083Sn Critical Current Density Surface at Zero Strain"
)
ax.view_init(elev=30, azim=45)  # Rotate the plot around the z-axis by 45 degrees
fig.colorbar(surf, shrink=0.5, aspect=10)

# Save the plot as a PNG file
plt.savefig("western_superconducting_Nb3Sn_zero_strain.png", dpi=1000)


# Convert the data to a format suitable for Plotly
fig_plotly = go.Figure(
    data=[
        go.Surface(z=j_scaling, x=temp_conductor, y=b_conductor, colorscale="Viridis")
    ]
)

# Update layout for better visualization
fig_plotly.update_layout(
    title="Western Superconducting Ltd. Nb\u2083Sn Critical Current Density Surface at Zero Strain",
    scene={
        "xaxis_title": "Temperature of superconductor (K)",
        "yaxis_title": "Magnetic field strength at superconductor (T)",
        "zaxis_title": "Critical current density (kA/mm²)",
    },
)
# Export the interactive plot to an HTML file
fig_plotly.write_html("western_superconducting_Nb3Sn_zero_strain.html")
