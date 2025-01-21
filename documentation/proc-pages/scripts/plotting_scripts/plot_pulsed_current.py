import matplotlib.pyplot as plt

# Define data points for each segment
x_red_blue = [0, 1, 2, 2.5, 5, 6]  # x-coordinates for the red and blue dashed lines

y_red = [0, 3, 2.5, 2.5, 1.5, 0]  # y-coordinates for the red dashed line

y_blue = [0, 4, -2, -2, -3, 0]  # y-coordinates for the blue dashed line

x_black = x_red_blue[1:]  # x-coordinates for the black solid line

y_black = [0, 3.5, 3.5, 3.5, 0]  # y-coordinates for the black solid line

# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Plot the PF coil current
ax.plot(x_red_blue, y_red, "r--", linewidth=2, label="PF coil")

# Plot the central solenoid current
ax.plot(x_red_blue, y_blue, "b--", linewidth=2, label="Central solenoid")

# Plot the plasma current
ax.plot(x_black, y_black, "black", linewidth=2, label="Plasma")

# Move the x-axis to 0 on the y-axis
ax.spines["bottom"].set_position("zero")

# Adjust the axes limits to allow space for labels
ax.set_xlim(0, 7.2)
ax.set_ylim(-3.5, 4.2)

# Annotate key points
ax.annotate(
    "Pulse Start",
    xy=(1, 0),
    xytext=(1, -0.75),
    arrowprops={"facecolor": "black", "arrowstyle": "->"},
    fontsize=10,
    ha="center",
)
ax.annotate(
    "Start of Flat-Top",
    xy=(2, 0),
    xytext=(2, 0.75),
    arrowprops={"facecolor": "black", "arrowstyle": "->"},
    fontsize=10,
    ha="center",
)
ax.annotate(
    "End of Flat-Top",
    xy=(5, 0),
    xytext=(5, -0.75),
    arrowprops={"facecolor": "black", "arrowstyle": "->"},
    fontsize=10,
    ha="center",
)

# Annotate intervals
ax.annotate(
    "",
    xy=(0, -3.2),
    xytext=(1, -3.2),
    xycoords="data",
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate("Precharge", xy=(0.5, -3.4), ha="center", va="center")
ax.annotate(
    "",
    xy=(1, -3.2),
    xytext=(2, -3.2),
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate(r"$I_{\text{P}}$ Ramp-Up", xy=(1.5, -3.4), ha="center", va="center")
ax.annotate(
    "",
    xy=(2, -3.2),
    xytext=(2.5, -3.2),
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate("Fusion\nRamp", xy=(2.25, -3.5), ha="center", va="center")
ax.annotate(
    "",
    xy=(2.5, -3.2),
    xytext=(5, -3.2),
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate("Burn", xy=(3.75, -3.4), ha="center", va="center")
ax.annotate(
    "",
    xy=(5, -3.2),
    xytext=(6, -3.2),
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate("Ramp Down", xy=(5.5, -3.4), ha="center", va="center")
ax.annotate(
    "",
    xy=(6, -3.2),
    xytext=(7, -3.2),
    arrowprops={"arrowstyle": "|-|", "shrinkA": 0, "shrinkB": 0},
)
ax.annotate("Between Pulse", xy=(6.5, -3.4), ha="center", va="center")

# Add axis labels
ax.set_xlabel("Time", fontsize=12)
ax.xaxis.set_label_coords(0.97, 0.5)
ax.set_ylabel("Current", fontsize=12)

# Remove axes ticks
ax.set_xticklabels([])
ax.set_yticklabels([])

# Add a title
ax.set_title("Current Profiles Over Time", fontsize=14)

# Add a legend
ax.legend()

# Add a grid for better readability
ax.grid(True, linestyle="--", alpha=0.6)

# Adjust layout and display the plot
plt.tight_layout()
plt.savefig("current_vs_time_plot.png")
