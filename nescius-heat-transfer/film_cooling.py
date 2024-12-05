import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
x = np.linspace(0, 2 * np.pi, 500)  # Generate 500 points from 0 to 2π
y = np.sin(x)  # Compute sine of each point

# Create the plot
plt.plot(x, y, label="Sine Wave")  # Plot x vs. y with a label

# Add labels and title
plt.xlabel("x (radians)")
plt.ylabel("sin(x)")
plt.title("Sine Wave Example")

# Add a legend
plt.legend()

# Display the plot
plt.show()
