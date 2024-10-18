import numpy as np
import pyvista as pv

# Parameters for the liquid region (inner cylinder)
liquid_diameter = 2.0  # Inner liquid diameter
length = 10.0          # Length of the pipe
num_theta = 12         # Number of radial sections (angular resolution)
num_z = 10             # Resolution along the z-axis (length of the pipe)
num_r_liquid = 5       # Resolution in radius for the liquid region

# Derived parameters
liquid_radius = liquid_diameter / 2

# Create a grid of z, r, and theta coordinates
z_co_ords = np.linspace(0, length, num_z + 1)  # Z-coordinates along the length
r_co_ords = np.linspace(0, liquid_radius, num_r_liquid + 1)  # Radial coordinates
theta_co_ords = np.linspace(0, 2 * np.pi, num_theta, endpoint=False)  # Theta coordinates (angles)

# Initialize an empty list to store the 3D points
points = []

# Loop through each combination of r, theta, and z to generate the points
for z in z_co_ords:
    for r in r_co_ords:
        for theta in theta_co_ords:
            # Convert cylindrical (r, theta, z) to Cartesian (x, y, z)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            points.append([x, y, z])

# Convert points to a NumPy array for further processing
points = np.array(points)

# Create a PyVista point cloud
cloud = pv.PolyData(points)

# Plot the points
plotter = pv.Plotter()
plotter.add_points(cloud, color="blue", point_size=5.0)
plotter.show()
