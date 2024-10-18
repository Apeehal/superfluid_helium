import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pyvista as pv


# Parameters
initial_radius = 1.0  # Start with this radius
height = 2.0          # Height of all cylinders
num_cylinders = 5     # Number of concentric cylinders
radius_step = 0.15    # Step to decrease the radius for each inner cylinder
theta_resolution = 8  # Divisions around circumference (theta)
z_resolution = 10     # Divisions along the height (z direction)
radial_resolution = 1 # Radial resolution stays 1 because we are manually layering

# Create a plotter
plotter = pv.Plotter()

# Create and plot points from concentric cylinders
for i in range(num_cylinders):
    radius = initial_radius - i * radius_step
    if radius <= 0:
        break  # Stop if radius becomes zero or negative

    # Create the cylinder mesh with the current radius
    cylinder = pv.CylinderStructured(
        radius=radius,
        height=height,
        theta_resolution=theta_resolution,
        z_resolution=z_resolution,
    )
    
    # Extract the points from the cylinder and add them as points to the plotter
    plotter.add_points(cylinder.points, color="red", point_size=10)


new_row = np.array([
    [-1.00000000e+00, 0, 0],
    [-7.77777778e-01, 0, 0],
    [-5.55555556e-01, 0, 0],
    [-3.33333333e-01, 0, 0],
    [-1.11111111e-01, 0, 0],
    [1.11111111e-01, 0, 0],
    [3.33333333e-01, 0, 0],
    [5.55555556e-01, 0, 0],
    [7.77777778e-01, 0, 0],
    [1.00000000e+00, 0, 0]
])




# If the number of columns in large_array is more than in new_row, extend new_row
if cylinder.points.shape[1] > new_row.shape[1]:
    # Create a new array with the same number of columns as large_array, filled with zeros in the last columns
    new_row_extended = np.zeros((new_row.shape[0], cylinder.points.shape[1]))
    new_row_extended[:, :new_row.shape[1]] = new_row  # Copy original data into the extended array
else:
    new_row_extended = new_row

# Concatenate the new_row with the large_array
combined_array = np.concatenate((cylinder.points, new_row_extended), axis=0)

# Print the result
plotter = pv.Plotter()



# Show the plot with concentric cylinder points and central axis point


# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(combined_array[:, 0], combined_array[:, 1], combined_array[:, 2], c='b', marker='o')  # c is the color, marker is the shape
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()