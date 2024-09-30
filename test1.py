import numpy as np
import pyvista as pv

# Pipe parameters
diameter = 2.0  # Diameter of the pipe
length = 10.0   # Length of the pipe
num_theta = 100  # Resolution in theta (number of points around the circumference)
num_z = 100     # Resolution along the z-axis (length of the pipe)
num_r = 50      # Resolution in radius (to fill the cylinder)

# Derived parameters
radius = diameter / 2

# Create arrays for the cylindrical coordinates
theta = np.linspace(0, 2 * np.pi, num_theta)  # Circumferential angles
z = np.linspace(0, length, num_z)             # Points along the length of the pipe
r = np.linspace(0, radius, num_r)             # Points from center to the edge

# Create 3D grids for r, theta, and z
r_grid, theta_grid, z_grid = np.meshgrid(r, theta, z)

# Convert cylindrical coordinates to Cartesian coordinates for 3D plotting
x_grid = r_grid * np.cos(theta_grid)
y_grid = r_grid * np.sin(theta_grid)

# Flatten the grids to create a list of points for the entire volume
x_flat = x_grid.ravel()
y_flat = y_grid.ravel()
z_flat = z_grid.ravel()

# Create a PolyData object to represent the solid volume of the cylinder
points = np.column_stack((x_flat, y_flat, z_flat))
mesh = pv.PolyData(points)

# ===========================
# Initialize Physical Properties
# ===========================

# Number of points in the mesh
num_points = points.shape[0]

# Initialize velocity vectors (3D vector for each point)
velocity = np.zeros((num_points, 3))  # (vx, vy, vz) for each point

# Initialize scalar fields (temperature, pressure, and density) for each point
temperature = np.zeros(num_points)  # Temperature at each point
pressure = np.zeros(num_points)     # Pressure at each point
density = np.zeros(num_points)      # Density at each point

# ===========================
# Visualization using PyVista
# ===========================

# Convert velocity, temperature, pressure, and density into a PyVista data structure
mesh["velocity_x"] = velocity[:, 0]  # Assign the x-component of velocity
mesh["velocity_y"] = velocity[:, 1]  # Assign the y-component of velocity
mesh["velocity_z"] = velocity[:, 2]  # Assign the z-component of velocity
mesh["temperature"] = temperature    # Assign temperature to each point
mesh["pressure"] = pressure          # Assign pressure to each point
mesh["density"] = density            # Assign density to each point

# Create a PyVista plotter to visualize the mesh and properties
plotter = pv.Plotter()
plotter.add_mesh(mesh, color='lightblue', point_size=2, render_points_as_spheres=True)
plotter.show()


