import numpy as np
import pyvista as pv
from pyvista import CellType


lambda_transition = 2.17


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

# Initialize a list to store the hexahedral cells
hexa_cells = []

# Define a helper function to convert 2D indices to 1D index in the points array
def idx(r_idx, theta_idx, z_idx):
    return z_idx * (num_r_liquid + 1) * num_theta + r_idx * num_theta + theta_idx

# Initialize an empty list to store cell labels
cell_labels = []

# Loop over the grid to define hexahedral cells
for z in range(num_z):
    for r in range(num_r_liquid):
        for theta in range(num_theta):
            # Get the indices of the eight points forming a hexahedral cell
            p1 = idx(r, theta, z)
            p2 = idx(r + 1, theta, z)
            p3 = idx(r + 1, (theta + 1) % num_theta, z)
            p4 = idx(r, (theta + 1) % num_theta, z)
            p5 = idx(r, theta, z + 1)
            p6 = idx(r + 1, theta, z + 1)
            p7 = idx(r + 1, (theta + 1) % num_theta, z + 1)
            p8 = idx(r, (theta + 1) % num_theta, z + 1)

            # Add the hexahedron (8 points) to the list
            hexa_cells.append([8, p1, p2, p3, p4, p5, p6, p7, p8])

            # Determine if the cell is on the outer surface or inner part
            if r == num_r_liquid - 1:  # Outer surface cells at max radius
                cell_labels.append(2)
            else:  # Inner cells
                cell_labels.append(1)

# Convert cells to a format suitable for PyVista UnstructuredGrid
hexa_cells = np.hstack(hexa_cells)  # Flatten the cell array

# Create the cell types array (all hexahedral cells in this case)
cell_types = np.full(len(hexa_cells) // 9, CellType.HEXAHEDRON, dtype=np.uint8)

# Create the UnstructuredGrid from points, cells, and cell types
mesh = pv.UnstructuredGrid(hexa_cells, cell_types, points)

# Add the cell labels as a scalar field to the grid
mesh.cell_data['Labels'] = np.array(cell_labels)



# ===========================
# Gravity Vector
# ===========================
#g = [9.8, np.pi, z]



# ===========================
# Compute Cell Volumes
# ===========================
volumes = mesh.compute_cell_sizes(volume=True)["Volume"]  # Compute the volume of each cell
mesh["cell_volume"] = volumes  # Attach the volumes to the mesh
mesh["body_label"] = np.array(cell_labels)  # Label cells with body identifiers (1 for liquid, 2 for solid wall)



# ===========================
# Assign Bodies
# ===========================

body_labels = mesh.cell_data["body_label"]
liquid_body = np.where(body_labels == 1)
solid_body = np.where(body_labels == 2)


# ===========================
# Initialize Temperature Field
# ===========================

# Initialize temperature array for the cells
temperature = np.zeros(mesh.n_cells)

# Assign 1.9K to the liquid region and 300K to the solid outer wall based on cell labels
temperature[mesh["body_label"] == 1] = 1.9    # Liquid region at 1.9K
temperature[mesh["body_label"] == 2] = 300.0  # Solid outer wall at 300K

# Assign temperature field to the mesh
mesh["temperature"] = temperature


# ===========================
# Initialize Presssure Field
# ===========================

pressure = np.zeros(mesh.n_cells)
pressure[mesh["body_label"] == 1] = 101325    # Liquid region at 1.9K
pressure[mesh["body_label"] == 2] = 101325    # Liquid region at 1.9K
mesh["pressure"] = pressure



# ===========================
# Initialize Density
# ===========================

normal_to_total_density = (temperature[liquid_body] / lambda_transition)**5.6

superfluid_to_total_density = (1-((temperature[liquid_body] / lambda_transition)**5.6))

normal_to_superfluid_density = normal_to_total_density/superfluid_to_total_density

#superfluid_density[mesh["body_label" == 1]] = normal_density[liquid_body][0]/(normal_to_total_density*(1/superfluid_to_total_density))


normal_density = np.zeros(mesh.n_cells)   # Density at each point

normal_density[mesh["body_label"] == 1] = 147.5 
mesh["normal_density"] = normal_density




superfluid_density = np.zeros(mesh.n_cells)   # Density at each point

superfluid_density_value = normal_density[liquid_body][0]/(normal_to_total_density*(1/superfluid_to_total_density))


superfluid_density[mesh["body_label"] == 1] = superfluid_density_value[0]
mesh["superfluid_density"] = superfluid_density



# ===========================
# Initialize Entropy
# ===========================

entropy = np.zeros(mesh.n_cells)
entropy[mesh["body_label"] == 1] = 1500    # Liquid region at 1.9K
mesh["pressure"] = pressure



# ===========================
# Initialize Velocity
# ===========================

# Initialize velocity vectors (3D vector for each point)
velocity = np.zeros((mesh.n_points, 3))  # (vx, vy, vz) for each point

# Initialize scalar fields (pressure and density) for each point

# Assign the physical fields to the mesh
mesh["velocity_x"] = velocity[:, 0]
mesh["velocity_y"] = velocity[:, 1]
mesh["velocity_z"] = velocity[:, 2]




"""
# ===========================
# COMPUTATION
# ===========================
"""
timestep = 0.1


# ===========================
# COMPUTE grad T
# ===========================


mesh_with_gradient_T = mesh.compute_derivative(scalars="temperature", gradient=True)
temperature_gradient = mesh_with_gradient_T.cell_data["gradient"]

t_gradient_r = temperature_gradient[:, 0]
t_gradient_theta = temperature_gradient[:, 1]
t_gradient_z = temperature_gradient[:, 2]

# Now you have the temperature gradient for each cell
#for i in range(mesh_with_gradient_T.n_cells):
#    print(f"Cell {i} temperature gradient: ({t_gradient_r[i]}, {t_gradient_theta[i]}, {t_gradient_z[i]})")

# ===========================
# COMPUTE grad P
# ===========================

mesh_with_gradient = mesh.compute_derivative(scalars="pressure", gradient=True)

pressure_gradient = mesh_with_gradient.cell_data["gradient"]

p_gradient_r = pressure_gradient[:, 0]
p_gradient_theta = pressure_gradient[:, 1]
p_gradient_z = pressure_gradient[:, 2]

# Now you have the temperature gradient for each cell

#for i in range(mesh_with_gradient.n_cells):
#    print(f"Cell {i} temperature gradient: ({p_gradient_r[i]}, {p_gradient_theta[i]}, {p_gradient_z[i]})")

# ===========================
# Superfluid Velocity
# ===========================



rho_vs_r = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vs_r[liquid_body] = (superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_r[liquid_body] - (superfluid_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))*p_gradient_r[liquid_body]) * timestep
mesh.cell_data['rho_vs_r'] = rho_vs_r

rho_vs_theta = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vs_theta[liquid_body] = (superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_theta[liquid_body] - (superfluid_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))*p_gradient_theta[liquid_body]) * timestep
mesh.cell_data['rho_vs_theta'] = rho_vs_theta

rho_vs_z = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vs_z[liquid_body] = (superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_z[liquid_body] - (superfluid_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))*p_gradient_z[liquid_body]) * timestep
mesh.cell_data['rho_vs_z'] = rho_vs_z


vs_r = np.zeros((mesh.n_cells)) 
vs_r[liquid_body] = rho_vs_r[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vs_r'] = vs_r

vs_theta = np.zeros((mesh.n_cells)) 
vs_theta[liquid_body] = rho_vs_theta[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vs_theta'] = vs_theta

vs_z = np.zeros((mesh.n_cells)) 
vs_z[liquid_body] = rho_vs_z[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vs_z'] = vs_z


vs_abs = np.zeros((mesh.n_cells))
vs_abs[liquid_body] = np.sqrt(vs_z[liquid_body]**2+vs_theta[liquid_body]**2+vs_r[liquid_body]**2)
mesh.cell_data['rho_vs_abs'] = vs_abs



# ===========================
# Normal Velocity
# ===========================


rho_vn_r = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vn_r[liquid_body] = (-superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_r[liquid_body] - (normal_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))    *p_gradient_r[liquid_body]) * timestep
mesh.cell_data['rho_vn_r'] = rho_vn_r

rho_vn_theta = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vn_theta[liquid_body] = (-superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_theta[liquid_body] - (normal_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))*p_gradient_theta[liquid_body]) * timestep
mesh.cell_data['rho_vn_theta'] = rho_vn_theta

rho_vn_z = np.zeros((mesh.n_cells))  # (vx, vy, vz) for each point
rho_vn_z[liquid_body] = (-superfluid_density[liquid_body] * entropy[liquid_body] * t_gradient_z[liquid_body] - (normal_density[liquid_body]/(superfluid_density[liquid_body]+normal_density[liquid_body]))*p_gradient_z[liquid_body]) * timestep
mesh.cell_data['rho_vn_z'] = rho_vn_z


vn_r = np.zeros((mesh.n_cells)) 
vn_r[liquid_body] = rho_vn_r[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vn_r'] = vn_r

vn_theta = np.zeros((mesh.n_cells)) 
vn_theta[liquid_body] = rho_vn_theta[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vn_theta'] = vn_theta

vn_z = np.zeros((mesh.n_cells)) 
vn_z[liquid_body] = rho_vn_z[liquid_body]/superfluid_density[liquid_body]
mesh.cell_data['vn_z'] = vn_z


vn_abs = np.zeros((mesh.n_cells))
vn_abs[liquid_body] = np.sqrt(vn_z[liquid_body]**2+vn_theta[liquid_body]**2+vn_r[liquid_body]**2)
mesh.cell_data['rho_vn_abs'] = vn_abs






# ===========================
# Visualization using PyVista
# ===========================

# Create a PyVista plotter to visualize the liquid and outer wall, colored by temperature
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars=vn_abs, show_edges=True, opacity=0.7)
plotter.show()

