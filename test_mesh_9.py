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

# Create an empty list to store the cells (quadrilateral faces)
faces = []

# Define a helper function to convert 2D indices to 1D index in the points array
def idx(r_idx, theta_idx, z_idx):
    return z_idx * (num_r_liquid + 1) * num_theta + r_idx * num_theta + theta_idx

# Loop over the grid to define faces (connect points)
for z in range(num_z):
    for r in range(num_r_liquid):
        for theta in range(num_theta):
            # Get the indices of the four points forming a quadrilateral face
            p1 = idx(r, theta, z)
            p2 = idx(r + 1, theta, z)
            p3 = idx(r + 1, (theta + 1) % num_theta, z)
            p4 = idx(r, (theta + 1) % num_theta, z)
            p5 = idx(r, theta, z + 1)
            p6 = idx(r + 1, theta, z + 1)
            p7 = idx(r + 1, (theta + 1) % num_theta, z + 1)
            p8 = idx(r, (theta + 1) % num_theta, z + 1)

            # Add the quadrilateral for this section of the mesh
            # Create faces using points in the z-direction, radial, and theta directions
            faces.append([4, p1, p2, p3, p4])  # Lower quad (in z-plane)
            faces.append([4, p5, p6, p7, p8])  # Upper quad (in z-plane)
            faces.append([4, p1, p2, p6, p5])  # Radial quad
            faces.append([4, p4, p3, p7, p8])  # Next radial quad

for i in range(num_r_liquid+1):
    for z in range(num_z):
        for theta in range(num_theta):
            # Close the inner radial faces (at r = 0)
            p1 = idx(0, theta, z)
            p2 = idx(0, (theta + 1) % num_theta, z)
            p3 = idx(0, theta, z + 1)
            p4 = idx(0, (theta + 1) % num_theta, z + 1)
            faces.append([3, p1, p3, p2])  # Triangle for the lower section
            faces.append([3, p2, p3, p4])  # Triangle for the upper section

            # Close the outer radial faces (at the last radial position)
            p1 = idx(i, theta, z)
            p2 = idx(i, (theta + 1) % num_theta, z)
            p3 = idx(i, theta, z + 1)
            p4 = idx(i, (theta + 1) % num_theta, z + 1)
            faces.append([4, p1, p2, p4, p3])  # Quadrilateral for the outer section#

# Create a PyVista mesh from the points and faces
faces = np.hstack(faces)  # Flatten the face array
mesh = pv.PolyData(points, faces)






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












# Plot the mesh with the radial sectioning
plotter = pv.Plotter()
plotter.add_mesh(mesh, color="lightblue", show_edges=True, opacity=0.7)
plotter.show()
