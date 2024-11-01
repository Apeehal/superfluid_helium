import numpy as np
import pyvista as pv

lambda_transition = 2.17

# Pipe parameters for the liquid region (inner cylinder)
liquid_diameter = 2.0  # Inner liquid diameter (same as before)
length = 10.0          # Length of the pipe
num_theta = 12         # Number of radial sections (angular resolution)
num_z = 10             # Resolution along the z-axis (length of the pipe)
num_r_liquid = 5       # Resolution in radius for the liquid region

# Outer solid wall parameters
wall_thickness = 0.5   # Thickness of the solid outer wall
num_r_wall = 2         # Resolution in radius for the solid outer wall

# Derived parameters
liquid_radius = liquid_diameter / 2
outer_wall_radius = liquid_radius + wall_thickness  # Outer wall extends beyond liquid

# Create arrays for the cylindrical coordinates
theta = np.linspace(0, 2 * np.pi, num_theta + 1)  # Circumferential angles (add 1 to close loop)
z = np.linspace(0, length, num_z + 1)             # Points along the length of the pipe (add 1 for edges)
r_liquid = np.linspace(0, liquid_radius, num_r_liquid + 1)  # Radial points for the liquid region
r_wall = np.linspace(liquid_radius, outer_wall_radius, num_r_wall + 1)  # Radial points for the solid wall

# Storage for the points and the cells (volumes)
points = []
cells = []
cell_body_labels = []

# ===========================
# Liquid Region (Body 1)
# ===========================

for i in range(num_r_liquid):
    for j in range(num_theta):
        for k in range(num_z):
            # Define the corner points of each wedge cell in the liquid region
            r_inner = r_liquid[i]
            r_outer = r_liquid[i+1]
            theta_start = theta[j]
            theta_end = theta[j+1]
            z_bottom = z[k]
            z_top = z[k+1]

            # Create the eight corner points of the wedge cell in Cartesian coordinates
            p0 = [r_inner * np.cos(theta_start), r_inner * np.sin(theta_start), z_bottom]  # Bottom-inner-start
            p1 = [r_outer * np.cos(theta_start), r_outer * np.sin(theta_start), z_bottom]  # Bottom-outer-start
            p2 = [r_outer * np.cos(theta_end), r_outer * np.sin(theta_end), z_bottom]      # Bottom-outer-end
            p3 = [r_inner * np.cos(theta_end), r_inner * np.sin(theta_end), z_bottom]      # Bottom-inner-end
            p4 = [r_inner * np.cos(theta_start), r_inner * np.sin(theta_start), z_top]     # Top-inner-start
            p5 = [r_outer * np.cos(theta_start), r_outer * np.sin(theta_start), z_top]     # Top-outer-start
            p6 = [r_outer * np.cos(theta_end), r_outer * np.sin(theta_end), z_top]         # Top-outer-end
            p7 = [r_inner * np.cos(theta_end), r_inner * np.sin(theta_end), z_top]         # Top-inner-end

            # Add the points to the points list
            point_indices = []
            for pt in [p0, p1, p2, p3, p4, p5, p6, p7]:
                if pt not in points:  # Avoid duplicating points
                    points.append(pt)
                point_indices.append(points.index(pt))

            # Define the cell (hexahedral wedge with 8 corner points)
            cells.append([8, point_indices[0], point_indices[1], point_indices[2], point_indices[3],
                              point_indices[4], point_indices[5], point_indices[6], point_indices[7]])
            cell_body_labels.append(1)  # Label this cell as part of the liquid region

# ===========================
# Solid Outer Wall (Body 2)
# ===========================

for i in range(num_r_wall):
    for j in range(num_theta):
        for k in range(num_z):
            # Define the corner points of each wedge cell in the solid outer wall
            r_inner = r_wall[i]
            r_outer = r_wall[i+1]
            theta_start = theta[j]
            theta_end = theta[j+1]
            z_bottom = z[k]
            z_top = z[k+1]

            # Create the eight corner points of the wedge cell in Cartesian coordinates
            p0 = [r_inner * np.cos(theta_start), r_inner * np.sin(theta_start), z_bottom]  # Bottom-inner-start
            p1 = [r_outer * np.cos(theta_start), r_outer * np.sin(theta_start), z_bottom]  # Bottom-outer-start
            p2 = [r_outer * np.cos(theta_end), r_outer * np.sin(theta_end), z_bottom]      # Bottom-outer-end
            p3 = [r_inner * np.cos(theta_end), r_inner * np.sin(theta_end), z_bottom]      # Bottom-inner-end
            p4 = [r_inner * np.cos(theta_start), r_inner * np.sin(theta_start), z_top]     # Top-inner-start
            p5 = [r_outer * np.cos(theta_start), r_outer * np.sin(theta_start), z_top]     # Top-outer-start
            p6 = [r_outer * np.cos(theta_end), r_outer * np.sin(theta_end), z_top]         # Top-outer-end
            p7 = [r_inner * np.cos(theta_end), r_inner * np.sin(theta_end), z_top]         # Top-inner-end

            # Add the points to the points list
            point_indices = []
            for pt in [p0, p1, p2, p3, p4, p5, p6, p7]:
                if pt not in points:  # Avoid duplicating points
                    points.append(pt)
                point_indices.append(points.index(pt))

            # Define the cell (hexahedral wedge with 8 corner points)
            cells.append([8, point_indices[0], point_indices[1], point_indices[2], point_indices[3],
                              point_indices[4], point_indices[5], point_indices[6], point_indices[7]])
            cell_body_labels.append(2)  # Label this cell as part of the solid outer wall

# Convert points and cells to NumPy arrays for PyVista
points = np.array(points)
cells = np.array(cells)

# Create the PyVista Unstructured Grid
cell_types = np.full(len(cells), pv.CellType.HEXAHEDRON, dtype=np.uint8)  # Define cell type as HEXAHEDRON
mesh = pv.UnstructuredGrid(cells, cell_types, points)


plotter = pv.Plotter()
plotter.add_mesh(mesh, show_edges=True, opacity=0.7)
plotter.show()
