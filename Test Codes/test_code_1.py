import numpy as np
import matplotlib.pyplot as plt

# Define the Node class with varying k and alpha based on superfluid density, entropy, etc.
class Node:
    def __init__(self, x, y, temperature=0.0, pressure=0.0, velocity=(0.0, 0.0),
                 density_superfluid=1.0, density_normal=1.0, entropy=1.0, viscosity=0.0,alpha = 0.0):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.pressure = pressure
        self.velocity = velocity
        self.density_superfluid = density_superfluid
        self.density_normal = density_normal
        self.entropy = entropy
        self.viscosity = viscosity
        self.k = self.calculate_k()  # Thermal conductivity (initial calculation)
        self.alpha = alpha  # Initialize alpha; could be dynamically computed if needed

    # Function to calculate k based on density, entropy, and viscosity
    def calculate_k(self):
        term = (num_nodes_y * 10**-3) * (self.density_superfluid + self.density_normal) * self.entropy
        k = (term ** 2) / (8 * self.viscosity)
        return k

# Define grid dimensions
num_nodes_x = 60
num_nodes_y = 40
dx = dy = 0.01  # Spatial step (assuming uniform grid)
dt = 0.01    # Time step size
time_steps = 1  # Number of time steps
q_top = 0.2   # Heat flux at the top boundary
q_bottom = 0.2 # Heat flux at the bottom boundary
lambda_transition = 2.17


# Create a grid of nodes with varying density_superfluid, density_normal, entropy, and viscosity
def initialize_node(x, y):
    # Example initialization: You can modify the logic to reflect realistic conditions
    temperature = 1.9
    density_normal = 147.5
    density_superfluid = (density_normal/((temperature/lambda_transition)**5.6)) - density_normal
    entropy = 756.95  # Entropy can be constant or vary
    viscosity = 2.5*10**-6
    return Node(x, y, temperature=temperature, density_normal=density_normal, density_superfluid=density_superfluid, entropy=entropy, viscosity=viscosity)

nodes = [[initialize_node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]



c_p = 5000
# Function to apply FDM and calculate next temperature values with varying k and Neumann BCs
def fdm_step(nodes, dt, dx, dy, q_top, q_bottom):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature, node.pressure, node.velocity, node.density_superfluid,
                       node.density_normal, node.entropy, node.viscosity) for x, node in enumerate(row)] for y, row in enumerate(nodes)]
    
    # Update interior nodes using FDM
    for y in range(num_nodes_y):
        for x in range(num_nodes_x):

            #central nodes
            if x != 0 and x != num_nodes_x-1 and y != 0 and y != num_nodes_y-1:
                T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
                T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
                k_local = nodes[y][x].calculate_k()   # Calculate k dynamically based on node properties
                print(nodes[y][x].temperature)
                alpha_local = k_local/((nodes[y][x].density_superfluid+nodes[y][x].density_normal)*c_p) # Use local alpha value (could also be dynamic if needed)
                new_nodes[y][x].temperature = nodes[y][x].temperature + alpha_local * dt * (T_xx + T_yy)

            #right central nodes
            if x == num_nodes_x - 1 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x-1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                k_local = nodes[y][x].calculate_k()   # Calculate k dynamically based on node properties
                alpha_local = k_local/((nodes[y][x].density_superfluid+nodes[y][x].density_normal)*c_p) # Use local alpha value (could also be dynamic if needed)
                new_nodes[y][x].temperature = nodes[y][x].temperature + alpha_local * dt * (T_avg)
            
            #left central nodes
            if x == 0 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x+1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                k_local = nodes[y][x].calculate_k()   # Calculate k dynamically based on node properties
                alpha_local = k_local/((nodes[y][x].density_superfluid+nodes[y][x].density_normal)*c_p) # Use local alpha value (could also be dynamic if needed)
                new_nodes[y][x].temperature = nodes[y][x].temperature + alpha_local * dt * (T_avg)



    ################################################################
    # Neumann Boundary Conditions with heat flux at top and bottom##
    ################################################################

    # Bottom boundary (specified heat flux q_bottom)
            if y == 0 and (x != 0 and x != num_nodes_x - 1):
                #T_interior = nodes[1][x].temperature  # First interior node above the boundary
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y+1][x].temperature) / 3
                k_local = nodes[0][x].calculate_k()   # Calculate k dynamically based on node properties
                new_nodes[0][x].temperature = T_avg + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

    # Top boundary (specified heat flux q_top)
            if y == num_nodes_y-1  and (x != 0 and x != num_nodes_x - 1):
                #T_interior = nodes[num_nodes_y - 2][x].temperature  # Last interior node below the boundary
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y-1][x].temperature) / 3
                k_local = nodes[num_nodes_y - 1][x].calculate_k()   # Calculate k dynamically based on node properties
                new_nodes[num_nodes_y - 1][x].temperature = T_avg + (q_top * dy) / k_local  # Apply Neumann BC for top boundary
    
    # Left and right boundaries (assuming no flux for simplicity, can be adjusted if needed)
    #for y in range(1, num_nodes_y-1):
    #    new_nodes[y][0].temperature = new_nodes[y][1].temperature  # Left boundary (no flux)
    #   new_nodes[y][num_nodes_x - 1].temperature = new_nodes[y][num_nodes_x - 2].temperature  # Right boundary (no flux)

            if x==0 and y==0:
                T_avg = (nodes[y][x+1].temperature + nodes[y+1][x].temperature)/2
                k_local = nodes[0][0].calculate_k() 
                new_nodes[0][x].temperature = T_avg + (q_bottom * dy) / k_local


            if x==0 and y==num_nodes_y-1:
                T_avg = (nodes[y-1][x].temperature + nodes[y][x+1].temperature)/2
                k_local = nodes[y][0].calculate_k() 
                new_nodes[y][0].temperature = T_avg + (q_bottom * dy) / k_local
     
            if x==num_nodes_x-1 and y==0:
                T_avg = (nodes[y][x-1].temperature + nodes[y+1][x].temperature)/2
                k_local = nodes[0][num_nodes_x-1].calculate_k() 
                new_nodes[0][num_nodes_x-1].temperature = T_avg + (q_top * dy) / k_local
 
            if x==num_nodes_x-1 and y==num_nodes_y-1:
                T_avg = (nodes[y-1][x].temperature + nodes[y][x-1].temperature)/2
                k_local = nodes[num_nodes_y-1][num_nodes_x-1].calculate_k() 
                new_nodes[num_nodes_y-1][num_nodes_x-1].temperature = T_avg + (q_top * dy) / k_local
 

    return new_nodes



    # Return the updated grid
    

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, dt, dx, dy, q_top, q_bottom)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution with Varying k and Neumann BC after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()
