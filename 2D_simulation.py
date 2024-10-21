import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Define the grid dimensions
num_nodes_x = 60  # Number of nodes in the x-direction
num_nodes_y = 40  # Number of nodes in the y-direction

# Generate the grid of nodes
grid = [(x, y) for y in range(num_nodes_y) for x in range(num_nodes_x)]

lambda_transition = 2.17
#viscosity = (3.5*10)**-6 #assume constant viscosity
viscosity = 2.5*10**-6
class Node:
    def __init__(self, x, y, heat_flux=(0.0,0.0), temperature=0.0, temperature_gradient=(0.0,0.0), pressure=0.0, pressure_gradient=(0.0,0.0), velocity_normal=(0.0, 0.0), velocity_superfluid=(0.0,0.0), density_normal=0.0, density_superfluid=0.0, entropy=0.0):
    
        self.x = x
        self.y = y
        self.heat_flux = heat_flux
        self.temperature = temperature
        self.temperature_gradient = temperature_gradient
        self.pressure = pressure
        self.pressure_gradient = pressure_gradient
        self.velocity_normal = velocity_normal
        self.velocity_superfluid = velocity_superfluid
        self.density_normal = density_normal
        self.density_superfluid = density_superfluid
        self.entropy = entropy


    def __repr__(self):
        return f"Node({self.x}, {self.y}, q={self.heat_flux}, T={self.temperature}, grad_T={self.temperature_gradient}, P={self.pressure}, grad_P={self.pressure_gradient}, v_n={self.velocity_normal}, v_s={self.velocity_superfluid}, rho_n={self.density_normal}, rho_s={self.density_superfluid}, s={self.entropy})"



# Create a grid of nodes
nodes = [[Node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]


############################
# Assign Initial Conditions
############################

initial_temperature = 1.9
for row in nodes:
    for node in row:
        node.temperature = initial_temperature


initial_pressure = 0.101325
for row in nodes:
    for node in row:
        node.pressure = initial_pressure


initial_entropy = 756.95
for row in nodes:
    for node in row:
        node.entropy = initial_entropy


initial_normal_density = 147.5
for row in nodes:
    for node in row:
        node.density_normal = initial_normal_density


initial_superfluid_density = (initial_normal_density/((initial_temperature/lambda_transition)**5.6)) - initial_normal_density
for row in nodes:
    for node in row:
        node.density_superfluid = initial_superfluid_density



initial_heat_flux = 0
for row in nodes:
    for node in row:
        node.heat_flux = initial_heat_flux


initial_normal_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_normal = initial_normal_density


initial_superfluid_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_superfluid = initial_superfluid_density


############################
# Assign boundary conditions
############################
for y in range(num_nodes_y):
    for x in range(num_nodes_x):
        if y == 0:  # Bottom boundary
            nodes[y][x].heat_flux = 0.2
            #nodes[y][x].temperature = 300

        elif y == num_nodes_y - 1:  # Top boundary
            nodes[y][x].heat_flux = 0.2
           
        if x == 0:  # Left boundary
            nodes[y][x].pressure = 101325



############################
############################
# Start Computation ########
############################
############################
"""
for y in range(num_nodes_y):
    for x in range(num_nodes_x):
        #Top and Bottom Nodes
        if y == 0 and (x != 0 and x != num_nodes_x - 1):
            print(nodes[y][x].heat_flux)

"""
############################
# Solve for Temperatures ###
############################
timestep = 1000


def temperature_averaging_step(nodes):

    new_nodes = [[Node(x, y, node.temperature) for x, node in enumerate(row)] for y, row in enumerate(nodes)]


    for y in range(num_nodes_y):
        for x in range(num_nodes_x):
            
            #Top and Bottom Nodes
            if y == 0 and (x != 0 and x != num_nodes_x - 1):
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y+1][x].temperature) / 3
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(3*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T


            if y == num_nodes_y-1  and (x != 0 and x != num_nodes_x - 1):
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y-1][x].temperature) / 3
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(3*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T


            #corner nodes
            if x==0 and y==0:
                T_avg = (nodes[y][x+1].temperature + nodes[y+1][x].temperature)/2
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(2*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T



            if x==0 and y==num_nodes_y-1:
                T_avg = (nodes[y-1][x].temperature + nodes[y][x+1].temperature)/2
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(2*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T


            if x==num_nodes_x-1 and y==0:
                T_avg = (nodes[y][x-1].temperature + nodes[y+1][x].temperature)/2
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(2*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T


            if x==num_nodes_x-1 and y==num_nodes_y-1:
                T_avg = (nodes[y-1][x].temperature + nodes[y][x-1].temperature)/2
                k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
                Q_component = nodes[y][x].heat_flux
                T = T_avg + ((Q_component*timestep))*(((1*10**-3)**2)/(2*k))
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T






            #CentralNodes
            if x != 0 and x != num_nodes_x-1 and y != 0 and y != num_nodes_y-1:  
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature)/4
                T = T_avg 
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T

            #right central nodes
            if x == num_nodes_x - 1 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x-1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                T = T_avg
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T

            #left central nodes
            if x == 0 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x+1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                T = T_avg 
                new_nodes[y][x].temperature = new_nodes[y][x].temperature + T

            new_nodes[y][x].density_normal = initial_normal_density
            new_nodes[y][x].density_superfluid = initial_superfluid_density
            new_nodes[y][x].entropy = initial_entropy

            if y == 0:  # Bottom boundary
                new_nodes[y][x].heat_flux = 0.2  

            elif y == num_nodes_y - 1:  # Top boundary
                new_nodes[y][x].heat_flux = 0.2  
                    
            else:

                new_nodes[y][x].heat_flux = 0

    return new_nodes

"""
    temperature_field = np.array([[node.temperature for node in row] for row in nodes])
    grad_T_y, grad_T_x = np.gradient(temperature_field)

    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)                    
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, grad_T_x, grad_T_y)
    plt.title('Temperature Gradient Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()


    pressure_field = np.array([[node.pressure for node in row] for row in nodes])

    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)
"""

"""
    grad_P_y, grad_P_x = np.gradient(pressure_field)

    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, grad_P_x, grad_P_y)
    plt.title('Pressure Gradient Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()


    #return new_nodes

"""







# Apply the temperature averaging step for multiple iterations
iterations = 1000
 # Number of times we want to average temperatures
for i in range(iterations):
    nodes = temperature_averaging_step(nodes)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature distribution
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar and labels
plt.colorbar(label='Temperature (K)')
plt.title(f'Temperature Distribution after {iterations} Averaging Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Show the plot
plt.show()










































"""
# Assuming delta_x and delta_y are grid spacings in x and y directions
delta_x = 1.0  # Assuming uniform grid spacing
delta_y = 1.0

# Compute the temperature gradient for interior nodes using central differences
for y in range(1, num_nodes_y - 1):
    for x in range(1, num_nodes_x - 1):
        # Central difference for the x-direction (dT/dx)
        grad_T_x = (nodes[y][x+1].temperature - nodes[y][x-1].temperature) / (2 * delta_x)
        
        # Central difference for the y-direction (dT/dy)
        grad_T_y = (nodes[y+1][x].temperature - nodes[y-1][x].temperature) / (2 * delta_y)

        # Store the temperature gradient as a tuple (dT/dx, dT/dy)
        nodes[y][x].temperature_gradient = (grad_T_x, grad_T_y)

# Handle boundary nodes (forward or backward differences)
# Bottom row (excluding corners)
for x in range(1, num_nodes_x - 1):
    grad_T_x = (nodes[0][x+1].temperature - nodes[0][x-1].temperature) / (2 * delta_x)
    grad_T_y = (nodes[1][x].temperature - nodes[0][x].temperature) / delta_y  # Forward difference

    # Store temperature gradient in the node
    nodes[0][x].temperature_gradient = (grad_T_x, grad_T_y)

# Top row (excluding corners)
for x in range(1, num_nodes_x - 1):
    grad_T_x = (nodes[num_nodes_y-1][x+1].temperature - nodes[num_nodes_y-1][x-1].temperature) / (2 * delta_x)
    grad_T_y = (nodes[num_nodes_y-1][x].temperature - nodes[num_nodes_y-2][x].temperature) / delta_y  # Backward difference

    # Store temperature gradient in the node
    nodes[num_nodes_y-1][x].temperature_gradient = (grad_T_x, grad_T_y)

# Left column (excluding corners)
for y in range(1, num_nodes_y - 1):
    grad_T_x = (nodes[y][1].temperature - nodes[y][0].temperature) / delta_x  # Forward difference
    grad_T_y = (nodes[y+1][0].temperature - nodes[y-1][0].temperature) / (2 * delta_y)

    # Store temperature gradient in the node
    nodes[y][0].temperature_gradient = (grad_T_x, grad_T_y)

# Right column (excluding corners)
for y in range(1, num_nodes_y - 1):
    grad_T_x = (nodes[y][num_nodes_x-1].temperature - nodes[y][num_nodes_x-2].temperature) / delta_x  # Backward difference
    grad_T_y = (nodes[y+1][num_nodes_x-1].temperature - nodes[y-1][num_nodes_x-1].temperature) / (2 * delta_y)

    # Store temperature gradient in the node
    nodes[y][num_nodes_x-1].temperature_gradient = (grad_T_x, grad_T_y)

# Handle corner nodes (e.g., using forward/backward differences)
# Bottom-left corner
grad_T_x = (nodes[0][1].temperature - nodes[0][0].temperature) / delta_x  # Forward x
grad_T_y = (nodes[1][0].temperature - nodes[0][0].temperature) / delta_y  # Forward y
nodes[0][0].temperature_gradient = (grad_T_x, grad_T_y)

# Bottom-right corner
grad_T_x = (nodes[0][num_nodes_x-1].temperature - nodes[0][num_nodes_x-2].temperature) / delta_x  # Backward x
grad_T_y = (nodes[1][num_nodes_x-1].temperature - nodes[0][num_nodes_x-1].temperature) / delta_y  # Forward y
nodes[0][num_nodes_x-1].temperature_gradient = (grad_T_x, grad_T_y)

# Top-left corner
grad_T_x = (nodes[num_nodes_y-1][1].temperature - nodes[num_nodes_y-1][0].temperature) / delta_x  # Forward x
grad_T_y = (nodes[num_nodes_y-1][0].temperature - nodes[num_nodes_y-2][0].temperature) / delta_y  # Backward y
nodes[num_nodes_y-1][0].temperature_gradient = (grad_T_x, grad_T_y)

# Top-right corner
grad_T_x = (nodes[num_nodes_y-1][num_nodes_x-1].temperature - nodes[num_nodes_y-1][num_nodes_x-2].temperature) / delta_x  # Backward x
grad_T_y = (nodes[num_nodes_y-1][num_nodes_x-1].temperature - nodes[num_nodes_y-2][num_nodes_x-1].temperature) / delta_y  # Backward y
nodes[num_nodes_y-1][num_nodes_x-1].temperature_gradient = (grad_T_x, grad_T_y)



# Extract x, y coordinates, and temperature gradient values for plotting
X = np.array([[node.x for node in row] for row in nodes])
Y = np.array([[node.y for node in row] for row in nodes])
grad_T_x = np.array([[node.temperature_gradient[0] for node in row] for row in nodes])
grad_T_y = np.array([[node.temperature_gradient[1] for node in row] for row in nodes])

# Create the quiver plot
plt.figure(figsize=(8, 6))
plt.quiver(X, Y, grad_T_x, grad_T_y, scale=5, color='blue')  # scale controls the arrow length
plt.title('Temperature Gradient (Vector Field)')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')
plt.gca().invert_yaxis()  # Invert y-axis for consistency with grid layout
plt.grid(True)
plt.show()
















# Function to apply FDM and calculate next temperature values
def fdm_step(nodes, alpha, dt, dx, dy):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature) for x, node in enumerate(row)] for y, row in enumerate(nodes)]
    
    for y in range(1, num_nodes_y - 1):  # Avoid boundary nodes
        for x in range(1, num_nodes_x - 1):  # Avoid boundary nodes
            T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
            T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
            new_nodes[y][x].temperature = nodes[y][x].temperature + alpha * dt * (T_xx + T_yy)
    
    # Return the updated grid
    return new_nodes

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, alpha, dt, dx, dy)


# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()









# Extract temperature values into a 2D array
heat_flux_data = np.array([[node.heat_flux for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(heat_flux_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Pressure (Pa)')

# Add labels and title
plt.title('Temperature Distribution with Global Initial Condition')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()












# Extract temperature values into a 2D array
pressure_data = np.array([[node.pressure for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(pressure_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Pressure (Pa)')

# Add labels and title
plt.title('Temperature Distribution with Global Initial Condition')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()
"""
