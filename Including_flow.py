import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd


# Constants
num_nodes_x = 60  # Number of nodes in the x direction
num_nodes_y = 40  # Number of nodes in the y direction
dx = 0.001  # Spatial step in the x direction
dy = 0.001  # Spatial step in the y direction
q_top = 0.2  # Heat flux at the top boundary
q_bottom = 0.2  # Heat flux at the bottom boundary
dt = 0.000000001  # Time step for FDM update
num_timesteps = 10  # Total number of timesteps

# Node class
class Node:
    def __init__(self, x, y, temperature=1.9, density_superfluid=162.9, density_normal=147.5, csv_file='helium_Property_Table.csv', entropy=756.95, viscosity=2.5*10**-6, pressure = 0.876,v_s=(0.0,0.0),v_n=(0.0,0.0),grad_P=(0.0,0.0),grad_T=(0.0,0.0)):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.density_superfluid = density_superfluid
        self.density_normal = density_normal
        self.entropy = entropy
        self.viscosity = viscosity
        self.pressure = pressure
        self.v_s = v_s
        self.v_n = v_n
        self.grad_T = grad_T
        self.grad_P = grad_P
    

        df = pd.read_csv(csv_file, header=None)
        
        # Assuming the CSV format is [temperature, density, entropy]
        self.temperatures = df[0].values  # First column
        self.densities = df[1].values      # Second column
        self.entropies = df[2].values       # Third column

        # Create the interpolation function
        self.density_interpolator = interp1d(self.temperatures, self.densities, kind='linear', fill_value='extrapolate')
        self.entropy_interpolator = interp1d(self.temperatures, self.entropies, kind='linear', fill_value='extrapolate')


    def update_properties(self):
        # Update properties based on the current temperature
        self.density_superfluid = self.calculate_density_superfluid()
        self.density_normal = self.calculate_density_normal()
        self.entropy = self.calculate_entropy()


    def calculate_density_normal(self):
        # Example: Simple relationship; modify based on your system
        return (self.density_interpolator(self.temperature))  # Adjust as needed

    def calculate_density_superfluid(self):
        rho_n = self.density_interpolator(self.temperature) # Example: Simple relationship; modify based on your system
        rho_s = (rho_n/((self.temperature/2.17)**5.6)) - rho_n
        return rho_s  # Adjust as needed


    def calculate_entropy(self):
        # Example: Simple relationship; modify based on your system
        return (self.entropy_interpolator(self.temperature))  # Adjust as needed

    def calculate_k(self):

        L = num_nodes_y * (dy**2)  # Length of the grid in meters
        total_density = self.density_superfluid + self.density_normal
        k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        return k

    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        total_density = self.density_superfluid + self.density_normal
        alpha = self.calculate_k() / (total_density * c_p)
        return alpha



def initialize_grid():
    nodes = []
    for y in range(num_nodes_y):
        row = []
        for x in range(num_nodes_x):
            # Set pressure to 1.0 atm for the leftmost column, otherwise set to 0.1 atm
            pressure = 101325 if x == 0 else 10132.5 #pressure differential
            node = Node(x, y, temperature=1.9, pressure=pressure,density_normal=147.5, density_superfluid=162.9, entropy=756.95,v_s=[0.0,0.0],v_n=[0.0,0.0],grad_P= (0, 0), grad_T= (0, 0))
            row.append(node)
        nodes.append(row)
    return nodes


# Apply Neumann boundary conditions to update boundary nodes only
def apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy):
    # Bottom boundary
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node below the boundary
        k_local = nodes[0][x].calculate_k()  # Use local k value for boundary node
        nodes[0][x].temperature = T_interior + (q_top * dy) / k_local  # Apply Neumann BC for top boundary

    # Top boundary
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # First interior node above the boundary
        k_local = nodes[num_nodes_y - 1][x].calculate_k()  # Use local k value for boundary node
        nodes[num_nodes_y - 1][x].temperature = T_interior + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

# Update interior nodes based on the FDM method
def update_interior_nodes(nodes, dt):
    new_temperatures = np.array([[node.temperature for node in row] for row in nodes])

    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx = (nodes[y][x + 1].temperature - 2 * nodes[y][x].temperature + nodes[y][x - 1].temperature) / (dx ** 2)
            T_yy = (nodes[y + 1][x].temperature - 2 * nodes[y][x].temperature + nodes[y - 1][x].temperature) / (dy ** 2)
            alpha = nodes[y][x].calculate_alpha()
            new_temperatures[y][x] += alpha * (T_xx + T_yy) * dt  # Update using alpha and FDM

    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            nodes[y][x].temperature = new_temperatures[y][x]

    for y in range(num_nodes_y):
        nodes[y][0].temperature = nodes[y][1].temperature  # Left boundary (no flux)
        nodes[y][num_nodes_x - 1].temperature = nodes[y][num_nodes_x - 2].temperature  # Right boundary (no flux)



def grad_T(nodes):
    # Create a 2D array of temperatures
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    
    # Calculate the gradients
    grad_T_y, grad_T_x = np.gradient(temperatures)
    
    # Assign gradients back to each node
    for y, row in enumerate(nodes):
        for x, node in enumerate(row):
            node.grad_T = (grad_T_x[y, x], grad_T_y[y, x])  # Update node's grad_T



def grad_P(nodes):
    # Create a 2D array of pressures
    pressures = np.array([[node.pressure for node in row] for row in nodes])
    
    # Calculate the gradients
    grad_P_y, grad_P_x = np.gradient(pressures)
    
    # Assign gradients back to each node
    for y, row in enumerate(nodes):
        for x, node in enumerate(row):
            node.grad_P = (grad_P_x[y, x], grad_P_y[y, x])  # Update node's grad_P

 


def v_s(nodes, dt):

    grad_T_x = np.array([[node.grad_T[0] for node in row] for row in nodes])
    grad_T_y = np.array([[node.grad_T[1] for node in row] for row in nodes])

    grad_P_x = np.array([[node.grad_P[0] for node in row] for row in nodes])
    grad_P_y = np.array([[node.grad_P[1] for node in row] for row in nodes])

    density_superfluid = np.array([[node.density_superfluid for node in row] for row in nodes])
    density_normal = np.array([[node.density_normal for node in row] for row in nodes])
    entropy = np.array([[node.entropy for node in row] for row in nodes])
    
    # Calculate total density for each node
    total_density = density_superfluid + density_normal

    v_s_x = ((density_superfluid * entropy * grad_T_x - (density_superfluid/total_density) * (grad_P_x) )*dt) / density_superfluid
    v_s_y = ((density_superfluid * entropy * grad_T_y - (density_superfluid/total_density) * (grad_P_y) )*dt) / density_superfluid #gravity vector has been excluded as it is assumed we are looking at a birdseye crossection
    
    for y in range(len(nodes)):
        for x in range(len(nodes[y])):
            nodes[y][x].v_s = (v_s_x[y, x], v_s_y[y, x])

    return [v_s_x, v_s_y]



def v_n(nodes, dt):

    viscosity=2.5*10**-6

    grad_T_x = np.array([[node.grad_T[0] for node in row] for row in nodes])
    grad_T_y = np.array([[node.grad_T[1] for node in row] for row in nodes])

    grad_P_x = np.array([[node.grad_P[0] for node in row] for row in nodes])
    grad_P_y = np.array([[node.grad_P[1] for node in row] for row in nodes])

    density_superfluid = np.array([[node.density_superfluid for node in row] for row in nodes])
    density_normal = np.array([[node.density_normal for node in row] for row in nodes])
    entropy = np.array([[node.entropy for node in row] for row in nodes])


    
    # Calculate total density for each node
    total_density = density_superfluid + density_normal




#################################################
###Calculating new normal velocity in x-direction
#################################################


    #Bottom Nodes:
    for x in range(2,num_nodes_x-1):
        laplacian_bottom = nodes[1][x].v_n[0] + nodes[0][x+1].v_n[0] + nodes[0][x-1].v_n[0] - 4* (nodes[0][x].v_n[0])# Access the y-component of v_n
        v_n_x_bottom = ((-density_superfluid[0][x]*entropy[0][x]*grad_T_x[0][x] - (density_normal[0][x]/total_density[0][x])*grad_P_x[0][x] + (viscosity*laplacian_bottom)/(dy**2))*dt)/density_normal[0][x]
        nodes[0][x].v_n[0] += v_n_x_bottom


    #Top Nodes    
    for x in range(2,num_nodes_x-1):
        laplacian_top = nodes[num_nodes_y-2][x].v_n[0] + nodes[num_nodes_y-1][x+1].v_n[0] + nodes[num_nodes_y-1][x-1].v_n[0] - 4* (nodes[num_nodes_y-1][x].v_n[0])# Access the y-component of v_n
        v_n_x_top = ((-density_superfluid[num_nodes_y-1][x]*entropy[num_nodes_y-1][x]*grad_T_x[num_nodes_y-1][x] - (density_normal[num_nodes_y-1][x]/total_density[num_nodes_y-1][x])*grad_P_x[num_nodes_y-1][x] + (viscosity*laplacian_top)/(dy**2))*dt)/density_normal[num_nodes_y-1][x]
        nodes[num_nodes_y-1][x].v_n[0] += v_n_x_top


    
    #Corner Nodes:
    for y in range(num_nodes_y):
        for x in range(num_nodes_x):
            
            #bottom left
            if x==0 and y==0:
                laplacian_bottom_left = (nodes[0][1].v_n[0] + nodes[1][0].v_n[0] - 4*nodes[0][0].v_n[0])
                v_n_x_bottom_left = ((-density_superfluid[0][0]*entropy[0][0]*grad_T_x[0][0] - (density_normal[0][0]/total_density[0][0])*grad_P_x[0][0] + (viscosity*laplacian_bottom_left)/(dy**2))*dt)/density_normal[0][0]
                nodes[0][0].v_n[0] += v_n_x_bottom_left
                print(laplacian_bottom_left)


            #top left
            if x==0 and y==num_nodes_y-1:
                laplacian_top_left = (nodes[num_nodes_y-2][0].v_n[0] + nodes[num_nodes_y-1][1].v_n[0] - 4*nodes[num_nodes_y-1][0].v_n[0])
                v_n_x_top_left = ((-density_superfluid[num_nodes_y-1][0]*entropy[0][num_nodes_y-1]*grad_T_x[num_nodes_y-1][0] - (density_normal[num_nodes_y-1][0]/total_density[num_nodes_y-1][0])*grad_P_x[num_nodes_y-1][0] + (viscosity*laplacian_top_left)/(dy**2))*dt)/density_normal[num_nodes_y-1][0]
                nodes[num_nodes_y-1][0].v_n[0] += v_n_x_top_left

            #bottom right
            if x==num_nodes_x-1 and y==0:
                laplacian_bottom_right = (nodes[0][num_nodes_x-2].v_n[0] + nodes[1][num_nodes_x-1].v_n[0] - 4*nodes[0][num_nodes_x-1].v_n[0])
                v_n_x_bottom_right = ((-density_superfluid[0][num_nodes_x-1]*entropy[0][num_nodes_x-1]*grad_T_x[0][num_nodes_x-1] - (density_normal[0][num_nodes_x-1]/total_density[0][num_nodes_x-1])*grad_P_x[0][num_nodes_x-1] + (viscosity*laplacian_bottom_right)/(dy**2))*dt)/density_normal[0][num_nodes_x-1]
                nodes[0][num_nodes_x-1].v_n[0] += v_n_x_bottom_right


            
            #top right
            if x==num_nodes_x-1 and y==num_nodes_y-1:
                laplacian_top_right = (nodes[num_nodes_y-1][num_nodes_x-2].v_n[0] + nodes[num_nodes_y-2][num_nodes_x-1].v_n[0] - 4*nodes[num_nodes_y-1][num_nodes_x-1].v_n[0])
                v_n_x_top_right = ((-density_superfluid[num_nodes_y-1][num_nodes_x-1]*entropy[num_nodes_y-1][num_nodes_x-1]*grad_T_x[num_nodes_y-1][num_nodes_x-1] - (density_normal[num_nodes_y-1][num_nodes_x-1]/total_density[num_nodes_y-1][num_nodes_x-1])*grad_P_x[num_nodes_y-1][num_nodes_x-1] + (viscosity*laplacian_top_right)/(dy**2))*dt)/density_normal[num_nodes_y-1][num_nodes_x-1]
                nodes[num_nodes_y-1][num_nodes_x-1].v_n[0] += v_n_x_top_right


    #Left Central Nodes:

    for y in range(1,num_nodes_y-1):
        laplacian_left = (nodes[y+1][0].v_n[0] + nodes[y-1][0].v_n[0] + nodes[y][1].v_n[0] - 4*nodes[y][0].v_n[0])
        v_n_x_left = ((-density_superfluid[y][0]*entropy[y][0]*grad_T_x[y][0] - (density_normal[y][0]/total_density[y][0])*grad_P_x[y][0] + (viscosity*laplacian_left)/(dy**2))*dt)/density_normal[y][0]
        nodes[y][0].v_n[0] += v_n_x_left


    #right Central Nodes:
    for y in range(1,num_nodes_y-1):
        laplacian_right = (nodes[y+1][num_nodes_x-1].v_n[0] + nodes[y-1][num_nodes_x-1].v_n[0] + nodes[y][num_nodes_x-1].v_n[0] - 4*nodes[y][num_nodes_x-1].v_n[0])
        v_n_x_right = ((-density_superfluid[y][num_nodes_x-1]*entropy[y][num_nodes_x-1]*grad_T_x[y][num_nodes_x-1] - (density_normal[y][num_nodes_x-1]/total_density[y][num_nodes_x-1])*grad_P_x[y][num_nodes_x-1] + (viscosity*laplacian_right)/(dy**2))*dt)/density_normal[y][num_nodes_x-1]
        nodes[y][num_nodes_x-1].v_n[0] += v_n_x_right


    
    #central nodes

    for y in range(2, num_nodes_y - 2):
        for x in range(2, num_nodes_x - 2):
            laplacian = (nodes[y+1][x].v_n[0]+nodes[y-1][x].v_n[0]+nodes[y][x+1].v_n[0]+nodes[y][x-1].v_n[0] - 4*nodes[y][x].v_n[0])
            v_n_x = ((-density_superfluid[y][x]*entropy[y][x]*grad_T_x[y][x] - (density_normal[y][x]/total_density[y][x])*grad_P_x[y][x] + (viscosity*laplacian)/(dy**2))*dt)/density_normal[y][x]
            nodes[y][x].v_n[0] += v_n_x




#################################################
###Calculating new normal velocity in y-direction
#################################################
    #Bottom Nodes:
    for x in range(2,num_nodes_x-2):
        laplacian_bottom = nodes[1][x].v_n[1] + nodes[0][x+1].v_n[1] + nodes[0][x-1].v_n[1] - 4* (nodes[0][x].v_n[1])# Access the y-component of v_n
        v_n_y_bottom = ((-density_superfluid[0][x]*entropy[0][x]*grad_T_y[0][x] - (density_normal[0][x]/total_density[0][x])*grad_P_y[0][x] + (viscosity*laplacian_bottom)/(dy**2))*dt)/density_normal[0][x]
        nodes[0][x].v_n[1] += v_n_y_bottom


    #Top Nodes    
    for x in range(2,num_nodes_x-2):
        laplacian_top = nodes[num_nodes_y-2][x].v_n[1] + nodes[num_nodes_y-1][x+1].v_n[1] + nodes[num_nodes_y-1][x-1].v_n[1] - 4* (nodes[num_nodes_y-1][x].v_n[1])# Access the y-component of v_n
        v_n_y_top = ((-density_superfluid[num_nodes_y-1][x]*entropy[num_nodes_y-1][x]*grad_T_y[num_nodes_y-1][x] - (density_normal[num_nodes_y-1][x]/total_density[num_nodes_y-1][x])*grad_P_y[num_nodes_y-1][x] + (viscosity*laplacian_top)/(dy**2))*dt)/density_normal[num_nodes_y-1][x]
        nodes[num_nodes_y-1][x].v_n[1] += v_n_y_top


    #Corner Nodes:
    for y in range(num_nodes_y):
        for x in range(num_nodes_x):
            
            #bottom left
            if x==0 and y==0:
                laplacian_bottom_left = (nodes[0][1].v_n[1] + nodes[1][0].v_n[1] - 4*nodes[0][0].v_n[1])
                v_n_y_bottom_left = ((-density_superfluid[0][0]*entropy[0][0]*grad_T_y[0][0] - (density_normal[0][0]/total_density[0][0])*grad_P_y[0][0] + (viscosity*laplacian_bottom_left)/(dy**2))*dt)/density_normal[0][0]
                nodes[0][0].v_n[1] += v_n_y_bottom_left


            #top left
            if x==0 and y==num_nodes_y-1:
                laplacian_top_left = (nodes[num_nodes_y-2][0].v_n[1] + nodes[num_nodes_y-1][1].v_n[1] - 4*nodes[num_nodes_y-1][1].v_n[1])
                v_n_y_top_left = ((-density_superfluid[num_nodes_y-1][0]*entropy[0][num_nodes_y-1]*grad_T_y[num_nodes_y-1][0] - (density_normal[num_nodes_y-1][0]/total_density[num_nodes_y-1][0])*grad_P_y[num_nodes_y-1][0] + (viscosity*laplacian_top_left)/(dy**2))*dt)/density_normal[num_nodes_y-1][0]
                nodes[num_nodes_y-1][0].v_n[1] += v_n_y_top_left

            #bottom right
            if x==num_nodes_x-1 and y==0:
                laplacian_bottom_right = (nodes[0][num_nodes_x-2].v_n[1] + nodes[1][num_nodes_x-1].v_n[1] - 4*nodes[0][num_nodes_x-1].v_n[1])
                v_n_y_bottom_right = ((-density_superfluid[0][num_nodes_x-1]*entropy[0][num_nodes_x-1]*grad_T_y[0][num_nodes_x-1] - (density_normal[0][num_nodes_x-1]/total_density[0][num_nodes_x-1])*grad_P_y[0][num_nodes_x-1] + (viscosity*laplacian_bottom_right)/(dy**2))*dt)/density_normal[0][num_nodes_x-1]
                nodes[0][num_nodes_x-1].v_n[1] += v_n_y_bottom_right


            
            #top right
            if x==num_nodes_x-1 and y==num_nodes_y-1:
                laplacian_top_right = (nodes[num_nodes_y-1][num_nodes_x-2].v_n[1] + nodes[num_nodes_y-2][num_nodes_x-1].v_n[1] - 4*nodes[num_nodes_y-1][num_nodes_x-1].v_n[1])
                v_n_y_top_right = ((-density_superfluid[num_nodes_y-1][num_nodes_x-1]*entropy[num_nodes_y-1][num_nodes_x-1]*grad_T_y[num_nodes_y-1][num_nodes_x-1] - (density_normal[num_nodes_y-1][num_nodes_x-1]/total_density[num_nodes_y-1][num_nodes_x-1])*grad_P_y[num_nodes_y-1][num_nodes_x-1] + (viscosity*laplacian_top_right)/(dy**2))*dt)/density_normal[num_nodes_y-1][num_nodes_x-1]
                nodes[num_nodes_y-1][num_nodes_x-1].v_n[1] += v_n_y_top_right











    #Left Central Nodes:

    for y in range(1,num_nodes_y-1):
        laplacian_left = (nodes[y+1][0].v_n[1] + nodes[y-1][0].v_n[1] + nodes[y][1].v_n[1] - 4*nodes[y][0].v_n[1])
        v_n_x_left = ((-density_superfluid[y][0]*entropy[y][0]*grad_T_y[y][0] - (density_normal[y][0]/total_density[y][0])*grad_P_y[y][0] + (viscosity*laplacian_left)/(dy**2))*dt)/density_normal[y][0]
        nodes[y][0].v_n[1] += v_n_x_left

    #right Central Nodes:
    for y in range(1,num_nodes_y-1):
        laplacian_right = (nodes[y+1][num_nodes_x-1].v_n[1] + nodes[y-1][num_nodes_x-1].v_n[1] + nodes[y][num_nodes_x-1].v_n[1] - 4*nodes[y][num_nodes_x-1].v_n[1])
        v_n_x_right = ((-density_superfluid[y][num_nodes_x-1]*entropy[y][num_nodes_x-1]*grad_T_y[y][num_nodes_x-1] - (density_normal[y][num_nodes_x-1]/total_density[y][num_nodes_x-1])*grad_P_y[y][num_nodes_x-1] + (viscosity*laplacian_right)/(dy**2))*dt)/density_normal[y][num_nodes_x-1]
        nodes[y][num_nodes_x-1].v_n[1] += v_n_x_right

    
    #central nodes

    for y in range(2, num_nodes_y - 2):
        for x in range(2, num_nodes_x - 2):
            laplacian = (nodes[y+1][x].v_n[1]+nodes[y-1][x].v_n[1]+nodes[y][x+1].v_n[1]+nodes[y][x-1].v_n[1] - 4*(nodes[y][x].v_n[1]))
            v_n_x = ((-density_superfluid[y][x]*entropy[y][x]*grad_T_x[y][x] - (density_normal[y][x]/total_density[y][x])*grad_P_x[y][x] + (viscosity*laplacian)/(dy**2))*dt)/density_normal[y][x]
            nodes[y][x].v_n[1] += v_n_x  # Ensure it's a scalar



def plot_grad_P(nodes):
    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)
    pressure = np.array([[node.pressure for node in row] for row in nodes])
    grad_P_y, grad_P_x = np.gradient(pressure)
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, grad_P_x, grad_P_y)
    plt.title('Pressure Gradient Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()

def plot_grad_T(nodes):
    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)
    temperature = np.array([[node.temperature for node in row] for row in nodes])
    grad_P_y, grad_P_x = np.gradient(temperature)
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, grad_P_x, grad_P_y)
    plt.title('Temperature Gradient Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()


def plot_v_s(nodes):
    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)
    v_s_x = np.array([[node.v_s[0] for node in row] for row in nodes])
    v_s_y = np.array([[node.v_s[1] for node in row] for row in nodes])
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, v_s_x,v_s_y)
    plt.title('Superfluid Velocity Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()




def plot_v_n(nodes):
    x_mesh = np.arange(num_nodes_x)
    y_mesh = np.arange(num_nodes_y)
    X, Y = np.meshgrid(x_mesh, y_mesh)
    v_n_x = np.array([[node.v_n[0] for node in row] for row in nodes])
    v_n_y = np.array([[node.v_n[1] for node in row] for row in nodes])
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, v_n_x,v_n_y)
    plt.title('Normal Velocity Field')
    plt.xlabel('X-axis (nodes)')
    plt.ylabel('Y-axis (nodes)')
    plt.gca().invert_yaxis()  # Invert y-axis to align with the grid
    plt.show()



# Visualization function
def visualize_results(nodes):
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    plt.imshow(temperatures, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Temperature (K)')
    plt.title(f'Temperature Distribution After {num_timesteps} Timesteps')
    plt.xlabel('X Nodes')
    plt.ylabel('Y Nodes')
    plt.show()


# Main simulation function
def run_simulation():
    nodes = initialize_grid()  # Initialize the grid

    # Iterate over the specified number of timesteps
    for t in range(num_timesteps):
        apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy)  # Update boundary nodes
        update_interior_nodes(nodes, dt)  # Update interior nodes based on FDM
        grad_P(nodes)
        grad_T(nodes)
        v_s(nodes, dt)
        v_n(nodes, dt)
  

        # Update properties for each node based on the new temperature after updating temperatures
        for row in nodes:
            for node in row:
                node.update_properties()

    plot_v_n(nodes)
    grad_T(nodes)
    grad_P(nodes)
    visualize_results(nodes)



run_simulation()

