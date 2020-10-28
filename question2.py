from toolbox import euler
from sys import argv
import numpy as np
import matplotlib.pyplot as mp

print("-"*40, "Part A", "-"*40)
print("The derrivation for the unitless equation is in a pdf called question2.pdf")
print("-"*40, "Part B", "-"*40)

if len(argv) != 3:
    print(f"Error! Proper usage requires:\n {argv[0]} <initial velocity[m/s]> <launch angle [degrees]>")
    quit()

#define the original constants for later use
g = 9.8
v_0 = float(argv[1]) #m/s
theta = float(argv[2])*np.pi/180
vx_0 = v_0*np.cos(theta) #m/s
vy_0 = v_0*np.sin(theta) #m/s
t_f = 1 #s
N = 100
dt = (t_f)/N
t = np.arange(0,1,dt)

#define the drag coefficient arrray tobe used in a for loop
b_array = np.array([0.0,0.1,0.2,0.3,0.4,0.5], float)

#Initialize the graph with a title
mp.title(f"Projectile Motion with Drag.\n Initial Velocity: {v_0} Launch angle: {theta}")
mp.ylabel(r"Y Position [m]")
mp.xlabel(r"X Position [m]")

#Create a for loop to iteate
for b in b_array:

    #define scales for the iteration's drag coefficent
    V_scale = (v_0)    # in units [L/T]
    T_scale = (v_0/g)    # in units [T]
    X_scale = T_scale*V_scale # in units [L]

    #create unitary initial values
    vx_0 = vx_0/V_scale
    vy_0 = vy_0/V_scale

    #create the unitary vector for this iteration
    xi = np.array([0.0,0.0,vx_0,vy_0], float)

    #create an array to store the x and y data
    x = []
    y = []

    #now for corpus
    for time in t:
        x.append(xi[0] * X_scale)
        y.append(xi[1] * X_scale)
        xi = euler(xi, dt, b)

    mp.plot(x,y, label=f"b = {b}")

mp.legend()
mp.show()





