##############################################################################################
#-------------------------------The Three Body Problem (3D)----------------------------------#
#---------------------------------By Geoffrey M. Wakeley-------------------------------------#
##############################################################################################
##############################################################################################
########-Python Packages-#####################################################################
##############################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter
from mpl_toolkits import mplot3d
from datetime import datetime
import os
start_time = datetime.now()
my_path = os.path.abspath('C:/Users/Geoffrey Wakeley/Documents/Central Principles Series/Software/Classical Mechanics Software/Newtonian Mechanics/The Three Body Problem')

##############################################################################################
########-Constants-###########################################################################
##############################################################################################

G                      = 6.674e-11
AU                     = 1.495978707e11
M_s                    = 1.989E30
initial                = -1.5  #float(input("Input Lower Domain: "))
final                  = 1.5   #float(input("Input Higher Domain: "))
m1                     = 2.2343 #float(input("Input Mass (First Object): "))
m2                     = 1.4 #float(input("Input Mass (Second Object): "))
m3                     = 2.5 #float(input("Input Mass (Third Object): "))
x1_0                   = 0 #float(input("Input Initial X-Direction Position (First Object): "))
y1_0                   = 0 #float(input("Input Initial Y-Direction Position (First Object): "))
z1_0                   = 0
x2_0                   = 1 #float(input("Input Initial X-Direction Position (Second Object): "))
y2_0                   = 0.5 #float(input("Input Initial Y-Direction Position (Second Object): "))
z2_0                   = 1
x3_0                   = -1 #float(input("Input Initial Y-Direction Position (Third Object): "))
y3_0                   = -1.2 #float(input("Input Initial Y-Direction Position (Third Object): "))
z3_0                   = -1
vx1_0                  = 0.1 #float(input("Input Initial X-Direction Velocity (First Object): "))
vy1_0                  = -0.15
vx2_0                  = 0.1 #float(input("Input Initial X-Direction Velocity (Second Object): "))
vy2_0                  = 0.14 #float(input("Input Initial Y-Direction Velocity (Second Object): "))
vx3_0                  = 0.1 #float(input("Input Initial X-Direction Velocity (Third Object): "))
vy3_0                  = 0.15 #float(input("Input Initial Y-Direction Velocity (Third Object): "))
vz1_0                  = 0.1
vz2_0                  = -0.14
vz3_0                  = 0.1

##############################################################################################
########-Program-#############################################################################
##############################################################################################

print('Program Initalisation Successful (Packages and Constants Loaded Successfully)')

##############################################################################################
#-------------The Three Body Problem---------------------------------------------------------#
##############################################################################################

def dSdt(t, S):
    x1,y1,z1,x2,z2,y2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3 = S
    r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    r13 = np.sqrt((x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2)
    r23 = np.sqrt((x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2)
    return [vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3,
            m2/r12**3 * (x2-x1) + m3/r13**3 * (x3-x1), 
            m2/r12**3 * (y2-y1) + m3/r13**3 * (y3-y1),
            m2/r12**3 * (z2-z1) + m3/r13**3 * (z3-z1), 
            m1/r12**3 * (x1-x2) + m3/r23**3 * (x3-x2), 
            m1/r12**3 * (y1-y2) + m3/r23**3 * (y3-y2),
            m1/r12**3 * (z1-z2) + m3/r23**3 * (z3-z2),
            m1/r13**3 * (x1-x3) + m2/r23**3 * (x2-x3),
            m1/r13**3 * (y1-y3) + m2/r23**3 * (y2-y3),
            m1/r13**3 * (z1-z3) + m2/r23**3 * (z2-z3)]


time = np.linspace(0, 20, 1000)

Solution = solve_ivp(dSdt, (0,20), y0=[x1_0,y1_0,z1_0,x2_0,y2_0,z2_0,x3_0,y3_0,z3_0,
                       vx1_0,vy1_0,vz1_0,vx2_0,vy2_0,vz2_0,vx3_0,vy3_0,vz3_0],
                method = 'DOP853', t_eval=time, rtol=1e-10, atol=1e-13)

t = Solution.t
x1 = Solution.y[0]
y1 = Solution.y[1]
z1 = Solution.y[2]
x2 = Solution.y[3]
y2 = Solution.y[4]
z2 = Solution.y[5]
x3 = Solution.y[6]
y3 = Solution.y[7]
z3 = Solution.y[8]

tt = 1/np.sqrt(G*M_s/(AU)**3) # seconds
tt = tt/(60*60*24*365.25)*np.diff(t)[0] # per time step (in years)

print('Three Body Problem Solved Successfully')

##############################################################################################
########-Plots and Animations-##############################################################
##############################################################################################

def animate(i):
    line.set_data([x1[i], x2[i], x3[i]],[y1[i], y2[i], y3[i]])
    line.set_3d_properties([z1[i],z2[i],z3[i]])
    text.set_text('Time = {:.1f} Years'.format(i*tt))

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111,projection='3d')
line, = ax.plot([], [], [],'ro',lw=3,markersize=6,)
text = ax.text(-10, 10, 12, '%s', size=15, zorder=1, color='#000000')
#text = plt.text(0,1.75,'Text',fontsize=20,backgroundcolor='white',ha='center')
ax.set_xlim3d(-10,10)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10,10)
ax.set_xlabel('$x \ [a.u]$',fontsize=15)
ax.set_ylabel('$y \ [a.u]$',fontsize=15)
ax.set_zlabel('$z \ [a.u]$',fontsize=15)
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
ani.save('C:/Users/Geoffrey Wakeley/Documents/Central Principles Series/Software/Classical Mechanics Software/Newtonian Mechanics/The Three Body Problem/Three_Body_Orbit_3D.gif',writer='pillow',fps=30)

##############################################################################################
########-The Computational Time-##############################################################
##############################################################################################

end_time = datetime.now()
time=end_time-start_time
times = time.total_seconds()
print('Time: {}'.format(times))
test = '{}\n'.format(times)
name = r"Timedata"
with open(os.path.join(my_path, "The_Three_Body_Problem_3D_Computational_Time.txt" ),
mode='a') as file:
    file.write(test)

print('Program Finalised Successfully')