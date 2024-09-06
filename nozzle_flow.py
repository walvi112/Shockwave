import math
import sympy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def check_choked(Pe, Po, gamma):
    """ Checking if a nozzle is chocked, inputting the exit pressure, reservoir pressure, heat capacity ratio """
    Me = (((Po / Pe) ** ((gamma - 1) / gamma) - 1) * (2 / (gamma - 1))) ** 0.5
    if Me >= 1:
        return Me
    else:
        return False

def nozzle_flow_properties(area, throat_area, gamma, p0=None, t0=None):
    """ Return the mach number, pressure and temperature at a location in the nozzle, inputting the area at that location, throat area of the
     nozzle, heat capacity ratio, reservoir pressure (optional) and reservoir temperature (optional)
     If the reservoir pressure and temperature are missing, return mach number only """
    class flow_properties:
        def __init__(self, mach):
            self.mach = mach
            if p0 and t0:
                self.temperature = t0 / (1 + ((gamma - 1) / 2) * mach ** 2)
                self.pressure = p0 / ((1 + ((gamma - 1) / 2) * mach ** 2) ** (gamma / (gamma - 1)))

    if area <= 0 or throat_area <= 0 or gamma <= 0:
        # messagebox.showinfo('Error', 'Please check the values')
        # Raise error if the arguments are not physical
        raise ValueError
    else:
        mach = sympy.symbols("M")
        func = area / throat_area - (
                (1 / mach) * ((2 + (gamma - 1) * (mach ** 2)) / (1 + gamma)) ** ((gamma + 1) / (2 * (gamma - 1))))
        func_deri = sympy.diff(func, mach)
        solutions = []
        for value in [0.01, 5]:
            iteration = 1
            while abs(func.subs(mach, value)) > 0.0001:
                if iteration >= 40:
                    # messagebox.showinfo(title="No Solution", message="Could not find any root after 20 iterations.")
                    print("No solution")
                    return
                else:
                    value = value - func.subs(mach, value) / func_deri.subs(mach, value)
                    iteration += 1
            solutions.append(flow_properties(value))
    return solutions

# Input the constant: exit pressure, reservoir pressure, reservoir temperature, heat capacity ratio
Pe = 16.8727 * 10 ** 3
P0 = 20.4 * 10 ** 6
T0 = 3500
gamma = 1.22

# Read SSME nozzle shape data from csv file into pandas dataframe
my_data = pd.read_csv('SSME_nozzle.csv', delimiter=',', names=["X", "+Z", "-Z"])
# Calculate the area along the X axis, the throat and exit area
my_data["AREA"] = math.pi * my_data["+Z"] ** 2
exit_area = my_data["AREA"].max()
throat_area = my_data["AREA"].min()
throat_index = my_data["AREA"].idxmin()
length_of_data = my_data.shape[0]


# Check if the nozzle is choked, continue calculating only if the nozzle is choked
Me = check_choked(Pe, P0, gamma)
if Me:
    # Creating columns for mach, temperatue and pressure ratio
    my_data[["Mach", "Temperature", "P/Po"]] = 0
    data_dict = my_data.to_dict()
    for index in range(length_of_data):
        # Calculating flow properties at each location in the dataset
        flow_properties = nozzle_flow_properties(area=data_dict["AREA"][index], throat_area=throat_area, gamma=gamma, p0=P0, t0=T0)
        if index < throat_index:
            data_dict["Mach"][index] = flow_properties[0].mach
            data_dict["P/Po"][index] = flow_properties[0].pressure / P0
            data_dict["Temperature"][index] = flow_properties[0].temperature
        else:
            data_dict["Mach"][index] = flow_properties[1].mach
            data_dict["P/Po"][index] = flow_properties[1].pressure / P0
            data_dict["Temperature"][index] = flow_properties[1].temperature
    my_data = pd.DataFrame(data_dict)
else:
    # print("Nozzle is not choked")
    raise Exception("The nozzle is not choked")


# Plotting
x,y = np.meshgrid(my_data.X,np.arange(my_data["-Z"].min(),my_data["+Z"].max(),1))
X_element = x.shape[0]
Y_element = x.shape[1]
temperature_matrix = np.zeros((x.shape[0],y.shape[1]))
for i in range(X_element):
    for j in range(Y_element):
        # Create temperature distribution matrix
        temperature_matrix[i][j] = my_data.Temperature[j]

plt.figure()
plot = plt.subplot(1,1,1)
ax1 = plt.gca()
ax2 = ax1.twinx()
ax3 = ax1.twinx()
orig_map=plt.cm.get_cmap('plasma')
reversed_map = orig_map.reversed()
# Plot temperature contour
plt.contourf(x,y,temperature_matrix, cmap = reversed_map)
plt.colorbar(location='right')
# Plot SSME nozzle shape
ax3.plot(my_data.X,my_data["+Z"],"k")
ax3.plot(my_data.X,my_data["-Z"],"k")
ax3.fill_between(my_data.X,my_data["-Z"],my_data["-Z"].min()-5,color='w')
ax3.fill_between(my_data.X,my_data["+Z"],my_data["+Z"].max()+5,color='w')

# Plot mach number along X axis
mach_plot = ax1.plot(my_data.X,my_data["Mach"],"r",label='Mach')
# Plot pressure along X axis
pressure_plot = ax2.plot(my_data.X,my_data["P/Po"],"b",label='P/Po')
ax1.set_xlabel("X [cm]")
ax1.set_ylabel("Mach")
ax2.set_ylabel("P/Po")
ax3.axis('off')
plt.xlim(my_data.X.min(),my_data.X.max())
lns = mach_plot + pressure_plot
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc=0)
ax1.set_zorder(ax3.get_zorder()+1)
ax1.set_frame_on(False)
ax2.set_zorder(ax3.get_zorder()+1)
plt.show()
