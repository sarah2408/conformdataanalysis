#!/usr/bin/env python2.7
'''This is a python script to plot the recorded data from the BWE Ltd. Conform 315i machine. It takes the raw  file from the machine and will give two plots relevant to understanding how the machine responded during the extrusion. Data analysed from this script is particularly useful when Conforming metal particulates
'''

import dateutil
import time
'''from Tkinter import *'''
import csv, sys
import numpy as np
import pandas as pd
import matplotlib as plt
from matplotlib.pyplot import *
from pylab import *
'''from StringIO import StringIO'''
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib import rc
rc('mathtext', default='regular')
import math
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.ticker import NullFormatter
import pickle
from mpl_toolkits.mplot3d import Axes3D                                      # importing libraries - not all are used
matplotlib.rcParams.update({'font.size': 32})




#--------Get filename from user and append .txt to the filename---------
print ("Please enter the filename without extension:(test.txt)")

file = raw_input()		#Get keyboard input from the terminal

if not file:
	file = 'test'		#Default to the test data if nothing is entered
	isfilenew = 'Yes'
else:
	print("Is the file in the new format?:")  #If a filename is entered then ask if it is in the old format - if it is then ask the user (further down the code) for the correct field headings for the plot
	isfilenew = raw_input()

file_ext = ".csv"			#Gives the data file file extension. Can be changed to csv.

datafile = file + file_ext	#Append .txt file extension to the filename






#--------Transfer data from the data file into a pandas dataframe---------

'''data = np.genfromtxt(datafile, dtype=None, delimiter='\t', names=True)''' '''initialise numpy array for holding data'''

data = pd.read_csv(datafile) # for reading .csv files. use "pd.read_table" for .tsv files


motor_amps = 'Conform Motor Current (Amps)' #initialising variable names for name-based column referencing in a dataframe
wheel_temp = 'Wheel Temperature (°C)'
abut_temp = 'Abutment 1 Temperature (°C)'
wheel_col = 'Conform Wheel Speed (rpm)'



IS08601=data['Timestamp']  #timestamp in original data is in ISO8601 datetime format

timeplaceholder=[]

for i in IS08601:                                 #for loop for parsing each entry in the timestamp series as a datetime object
    dateformat=dateutil.parser.parse(i)          #
    timeplaceholder.append(dateformat)
   

data['Timestamp']=timeplaceholder
data['delta_t'] = ( data['Timestamp'] - data['Timestamp'].iloc[0])  #converting datetime object to elapsed seconds
data['elapsed_seconds'] = data['delta_t'].dt.total_seconds()





#If the datafile is in an old format style ask the user for the heading names so that they can still be plotted without having to go back into the source code. This can be modified in the future to change the field headings automatically depending on if isfilenew = Yes or isfilenew = No.
if isfilenew == 'No' or 'no' or 'n' or 'N':
	# If the data file is old (ie. isfilenew = No) ask the user to enter the field headings. If left blank the headings will default back to what the new headings should be.
	if not wheel_col:
		print('Please enter the name of the wheel speed column heading (default: ConformWheelSpeedAvgA"Date/Time\n')
		wheel_col = input()
	else:
		wheel_col = 'Conform Wheel Speed (rpm)'
		print("Using heading for wheel speed: \t%s" % wheel_col)

	if not motor_amps:
		print('\nPlease enter the name of the motor amps column heading (default: ConformMotorAmpsAvgDateTime\n')
		motor_amps = input()
	else:
		motor_amps = 'Conform Motor Current (Amps)'
		print("Using heading for motor amps: \t%s" % motor_amps)

	if not wheel_temp:
		print('\nPlease enter the name of the wheel temperature column heading (default: ConformWheelTempDateTime\n')
		wheel_temp = input()
	else:
		wheel_temp = 'Wheel Temperature (°C)'
		print("Using heading for wheel temperature: \t%s" % wheel_temp)

	if not abut_temp:
		print('\nPlease enter the name of the abutment temperature column heading (default: ConformAbutTemp1DateTime\n')
		abut_temp = input()
	else:
		abut_temp = 'Abutment 1 Temperature (°C)'
		print("Using heading for abutment temperature: \t%s" % abut_temp)
        
        
data = data.rename(columns={ motor_amps : 'motor_current_amps', wheel_temp : 'wheel_temp_c', abut_temp : 'abut_temp_c', wheel_col : 'wheel_speed', 'elapsed_seconds' : 'elapsed_seconds'}) #renaming column headings to be consistent


#---------Calculate Abutment Stress--------
def wheel_motor(x,y):
	try:
		return x/y
	except RuntimeWarning:
		return 0

motor_voltage = 415.0	#units [V]
motor_eff = 0.8			#units [fractional]
wheel_radius = 0.15125 	#units [m]
abutment_area = 265.0	#units [mm^2]
wheel_speed_rads = data['wheel_speed']*2*math.pi/60

Wheel_Torque = motor_voltage * motor_eff * data['motor_current_amps']/wheel_speed_rads

Abut_Stress = Wheel_Torque/(wheel_radius * abutment_area)
#print Abut_Stress

outfile=data.filter(['Timestamp','elapsed_seconds','motor_current_amps', 'wheel_temp_c', 'abut_temp_c', 'wheel_speed'], axis=1) #saves an output file with new format data.
outfile['abut_stress']=Abut_Stress
outfile.to_csv('outfile.csv', index=False)

#---------------Calculate the grip length---------------
# power law function

a = 5.8078
n = 1.1692
fill_height = (Abut_Stress/a)**(1.0/n)

#grip_length = wheel_radius*asin(fill_height/wheel_radius)

#---------------Print the basic equipment graph----------------

fig = plt.figure()

fig.subplots_adjust(left=0.13, bottom=0.12, right=0.91, top=0.9, wspace=0.3, hspace=0.5)

ax = fig.add_subplot(111)

#colour = data['ConformAbutTemp1DateTime']

#ax.scatter(data['DateTime'],Abut_Stress,label = 'Abutment Stress',s=1, c='b', color = 'k')
ax.scatter(data['elapsed_seconds'], Abut_Stress,label = 'Abutment Stress Calculated', linewidth = 0, marker='o',s=5, c='k')
#ax.plot(data['DateTime'],Abut_Stress,label = 'Abutment Stress',linestyle = "-", color = 'k' )

ax.plot(data['elapsed_seconds'],data['wheel_temp_c'],label = 'Wheel Temperature',linestyle = "-", color = 'b' )
ax.plot(data['elapsed_seconds'],data['abut_temp_c'],label = 'Abutment Temperature',linestyle = "-", color = 'g' )
ax2 = ax.twinx()
ax2.plot(data['elapsed_seconds'],data['wheel_speed'],label = 'Wheel Speed',linestyle = "-", color = 'r' )
#ax.plot(data['DateTime'],fill_height*10,label = 'Fill Height x 10(mm)',linestyle = "-", color = 'c' )



#-------------Get legend names to be plotted in the same legend box-------

h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
lg = ax.legend(h1+h2, l1+l2, bbox_to_anchor=(0., 1.01, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.,fontsize=32,)
lg.get_frame().set_linewidth(0.5)

#----------------------Settings for the first figure window---------------

ax.set_ylim(0,1000)
ax2.set_ylim(0,12)
#ax.set_xlim(0,len(data['DateTime']))
ax.set_xlim(0,2000)
ax.set_title(r"Equipment Data")
ax.set_xlabel(r"Time (s)")
ax.set_ylabel(r"Abutment Stress (MPa)    Temperature($^\circ$C)")
ax2.set_ylabel(r"Wheel Speed (RPM)")

ax.get_xaxis().tick_bottom()
ax2.get_xaxis().tick_bottom()


#-------------------------Second Plot-------------------------------
#--------Fit the exponential function to the data------------

def func(x,A,B):
	return A*exp(-B*x)/x

par_guess = (50.0,3000.0)

#popt, pcov = curve_fit(func,data['wheel_speed'],Abut_Stress,p0=par_guess) #issues with NaNs or infs

#print popt
#print pcov


#---------------Wheel Speed vs Abutment Stress Graph Plot----------------

fig = plt.figure()
ax = plt.subplot2grid((2,2),(0,0),rowspan=2,colspan=2)


colour = data['abut_temp_c']

plot2 = ax.scatter(data['wheel_speed'], Abut_Stress, linewidth = 0, marker='o',s=10, c=colour)

fig.colorbar(plot2,cmap=cm.CMRmap, norm=Normalize(clip=False),label = 'Abutment Temperature')


#ax.plot(xdata,(par_guess[0]/xdata)*exp(-par_guess[1]*xdata))

#---------Settings for the second figure window---------------

ax.set_ylim(0,2000)
ax.set_xlim(0,12)
ax.set_title(r"Powder Fill Response")
ax.set_xlabel(r"Wheel Speed (RPM)")
ax.set_ylabel(r"Abutment Stress (MPa)")

#---------Plot grip length vs powder material graph---------------

fig = plt.figure()
ax = plt.subplot2grid((2,2), (0,0), rowspan=2,colspan=2)

#-----Experimental Data------
# Has the form (fill height (mm), abutment stress (MPa))

#Al1100 = [(7.0,3.0),(13.0,18.0),(18.0,32.0),(24.0,43.0),(30.0,50.0),(40.0,70.0)]
CP_Ti2 = [(7.0,30.0),(13.0,95.0),(18.0,185.0),(24.0,235.0),(30.0,310.0),(40.0,465.0),(105.0,1215.0)]
CP_Ti = [(7.0,30.0),(13.0,95.0),(18.0,185.0),(24.0,235.0),(30.0,310.0),(40.0,465.0)]
CP_Ti_abutment = [(18.0,557.0),(30.0,854.0),(40.0,1152.0)]
#Ti_64 = [(18.0,285.0),(30.0,500.0),(40.0,725.0)]
#Ti_5553 = [(18.0,330.0),(30.0,560.0),(40.0,820.0)]

#x1,y1 = zip(*Al1100)
x2,y2 = zip(*CP_Ti)
x21,y21 = zip(*CP_Ti2)
x22,y22 = zip(*CP_Ti_abutment)
#x3,y3 = zip(*Ti_64)
#x4,y4 = zip(*Ti_5553)

fill = arange(0,190,1)

grip = range(0,190,1)
print grip

def ex_fit(fill,a,n):
	return a*fill**n
weights = (1.0,9.0)

def lin_fit(fill,m):
	return m*fill

p_guess=(5.0,1.0)
#p_opt1,cov1 = curve_fit(ex_fit,x1,y1,p0=p_guess)
p_opt2,cov2 = curve_fit(ex_fit,x2,y2,p0=p_guess)
m_22,c_22 = polyfit(x22,y22,1)
#p_opt3,cov3 = curve_fit(ex_fit,x3,y3,p0=p_guess)
#p_opt4,cov4 = curve_fit(ex_fit,x4,y4,p0=p_guess)

print m_22
print c_22

die_stress = m_22*fill + c_22


#ax.scatter(x1,y1,marker="x",color="r",label="Al1100")
#ax.plot(grip,p_opt1[0]*grip**p_opt1[1],"r-",label = "Al1100 fit")

ax.scatter(x2,y2,marker="x",color="r",label="CP-Ti Torque Abutment Stress (20$^\circ$C)")
ax.scatter(x21,y21,marker="x",color="r")
ax.scatter(x22,y22,marker="^",color="b",label="DEFORM Abutment Die Stress (20$^\circ$C)")

ax.plot(grip,p_opt2[0]*grip**p_opt2[1],"r--",label = "Torque fit (20$^\circ$C)")
ax.plot(grip,die_stress,"b--",label = "Die Stress fit (20$^\circ$C)")

#ax.annotate(r'Fitted with $y=9.43x^{1.02}$', xy=(20,1800))

#ax.scatter(x3,y3,marker="x",color="b",label="Ti-6-4")
#ax.plot(grip,p_opt3[0]*grip**p_opt3[1],"b-",label = "Ti-6-4 fit")

#ax.scatter(x4,y4,marker="x",color="c",label="Ti-5553")
#ax.plot(grip,p_opt4[0]*grip**p_opt4[1],"c-",label = "Ti-5553 fit")


ax.set_ylim(0,1500)
ax.set_xlim(0,80)
ax.set_title(r"CP-Ti Truncated Workpiece DEFORM Simulations")
ax.set_xlabel(r"Grip Length")
ax.set_ylabel(r"Abutment Stress (MPa)")



h1, l1 = ax.get_legend_handles_labels()
ax.legend(h1, l1, loc="best")

'''
#---------Parameter a vs yield stress---------------
data_a = [p_opt1[0],p_opt2[0],p_opt3[0],p_opt4[0]]
data_n = [p_opt1[1],p_opt2[1],p_opt3[1],p_opt4[1]]
stress_yield = [70.0,380.0,720.0,1100.0]

ax2 = plt.subplot2grid((2,2), (0, 1))
ax2.set_title("A")
ax2.set_xlabel("Yield Stress (MPa)")
ax2.set_ylabel("A")
ax2.scatter(stress_yield,data_a,marker="x",color="k",label="Parameter a")

#---------Parameter n vs yield stress---------------
ax3 = plt.subplot2grid((2,2), (1, 1))
ax3.set_title("n")
ax3.set_xlabel("Yield Stress (MPa)")
ax3.set_ylabel("n")
ax3.scatter(stress_yield,data_n,marker="x",color="k",label="Parameter n")
'''
plt.show()
input()

#call(["clear"])	#remove crud from the terminal window after plotting has finished
