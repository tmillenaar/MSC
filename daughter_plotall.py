#Importing relevant packages:
import numpy as np			#I use this for storing the data from textfiles in array
import os				#I use this for terminal commands (not very relevant in this 'light version')
import matplotlib.pyplot as plt		#Is the bread and butter of the script (used for plotting)
import math

#+++++++++++++#
#+  Modify:  +#
#+++++++++++++#

startyear = 1
start_month = 3 #[0 = jan,..., 11 = dec] the fortran code starts with april
endyear = 10
step = 1 #skip months in output data, 1 is every month, 2 is every other month etc.

Plotname = 'advDiv 2 basins plot' #Input name of savefile as string
plot_legend = True   #For many plotted lines, legend may be too long to be usefull

share_x_axis = True  #Recommended unless the values differ a lot between the two basins
share_y_axis = True  #Highly recommended


#
#---------------------------------
#        End of modifyable part






#+++++++++++++++#
#+  Plotting:  +#
#+++++++++++++++#

#Note, 12 steps (n) in h_n is 1 year

#h consists of:
plotvariables = [  'depth',  'Temperature',  'Salinity',  'Density'  ]
plot_units =    [  '[m]',  '[degC]',  '[kg/m^3]',  '[kg/m^3]'  ]
#h consists of: (    0  ,      1  ,      2      ,      3       )

rgb = list(range(0,4))
month = ['jan','feb','mar','apr','may','jun','jul','aug','sep','okt','nov','dec']
#month = ['apr','may','jun','jul','aug','sep','okt','nov','dec','jan','feb','mar']
#	   0     1     2     3     4     5     6     7     8     9    10    11
month_nr_original = 3+11-12 #First month in dataset (april in this case)

startyear += (start_month-3.)/12.
if ((start_month-3.)/12. < 0):
  ending = math.ceil(startyear)
else:
  ending = math.floor(startyear)
if (startyear<0):
  startyear = 0
  start_month = 3
  print('note, the run starts in april year 0')
  
if not os.path.exists('pics'):
    os.makedirs('pics')

#data = [np.loadtxt('h_'+str(i)+'.xy') for i in range (10*12,15*12)]
columnname = ['Main column', 'Secundary column']

#Make Temp, Salinity and Density plots
#---------------------------------
if (share_x_axis == True):
  share_x = 'col'
else:
  share_x = 'none'
  
if (share_y_axis == True):
  share_y = 'row'
else:
  share_y = 'none'
    
year = 0	
	
if (plot_legend == True):
    ny, nx = 12, int((endyear-ending))+1
    r, g, b = [(np.random.random(ny*nx)).reshape((ny, nx)) for i in range(3)]
    for i in range(12):
        for j in range(int( np.size(r)/12 )):
            r[i,j] = 1
            b[i,j] = 1
            g[i,j] = 1
            c = (np.dstack([r,g,b]))	
	
f, axarr=plt.subplots(1, 3, sharey=share_y, sharex=share_x)

for j in range(3): 					#Temp, sal and density loop
	month_nr = month_nr_original
	for i in range(int(startyear*12+1),int(endyear*12+1),step):   #timeloop
	    legend_label = 'year '+str(int(((i+2)/12)))+', '+str(month[month_nr])
	    data = np.loadtxt('Secundary_basin_output_'+str(i)+'.xy')
    # get rgb color based on month and year
	    rgb[0] = 1-0.75*float(i-startyear*12)/(endyear*12-startyear*12) #- 0.1*(month_nr)/12.
	    rgb[1] = (month_nr)/12. #- 0.005*float(i-startyear*12)/(endyear*12-startyear*12) #min division=12 (maybe 30 if step is 6)
	    if (rgb[1] < 0):
		    rgb[1] = 0
	    rgb[2] = 1-rgb[1]
	    rgb[3] = 0.6+0.4*rgb[1]#-0.4*(month_nr)/12
	    if (j==0):
		    r[month_nr,year] = rgb[0]
		    g[month_nr,year] = rgb[1]
		    b[month_nr,year] = rgb[2]  
    # Plot data of timestep i
	    axarr[j].plot(data[:,j+1],data[:,0], label=legend_label, color=(rgb[0],rgb[1],rgb[2],rgb[3]))
    # keep track of month of current iteration
	    month_nr += step
	    if (month_nr > 11):
		    month_nr = month_nr-12
		    if (j==0):
		        year+=1

	axarr[j].invert_yaxis()
	#plt.gca().invert_yaxis() ##Invert y if not using axarr but plt
	axarr[j].set_xlabel(plotvariables[j+1]+' '+plot_units[j+1])
	axarr[0].set_ylabel('Depth [m]')#+columnname[k])
    # axarr[1].set_title('Main column')
    # axarr[1].set_title('Secundary column')
	# if plot_legend == True:
		# axarr[2].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
	##rotate density values
	# plt.sca(axarr[2])
	# plt.xticks(rotation='20')
	plt.subplots_adjust(left=None, bottom=0.15, right=0.74, top=None,
                wspace=None, hspace=None)
plt.savefig('pics/'+Plotname+'.png', format='png', figsize=(11, 7), dpi=300, bbox_inches='tight')   



#Mixing plot (imld)
#---------------------------------
f, axarr2=plt.subplots(2, 1, sharex='col')
for k in range (2):
  #plt.figure(7+k)

	if ( (os.stat("main_mixing.xy").st_size != 0) or (os.stat("secundary_mixing.xy").st_size != 0) ):
	    counter = 0
	    imld = []
	    counters = []
	    if (k==0):
	        if os.stat("main_mixing.xy").st_size != 0:
	          tseries = np.loadtxt('main_mixing.xy')
	    elif (k==1):
	        if os.stat("secundary_mixing.xy").st_size != 0 :
	          tseries = np.loadtxt('secundary_mixing.xy')
	    for i in range(len(tseries[:,1])):
	        if (i > startyear*12+1 and i < endyear*12+1):
	          imld.append(tseries[i,1])
	          counters.append(counter)
	          counter += 1
	    years = [int(startyear+(3./12.)+(counters[item]+(3./12.))/12.) for item in counters]
	    #plt2 = plt.plot(years[:],imld[:])
	    axarr2[k].invert_yaxis()
	    axarr2[1].set_xlabel('Time [yr]')
	    axarr2[0].set_ylabel('Depth of unstable column [m]                                            ')
	    if (k==0):
	        axarr2[k].set_title('Main Column')
	    elif (k==1):
	        axarr2[k].set_title('Secundary column')
	    axarr2[k].plot(tseries[:,0]/12.,tseries[:,1], 's', ms=0.5)
		
#Plotting legend:
if (plot_legend == True):
    plt.figure()
    plt.xlabel('Year')
    c = (np.dstack([r,g,b]))
    plt.imshow(c, interpolation='nearest', extent=[math.floor(startyear+3./12.)-0.5,math.ceil(endyear)+0.5,12.5,0.5])
    ax = plt.axes()
    ax.set_title('Legend')
    ax.set_yticks(range(1,len(month)+1))
    #ax.set_xticks(range(int(startyear),endyear+1))
    ax.set_yticklabels(month)
    plt.savefig('pics/'+'Legend.png', format='png', dpi=300, bbox_inches='tight')   
plt.savefig('pics/'+Plotname+str(' Mixing')+'.png', format='png', figsize=(11, 7), dpi=300, bbox_inches='tight')   
plt.show()




