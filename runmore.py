import os
import numpy as np
import matplotlib.pyplot as plt
import math


#To be modified:

startyear = 10        ## Starts at year 0, month 3
start_month = 3        ## [0=jan, 11=dec]
endyear = 12

plot_legend = True
run_program = True     #if False, tries to plot with existing data

varname = 'Salt_diff_amp = ' #later on, adds the value of 'values[i]' behind it as well


# # For the following line: put in a CSV string and make sure "MYVALUE" appears in setup.txt as value for the approriate variable
# # 'Values' may not be empty
#values = [0.05,0.1,0.5,1.0,2,3,4,5] 
values = [1]




#End of modificationable part 
#---------------------------------------------
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

# Modifying setup file and running:
for p in range(0,len(values)):
  if (run_program == True):
    #Remove files of last run
    os.system('rm *.xy')

    # Read in the file
    with open('setup.txt', 'r') as file :
      filedata = file.read()
    
    # Replace the target string
    filedata = filedata.replace( 'MYVALUE', str(values[p]) )
    #filedata = filedata.replace( 'MYVALUE2', str(values2[p]) )
    filedata = filedata.replace( 'tmax_in_years', str(endyear) ) #Set max nr of years in setup file
    
    # Write the file out again
    with open('setup_data.txt', 'w') as file:
      file.write(filedata)
      
    print ('')
    print ('Starting run  ' + str(p+1)+'/'+str(len(values)))
    print ('___________________________________________')
    os.system('gfortran advDiv_2_basins_DD_original.f90')
    # os.system('./a.out')
    os.system('a.exe')
    
  month  = ['jan','feb','mar','apr','may','jun','jul','aug','sep','okt','nov','dec']
  #month = ['apr','may','jun','jul','aug','sep','okt','nov','dec','jan','feb','mar']
  #          0     1     2     3     4     5     6     7     8     9    10    11
  month_nr = int(start_month)
  if (month_nr > 11):
    month_nr -= 12
  elif (month_nr < 0):
    month_nr += 12
    

  ## Modifying gnuplot script for appropriate filenames
  #with open('heatmap_original', 'r') as file :
    #script = file.read()
  
  #script = script.replace( 'Mainoutput', '"pics/'+str(varname)+'_'+str(values[p])+'MAINmixing.png"')
  #script = script.replace( 'Secundaryoutput', '"pics/'+str(varname)+'_'+str(values[p])+'SECmixing.png"')
  #with open('heatmap', 'w') as file:
    #file.write(script)




  #!!!!!!!!!!!!!!#
  #!  Plotting  !#
  #!!!!!!!!!!!!!!#

  print ('')
  print ('Plotting...')
  print ('')

  #12 steps (n) in h_n is 1 year

  #h consists of:
  plotvariables = [  'depth',  'Temperature',  'Salinity',  'Density'  ]
  plot_units = [  '[m]',  '[degC]',  '[kg/m^3]',  '[kg/m^3]'  ]
  #h consists of: (     0     ,    1  ,   2  ,      3       )
  plotvar = 2

  rgb = list(range(0,4))
  #data = [np.loadtxt('h_'+str(i)+'.xy') for i in range (10*12,15*12)]
  columnname = ['Main column', 'Secundary column']

  step = 1 #skip months in output data, 1 is every, 2 is every other etc


  year = 0
  #a = np.random.random((12, 11))
  #a[0,0]=0
  #a[1,0]=0
  #a[2,0]=0
  
  if (plot_legend == True):
    ny, nx = 12, int((endyear-ending))+1
    r, g, b = [(np.random.random(ny*nx)).reshape((ny, nx)) for i in range(3)]
    for i in range(12):
        for j in range(int( np.size(r)/12 )):
            r[i,j] = 1
            b[i,j] = 1
            g[i,j] = 1
            c = (np.dstack([r,g,b]))
      
  f, axarr=plt.subplots(2, 3, sharey='row', sharex='col')
  for k in range (2):
    for j in range(3):
      month_nr = start_month
      for i in range(int(startyear*12),int(endyear*12)+1,step):
          #print(a[:,1])
          legend_label = 'year '+str(int(((i+2)/12)))+', '+str(month[month_nr])
          if (k==0):
            data = np.loadtxt('Main_basin_output_'+str(i)+'.xy')
          elif (k==1):
            data = np.loadtxt('Secundary_basin_output_'+str(i)+'.xy')
      ##get rgb color based on month and year
      ##Color palette vibrant:
          rgb[0] = 1-0.75*float(i-startyear*12)/(endyear*12-startyear*12) #- 0.1*(month_nr)/12.
          rgb[1] = (month_nr)/12. #- 0.005*float(i-startyear*12)/(endyear*12-startyear*12) #min division=12 (maybe 30 if step is 6)
          if (rgb[1] < 0):
              rgb[1] = 0
          rgb[2] = 1-rgb[1]
          rgb[3] = 0.6+0.4*rgb[1]#-0.4*(month_nr)/12
      ##Color palette brown-green:
          #rgb[0] = 0.1+0.4*(month_nr)/12.+0.55*float(i-startyear*12)/(endyear*12-startyear*12)
          #rgb[1] = 0.2-0.15*(month_nr)/12.+0.6*float(i-startyear*12)/(endyear*12-startyear*12)#(month_nr)/12. #- 0.005*float(i-startyear*12)/(endyear*12-startyear*12) #min division=12 (maybe 30 if step is 6)
          #if (rgb[1] < 0):
              #rgb[1] = 0
          #rgb[2] = 0.0+0.0*(month_nr)/12.+0.0*(i-startyear*12)/(endyear*12-startyear*12)
          #rgb[3] = 1#0.6+0.4*rgb[1]#-0.4*(month_nr)/12.
          if (j==0 and k == 0):
            r[month_nr,year] = rgb[0]
            g[month_nr,year] = rgb[1]
            b[month_nr,year] = rgb[2]     
      ## Plot data of timestep i
          # data2 = range(0,len(data[:,0]))
          # data3 = range(0,len(data[:,0]))
          # for qp in range(0,len(data[:,0])):
            # data2[qp] = data[qp,0]
            # data3[qp] = data[qp,j+1]
          #plt1 = plt.plot(data[:,j+1],data[:,0], label=legend_label, color=(rgb[0],rgb[1],rgb[2],rgb[3]))
          axarr[k, j].plot(data[:,j+1],data[:,0], linewidth=1.0, color=(rgb[0],rgb[1],rgb[2],rgb[3]))
      ## keep track of month of current iteration
          month_nr += step
          if (month_nr > 11):
              month_nr = month_nr-12
              if (j==0 and k == 0):
                year+=1
      if (j<2):
        axarr[k,j].tick_params(axis='x', pad=11)#lower x labels to match the tilted density label
      #axarr[0,2].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
      axarr[k,j].invert_yaxis()
      #plt.gca().invert_yaxis()
      axarr[1,j].set_xlabel(plotvariables[j+1]+' '+plot_units[j+1])
      axarr[k,0].set_ylabel('Depth [m]')#+columnname[k])
      axarr[0,1].set_title('Deep basin')
      axarr[1,1].set_title('Marginal basin', y=1.0)
      plt.subplots_adjust(left=None, bottom=0.15, right=0.74, top=None,
          wspace=None, hspace=0.3)
    ##rotate density values
    plt.sca(axarr[k, 2])
    plt.xticks(rotation='20')
  plt.savefig('pics/'+str(varname)+'_'+str(values[p])+'.png', format='png', figsize=(11, 7), dpi=300, bbox_inches='tight')   
  #plotting mixing column in seperate figure
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
          axarr2[k].set_title('Deep basin')
        elif (k==1):
          axarr2[k].set_title('Marginal basin', y=1.0)
        axarr2[k].plot(tseries[:,0]/12.,tseries[:,1], linewidth=1.0, color=(69./255.,184./255.,222./255.))
        axarr2[k].plot(tseries[:,0]/12.,tseries[:,2], linewidth=1.0, color=(158./255.,61./255.,0))

  plt.savefig('pics/'+str(varname)+'_'+str(values[p])+'mixing.png', format='png', figsize=(11, 7), dpi=300, bbox_inches='tight')   

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
#os.system( 'gnuplot heatmap' )  
print ('')
print ('Done! Enjoy the plots in folder /pics :)')
print ('')
if ( len(values) == 1 ):
  plt.show()
