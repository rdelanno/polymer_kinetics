'''
# ----Avrami parameter determination 1.0 - R. Delannoy 2024--------------------
# -----------------------------------------------------------------------------

# This program allows for the determination of the Avrami parameter from
# isothermal DSC curves. 

# Files to process have to be located in the same folder and formatted like the 
# following example. To do so, save each curve in the software with the 
# 'Export in Tabular Format to Text File' tool found in the Import/Export drop
# panel

#          Index             t      Heatflow            Tr
#            [#]             [s]             [mW]          [°C]
#             68  6.80000e+001 -7.03360e-001  6.50000e+001
#             69  6.90000e+001 -7.03594e-001  6.50000e+001
#             70  7.00000e+001 -7.05011e-001  6.50000e+001
#             71  7.10000e+001 -7.06709e-001  6.50000e+001
#             ..  ..           ..             ..
#           1798  1.79800e+003 -5.45504e-001  6.50000e+001
#           1799  1.79900e+003 -4.98015e-001  6.50000e+001
# ^&]2[RD-PHBV-MELT-ISOT-65C-260124            

# Induction time has to be cut beforehand (within the DSC software or manually)
'''

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import stats
    
# --------------------------Main parameters------------------------------------
# -----------------------------------------------------------------------------

maxtime = 3000.0                    # max time (s) for general graph
[minheat,maxheat] = [-3.,0.1]       # min,max heatflow (mW) for general graph

baseline_hw = 200           # half width of window for baseline determination
target_R_sq = 0.995         # target R² for linear portion determination

# /output/ folder has to be created for the following to operate

outgraph = 'ye'     # 'yes' if you want to export the graphs
outresults = 'ye'   # 'yes' if you want to export results in a txt file
outcurves = 'ye'    # 'yes' is you want to export the curves data

# Be mindful that exporting will take more time for the script to run. A good
# practice would be to try the script first, control, and then export.

# --------------------------Functions------------------------------------------
# -----------------------------------------------------------------------------

# Finding a baseline by defining it as the 'flattest' region of the curve------
def subtract_flat_baseline(x,y,hw):
    mean2med = abs(max(y)-min(y))       #define initial mean to median and
    height = abs(max(y)-min(y))         #height values as largest span of curve
                                        
    baseline = y[0]                     #define initial baseline
    baseline_x = x[0]                   #define initial baseline center
    
    for i in range(len(x)-2*hw):        #exclude edges
        row = i + hw
        
        box =[]                         #create local data box of half width hw
        for j in range(2*hw):
            box.append(y[row-hw+j]) 
            
        mean2med_i = abs(np.mean(box) - y[row])     #dist. from mean to median
        height_i = abs(max(box)-min(box))           #height of box
        
        
        if mean2med_i < mean2med:       #if mean2med and height of box smallest
           if height_i < height:        #one found yet, then define it as new
               mean2med = mean2med_i    #baseline
               height = height_i
               baseline_x = x[row]
               baseline = y[row]
               
    return(x,y-baseline,baseline_x,baseline)

# -------------------General graph and output files----------------------------
# -----------------------------------------------------------------------------

# Figure showing both raw and baseline-corrected data--------------------------

figdata, axs = plt.subplots(1,2, sharex=True, sharey=True)

axs[0].set_title('Raw data')
axs[0].set_xlabel('Time ($s$)')
axs[0].set_ylabel('Heatflow ($mW$)')
axs[0].set_xlim(0.,maxtime)
axs[0].set_ylim(minheat,maxheat)

axs[1].set_title('Baseline corrected data')
axs[1].set_xlabel('Time ($s$)')

figcorr, axcorr = plt.subplots()

axcorr.set_xlabel('Time ($s$)')
axcorr.set_ylabel('Heatflow ($mW$)')
axcorr.set_xlim(0.,maxtime)
axs[0].set_ylim(minheat,maxheat)

# Figure showing the Avrami results--------------------------------------------

fig_avrami, ax_avrami = plt.subplots()
ax_avrami.set_xlabel('ln(time)')
ax_avrami.set_ylabel('ln(-ln(1-X_r))')

fig_param, ax_param = plt.subplots(1,2)

ax_param[0].set_xlabel('Temperature (°C)')
ax_param[0].set_ylabel('n_Avrami')

ax_param[1].set_xlabel('Temperature (°C)')
ax_param[1].set_ylabel('k_Avrami')


# Initialise general output results files--------------------------------------
 
if outresults == 'yes':
    results = np.array(['Filename','Temperature','Baseline_center','Baseline',
                        'Area','n_Avrami','k_Avrami','R²'])
    results = np.vstack((results,['-','°C','s','mW','mJ','-','-','-']))

# --------------------------Main-----------------------------------------------
# -----------------------------------------------------------------------------
#Create an output folder in the directory if there is none---------------------

if outresults == 'yes' or outgraph == 'yes' or outcurves == 'yes':
    if not 'output' in os.listdir():
        os.mkdir('output')

# Find txt files in the folder and create an array of file names---------------
list_txt = []
for file in os.listdir():
    if '.txt' in file:
        list_txt.append(file)

nb_files = len(list_txt)
             
print(nb_files,' files to process:',list_txt)

# Analyse each file one by one-------------------------------------------------

for file in list_txt:
    
    print('File name:',file,'\n')
    
    file_label = file.split('.txt')[0]          #label used for graph legend
    
# Import data from txt file----------------------------------------------------
    data = np.genfromtxt(file,skip_header=2,skip_footer=5,usecols=(1,2),unpack=True)
                                                #data from 2nd and 3rd columns  
                                                #is imported from the txt file
    [x, y] = [data[0]-data[0,0],data[1]]        
    datasize = len(x)                     
    
    axs[0].plot(x,y, label=file_label)          #plot raw data
    axs[0].legend(frameon=False)   
    
    data = np.genfromtxt(file,skip_header=2,skip_footer=3+datasize,
                         usecols=(3),unpack=True)
    temp = data[0]
        
# Find and subtract baseline---------------------------------------------------
    [x,corr_y,baseline_x,baseline] = subtract_flat_baseline(x,y,baseline_hw)
                                                #use defined function to find
                                                #flat baseline
    
    axs[1].plot(x,corr_y, label=file_label)     #plot baseline-corrected data
    axs[1].legend(frameon=False)
    
    axcorr.plot(x,corr_y, label=file_label)     #plot baseline-corrected data
    axcorr.legend(frameon=False)
    
# Integrate data with trapezium rule and calculate cumulative area-------------
#Y = (1/2)* (y(n) +y(n-1)) / (x(n)-x(n-1))    

    area = [] 
    for i in range(datasize):
        if i == 0:
            area.append(0)
        else:    
            area.append(area[i-1] + (1.0/2.0)*(corr_y[i]+corr_y[i-1])/(x[i]-x[i-1]))

# Calculate relative crystallinity over time-----------------------------------
    rel_cryst = area/area[-1]

# Find the linear part of the relative crystallinity curve---------------------

# Expected curve is sigmoidal from 0 to 1. The most linear non-zero part of such 
# curve should be found where the slope is highest. As such, 1st-order diff is 
# used to define a box of data. Limits of the box are defined by a threshold 
# value, defined as the percentage of the max value of the 1st-order diff. 
# Linear regression is then operated on this box. This loops by increasing 
# threshold until the resulting R-squared reaches its target value.

    diff1_rel_cryst = np.diff(rel_cryst,n=1)            #1st-order differential
    maxdiff = max(diff1_rel_cryst)

    R_squared = 0.0                                     #initialise R-squared
    threshold = 0.4                                     #initialise threshold
    
    while R_squared < target_R_sq:    #loop repeat until target R-squared is reached
        
        [x_cut,y_cut] = [[],[]]
        
        for i in range(len(diff1_rel_cryst)):
            if diff1_rel_cryst[i] > threshold * maxdiff: 
                x_cut.append(x[i+1])
                y_cut.append(rel_cryst[i+1])
                
        slope_r_c, intercept_r_c, r, p, std_err = stats.linregress(x_cut, y_cut)
                                                            #linear reg of box
        R_squared = r**2
        print ('R-squared:',R_squared)
        
        threshold = threshold*1.05 #slightly increase threshold after each loop
    
# Calculate the Avrami parameters----------------------------------------------

    ones = np.ones(len(x_cut))
    avrami_y = np.log(-np.log(ones-y_cut))              #ln(-ln(1-rel_cryst))
     
    n_avrami, ln_k_avrami, r_avrami, p, std_err = stats.linregress(np.log(x_cut),avrami_y)


    ax_avrami.plot(np.log(x_cut),avrami_y,label=file_label+' n='+str(n_avrami)[:4])
    ax_avrami.legend(frameon=False,fontsize='small')
    
    ax_param[0].plot(temp,n_avrami,'o')
    ax_param[1].annotate(str(n_avrami)[:4],xy=(temp,n_avrami+0.2),ha='center',size=8)
    
    ax_param[1].plot(temp,ln_k_avrami,'o')
    ax_param[1].annotate(str(ln_k_avrami)[:5],xy=(temp,ln_k_avrami+0.2),ha='center',size=8)
    
# Define local graphs----------------------------------------------------------    
# These graphs are used to control baseline and linear portion have been found

    figlocal, axlocal = plt.subplots(1,2, sharex=True)
    
    axlocal[0].set_title('Baseline correction')
    axlocal[0].set_xlabel('Time ($s$)')
    axlocal[0].set_xlim(0.,max(x))
    axlocal[0].set_ylim(-2.0,0.1)
    
    axlocal[0].axvspan(baseline_x-baseline_hw,baseline_x+baseline_hw,
                       color="y",alpha=0.3)
    axlocal[0].plot(x,y,label=file_label)
    axlocal[0].plot(x,corr_y)
    axlocal[0].plot([0.,max(x)],[0.,0.],'--',color = "k")
    axlocal[0].legend(frameon=False)
    
    axlocal[1].set_title('Relative crystallinity')
    axlocal[1].set_xlabel('Time ($s$)')
    axlocal[1].set_xlim(0.,max(x))
    axlocal[1].set_ylim(-0.2,1.1)
             
    axlocal[1].axvspan(x_cut[0],x_cut[-1], color="y", alpha=0.3)
    axlocal[1].plot(x, rel_cryst, label='rel_cryst')
    axlocal[1].plot(x[1:],100*diff1_rel_cryst, label='100*diff')
    axlocal[1].plot(x,slope_r_c*np.array(x) + intercept_r_c,
                    '--',color = "k")
    axlocal[1].legend(frameon=False)
    
    if outgraph == 'yes':
        figlocal.savefig('output/'+file_label+'_control.png',
                         dpi=300,bbox_inches='tight')



# Fill out output arrays and export local txt files----------------------------
    
    if outresults == 'yes':
        results = np.vstack((results,
                             [file_label,temp,str(baseline_x),str(baseline),
                              str(max(area)),str(n_avrami),
                              str(ln_k_avrami),str(r_avrami**2)]))
        
    if outcurves == 'yes':
        curves = ['Time','Heatflow','Cumulative area',
                  'Relative crystallinity']
        
        curves = np.vstack ((curves,['s','mW','mJ','-']))
    
        for i in range(datasize):
            curves = np.vstack ((curves,[x[i],corr_y[i],area[i],rel_cryst[i]]))
        
        np.savetxt('output/'+file_label+'_output.txt',curves,fmt='%s')                                                    
                                                        
# Print useful informations----------------------------------------------------
    print('\n','baseline:',baseline_x, baseline,'\n',
          'area (in J):',str(abs(area[-1]))[:6],'\n',
          'n:',str(n_avrami)[:4],'k:',str(ln_k_avrami),
          'R²:',str(r_avrami**2)[:7])    
    print('---------------------')

# Export output txt files and graphs-------------------------------------------
 
if outgraph == 'yes':
    fig_avrami.savefig('output/avrami.png',dpi=300,bbox_inches='tight')

if outresults == 'yes':
    np.savetxt('output/results.txt', results,fmt='%s')
