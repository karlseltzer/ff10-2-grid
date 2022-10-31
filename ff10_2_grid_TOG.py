from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import time

startTime = datetime.now()

####################################################################################################
### ff10_2_grid_TOG.py: A python script that generates CMAQ-ready emissions from a ff10 file.
### Note that this is a very basic script: it only works for VOC-only FF10 files, only generates
### monthly files (i.e., all days in a month have the same emissions; no day-of-week variation),
### and does not allow spatial allocation, temporal allocation, or speciation that varies by
### county or state for an individual SCC.
####################################################################################################

####################################################################################################
### User Input
### Select modeling domain (only 12US1 currently available):
DOMAIN     = '12US1'
### Select photochemical modeling mechanism (CB6R3_AE8, CRACMMv0.3)
MECHANISM  = 'CB6R3_AE8'
### Name of ff10 file to grid:
FF10       = 'VCPy_Pesticides_SmokeFlatFile_EQUATESv2_wPSS_2018.csv'
### Emissions year:
YEAR       = '2018'
### Emissions source (used for output file name):
SOURCE     = 'VCPy_Pesticides'
####################################################################################################

####################################################################################################
### Import data
### Timezone offsets:
f1             = Dataset('./input/timezones_'+DOMAIN+'.nc', 'r')
TIMEZONES      = f1.variables['TZOFFSET'][0,0,:,:]      # Row, Col
num_rows       = f1.NROWS
num_cols       = f1.NCOLS
f1.close()
### FF10 emissions file:
FF10           = np.genfromtxt('./input/'+FF10,delimiter=",",usecols=(1,5,8))#,skip_header=2)  # FIPS, SCC, ann_value
FF10           = FF10[1:,:]
### Generate array of unique SCCs from FF10 file:
SCCs           = np.unique(FF10[:,1])
### Temporal allocation files:
DAILY          = np.genfromtxt('./input/daily_allocation_profiles.csv',delimiter=",",skip_header=1)
DAILY[:,1:]    = DAILY[:,1:] / (1/24)       # Normalize to 1
MONTHLY        = np.genfromtxt('./input/monthly_allocation_profiles.csv',delimiter=",",skip_header=1)
MONTHLY[:,1:]  = MONTHLY[:,1:]
### Allocation information:
ALLOCATION     = np.genfromtxt('./input/allocation_assignments.csv',delimiter=",",skip_header=1)
#### GSPRO:
#GSPRO          = np.genfromtxt('./input/gspro.'+MECHANISM+'_criteria.20220322.csv',delimiter=",",usecols=(0,2,3),skip_header=1)
#MODELSPECS     = np.genfromtxt('./input/gspro.'+MECHANISM+'_criteria.20220322.csv',delimiter=",",dtype='str',usecols=(1),skip_header=1)
#### Generate arry of unique mechanism surrogates from GSPRO:
#SURROGATES     = np.unique(MODELSPECS[:])
#### GSCNV:
#GSCNV          = np.genfromtxt('./input/gscnv.'+MECHANISM+'_criteria.20220322.csv',usecols=(2,3),delimiter=",",skip_header=1)
### GSPRO:
GSPRO          = np.genfromtxt('./input/gspro.'+MECHANISM+'_CRITERIA_VOC.CMAQ.2022-10-26.txt',usecols=(0,3,4))
MODELSPECS     = np.genfromtxt('./input/gspro.'+MECHANISM+'_CRITERIA_VOC.CMAQ.2022-10-26.txt',dtype='str',usecols=(2))
### Generate arry of unique mechanism surrogates from GSPRO:
SURROGATES     = np.unique(MODELSPECS[:])
### GSCNV:
GSCNV          = np.genfromtxt('./input/gscnv.'+MECHANISM+'_CRITERIA_VOC.CMAQ.2022-10-26.txt',usecols=(2,3))
####################################################################################################

####################################################################################################
### Number of days accounts for leap years
if YEAR == '2000' or YEAR == '2004' or YEAR == '2008' or YEAR == '2012' or YEAR == '2016' or YEAR == '2020':
    daysinmonth = [31,29,31,30,31,30,31,31,30,31,30,31]
    julian = [1,32,61,91,121,153,183,214,245,275,306,336]
else:
    daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
    julian = [1,32,60,91,121,152,182,213,244,274,305,335]
MONTHS     = ['01','02','03','04','05','06','07','08','09','10','11','12']
####################################################################################################

### Add code here to check that all assigned information in ALLOCATION is available

################################################################################################
### Generate daily profile w/ timezone offsets
bulk_timezone_offset  = np.zeros((len(SURROGATES),25,num_rows,num_cols))
temp_offset           = np.zeros((24*2))

for i in range(len(SCCs)):
    target_daily     = ALLOCATION[ALLOCATION[:,0]==int(SCCs[i])]
    target_daily     = DAILY[DAILY[:,0]==target_daily[0,2]]
    temp_offset[:24] = target_daily[0,1:]
    temp_offset[24:] = target_daily[0,1:]
    
    for j in range(num_rows):
        for k in range(num_cols):
            offset = TIMEZONES[j,k]
            if np.isnan(offset): offset = 0
            else:
                bulk_timezone_offset[i,0:24,j,k] = temp_offset[24-offset:24+(24-offset)]
                bulk_timezone_offset[i,-1,j,k]   = bulk_timezone_offset[i,0,j,k]
################################################################################################

####################################################################################################
#### Initiate final array:
final_array = np.zeros((12,len(SURROGATES),25,1,num_rows,num_cols)) # surrogates, timesteps, layer, row, col
####################################################################################################

####################################################################################################
### Loop through FF10 file and grid:
for i in range(len(FF10)):

    data_temp      = ALLOCATION[ALLOCATION[:,0]==FF10[i,1]]
    target_spec    = GSPRO[GSPRO[:,0]==data_temp[0,1]]
    target_surro   = MODELSPECS[GSPRO[:,0]==data_temp[0,1]]
    target_monthly = MONTHLY[MONTHLY[:,0]==data_temp[0,3]]
    target_monthly = target_monthly[0,1:]
    target_spatial = data_temp[0,4]
    TOG2VOC        = GSCNV[GSCNV[:,0]==data_temp[0,1]]
    TOG            = FF10[i,2] * TOG2VOC[0,1] * 907185    # tons per year --> grams per year
    
    ################################################################################################    
    ### Pull desired timezone_offset from bulk_timezeone_offset:
    for j in range(len(SCCs)):
        if SCCs[j]==FF10[i,1]:
            timezone_offset = bulk_timezone_offset[j,:,:,:]
        else: pass
    ################################################################################################

    ################################################################################################    
    ### Translate % monthly allocation to per second values:
    for j in range(len(target_monthly)):
        target_monthly[j] = target_monthly[j] / (daysinmonth[j] * 86400)
    ################################################################################################

    ################################################################################################    
    ### Generate state/county FIPS:
    if FF10[i,0] / 10000 < 1.0:
        statefips = '0'+str(int(FF10[i,0]/1000))
        temp_cty  = int(FF10[i,0] - int(FF10[i,0]/1000) * 1000)
        if temp_cty / 10 < 1.0:
            countyfips = '00'+str(temp_cty)
        elif temp_cty / 100 < 1.0 and temp_cty / 10 > 1.0:
            countyfips = '0'+str(temp_cty)
        else:
            countyfips = str(temp_cty)
    else:
        statefips = str(int(FF10[i,0]/1000))
        temp_cty  = int(FF10[i,0] - int(FF10[i,0]/1000) * 1000)
        if temp_cty / 10 < 1.0:
            countyfips = '00'+str(temp_cty)
        elif temp_cty / 100 < 1.0 and temp_cty / 10 > 1.0:
            countyfips = '0'+str(temp_cty)
        else:
            countyfips = str(temp_cty)
    ################################################################################################

    ################################################################################################    
    ### Get sub-county spatial allocation:
    if int(statefips) <= 30:                  # UNCOMMENT WHEN RUNNING CONUS
        f1    = Dataset('/work/MOD3DEV/kseltzer/gridmasks/12US1_gridmasks/GRIDMASK_12US1_COUNTY_'+str(int(target_spatial))+'_GROUP1.nc', 'r')
    else:
        f1    = Dataset('/work/MOD3DEV/kseltzer/gridmasks/12US1_gridmasks/GRIDMASK_12US1_COUNTY_'+str(int(target_spatial))+'_GROUP2.nc', 'r')
    try:
        SPATIAL   = f1.variables['POP_FIPS_'+statefips+countyfips][0,0,:,:]
    except KeyError:
        print('Spatial variable not found for: '+statefips+countyfips+'; SCC: ',int(FF10[i,1]))
        continue
    f1.close()
    ################################################################################################

    ################################################################################################
    #### Initiate temp array:
    temp_array = np.zeros((len(SURROGATES),25,1,num_rows,num_cols)) # surrogates, timesteps, layer, row, col
    ### Add emissions to final_array:
    for j in range(len(target_surro)):
        for k in range(len(SURROGATES)):
            if target_surro[j] == SURROGATES[k]:
                temp_array[k,:,0,:,:] += TOG * SPATIAL[:,:] * target_spec[j,1] * timezone_offset[:,:,:] / target_spec[j,2]
            else: pass
    for j in range(len(target_monthly)):
        final_array[j,:,:,:,:,:] += temp_array[:,:,:,:,:] * target_monthly[j]
    ################################################################################################

for i in range(len(daysinmonth)):

    ################################################################################################
    ### Generate the TFLAG variables
    tflag_array = np.zeros((25,len(SURROGATES),2))
    
    tflag_array[:-1,:,0] = int(str(int(YEAR)*1000+julian[i]))
    tflag_array[-1,:,0]  = int(str(int(YEAR)*1000+julian[i]+1))

    for j in range(24):
        tflag_array[j,:,1] = int(j*10000)
    tflag_array[-1,:,1]    = int(0)
    ################################################################################################
   
    ################################################################################################
    f1 = Dataset('./output/'+SOURCE+'_'+MECHANISM+'_'+DOMAIN+'_'+YEAR+MONTHS[i]+'.nc','w',format='NETCDF4_CLASSIC')
        
    TSTEP                    = f1.createDimension('TSTEP',25)
    DATETIME                 = f1.createDimension('DATETIME',2)
    LAY                      = f1.createDimension('LAY',1)
    VAR                      = f1.createDimension('VAR',len(SURROGATES))
    ROW                      = f1.createDimension('ROW',num_rows)
    COL                      = f1.createDimension('COL',num_cols)
    
    f1.createVariable('TFLAG',np.int32,('TSTEP','VAR','DATETIME'))

    f1.variables['TFLAG'].units     = '<YYYYDDD,HHMMSS>'
    f1.variables['TFLAG'].long_nam  = 'TFLAG'
    f1.variables['TFLAG'].var_desc  = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'
    f1.variables['TFLAG'][:]        = tflag_array[:]

    for j in range(len(SURROGATES)):
        spaces = 16 - len(SURROGATES[j])
        f1.createVariable(SURROGATES[j],np.float32,('TSTEP','LAY','ROW','COL'))
        f1.variables[SURROGATES[j]].units     = 'moles/s         '
        f1.variables[SURROGATES[j]].long_name = SURROGATES[j]+' '*spaces
        if j == 0:
            VARLIST_str = SURROGATES[j]+' '*spaces
        else:
            VARLIST_str = VARLIST_str+(SURROGATES[j]+' '*spaces)
        spaces = 66 - len(SURROGATES[j])
        f1.variables[SURROGATES[j]].var_desc = 'Model species '+SURROGATES[j]+' '*spaces
        f1.variables[SURROGATES[j]][:]      = final_array[i,j,:,:,:,:]

    # Global Attributes
    f1.IOAPI_VERSION         = 'Compatible with I/O API v3.2'+(' '*52)
    f1.EXEC_ID               = SOURCE+' '+YEAR+' '+DOMAIN+(' '*(78-len(SOURCE+YEAR+DOMAIN)))
    f1.FTYPE                 = 1
    f1.CDATE                 = int(time.strftime('%Y%m%d'))
    f1.CTIME                 = int(time.strftime('%H%M%S'))
    f1.WDATE                 = int(time.strftime('%Y%m%d'))
    f1.WTIME                 = int(time.strftime('%H%M%S'))
    f1.SDATE                 = int(YEAR)*1000+julian[i]
    f1.STIME                 = 0
    f1.TSTEP                 = 10000
    f1.NTHIK                 = 1
    f1.NCOLS                 = num_cols
    f1.NROWS                 = num_rows
    f1.NLAYS                 = 1
    f1.NVARS                 = len(SURROGATES)
    f1.GDTYP                 = 2
    f1.P_ALP                 = 33.
    f1.P_BET                 = 45.
    f1.P_GAM                 = -97.
    f1.XCENT                 = -97.
    f1.YCENT                 = 40.
    f1.XORIG                 = -2556000.
    f1.YORIG                 = -1728000.
    f1.XCELL                 = 12000.
    f1.YCELL                 = 12000.
    f1.VGTYP                 = -9999
    f1.VGTOP                 = 5000.
    f1.VGLVLS                = 1., 0.
    f1.GDNAM                 = DOMAIN+'_'+str(num_cols)+'x'+str(num_rows)+(' '*(14-len(DOMAIN+str(num_cols)+str(num_rows))))
    f1.UPNAM                 = 'MRGGRID         '
    f1.VARLIST               = VARLIST_str
    f1.FILEDESC              = SOURCE+' Emissions, '+DOMAIN+', '+YEAR+(' '*(66-len(SOURCE+DOMAIN+YEAR)))
    f1.HISTORY               = 'Created by Karl Seltzer - seltzer.karl@epa.gov'+(' '*34)

    f1.close()
    ################################################################################################
        
print("Time to generate speciated, gridded emissions: ",datetime.now() - startTime)