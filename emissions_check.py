from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import time

####################################################################################################
### User Input
### Select modeling domain (only 12US1 currently available):
DOMAIN     = '12US1'
### Select photochemical modeling mechanism (CB6r3_ae8, CRACMMv0.3)
MECHANISM  = 'CB6r3_ae8'
### Name of ff10 file to grid:
FF10       = 'AsphaltPaving_SmokeFlatFile_2018.csv'
### Emissions year:
YEAR       = '2018'
### Emissions source (used for output file name):
SOURCE     = 'AsphaltPaving'
####################################################################################################

####################################################################################################
### Import data
### Timezone offsets:
f1             = Dataset('./input/timezones_'+DOMAIN+'.nc', 'r')
num_rows       = f1.NROWS
num_cols       = f1.NCOLS
f1.close()
### FF10 emissions file:
FF10           = np.genfromtxt('./input/'+FF10,delimiter=",",usecols=(8))  # ann_value
FF10           = FF10[1:]
### GSPRO:
MODELSPECS     = np.genfromtxt('./input/gspro.'+MECHANISM+'_criteria.20220322.csv',delimiter=",",dtype='str',usecols=(1,3),skip_header=1)
### Generate arry of unique mechanism surrogates from GSPRO:
SURROGATES     = np.unique(MODELSPECS[:,0])
MWs            = np.zeros((len(SURROGATES)))
for i in range(len(SURROGATES)):
    for j in range(len(MODELSPECS)):
        if SURROGATES[i] == MODELSPECS[j,0]:
            MWs[i] = float(MODELSPECS[j,1])
            break
        else: pass
####################################################################################################

####################################################################################################
### Number of days accounts for leap years
if YEAR == '2000' or YEAR == '2004' or YEAR == '2008' or YEAR == '2012' or YEAR == '2016' or YEAR == '2020':
    daysinmonth = [31,29,31,30,31,30,31,31,30,31,30,31]
else:
    daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
MONTHS     = ['01','02','03','04','05','06','07','08','09','10','11','12']
####################################################################################################

####################################################################################################
#### Initiate final array:
final_array = np.zeros((12,len(SURROGATES),num_rows,num_cols)) # month, surrogates, row, col
####################################################################################################

for i in range(len(daysinmonth)):
    
    ################################################################################################    
    ### Get emissions file:
    f1        = Dataset('./output/'+SOURCE+'_'+MECHANISM+'_'+DOMAIN+'_'+YEAR+MONTHS[i]+'.nc', 'r')
    for j in range(len(SURROGATES)):
        final_array[i,j,:,:] = np.average(f1.variables[SURROGATES[j]][0:24,0,:,:],axis=0) * 86400 * \
                               daysinmonth[i] * MWs[j] * 1.10231e-6  # mol/sec --> ton/month
    f1.close()
    ################################################################################################

print("Mass of emissions in FF10 file: ",np.sum(FF10))
print("Mass of emissions in gridded files: ",np.sum(final_array[:,:,:,:]))