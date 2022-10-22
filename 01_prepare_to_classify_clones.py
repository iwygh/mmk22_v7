#%%
# This file finds objects in the SBDB and MPC databases that might be Plutinos
# and prepares for them to be classified on the Puma cluster. After running
# this file, simply upload the entire folder to Puma and run the command
# 'sbatch slurm_classify_clones.slurm'.
import aa_utilities as ut
import pandas as pd
# This is the date we downloaded the SBDB and MPC databases. It is the epoch
# that will be used to compute population-level statistics and make population-
# level figures after we have identified the Plutino sample.
year = 2022
month = 10
day = 12
hour = 0
minute = 0
second = 0
Nclones = 300
# Save settings for epoch and number of clones.
dictionary = {'year':[year],'month':[month],'day':[day],'hour':[hour],'minute':[minute],\
              'second':[second],'Nclones':[Nclones]}
df_out = pd.DataFrame(dictionary)
df_out.to_csv('epoch_settings.csv',index=False)
# Make a string of the date, for use in filenames.
datestr = ut.datestr([year,month,day,hour,minute,second],'short')
# We will need the Julian date to look up orbital elements at the epoch
# chosen for making population-level figures and statistics.
JD = ut.JDfun(year,month,day,hour,minute,second)
# Retrieve all the information we really need from the MPCORB file.
ut.dat2csv(datestr)
# Reduce entire SBDB catalog to a smaller population as described in README.docx.
# Each object that remains will be checked to see if it is Resonant or not.
Nstart,Njobs = ut.sbdb_reduce(JD,datestr)
print('starting number of objects',Nstart)
print('objects in amin,amax range',Njobs)
# Generate clones as described in README.docx.
ut.generate_clones(JD,datestr,Nclones)
# Generate as many instances of classify_clones_#.py as there are objects to
# test. Each instance will classify a separate object on the cluster.
template_file = '01a_classify_clones_template.py'
with open(template_file,'r') as file:
    template_data = file.read()
for i in range(Njobs):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'classify_clones_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+str(int(i+1)))
    with open(outfile,'w') as file:
        file.write(filedata)
# Generate a batch submission file for the cluster to classify the clones with.
# On the cluster, run "sbatch slurm_classify_clones.slurm".
template_file = '01b_slurm_classify_clones_template.slurm'
with open(template_file,'r') as file:
    template_data = file.read()
filedata = template_data
filedata = filedata.replace('NJOBS',str(Njobs))
outfile = 'slurm_classify_clones.slurm'
with open(outfile,'w') as file:
    file.write(filedata)
