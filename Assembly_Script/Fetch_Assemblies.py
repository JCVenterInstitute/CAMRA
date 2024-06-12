#!/usr/bin/env python
# coding: utf-8

# In[46]:
import os
import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime
'''
USE:
    python3 From_BVBRC-JSON_getFiles.py <BVBRC-ASSEMBLY LOCATION> <GBUCKET-ASEMBLY-LOCATION> <OUTPUT LOCATION>

EXAMPLE: 
    From_BVBRC-JSON_GetFIles.py /Dany@bvbrc/home/Assembly_BatchB gs://fc-1bde9971-12a8-416d-a2ff-e314c7ade/assemblies-test/ASSEMBLY/ /Users/dmatute/Assemblies/IFAIN_Retro_BatchB

PURPOSE: 
    Once the assemblies are done, statuses can be check on BV-BRC, this script copies the assemblies locally and to the designated directory in the google bucket. 
'''



# In[13]:





# In[58]:


# Get the location of the assemblies in the BVBRC to get the JSON files which are stored locally
if len(sys.argv) != 4:
    sys.exit(1)

# Input the location of the assembly in BVBRC
assembly_BVBRC = sys.argv[1] # /DanyMB@bvbrc/home/IFAIN_Retro_BatchB
assembly_gbucket = sys.argv[2] # gs://fc-1bde9971-12a8-416d-a2ff-e314c7ade234/CAMRA-IFAIN/Retrospective/Assembly
assembly_outputdir = sys.argv[3] #/Users/dmatute/Documents/CAMRA/Bioinformatics/IFAIN
now = datetime.now()



# In[59]:


if assembly_outputdir[-1] != "/": assembly_outputdir = assembly_outputdir+"/"
sample_dir = assembly_outputdir + os.path.basename(assembly_BVBRC) +"_" + now.strftime("%d%m%y-%H%M") # /Users/dmatute/Assemblies/Assembly_BatchB_04022025-1415
os.mkdir(sample_dir)
assembly_BVBRC = assembly_BVBRC+"/Assembly"
if assembly_BVBRC[:2] != "ws:":
    assembly_BVBRC ="ws:"+assembly_BVBRC
print("> BVBRC LOC = " , assembly_BVBRC, " SAVING TO = ", sample_dir)
BVBRC_fetch_command = ['p3-cp', '-r', assembly_BVBRC, sample_dir] #/Dany@bvbrc/home/Assembly_BatchB/Assembly, /Users/dmatute/Assemblies/IFAIN_Retro_BatchB/Assembly_Batch_B04022025-1415
subprocess.run(BVBRC_fetch_command)


# In[63]:


sample_dir = sample_dir+"/Assembly"
files = os.listdir(sample_dir)

#add .json to jason files
for filename in files:
    if not os.path.isdir(sample_dir+"/"+filename) and not filename.startswith(".") and not filename.endswith(".json"):  # Exclude hidden files, and json files
        old_path = os.path.join(sample_dir, filename)
        new_path = os.path.join(sample_dir, filename + ".json")
        os.rename(old_path, new_path)

#unhide the assembly directories
for filename in files:
    if filename!= "." or filename!= ".." or filename!= ".DS_Store" and filename.isdir:  # Exclude hidden files
        old_path = os.path.join(sample_dir, filename)
        new_path = os.path.join(sample_dir, filename[1:])
        os.rename(old_path, new_path)

#Move the json to their assembly directory
files = os.listdir(sample_dir)
for filename in files:
    if filename.endswith(".json") :
        name = filename.split(".")[0]
        mv_command = ["mv", sample_dir+"/"+filename, sample_dir+"/"+name+"/"+filename]
        subprocess.run(mv_command)



# In[ ]:


# parse_command = ["python", "Parse_NewBatch_to_AssemblyMetrics.py", sample_dir]
# subprocess.run(parse_command)


# In[ ]:


gsutil_command =['gsutil', '-m', 'cp', '-r', sample_dir+ '/', assembly_gbucket +'/']
subprocess.run(gsutil_command)

