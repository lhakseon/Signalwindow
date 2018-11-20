# This python script is originated from local_submit.py of LQanalzer.
# LQanalyzer was made by John Almond

timeWait=1

import os, getpass, sys
from batch_script import *

def now():
    from datetime import datetime
    return str(datetime.now().month) + str(datetime.now().day)+ str(datetime.now().hour)+str(datetime.now().minute)

sample = "SE_PU200"

isfile = os.path.isfile
join = os.path.join
#number_of_files = sum(1 for item in os.listdir(InputDir) if isfile(join(InputDir, item))) # count total # of input files to processe
number_of_files = 267

print "Job has " + str(number_of_files) + " files to process:"

number_of_cores = 100

# determine the # of files per job
nfilesperjobs= 0
for i in range(1,number_of_files):
    if not i%number_of_cores:
        nfilesperjobs+=1

if number_of_cores == 1:
    nfilesperjobs = number_of_files

if nfilesperjobs == 0:
    nfilesperjobs=1
 
files_torun = (nfilesperjobs*number_of_cores)
remainder = number_of_files - (nfilesperjobs*number_of_cores)

print "Each job will process  " + str(nfilesperjobs) + "/" + str(nfilesperjobs+1) + " files"

###################################################
## counters
###################################################
nfiles=0
count=1
total_nsamples=0
filesprocessed=0
nfiles_file=0
n_remainder_files=0
check_array = []

###################################################
# Setup work area on var tmp
###################################################
workspace = "/u/user/syrian14/NewSW-Rate/Trackisolation/CMSSW9/PixTRK-Eff/ForKNU-eff/work/"
if not (os.path.exists(workspace)):
        os.system("mkdir " + workspace)
out_end=sample
output=workspace + sample + "_" + now() + "/"
outputdir= output+ "output/"
outputdir_tmp= output+ "output_tmp/"
if not (os.path.exists(output)):
    os.system("mkdir " + output)
    print "Making tmp working directory to run Job  : " + output

if(os.path.exists(outputdir)):
    number_of_outputfiles = sum(1 for item in os.listdir(outputdir) if isfile(join(outputdir, item)))
    if  not number_of_outputfiles ==0:
       os.system("rm " + outputdir + "/*.root")
       print "Emptying output directory as this should be empty for new job"

if not (os.path.exists(outputdir)):
    os.system("mkdir " + outputdir)
    os.system("mkdir " + outputdir + "Tree")
    os.system("mkdir " + outputdir_tmp)


###################################################
## Make subjob directories
###################################################
printedworkdir =  output + "Job_[" + str(1) + "-" + str(number_of_cores) + "]/"
for i in range(1,number_of_cores+1):
    workdir =  output + "Job_" + str(i) + "/"
    if not (os.path.exists(workdir)):
            os.system("mkdir " + workdir)
            cpscript = "cp ./test.C ./test.h ./x_test.C " + output+ "Job_" + str(i) + "/."
            os.system(cpscript)
            set_each_path = "find " + output + "Job_" + str(i) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "txt_to_path/" + sample + "_%s" % (i) + ".txt" + "/g {} +"
            set_tree_temp_path = "find " + output + "Job_" + str(i) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (i) + ".root" + "/g {} +"
            set_tree_output_path = "find " + output + "Job_" + str(i) + "/." + " -name 'x_test.C' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (i) + ".root" + "/g {} +"
            os.system(set_each_path)
            os.system(set_tree_temp_path)
            os.system(set_tree_output_path)

            batchcfgfile = output+ "Job_" + str(i) + "/" +  "run.sh"
            batchconfigfile=open(batchcfgfile,'w')
            batchconfigfile.write(makeBatchConfigFile(output + "Job_" + str(i) + "/"))
            batchconfigfile.close()
            

            if i==1:
               print "making sub work directories " + printedworkdir


fr = open('./inputlist_my.txt', 'r')
outsamplename = sample

# find . -name 'test.h' -type f -exec sed -i s/PU140_2.txt/\"PU140_2.txt\"/g {} +
### specify the location of the macro for the subjob     
printedrunscript = output+ "Job_[1-" + str(number_of_cores)  + "]/runJob_[1-" + str(number_of_cores)  + "].C"

for line in fr:
    if nfiles < files_torun:
        if nfiles == 0 :
            #cpscript = "cp ./test.C ./test.h ./x_test.C ./run.sh " + output+ "Job_" + str(count) + "/."
            #os.system(cpscript)
            #set_each_path = "find " + output + "Job_" + str(count) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "txt_to_path/" + sample + "_%s" % (count) + ".txt" + "/g {} +"
            #set_tree_temp_path = "find " + output + "Job_" + str(count) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (count) + ".root" + "/g {} +"
            #set_tree_output_path = "find " + output + "Job_" + str(count) + "/." + " -name 'x_test.C' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (count) + ".root" + "/g {} +"
            #os.system(set_each_path)
            #os.system(set_tree_temp_path)
            #os.system(set_tree_output_path)
            filelist = output+ "Job_" + str(count) + "/" + sample + "_%s" % (count) + ".txt"
            fwrite = open(filelist, 'w')
            print "Making file : " + printedrunscript
            fwrite.write(line)
            filesprocessed+=1
            nfiles_file+=1
            nfiles+=1
            if files_torun == 1:
                fwrite.close()
            continue

        #End of file
        if not nfiles % nfilesperjobs:
            if not nfiles == number_of_files :
                # set counters
                nfiles_file=0
                count+=1
                # close files
                fwrite.close()
                ### Make next set of scripts
                #cpscript = "cp ./test.C ./test.h ./x_test.C ./run.sh " + output+ "Job_" + str(count) + "/."
                #os.system(cpscript)
                #set_each_path = "find " + output + "Job_" + str(count) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "txt_to_path/" + sample + "_%s" % (count) + ".txt" + "/g {} +"
                #set_tree_temp_path = "find " + output + "Job_" + str(count) + "/." + " -name 'test.h' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (count) + ".root" + "/g {} +"
                #set_tree_output_path = "find " + output + "Job_" + str(count) + "/." + " -name 'x_test.C' -type f -exec sed -i s/" + "Tree_output/" + "Tree_" + sample + "_%s" % (count) + ".root" + "/g {} +"
                #os.system(set_each_path)
                #os.system(set_tree_temp_path)
                #os.system(set_tree_output_path)
                filelist = output+ "Job_" + str(count) + "/" + sample + "_%s" % (count) + ".txt"
                fwrite = open(filelist, 'w')
                fwrite.write(line)
                filesprocessed+=1
                nfiles_file+=1
            else:
                fwrite.write(line)
                filesprocessed+=1
                nfiles_file+=1
                print "File " + filelist + " contains " + str(nfiles_file) + " files"
        else:
            fwrite.write(line)
            filesprocessed+=1
            nfiles_file+=1

        if nfiles == number_of_files-1 :
            fwrite.close()

    else:
        n_remainder_files+=1
        filelist = output+ "Job_" + str(n_remainder_files) + "/" + sample + "_%s" % (n_remainder_files) + ".txt"
        fwrite = open(filelist, 'a')
        fwrite.write(line)
        filesprocessed+=1
        fwrite.close()
    nfiles+=1
fr.close()


#################################################################### 
### Check Final input files have no duplicates
#################################################################### 
no_duplicate=False
for check in range(1, number_of_cores+1):
    filelist = output+ "Job_" + str(check) + "/" + sample + "_%s" % (check) + ".txt"
    fcheck = open(filelist, 'r')
    nsamples=0
    for line in fcheck:
        nsamples+=1
        total_nsamples+=1
        no_duplicate= True
        for s in check_array:
            if s == line :
                print "DUPLICATE file : " + s
                no_duplicate=False
                sys.exit()
        check_array.append(line)
    print  "File " + filelist + " contains " + str(nsamples) + " files"
    fcheck.close()
print "Total Number of input files = " + str(total_nsamples)

if no_duplicate:
    print "Checking for duplicates: NONE found"
else:
     print "Checking for duplicates: Duplicate files found. Check script "

print "Total number of files processed = " + str(filesprocessed)


###################################################
### Run each .C file in background
###################################################
import thread,time
start_time = time.time()


wait_sub = 1
if number_of_cores < 10:
    wait_sub = 5

print "Running Pixel track analyzer jobs for: " + getpass.getuser()
for i in range(1,number_of_cores+1):
    script_path = output+ "Job_" + str(i)
    runcommand = "qsub -q cms -V run.sh"
    if number_of_cores == 1:
        print "Running single job " 
        runcommand = "qsub -q cms -V run.sh"
        os.chdir(script_path)
        os.system(runcommand)
        os.chdir(workspace)
    else:
        if i==1:
            print "Running " + script_path
        elif i== number_of_cores:
            print "Running " + script_path
        elif i==2:
            print "......"
        os.chdir(script_path)
        os.system(runcommand)

