import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl
import sys

if len(sys.argv)>1:
    testfolder = sys.argv[1]
else:
    testfolder = os.getcwd()

#-3
print(testfolder)
if os.path.isdir(testfolder):  #CONDITION: Is it a folder? If yes go on
    #print(testfolder)
    count=0
    errorfiles=[]
    for file in os.listdir(testfolder): #CONDITION: Is there more than 1 error files?
        if file.startswith("ErrorL1"):
            count=count+1
            errorfiles.append(file)
    if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
        errorfiles.sort() #Order the files
        #Check where you are
        print("You are in the folder "+testfolder+" where you have "+str(count)+" L1 errors")
        #Defining the errors and mesh parameters vectors
        numberofnodes = np.zeros(count)
        errors = np.zeros((count,3))
        #Opening the file to write the rates of convergence
        fid = open(testfolder+"/convergence.tex",'w')
        print("  N     L1 error H     Order H     L1 Error Hu     Order Hu     L1 Error Hv     Order Hv")
        fid.write("  N   & L1 error H  &  Order H & L1 Error Hu  &  Order Hu & L1 Error Hv  &  Order Hv\\ \n")  
        for indi in range(count): #Process the files
            try:
                fileID = open(testfolder+"/"+errorfiles[indi],'r')
                numberofnodes[indi] = int(errorfiles[indi][-8:-4])
                for indl, line in enumerate(fileID): # iterate over each line
                    if indl==0:
                        h, hu, hv = line.split() # split it by whitespace
                        errors[indi,0] = float(h) 
                        errors[indi,1] = float(hu)
                        errors[indi,2] = float(hv) 
                fileID.close()    
            except:
                errors[indi,:]=float("NaN")
            if indi>0:
                order = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(numberofnodes[indi]/numberofnodes[indi-1])) for j in range(3)])
            else:
                order = np.zeros(3)
            print(int(numberofnodes[indi]),"     ",format(errors[indi,0], '.3e'),"     ", format(order[0], '.3f'),"     ",format(errors[indi,1], '.3e'),"     ", format(order[1], '.3f'),"     ",format(errors[indi,2], '.3e'),"     ", format(order[2], '.3f'))
            fid.write(" "+str(int(numberofnodes[indi]))+"  &   "+format(errors[indi,0], '.3e')+"  &  "+format(order[0], '.3f')+"  &  "+format(errors[indi,1], '.3e')+" & "+format(order[1], '.3f')+"  &  "+format(errors[indi,2], '.3e')+" & "+format(order[2], '.3f')+" \\\\ \n")  
        fid.close()
        #Plot
        fig=pl.figure()
        pl.loglog(numberofnodes,errors[:,0],"-*",linewidth=1.5,label='Error H')
        pl.grid()
        pl.loglog(numberofnodes,errors[:,1],"-+",linewidth=1.5,label='Error Hu')
        pl.loglog(numberofnodes,errors[:,2],"-o",linewidth=1.5,label='Error Hv')
        # pl.loglog(numberofnodes,errors[2,0]/numberofnodes[2]**(-1)*numberofnodes**(-1),"--",linewidth=1,label="1st order")
        # pl.loglog(numberofnodes,errors[2,0]/numberofnodes[2]**(-2)*numberofnodes**(-2),"--",linewidth=1,label="2nd order")
        # pl.loglog(numberofnodes,errors[2,0]/numberofnodes[2]**(-3)*numberofnodes**(-3),"--",linewidth=1,label="3rd order")
        # pl.loglog(numberofnodes,errors[2,0]/numberofnodes[2]**(-4)*numberofnodes**(-4),"--",linewidth=1,label="4th order")
        pl.loglog(numberofnodes[2:],1.1*errors[-1,1]/numberofnodes[-1]**(-4)*numberofnodes[2:]**(-4),"--",linewidth=1,label="4th order")
        pl.loglog(numberofnodes[2:],1.1*errors[-1,1]/numberofnodes[-1]**(-5)*numberofnodes[2:]**(-5),"--",linewidth=1,label="5th order")
        pl.legend(loc='lower left')
        params = {'mathtext.default': 'regular' }   
        pl.xlabel("$N_{elements}$")
        pl.ylabel("$L^1$ error")
        pl.savefig(testfolder+"/convergence.pdf", format="pdf", bbox_inches="tight")
        pl.savefig(testfolder+"/convergence.png",dpi=600)
        pl.show(block=False)
