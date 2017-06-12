import os
import subprocess

thresh = 60000
string = "interesting/"
fp = 0
fn = 0

for userNum in range(1, 25):


    ### Verify faces, i is index of interesting face, userNum is index of user ###
    for i in range(1, 25):
        if (i<10):
            if(userNum<10):
                cmd = "../code/main --verifyface interesting/0" + str(i) +  ".tga neutral.user neutral/0" + str(userNum) + " neutral.face " + str(thresh)
                string = string + "0" + str(i) + ".tga"
            else:
                cmd = "../code/main --verifyface interesting/0" + str(i) +  ".tga neutral.user neutral/" + str(userNum) + " neutral.face " + str(thresh)
                string = string + "0" + str(i) + ".tga"

        else:
            if(userNum<10): 
                cmd = "../code/main --verifyface interesting/" + str(i) +  ".tga neutral.user neutral/0" + str(userNum) + " neutral.face " + str(thresh)
                string = string + str(i) + ".tga"
            else:
                cmd = "../code/main --verifyface interesting/" + str(i) +  ".tga neutral.user neutral/" + str(userNum) + " neutral.face " + str(thresh)
                string = string + str(i) + ".tga"
        
        os.system(cmd)
        proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
        output=proc.communicate()[0]
        
        if ((output == "true\n") and (userNum != i)):
            fp = fp+1.0
        if ((output == "false\n") and (userNum == i)):
            fn = fn+1.0

fp = float(fp/576)
fn = float(fn/576)
print "fp = " + str(float(fp*100)) + "%"  
print "fn = " + str(float(fn*100)) + "%"    
