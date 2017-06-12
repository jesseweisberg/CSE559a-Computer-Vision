import os
import subprocess

#numEigenfaces is j

### Recognize all faces in interesting folder ###
for j in range(1, 22, 2):
    #compute eigenfaces (from directory faces)
    cmd = "../code/main --eigenfaces " + str(j) + " 30 30 neutral/list.txt neutral.face"
    os.system(cmd)

    #creating user base
    cmd = "../code/main --constructuserbase neutral.face neutral/list.txt neutral.user"
    os.system(cmd)

    for i in range(1, 25):
        if (i<10):
            cmd = "../code/main --recognizeface interesting/0" + str(i) +  ".tga neutral.user neutral.face " + str(j)

        else:
            cmd = "../code/main --recognizeface interesting/" + str(i) +  ".tga neutral.user neutral.face " + str(j)
        #cmd = "../../PanoramaSampleMac sphrWarp IMG_" + str(i) + ".tga IMG_" + str(i) + "_warp.tga 215.466 0.0 0.0"
        os.system(cmd)