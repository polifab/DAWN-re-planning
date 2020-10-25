#! usr/bin/python

import os

os.chdir("C:/Users/User/Documents/GitHub/RAS_DAWN/Lavoro Vollaro/object_vectors")

dawn = open("dawn_rv.m", "x")

with open('dawn_vector.txt', 'rt') as myfile, open("dawn_rv.m", "w") as dawn:
    dawn.write("function data = dawn_rv() \n")
    dawn.write("data = [ \n")
    for myline in myfile:
        line = myline.strip('\n')
        if len(line) >= 2:
            if line[1] == 'X':
                line = line.replace("X =", "")
                line = line.replace("Y =", "")
                line = line.replace("Z =", "")
                dawn.write(line)
            elif line[1] == 'V':
                line = line.replace("VX=", "")
                line = line.replace("VY=", "")
                line = line.replace("VZ=", "")
                dawn.write(line)
                dawn.write('\n')
    dawn.write('];')
