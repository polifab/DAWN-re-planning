#! usr/bin/python

import os

os.chdir("FILES DIRECTORY")

ceres = open("ceres_rv.m", "x")
vesta = open("vesta_rv.m", "x")

with open('ceres_vector.txt', 'rt') as myfile, open("ceres_rv.m", "w") as ceres:
    ceres.write("function data = ceres_rv() \n")
    ceres.write("data = [ \n")
    for myline in myfile:
        line = myline.strip('\n')
        if len(line) >= 2:
            if line[1] == 'X':
                line = line.replace("X =", "")
                line = line.replace("Y =", "")
                line = line.replace("Z =", "")
                ceres.write(line)
            elif line[1] == 'V':
                line = line.replace("VX=", "")
                line = line.replace("VY=", "")
                line = line.replace("VZ=", "")
                ceres.write(line)
                ceres.write('\n')
    ceres.write('];')

with open('vesta_vector.txt', 'rt') as myfile, open("vesta_rv.m", "w") as vesta:
    vesta.write("function data = vesta_rv() \n")
    vesta.write("data = [ \n")
    for myline in myfile:
        line = myline.strip('\n')
        if len(line) >= 2:
            if line[1] == 'X':
                line = line.replace("X =", "")
                line = line.replace("Y =", "")
                line = line.replace("Z =", "")
                vesta.write(line)
            elif line[1] == 'V':
                line = line.replace("VX=", "")
                line = line.replace("VY=", "")
                line = line.replace("VZ=", "")
                vesta.write(line)
                vesta.write('\n')
    vesta.write('];')
