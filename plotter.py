import numpy as np
import matplotlib.pyplot as plt
title = ""
switchs = True
with open('points.txt', 'r') as reader:
     line = reader.readline()
     while line != '':  
         if(switchs):
             title=line
             switchs = False
             line = reader.readline()
         points = line.split(' ')
         plt.scatter(float(points[0]),float(points[1]))
         line = reader.readline()

plt.title(title)
plt.show()