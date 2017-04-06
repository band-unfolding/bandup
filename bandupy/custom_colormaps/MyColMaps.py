# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 15:08:58 2014

@author: manni
"""

import json
import matplotlib.cm as cm

flame = {'red':   ((0.0,  0.0, 0.0),
                   (0.337,  0.0, 0.0),
                   (0.682,  1.0, 1.0),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.102, 0.082, 0.082),
                   (0.2, 0.341, 0.341),
                   (0.339, 1.0, 1.0),
                   (0.669, 0.0, 0.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.102, 0.082, 0.082),
                   (0.2, 0.583, 0.583),
                   (0.339, 1.0, 1.0),
                   (0.669, 0.0, 0.0),
                   (1.0,  0.0, 0.0))}


def SaveColMap(ColMap, Name):
    #w = csv.writer(open(Name + ".cmap", "w"))
    #for key, val in ColMap.items():
    #    w.writerow([key, val])
    file = open(Name + ".cmap", 'w')
    json.dump(ColMap, file)
       

def LoadColMap(Filename, Name):
    #dict = {}
    #for key, val in csv.reader(open(Filename)):
    #    dict[key] = val
    file = open(Filename)
    dict = json.load(file)
    return cm.colors.LinearSegmentedColormap(Name,dict,2048)        
       
       
#SaveColMap(flame, "Flame")    
