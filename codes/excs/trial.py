#!/usr/bin/env python3
import pyten as ptn
import sys,os
pp=os.path.dirname(os.path.abspath(__file__))
pp = os.path.dirname(pp)
sys.path.append(pp)          #used this, be careful
from src import helpers, utils, DMRG 

fold = "/project/th-scratch/m/Mert.Kurttutan/trial.txt"

f = open(fold, 'a')
writer = csv.writer(f, delimiter=',')
writer.writerow([pp])
f.close()  