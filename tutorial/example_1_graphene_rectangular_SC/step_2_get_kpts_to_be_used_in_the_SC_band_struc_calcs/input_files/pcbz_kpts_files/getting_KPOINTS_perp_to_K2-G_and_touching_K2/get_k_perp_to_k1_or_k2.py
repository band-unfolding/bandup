#! /usr/bin/env python
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--out", help="")
parser.add_argument("--kpt", type=int, choices=[1,2], help="")
parser.add_argument("--a0", type=float, help="")
parser.add_argument("--alpha", type=float, help="")
args = parser.parse_args()
a0=args.a0
alpha=args.alpha

fourpi=4.0*np.pi
twopi=2.0*np.pi
sqrt3 = np.sqrt(3.0)

b1 = np.array([(2.0/sqrt3)*comp for comp in [0.5,-sqrt3/2.0,0.0]])
b2 = np.array([(2.0/sqrt3)*comp for comp in [0.5,+sqrt3/2.0,0.0]])
G_K1 = b1/3.0 + 2.0*b2/3.0 
G_K2 = -b1/3.0 + b2/3.0
v_perp_G_K1 = np.array([0.5,-sqrt3/2.0,0.0]) 
v_perp_G_K2 = np.array([1.0,0.0,0.0]) 


def perp_to_K1_G(alpha,a0):
    start = G_K1 - alpha*v_perp_G_K1*(a0/twopi)
    end = G_K1 + alpha*v_perp_G_K1*(a0/twopi)
    return start, end

def perp_to_K2_G(alpha,a0):
    start = G_K2 - alpha*v_perp_G_K2*(a0/twopi)
    end = G_K2 + alpha*v_perp_G_K2*(a0/twopi)
    return start, end

if args.kpt == 1:
    kpt = "K1"
    start, end  = perp_to_K1_G(alpha,a0)
if args.kpt == 2:
    kpt = "K2"
    start, end  = perp_to_K2_G(alpha,a0)

with open(args.out,'w') as f:
    header = str(a0) + " !Line perpend. to " + kpt + "-G and touching " + kpt
    f.write(header+"\n")
    f.write("25"+"\n")
    f.write("Line-mode   "+str(-alpha)+"\n")
    f.write("Cartesian"+"\n")
    str_start = str(start[0]) + " " + str(start[1]) + " " + str(start[2]) + "    1 ! " + str(-alpha)
    f.write(str_start+"\n")
    str_end = str(end[0]) + " " + str(end[1]) + " " + str(end[2]) + "    1 ! " + str(alpha)
    f.write(str_end+"\n")


