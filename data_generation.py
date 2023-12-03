# import libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import pandas as pd
import scipy.integrate as integrate
import math
import csv
from scipy.interpolate import interp1d

# import classes from TOVsolver (credit: https://github.com/amotornenko/TOVsolver)
from tov import *
from constants import *

# constants
# values from pQCD and χEFT
muPQCD, nPQCD, pPQCD = 2.6, 6.47, 3823*10**-3 # pQCD X=2
ePQCD = -pPQCD+nPQCD*muPQCD
muL, nL, pL = 0.978, 0.176, 3.542*10**-3 # stiff causality constraints
eL = -pL+nL*muL
n=25 # number of values for each varied property
filename='plot_data_25'

# calculates cs2min given high-density limit
def cs2_min(muH,nH,pH,eH):
    if (muH < muL): return -1 # exclude values where calculated muH from generated values is unphysical
    if (pH < pL + 0.1): return -1 # exclude pressure values too close or less than lower value 
    if (eH < 0): return -1 # exclude unphysical energy density values
    if (pH > muH*nH/2): return -1 # exclude pressure difference cannot be reached
    if (eH-eL < pH-pL): return -1 # exclude cs2>1
    
    n = 200
    cs2lim_arr = np.arange(0.01,1,(1-0.01)/n)
    if ((pH-pL)/(eH-eL) >= math.log(muH/muL,10)/math.log(nH/nL,10)):
        delta_p_arr = [integrate.quad(lambda mu: nH*(mu/muH)**(1/cs2lim),muL,muH)[0] for cs2lim in cs2lim_arr]
        if (max(delta_p_arr)<pH-pL): return -1 # if max allowed area is too small for pressure change
    else: 
        delta_p_arr = [integrate.quad(lambda mu: nL*(mu/muL)**(1/cs2lim),muL,muH)[0] for cs2lim in cs2lim_arr]
        if (min(delta_p_arr)>pH-pL): 
            return -1 # if min allowed area is too large for pressure change

    cs2lim_interpolation = interp1d(delta_p_arr, cs2lim_arr)
    return cs2lim_interpolation(pH-pL)

# Kurkela & Komoltsev (2022) interpolation method
def kk2022interpolation(cs2lim, muH, nH, pH, eH, isMin=False):
    if (muH < muL): return # exclude values where calculated muH from generated values is unphysical
    if (pH < pL + 0.1): return # exclude pressure values too close or less than lower value 
    if (eH < 0): return # exclude unphysical energy density values
    if (pH > muH*nH/2): return # exclude pressure difference cannot be reached
    if (eH-eL < pH-pL): return # exclude cs2>1

    if (nL*(muH/muL)**(1/cs2lim)>nH): return # min causality line is increases above upper bound
    if (nH*(muL/muH)**(1/cs2lim)<nL): return # max causality line is decreases below lower bound

    mu = np.arange(muL,muH,(muH-muL)/100)
    muc = (muH**(1/cs2lim)*(cs2lim*(muL*nL-muH*nH+pH-pL)+pH-pL))/(cs2lim*(nL*(muH/muL)**(1/cs2lim)-nH))
    if (muc < 0): return # exclude xintersection being less than 0
    muc = muc**(cs2lim/(cs2lim+1))
    mu1 = mu[mu<=muc]
    mu2 = mu[mu>muc]
    if(muc <= muL or muc >= muH): return # exclude possibility that intersection is not in mu range

    n_min = np.concatenate([nL*(mu1/muL)**(1/cs2lim),(mu2/muL)**(1/cs2lim)*(cs2lim*(mu2*nH*(mu2/muH)**(1/cs2lim)-muH*nH+(pH-pL))+(pH-pL))/(cs2lim*(mu2*(mu2/muL)**(1/cs2lim)-muL))])
    n_max = np.concatenate([(mu1/muH)**(1/cs2lim)*(cs2lim*(mu1*nL*(mu1/muL)**(1/cs2lim)-muL*nL-(pH-pL))-(pH-pL))/(cs2lim*(mu1*(mu1/muH)**(1/cs2lim)-muH)),nH*(mu2/muH)**(1/cs2lim)])
    nc = n_max[0]*(mu/muL)**(1/cs2lim)

    p_min = 10**3*(np.concatenate([pL+cs2lim/(1+cs2lim)*(mu-muL*(muL/mu)**(1/cs2lim))*n_min]))
    p_max = 10**3*(np.concatenate([pL+cs2lim/(1+cs2lim)*(mu[n_min<=nc]-muL*(muL/mu[n_min<=nc])**(1/cs2lim))*n_min[n_min<=nc],pH+cs2lim/(1+cs2lim)*(mu[n_min>nc]-muH*(muH/mu[n_min>nc])**(1/cs2lim))*n_min[n_min>nc]]))
    
    e_min = np.concatenate([-p_min+10**3*mu*n_max])
    e_max = np.concatenate([-p_max+10**3*mu*n_min])
    if ((e_min[-1] > 1.1*eH*10**3) or (e_min[-1] < 0.9*eH*10**3)): return # excludes not reaching high-density energy limit with 1% error
    if ((p_max[-1] < 0.9*pH*10**3) or (p_max[-1] < 0.9*pH*10**3)): return # excludes not reaching high-density pressure limit with 1% error

    e_min = np.concatenate([[10**3*eL+0.01],e_min])
    e_max = np.concatenate([e_max,[10**3*eH+0.01]])
    p_min = np.concatenate([[10**3*pL+0.01],p_min])
    p_max = np.concatenate([p_max,[10**3*pH+0.01]])

    if (isMin):
        if ((pH-pL)/(eH-eL) >= math.log(muH/muL,10)/math.log(nH/nL,10)): return [[e_min,p_min]]
        else: return [[e_max,p_max]]
    else: return [[e_min,p_min],[e_max,p_max]]

# loads crust EOS data
crust_eos = pd.read_csv("data/crust_to_XEFT.csv")
crust_e_arr = np.array([eval(e) for e in crust_eos['E[MeV/fm**3]']])
crust_p_arr = np.array([eval(p) for p in crust_eos['P[MeV/fm**3]']])    

# maps EOS from p-ε space to M-R space
def pe_to_mr(e_arr,p_arr):
    e_eos = np.concatenate((crust_e_arr[crust_p_arr < min(p_arr)], e_arr))
    p_eos = np.concatenate((crust_p_arr[crust_p_arr < min(p_arr)], p_arr))
    tov = TOV(e_eos,p_eos,add_crust=False,plot_eos=False)
    crust_R_arr, crust_M_arr = [], []
    R_arr, M_arr = [], []
    
    for e in np.logspace(-1,2,100):
        R, M, prof = tov.solve(e)
        crust_R_arr.append(R)
        crust_M_arr.append(M)

    for e in e_arr:
        R, M, prof = tov.solve(e)
        R_arr.append(R)
        M_arr.append(M)

    return [[crust_R_arr,crust_M_arr],[R_arr,M_arr]]

# varied properties of the high-density limit
pH_arr = np.logspace(math.log(10**3*pL,10),4,n)
eH_arr = np.logspace(math.log(1600,10),4.5,n)
dn = (0.16*50-0.16*10)/n
nH_arr = np.arange(0.16*10,0.16*50+dn,dn)
muH_arr = (eH_arr[:,np.newaxis,np.newaxis] + pH_arr[np.newaxis,:,np.newaxis])/nH_arr[np.newaxis,np.newaxis,:]

# writes multidimensional array into 2D csv file
# first index: EOS number
# second index: 0 gives original EOS number and 1 gives the EOS (original EOS number of -1 is the pQCD-χEFT interpolation)
# third index: cs2lim (0 is cs2lim = 1, 1 is cs2lim = 1/3, 2 is cs2lim = cs2min)
# fourth index: min (0) or max (1) boundary
# fifth index: M-R crust (0), M-R EOS (1), or p-ε EOS (2)
# sixth index: x-axis (0) or y-axis (1)
# seventh index: all the data values
with open(filename+'.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    eos_num = 0

    # writes pQCD-χEFT interpolation of cs2lim = 1 into csv 
    interpolation = kk2022interpolation(1,muPQCD,nPQCD,pPQCD,ePQCD)
    temp_mr_pe_data = []
    idx_mmax_list = []
    if interpolation is not None: 
        for i in range(2):
            pe = interpolation[i]
            temp_pe_to_mr = pe_to_mr(pe[0],pe[1])
            writer.writerow([-1,0,i,0,0]+temp_pe_to_mr[0][0])
            writer.writerow([-1,0,i,0,1]+temp_pe_to_mr[0][1])
            writer.writerow([-1,0,i,1,0]+temp_pe_to_mr[1][0])
            writer.writerow([-1,0,i,1,1]+temp_pe_to_mr[1][1])
            writer.writerow([-1,0,i,2,0]+pe[0].tolist())
            writer.writerow([-1,0,i,2,1]+pe[1].tolist())

    # writes pQCD-χEFT interpolation of cs2lim = 1/3 into csv 
    interpolation = kk2022interpolation(1./3,muPQCD,nPQCD,pPQCD,ePQCD)
    temp_mr_pe_data = []
    idx_mmax_list = []
    if interpolation is not None: 
        for i in range(2):
            pe = interpolation[i]
            temp_pe_to_mr = pe_to_mr(pe[0],pe[1])
            writer.writerow([-1,1,i,0,0]+temp_pe_to_mr[0][0])
            writer.writerow([-1,1,i,0,1]+temp_pe_to_mr[0][1])
            writer.writerow([-1,1,i,1,0]+temp_pe_to_mr[1][0])
            writer.writerow([-1,1,i,1,1]+temp_pe_to_mr[1][1])
            writer.writerow([-1,1,i,2,0]+pe[0].tolist())
            writer.writerow([-1,1,i,2,1]+pe[1].tolist())

    # writes pQCD-χEFT interpolation of cs2min into csv 
    interpolation = kk2022interpolation(cs2_min(muPQCD,nPQCD,pPQCD,ePQCD),muPQCD,nPQCD,pPQCD,ePQCD,isMin=True)
    temp_mr_pe_data = []
    idx_mmax_list = []
    if interpolation is not None:  
        temp_pe_to_mr = pe_to_mr(interpolation[0][0],interpolation[0][1])
        writer.writerow([-1,2,0,0,0]+temp_pe_to_mr[0][0])
        writer.writerow([-1,2,0,0,1]+temp_pe_to_mr[0][1])
        writer.writerow([-1,2,0,1,0]+temp_pe_to_mr[1][0])
        writer.writerow([-1,2,0,1,1]+temp_pe_to_mr[1][1])
        writer.writerow([-1,2,0,2,0]+interpolation[0][0].tolist())
        writer.writerow([-1,2,0,2,1]+interpolation[0][1].tolist())

    # writes interpolation with the varying high-density limits into csv
    for idx_e in range(n):
        for idx_p in range(n):
            for idx_n in range(n):
                muH, nH, pH, eH = muH_arr[idx_e][idx_p][idx_n]/10**3, nH_arr[idx_n], pH_arr[idx_p]/10**3, eH_arr[idx_e]/10**3

                # cs2 limit of 1
                interpolation = kk2022interpolation(1,muH,nH,pH,eH)
                temp_mr_pe_data = []
                idx_mmax_list = []
                if interpolation is not None: 
                    for i in range(2):
                        pe = interpolation[i]
                        temp_pe_to_mr = pe_to_mr(pe[0],pe[1])
                        writer.writerow([eos_num,0,i,0,0]+temp_pe_to_mr[0][0])
                        writer.writerow([eos_num,0,i,0,1]+temp_pe_to_mr[0][1])
                        writer.writerow([eos_num,0,i,1,0]+temp_pe_to_mr[1][0])
                        writer.writerow([eos_num,0,i,1,1]+temp_pe_to_mr[1][1])
                        writer.writerow([eos_num,0,i,2,0]+pe[0].tolist())
                        writer.writerow([eos_num,0,i,2,1]+pe[1].tolist())

                # cs2 limit of 1/3
                interpolation = kk2022interpolation(1./3,muH,nH,pH,eH)
                temp_mr_pe_data = []
                idx_mmax_list = []
                if interpolation is not None: 
                    for i in range(2):
                        pe = interpolation[i]
                        temp_pe_to_mr = pe_to_mr(pe[0],pe[1])
                        writer.writerow([eos_num,1,i,0,0]+temp_pe_to_mr[0][0])
                        writer.writerow([eos_num,1,i,0,1]+temp_pe_to_mr[0][1])
                        writer.writerow([eos_num,1,i,1,0]+temp_pe_to_mr[1][0])
                        writer.writerow([eos_num,1,i,1,1]+temp_pe_to_mr[1][1])
                        writer.writerow([eos_num,1,i,2,0]+pe[0].tolist())
                        writer.writerow([eos_num,1,i,2,1]+pe[1].tolist())

                # cs2 minimum limit
                cs2lim = cs2_min(muH,nH,pH,eH) # minimum cs2lim
                if (cs2lim != -1): interpolation = kk2022interpolation(cs2lim,muH,nH,pH,eH,isMin=True)
                else: interpolation = None
                temp_mr_pe_data = []
                idx_mmax_list = []
                if interpolation is not None:  
                    temp_pe_to_mr = pe_to_mr(interpolation[0][0],interpolation[0][1])
                    writer.writerow([eos_num,2,0,0,0]+temp_pe_to_mr[0][0])
                    writer.writerow([eos_num,2,0,0,1]+temp_pe_to_mr[0][1])
                    writer.writerow([eos_num,2,0,1,0]+temp_pe_to_mr[1][0])
                    writer.writerow([eos_num,2,0,1,1]+temp_pe_to_mr[1][1])
                    writer.writerow([eos_num,2,0,2,0]+interpolation[0][0].tolist())
                    writer.writerow([eos_num,2,0,2,1]+interpolation[0][1].tolist())

                # prints loop info (as a progress tracker)
                print("loop info",eos_num,idx_e,idx_p,idx_n,muH,nH,pH,eH)
                eos_num += 1