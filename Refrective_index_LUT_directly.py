import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib.patches as mpatches

# Below ranges come from Boucher's code that was used to generate ECHAM-HAM LUTs
dis_1_min = 0.001
dis_1_max = 25.
m_dis_1   = 101
ni_1_min  = 1.e-9
ni_1_max  = 1.00
m_ni_1    = 201
nr_1_min  = 1.33
nr_1_max  = 2.00
m_nr_1    = 101

dis_2_min = 0.4
dis_2_max = 40.
m_dis_2   = 101
ni_2_min  = 1.e-9
ni_2_max  = 1.00
m_ni_2    = 201
nr_2_min  = 1.33
nr_2_max  = 2.00
m_nr_2    = 101

# diss = 2 pi r_m/lamba, where r_m is the mean geometric radius (Grainger doc)
diss_1 = np.logspace(np.log10(dis_1_min), np.log10(dis_1_max),m_dis_1,endpoint=True)
nis_1  = np.logspace(np.log10(ni_1_min),np.log10(ni_1_max),m_ni_1,endpoint=True)
nrs_1  = np.linspace(nr_1_min,nr_1_max,m_nr_1,endpoint=True)

diss_2 = np.logspace(np.log10(dis_2_min), np.log10(dis_2_max),m_dis_2,endpoint=True)
nis_2  = np.logspace(np.log10(ni_2_min),np.log10(ni_2_max),m_ni_2,endpoint=True)
nrs_2  = np.linspace(nr_2_min,nr_2_max,m_nr_2,endpoint=True)

dataset_LUT = Dataset("./scaterring//lut_optical_properties_M7.nc")
dataset_LE_02=Dataset('./scaterring/TNO_AOP_LUT_02.nc')
dataset_LE_01=Dataset('./scaterring/TNO_AOP_LUT_01.nc')
#===Range from LUT LE information====
X_LE=np.zeros(100)
X_LE[1:]=np.logspace(-4,2,99)
MR_LE=np.linspace(1.12,1.9,40)
MI_LE=np.array([1.00000e-09,1.00000e-05,1.00000e-04,1.00000e-03,1.00000e-02,1.58489e-02,2.51189e-02,3.98107e-02,6.30957e-02,1.00000e-01,1.58489e-01,2.51189e-01,3.98107e-01,6.30957e-01,1.00000e+00])

def compareLUT(rm, ni, nr):
    
    sigma_a1 = RetrieveFromLUT(2.*np.pi*rm/0.44, ni, nr, 'sigma_1')*0.44**2
    sigma_a2 = RetrieveFromLUT(2.*np.pi*rm/0.55, ni, nr, 'sigma_1')*0.55**2
    sigma_a3 = RetrieveFromLUT(2.*np.pi*rm/0.86, ni, nr, 'sigma_1')*0.86**2
    sigma_b1 = RetrieveFromLUT(2.*np.pi*rm/0.44, ni, nr, 'sigma_2')*0.44**2
    sigma_b2 = RetrieveFromLUT(2.*np.pi*rm/0.55, ni, nr, 'sigma_2')*0.55**2
    sigma_b3 = RetrieveFromLUT(2.*np.pi*rm/0.86, ni, nr, 'sigma_2')*0.86**2

    print("Extinctions from LUT1 and LUT2 are {:4.2f} and {:4.2f}".format(sigma_a2,sigma_b2))
    print("AEs from LUT1 and LUT2 are {:4.2f} and {:4.2f}".format(-np.log10(sigma_a1/sigma_a3)/np.log10(0.44/0.86),-np.log10(sigma_b1/sigma_b3)/np.log10(0.44/0.86)))
    return
    
    
def RetrieveSigma(rm, wavel, ni, nr, LUTvar, LUT_model='ECHAM-HAM'): 
    
    if LUT_model== 'ECHAM-HAM':
        sigma = RetrieveFromLUT(2.*np.pi*rm/wavel, ni, nr, LUTvar, LUT_model)*wavel**2*1.e-12
    elif LUT_model=='LE':
        wave=wavel*1e-6
        sigma = RetrieveFromLUT(2.*np.pi*rm/wavel, ni, nr, LUTvar, LUT_model)*wave**2
    return sigma # [m2]
    
def calcAE(sigma_a, wavel_a, sigma_b, wavel_b):
    
    AE = -np.log10(sigma_a/sigma_b)/np.log10(wavel_a/wavel_b)
    
    return AE

def RetrieveMass(rm, dens, LUTvar):

    if LUTvar == 'sigma_1': # take into account units of Rm [mum] and dens [kg/m3] as well as size distribution
        mass  = 4.0/3.0*np.pi*rm**3*np.exp(9./2.*np.log(1.59)**2)*1.e-18*dens # [kg]
    if LUTvar == 'sigma_2':
        mass  = 4.0/3.0*np.pi*rm**3*np.exp(9./2.*np.log(2.00)**2)*1.e-18*dens # [kg]

    return mass # [kg]

def RetrieveCrossSection(rm, LUTvar, LUT_model='ECHAM-HAM'):

    if LUT_model== 'ECHAM-HAM':
        # rm is mean geometric size    
        if LUTvar == 'sigma_1': # take into account units of Rm [mum] as well as size distribution
            crosssection  = np.pi*rm**2*np.exp(2.*np.log(1.59)**2)*1.e-12 # [m2]
        if LUTvar == 'sigma_2':
            crosssection  = np.pi*rm**2*np.exp(2.*np.log(2.00)**2)*1.e-12 # [m2]
    elif LUT_model== 'LE':
        if LUTvar == 'sigma_1': # take into account units of Rm [mum] as well as size distribution
            crosssection  = np.pi*rm**2*np.exp(2.*np.log(1.59)**2)*1.e-12 # [m2]
        if LUTvar == 'sigma_2':
            crosssection  = np.pi*rm**2*np.exp(2.*np.log(2.00)**2)*1.e-12 # [m2]
        
    return crosssection # [m2]


def calcMEC(sigma, mass):
     
    return sigma/mass # [m2/kg]


def calcReff(Rm_1,Rm_2,mode): # checked vs Grainger on 202022/02/23
    
    # Rm is mean geometric radius, and Reff is effective particle size
    mom2 = Rm_1**2*np.exp(2.*np.log(1.59)**2)*(1.-mode)+Rm_2**2*np.exp(2.*np.log(2.00)**2)*mode
    mom3 = Rm_1**3*np.exp(9./2.*np.log(1.59)**2)*(1.-mode)+Rm_2**3*np.exp(9./2.*np.log(2.00)**2)*mode

    return mom3/mom2

def RetrieveFromLUT(dis, ni, nr, LUTvar, LUT_model):
    
    if LUT_model=='ECHAM-HAM':
        # Simple retrieval scheme based on NN
        if LUTvar == 'sigma_1':
            idis = (m_dis_1-1)*np.log10(dis/dis_1_min)/np.log10(dis_1_max/dis_1_min)
            ini  = (m_ni_1-1)*np.log10(ni/ni_1_min)/np.log10(ni_1_max/ni_1_min)
            inr  = (m_nr_1-1)*(nr-nr_1_min)/(nr_1_max-nr_1_min)
        
            if (idis < 0 or idis > m_dis_1-1):
                print("[ERROR] from RetrieveFromLut1: idis {} out of range.".format(idis))
            if (ini < 0 or ini > m_ni_1-1):
                print("[ERROR] from RetrieveFromLut1: ini {} out of range.".format(ini))    
            if (inr < 0 or inr > m_nr_1-1):
                print("[ERROR] from RetrieveFromLut1: inr {} out of range.".format(inr))
            
        if LUTvar == 'sigma_2':
            idis = np.round( (m_dis_2-1)*np.log10(dis/dis_2_min)/np.log10(dis_2_max/dis_2_min) )
            ini  = np.round( (m_ni_2-1)*np.log10(ni/ni_2_min)/np.log10(ni_2_max/ni_2_min) )
            inr  = np.round( (m_nr_2-1)*(nr-nr_2_min)/(nr_2_max-nr_2_min) )
        
            if (idis < 0 or idis > m_dis_2-1):
                print("[ERROR] from RetrieveFromLut2: idis {} out of range.".format(idis))
            if (ini < 0 or ini > m_ni_2-1):
                print("[ERROR] from RetrieveFromLut2: ini {} out of range.".format(ini))    
            if (inr < 0 or inr > m_nr_2-1):
                print("[ERROR] from RetrieveFromLut2: inr {} out of range.".format(inr))
                
        LUT = dataset_LUT.variables[LUTvar]
        # Need more sophisticated interpolation routines
        retrieved = LUT[ np.round(idis), np.round(ini),np.round(inr)] 
    elif LUT_model=='LE':
        # Simple retrieval scheme based on NN
        idis =(np.abs(X_LE-dis)).argmin()
        ini  = (np.abs(MI_LE-ni)).argmin()
        inr  = (np.abs(MR_LE-nr)).argmin()
        if (dis < X_LE[0] or dis > X_LE[-1]):
                print("[ERROR] from RetrieveFromLut: dis {} out of range.".format(dis))
                return
        if (ni < MI_LE[0] or ni > MI_LE[-1]):
                print("[ERROR] from RetrieveFromLut: ni {} out of range.".format(ni))
                return
        if (nr < MR_LE[0] or nr > MR_LE[-1]):
                print("[ERROR] from RetrieveFromLut: nr {} out of range.".format(nr))
                return
        if LUTvar == 'sigma_1':
            LUTvar_LE='ext_159'
            
        if LUTvar == 'sigma_2':
            LUTvar_LE='ext_200'
           
        LUT = dataset_LE_02.variables[LUTvar_LE]
        print('LE X: ',dis,'idis: ', idis, 'X_LE[idis]: ',X_LE[idis])
        print('LE nr: ',nr,'inr: ', inr, 'MR_LE[inr]: ',MR_LE[inr])
        print('LE ni: ',ni,'ini: ', ini, 'MI_LE[idni]: ',MI_LE[ini])
        retrieved = LUT[ ini, inr,idis]
    
    
    return retrieved

# Verify units of extinction by plotting efficiency

#wavel = 0.55
wavel=0.7016

size_701_LE=np.array([ 1.453,  3.084, 13.27 , 22.31 , 35.3  ])

#diss_1=size_701_LE[0:-1]
#diss_2=size_701_LE

Qeff_1 = np.zeros((diss_1.size))
Qeff_2 = np.zeros((diss_2.size))
Qeff_1_LE = np.zeros((diss_1.size))
Qeff_2_LE = np.zeros((diss_2.size))

Cext_1_LE = np.zeros((diss_1.size))
Cext_2_LE = np.zeros((diss_2.size))

#==Dust===
nr=1.517
ni=0.0009436

#==EC===
#nr=1.750
#ni=0.4442

test_le_lut=pd.read_csv('/Users/santiago/Documents/LE_outputs/Test_LUT_Dust.txt',delim_whitespace=True,index_col=0)
#test_le_lut=pd.read_csv('/Users/santiago/Documents/LE_outputs/Test_LUT_EC.txt',delim_whitespace=True,index_col=0)
test_le_lut['Qeff_1']=np.zeros(len(test_le_lut['ext_159']))
test_le_lut['Qeff_2']=np.zeros(len(test_le_lut['ext_200']))
#for idiss,diss in enumerate(size_701_LE[0:-2]):
for idiss,diss in enumerate(diss_1):    
    Qeff_1[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_1','ECHAM-HAM')/ \
                      RetrieveCrossSection(wavel*diss/2/np.pi, 'sigma_1')
                      
    Qeff_1_LE[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_1','LE')/ \
                      RetrieveCrossSection(wavel*diss/2/np.pi, 'sigma_1','LE')
    Cext_1_LE[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_1','LE')
    
    test_le_lut['Qeff_1'].iloc[idiss]=test_le_lut['ext_159'].iloc[idiss]*wavel**2*1.e-12/RetrieveCrossSection(wavel*test_le_lut.index[idiss]/2/np.pi, 'sigma_1')

#for idiss,diss in enumerate(size_701_LE):    
for idiss,diss in enumerate(diss_2):  
    Qeff_2[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_2','ECHAM-HAM')/ \
                      RetrieveCrossSection(wavel*diss/2/np.pi, 'sigma_2')
    Qeff_2_LE[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_2','LE')/ \
                      RetrieveCrossSection(wavel*diss/2/np.pi, 'sigma_2','LE')
    Cext_2_LE[idiss] = RetrieveSigma(wavel*diss/2/np.pi, wavel, ni, nr, 'sigma_2','LE')
    test_le_lut['Qeff_2'].iloc[idiss]=test_le_lut['ext_200'].iloc[idiss]*wavel**2*1.e-12/RetrieveCrossSection(wavel*test_le_lut.index[idiss]/2/np.pi, 'sigma_2')


fig, axs = plt.subplots(1,2,figsize=(12,5))

axs[0].plot(test_le_lut.index,test_le_lut['Qeff_1'],label='LE',linewidth=3)
axs[0].plot(diss_1,Qeff_1[:],label='ECHAM-HAM',linewidth=3)
axs[0].plot(axs[0].get_xlim(),[2,2],':')
axs[0].set_xscale('log')
axs[0].set_ylim(0,4)


axs[1].plot(test_le_lut.index,test_le_lut['Qeff_2'],label='LE',linewidth=3)
axs[1].plot(diss_2,Qeff_2[:],label='ECHAM-HAM',linewidth=3)
axs[1].plot(axs[1].get_xlim(),[2,2],':')
axs[1].set_xscale('log')
axs[1].set_ylim(0,4)

axs[0].set_xlabel('x []')
axs[0].set_ylabel('Qeff []')
axs[0].set_title('Dust LUT 1')

axs[1].set_xlabel('x []')
axs[1].set_ylabel('Qeff []')
axs[1].set_title('Dust LUT 2')
axs[1].legend()
fig.savefig('./Figures/Comparison_Qeff_ECHAM-HAM_LE_Dust.png',dpi=300)