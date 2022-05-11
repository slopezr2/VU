import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib.patches as mpatches

dis_1_min = 0.001
dis_1_max = 25.
m_dis_1   = 100
ni_1_min  = 1.e-9
ni_1_max  = 1.00
m_ni_1    = 200
nr_1_min  = 1.33
nr_1_max  = 2.00
m_nr_1    = 100

dis_2_min = 0.4
dis_2_max = 40.
m_dis_2   = 100
ni_2_min  = 1.e-9
ni_2_max  = 1.00
m_ni_2    = 200
nr_2_min  = 1.33
nr_2_max  = 2.00
m_nr_2    = 100

dataset_LUT = Dataset("./scaterring//lut_optical_properties_M7.nc")

def compareLUT(reff, ni, nr):
    
    sigma_a1 = RetrieveFromLUT(2.*np.pi*reff/0.44, ni, nr, 'sigma_1')*0.44**2
    sigma_a2 = RetrieveFromLUT(2.*np.pi*reff/0.55, ni, nr, 'sigma_1')*0.55**2
    sigma_a3 = RetrieveFromLUT(2.*np.pi*reff/0.86, ni, nr, 'sigma_1')*0.86**2
    sigma_b1 = RetrieveFromLUT(2.*np.pi*reff/0.44, ni, nr, 'sigma_2')*0.44**2
    sigma_b2 = RetrieveFromLUT(2.*np.pi*reff/0.55, ni, nr, 'sigma_2')*0.55**2
    sigma_b3 = RetrieveFromLUT(2.*np.pi*reff/0.86, ni, nr, 'sigma_2')*0.86**2

    print("Extinctions from LUT1 and LUT2 are {:4.2f} and {:4.2f}".format(sigma_a2,sigma_b2))
    print("AEs from LUT1 and LUT2 are {:4.2f} and {:4.2f}".format(-np.log10(sigma_a1/sigma_a3)/np.log10(0.44/0.86),-np.log10(sigma_b1/sigma_b3)/np.log10(0.44/0.86)))
    return
    
    
def RetrieveSigma(reff, wavel, ni, nr, LUTvar):
    
    sigma = RetrieveFromLUT(2.*np.pi*reff/wavel, ni, nr, LUTvar)*wavel**2*1.e-12
    
    return sigma # [m2]
    
def calcAE(sigma_a, wavel_a, sigma_b, wavel_b):
    
    AE = -np.log10(sigma_a/sigma_b)/np.log10(wavel_a/wavel_b)
    
    return AE

def RetrieveMass(reff, dens, LUTvar):

    if LUTvar == 'sigma_1': # take into account units of Reff [mum] and dens [kg/m3] as well as size distribution
        mass  = 4.0/3.0*np.pi*reff**3*np.exp(-3.*np.log(1.59)**2)*1.e-18*dens # [kg]
    if LUTvar == 'sigma_2':
        mass  = 4.0/3.0*np.pi*reff**3*np.exp(-3.*np.log(2.00)**2)*1.e-18*dens # [kg]

    return mass # [kg/m3]

def calcMEC(sigma, mass):
     
    return sigma/mass # [m2/kg]


def calcReff(Reff_1,Reff_2,mode):
    
    mom2 = Reff_1**2*np.exp(-3.*np.log(1.59)**2)*(1.-mode)+Reff_2**2*np.exp(-3.*np.log(2.00)**2)*mode
    mom3 = Reff_1**3*np.exp(-3.*np.log(1.59)**2)*(1.-mode)+Reff_2**3*np.exp(-3.*np.log(2.00)**2)*mode

    return mom3/mom2

def RetrieveFromLUT(dis, ni, nr, LUTvar):
    
    # Simple retrieval scheme based on NN
    if LUTvar == 'sigma_1':
        idis = (m_dis_1-1)*np.log10(dis/dis_1_min)/np.log10(dis_1_max/dis_1_min)
        ini  = (m_ni_1-1)*np.log10(ni/ni_1_min)/np.log10(ni_1_max/ni_1_min)
        inr  = (m_nr_1-1)*(nr-nr_1_min)/(nr_1_max-nr_1_min)
    
        if (idis < 0 or idis > m_dis_1-1):
            print("[ERROR] from RetrieveFromLut: idis {} out of range.".format(idis))
        if (ini < 0 or ini > m_ni_1-1):
            print("[ERROR] from RetrieveFromLut: ini {} out of range.".format(ini))    
        if (inr < 0 or inr > m_nr_1-1):
            print("[ERROR] from RetrieveFromLut: inr {} out of range.".format(inr))
        
    if LUTvar == 'sigma_2':
        idis = np.round( (m_dis_2-1)*np.log10(dis/dis_2_min)/np.log10(dis_2_max/dis_2_min) )
        ini  = np.round( (m_ni_2-1)*np.log10(ni/ni_2_min)/np.log10(ni_2_max/ni_2_min) )
        inr  = np.round( (m_nr_2-1)*(nr-nr_2_min)/(nr_2_max-nr_2_min) )
    
        if (idis < 0 or idis > m_dis_2-1):
            print("[ERROR] from RetrieveFromLut: idis {} out of range.".format(idis))
        if (ini < 0 or ini > m_ni_2-1):
            print("[ERROR] from RetrieveFromLut: ini {} out of range.".format(ini))    
        if (inr < 0 or inr > m_nr_2-1):
            print("[ERROR] from RetrieveFromLut: inr {} out of range.".format(inr))
            
    LUT = dataset_LUT.variables[LUTvar]
    
# Need more sophisticated interpolation routines
    retrieved = LUT[np.round(idis), np.round(ini), np.round(inr)] 
    
    return retrieved




refractive=pd.read_csv('/Users/santiago/Documents/LE_outputs/2008_Complete/aerosol-extinction_RH_V3.txt',delim_whitespace=True,index_col=None)
#refractive.rename(columns={'8.0':'Parameter'},inplace=True)
refractive.rename(columns={'0.0':'RH'},inplace=True)
rh_aux=np.zeros(len(refractive))
aux=refractive['Parameter'][:]
for i in range(len(rh_aux)):
    if aux[i]>=0:
        #rerh_aux[i]=aux[i]
        refractive['Parameter'][i]='c_ext'
    else:
        #rh_aux[i]=aux[i+aux[i]]
        if aux[i]==-1:
            refractive['Parameter'][i]='sizeparam'
        elif aux[i]==-5:
            refractive['Parameter'][i]='real_refractive'
        elif aux[i]==-6:
            refractive['Parameter'][i]='imaginary_refractive'
        elif aux[i]==-4:
            refractive['Parameter'][i]='extf_spec'  
        elif aux[i]==-8:
            refractive['Parameter'][i]='ext_lut'        
        # elif aux[i]==-5:
        #     refractive['Parameter'][i]='real_refr_spec'  
        # elif aux[i]==-6:
        #     refractive['Parameter'][i]='imagin_refr_spec'      

refractive.set_index(['aerosol', 'RH', 'Parameter'],inplace=True)
refractive.columns=[float(f) for f in refractive.columns]

# ==Graph specie optical parameter===
species={}
species[0]=['dust_ff ','dust_f  ','dust_ccc','dust_cc ','dust_c  ']
species[1]=['na_ff ','na_f  ','na_ccc','na_cc ','na_c  ']
species[2]=['no3a_c  ','no3a_f  ','so4a_f  ','nh4a_f  ']
species[3]=['pom_c  ','pom_f  ','ppm_c  ','ppm_f  ','ec_c  ','ec_f  ']

#species={}
#species[0]=['dust_c']
#species[1]=['so4a_c']

# rhs=[0,20,40,60,80,100] #0,20,40,60,80,100
rhs=[0]
parameters=['sizeparam','c_ext','real_refractive','imaginary_refractive']
# parameters=['c_ext']
# for rh in rhs:
#     for parameter in parameters:
#         for group in range(len(species)):
            
#             fig, ax = plt.subplots()
#             for specie in species[group]:
#                 refractive.loc[(specie,rh,parameter),:].plot(ax=ax,linewidth=2,markersize=5,marker='s')
#             ax.set_xscale('log')
#             ax.set_yscale('log')
#             ax.legend(species[group],loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
#             ax.set_xlabel('Short wave bands $\mu$m',fontsize=14)
#             if parameter=='sizeparam':
#                 ax.set_ylabel('Size parameter',fontsize=14)
#             elif parameter=='c_ext':
#                 ax.set_ylabel('Extinction Coefficient',fontsize=14)
#             elif parameter=='real_refractive':
#                 ax.set_ylabel('Real refractive Index',fontsize=14)
#             elif parameter=='imaginary_refractive':
#                 ax.set_ylabel('Imaginary refractive Index',fontsize=14)
            
#             plt.title(parameter+' '+specie+'_rh '+str(rh))
#             plt.savefig('./Figures/'+parameter+'_'+str(group)+'_rh_'+str(rh)+'.png',format='png', dpi=200,bbox_inches = "tight")
#             plt.show()
        


#==For ff and f, width is 1.59 (LUT1), for ccc, cc, and c, width is 2.0 (LUT2)
#aer_ff , 0.1623, 1.59
#aer_f  , 0.3444, 1.59
#aer_ccc, 1.4816, 2.00
#aer_cc , 2.4908, 2.00
#aer_c  , 3.9415, 2.00


species={}
species[0]=['dust_ff ','dust_f  ','dust_ccc','dust_cc ','dust_c  ']
species[1]=['na_ff ','na_f  ','na_ccc','na_cc ', 'na_c  ']
species[2]=['so4a_f  ','so4a_c  ']
species[3]=['no3a_f  ','no3a_c  ']
species[4]=['nh4a_f  ','nh4a_c  ']
species[5]=['ppm_f  ','ppm_c  ']
species[6]=['ec_f  ','ec_c  ']

ini_wave=0
end_wave=10
ispecies=0


wavelenghts=refractive.loc[(species[0][0].strip(),0.0,'ext_lut'),:].index[ini_wave:end_wave]
ext_LUT=pd.DataFrame()


size_parameter=np.zeros((len(wavelenghts)))

ratio=pd.DataFrame()  
radii=np.zeros((len(species[ispecies])))
j=0            
for specie in species[ispecies]:
    if 'f' in specie[-3:]:
        sigma='sigma_1'
    else:
        sigma='sigma_2'
    i=0
    ext_LUT[specie.strip()]=np.zeros((len(wavelenghts)))
    ext_LE=np.zeros((len(wavelenghts)))
    # if j==0:
    #     
    for wavelength in wavelenghts:
      size_parameter[i]=refractive.loc[(specie.strip(),0.0,'sizeparam'),wavelength]
      if i==0:
          radii[j]=size_parameter[i]*wavelength/(2*np.pi)
      ni=refractive.loc[(specie.strip(),0.0,'imaginary_refractive'),wavelength]
      nr=refractive.loc[(specie.strip(),0.0,'real_refractive'),wavelength]  
      aux=RetrieveFromLUT(size_parameter[i], ni, nr, sigma)*wavelength**2
      ext_LUT[specie.strip()][i] = aux
      ext_LE[i]=refractive.loc[(specie.strip(),0.0,'ext_lut'),wavelength]*wavelength**2
      i=i+1
    j=j+1
    
    
    ratio[specie.strip()]=ext_LE/ext_LUT[specie.strip()].values
    ratio[specie.strip()+'_higher']=''
    ratio[specie.strip()+'_higher'][ratio[specie.strip()]>=1]='LE'
    ratio[specie.strip()+'_higher'][ratio[specie.strip()]<1]='LUT'
    ratio[specie.strip()][ratio[specie.strip()]<1]=1/ratio[specie.strip()]

ext_LUT.index=wavelenghts
columns=[f.strip() for f in species[ispecies]]
columns_higher=[f+'_higher' for f in columns]
radii=[round(f,4) for f in radii]
fig2,ax2= plt.subplots()

bp_dict = ratio.boxplot(return_type='both',
    patch_artist = True,
)

colors =['C1' if (f=='LUT') else 'C0' for f in ratio[columns_higher].mode().values[0]]
for i,box in enumerate(bp_dict[1]['boxes']):
    box.set_facecolor(colors[i])
    box.set_color(colors[i])
for i,whis in enumerate(bp_dict[1]['whiskers']):
    whis.set_color(colors[round(np.floor(i/2))])

for i,median in enumerate(bp_dict[1]['medians']):
    median.set_color(colors[i])    
for i,whis in enumerate(bp_dict[1]['caps']):
    whis.set_color(colors[round(np.floor(i/2))])

#ax2.set_yscale('log')
ax2.set_xticklabels(radii)
ax2.set_xlabel('Mode radii')
ax2.set_ylabel('Extinction ratio')
LE_patch = mpatches.Patch(color='C0', label='LE>LUT')
LUT_patch = mpatches.Patch(color='C1', label='LUT>LE')
#plt.legend(handles=[LE_patch, LUT_patch])

specie_name=species[ispecies][0][0:-4]
fig2.savefig('./Figures/Ratio_LE_LUT_'+specie_name+'.png',dpi=300,bbox_inches='tight')
plt.show()


#===Graph or Extinction efficiency====
sigmas=[1.59,1.59,2,2,2]
c_ext_701_LE=refractive.loc[(columns,0.0,'ext_lut'),0.7016].values*0.7016**2
c_ext_701_LUT=ext_LUT.loc[0.7016,:].values
size_701_LE=refractive.loc[(columns,0.0,'sizeparam'),0.7016].values

Qext_701_LE=[ c/(np.pi*r*r*np.exp(2*np.log(sigma)**2 )) for c,r,sigma in zip(c_ext_701_LE,radii,sigmas)]
Qext_701_LUT=[ c/(np.pi*r*r*np.exp(2*np.log(sigma)**2 )) for c,r,sigma in zip(c_ext_701_LUT,radii,sigmas)]

fig3,ax3= plt.subplots()
ax3.plot(size_701_LE,Qext_701_LE,'s-',color='C0',label='LE',linewidth=3)
ax3.plot(size_701_LE,Qext_701_LUT,'s-',color='C1',label='ECHAM-HAM',linewidth=3)
ax3.plot(np.linspace(0,40,100),2*np.ones(100),'--',color='k',label=None,linewidth=1)
ax3.set_xlim([0,37])
ax3.set_xlabel('Size parameter')
ax3.legend()
ax3.set_title('Extiction efficiency Dust')
fig3.savefig('./Figures/Qext_Na_701nm.png',dpi=300,bbox_inches='tight')
