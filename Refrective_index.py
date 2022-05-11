import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

refractive=pd.read_csv('/Users/santiago/Documents/LE_outputs/2008_Complete/aerosol-extinction_RH_V2.txt',delim_whitespace=True,index_col=None)
#refractive.rename(columns={'8.0':'Parameter'},inplace=True)
refractive.rename(columns={'0.0':'RH'},inplace=True)
rh_aux=np.zeros(len(refractive))
aux=refractive['Parameter'][:]
for i in range(len(rh_aux)):
    if aux[i]>=0:
        rh_aux[i]=aux[i]
        refractive['Parameter'][i]='c_ext'
    else:
        rh_aux[i]=aux[i+aux[i]]
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
species[0]=['dust_ff','dust_f','dust_ccc','dust_cc','dust_c']
species[1]=['na_ff','na_f','na_ccc','na_cc','na_c']
species[2]=['no3a_c','no3a_f','so4a_f','nh4a_f']
species[3]=['pom_c','pom_f','ppm_c','ppm_f','ec_c','ec_f']

#species={}
#species[0]=['dust_c']
#species[1]=['so4a_c']

# rhs=[0,20,40,60,80,100] #0,20,40,60,80,100
rhs=[0]
parameters=['sizeparam','c_ext','real_refractive','imaginary_refractive']
# parameters=['c_ext']
for rh in rhs:
    for parameter in parameters:
        for group in range(len(species)):
            
            fig, ax = plt.subplots()
            for specie in species[group]:
                refractive.loc[(specie,rh,parameter),:].plot(ax=ax,linewidth=2,markersize=5,marker='s')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(species[group],loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
            ax.set_xlabel('Short wave bands $\mu$m',fontsize=14)
            if parameter=='sizeparam':
                ax.set_ylabel('Size parameter',fontsize=14)
            elif parameter=='c_ext':
                ax.set_ylabel('Extinction Coefficient',fontsize=14)
            elif parameter=='real_refractive':
                ax.set_ylabel('Real refractive Index',fontsize=14)
            elif parameter=='imaginary_refractive':
                ax.set_ylabel('Imaginary refractive Index',fontsize=14)
            
            plt.title(parameter+' '+specie+'_rh '+str(rh))
            plt.savefig('./Figures/'+parameter+'_'+str(group)+'_rh_'+str(rh)+'.png',format='png', dpi=200,bbox_inches = "tight")
            plt.show()
        
