import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import os
import re
import matplotlib.backends.backend_pdf


points = pd.DataFrame() #empty data frame to store values
points_cross_section = pd.DataFrame()
points_cross_section_repulsive = pd.DataFrame()

for root, dirs, files in os.walk("phase shift data/attractive/"):
    for file in files:
        print("Loading "+str(file))
        temp = re.findall(r'\d+', str(file))
        kappa=float(temp[0])
        beta = [float(temp[1])+float(temp[2])*10**(-1*len(temp[2])) if len(temp)==3 else float(temp[1])][0]
        #print('kappa='+str(kappa)+' beta='+str(beta))
        delta_llist1 = pd.read_csv("phase shift data/attractive/"+str(file))
        delta_llist = delta_llist1.loc[(delta_llist1['Converged']==True)]
        #print(delta_llist)
        deltal_differ_list = [lval*(lval+1)/(2*lval+1)*np.sin(delta_llist['delta_L'].iloc[lval+2]- delta_llist['delta_L'].iloc[lval])**2 for lval in range(len(delta_llist.index)-2)]
        deltal_differ_list_ident = [lval*(lval+1)/(2*lval+1)*np.sin(delta_llist['delta_L'].iloc[lval+2]- delta_llist['delta_L'].iloc[lval])**2 if lval%2!=0 else 0 for lval in range(len(delta_llist.index)-2)]
        sigmav_nonident_sum = np.sum(deltal_differ_list)
        sigmav_ident_sum = np.sum(deltal_differ_list_ident)
        difference = (sigmav_nonident_sum - sigmav_ident_sum)/sigmav_nonident_sum
        
        #print(sigmav_nonident_sum)
        #print(sigmav_ident_sum)

        points_cs = pd.DataFrame({'beta':np.array([beta]), 'kappa':np.array([kappa]), 'sigmav_nonident':np.array([sigmav_nonident_sum]), 'sigmav_ident':np.array([sigmav_ident_sum])})
        points_cs['difference'] = np.array([difference])
        #print(points_cs)
        points_cross_section=pd.concat([points_cross_section, points_cs], ignore_index=True)

for root, dirs, files in os.walk("phase shift data/repulsive/"):
    for file in files:
        print("Loading "+str(file))
        temp = re.findall(r'\d+', str(file))
        kappa=float(temp[0])
        beta = [float(temp[1])+float(temp[2])*10**(-1*len(temp[2])) if len(temp)==3 else float(temp[1])][0]
        #print('kappa='+str(kappa)+' beta='+str(beta))
        delta_llist1 = pd.read_csv("phase shift data/attractive/"+str(file))
        delta_llist = delta_llist1.loc[delta_llist1['Converged']==True]
        deltal_differ_list = [lval*(lval+1)/(2*lval+1)*np.sin(delta_llist['delta_L'].iloc[lval+2]- delta_llist['delta_L'].iloc[lval])**2 for lval in range(len(delta_llist.index)-2)]
        deltal_differ_list_ident = [lval*(lval+1)/(2*lval+1)*np.sin(delta_llist['delta_L'].iloc[lval+2]- delta_llist['delta_L'].iloc[lval])**2 if lval%2!=0 else 0 for lval in range(len(delta_llist.index)-2)]
        sigmav_nonident_sum = np.sum(deltal_differ_list)
        sigmav_ident_sum = np.sum(deltal_differ_list_ident)
        difference = (sigmav_nonident_sum - sigmav_ident_sum)/sigmav_nonident_sum
        
        #print(sigmav_nonident_sum)
        #print(sigmav_ident_sum)

        points_cs_rep = pd.DataFrame({'beta':np.array([beta]), 'kappa':np.array([kappa]), 'sigmav_nonident':np.array([sigmav_nonident_sum]), 'sigmav_ident':np.array([sigmav_ident_sum])})
        points_cs_rep['difference'] = np.array([difference])
        points_cross_section_repulsive=pd.concat([points_cross_section_repulsive, points_cs_rep], ignore_index=True)

pdf = matplotlib.backends.backend_pdf.PdfPages("sigmav_oddL_convergingdeltaL.pdf")
fig , axs = plt.subplots(2,2)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==1)],points_cross_section['sigmav_nonident'].loc[(points_cross_section['kappa']==1)], label='Non-Identical', s=7)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==1)],points_cross_section['sigmav_ident'].loc[(points_cross_section['kappa']==1)], label='Identical', s=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
#axs[0,0].set_xlabel(r'$\beta$')
axs[0,0].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,0].set_title('Attractive, $\kappa =1$')
#axs[0,0].legend()




axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==1)],points_cross_section_repulsive['sigmav_nonident'].loc[(points_cross_section_repulsive['kappa']==1)], label='Non-Identical', s=7)
axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==1)],points_cross_section_repulsive['sigmav_ident'].loc[(points_cross_section_repulsive['kappa']==1)], label='Identical', s=7)
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
#axs[0,1].set_xlabel(r'$\beta$')
#axs[0,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,1].set_title('Repuslive, $\kappa =1$', fontsize=10)
#axs[0,1].legend()

axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==2)],points_cross_section['sigmav_nonident'].loc[(points_cross_section['kappa']==2)], label='Non-Identical', s=7)
axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==2)],points_cross_section['sigmav_ident'].loc[(points_cross_section['kappa']==2)], label='Identical', s=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r'$\beta$')
axs[1,0].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,0].set_title('$\kappa=2$', fontsize=10)
#axs[1,0].legend()

axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==2)],points_cross_section_repulsive['sigmav_nonident'].loc[(points_cross_section_repulsive['kappa']==2)], label='Non-Identical', s=7)
axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==2)],points_cross_section_repulsive['sigmav_ident'].loc[(points_cross_section_repulsive['kappa']==2)], label='Identical', s=7)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlabel(r'$\beta$')
#axs[1,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,1].set_title('$\kappa=2$', fontsize=10)
axs[1,1].legend()
plt.tight_layout(pad=1.0)
pdf.savefig()

fig , axs = plt.subplots(2,2)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==5)],points_cross_section['sigmav_nonident'].loc[(points_cross_section['kappa']==5)], label='Non-Identical', s=7)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==5)],points_cross_section['sigmav_ident'].loc[(points_cross_section['kappa']==5)], label='Identical', s=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
#axs[0,0].set_xlabel(r'$\beta$')
axs[0,0].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,0].set_title('Attractive, $\kappa =5$', fontsize=10)
#axs[0,0].legend()


axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==5)],points_cross_section_repulsive['sigmav_nonident'].loc[(points_cross_section_repulsive['kappa']==5)], label='Non-Identical', s=7)
axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==5)],points_cross_section_repulsive['sigmav_ident'].loc[(points_cross_section_repulsive['kappa']==5)], label='Identical', s=7)
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
#axs[0,1].set_xlabel(r'$\beta$')
#axs[0,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,1].set_title('Repuslive, $\kappa =5$', fontsize=10)
#axs[0,1].legend()

axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==20)],points_cross_section['sigmav_nonident'].loc[(points_cross_section['kappa']==20)], label='Non-Identical', s=7)
axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==20)],points_cross_section['sigmav_ident'].loc[(points_cross_section['kappa']==20)], label='Identical', s=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r'$\beta$')
axs[1,0].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,0].set_title('$\kappa=20$', fontsize=10)
#axs[1,0].legend()

axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==20)],points_cross_section_repulsive['sigmav_nonident'].loc[(points_cross_section_repulsive['kappa']==20)], label='Non-Identical', s=7)
axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==20)],points_cross_section_repulsive['sigmav_ident'].loc[(points_cross_section_repulsive['kappa']==20)], label='Identical', s=7)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlabel(r'$\beta$')
#axs[1,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,1].set_title('$\kappa=20$', fontsize=10)
axs[1,1].legend()
plt.tight_layout(pad=1.0)
pdf.savefig()


fig , axs = plt.subplots(2,2)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==1)],points_cross_section['difference'].loc[(points_cross_section['kappa']==1)], label='Non-Identical', s=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
#axs[0,0].set_xlabel(r'$\beta$')
axs[0,0].set_ylabel(r'$\frac{\sigma_v^\mathrm{non-ident} - \sigma_v^\mathrm{ident}}{\sigma_v^\mathrm{non-ident}}$')
axs[0,0].set_title('Attractive, $\kappa =1$', fontsize=10)
#axs[0,0].legend()




axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==1)],points_cross_section_repulsive['difference'].loc[(points_cross_section_repulsive['kappa']==1)], label='Non-Identical', s=7)

axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
#axs[0,1].set_xlabel(r'$\beta$')
#axs[0,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,1].set_title('Repuslive, $\kappa =1$', fontsize=10)
#axs[0,1].legend()

axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==2)],points_cross_section['difference'].loc[(points_cross_section['kappa']==2)], label='Non-Identical', s=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r'$\beta$')
axs[1,0].set_ylabel(r'$\frac{\sigma_v^\mathrm{non-ident} - \sigma_v^\mathrm{ident}}{\sigma_v^\mathrm{non-ident}}$')
axs[1,0].set_title('$\kappa=2$', fontsize=10)
#axs[1,0].legend()

axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==2)],points_cross_section_repulsive['difference'].loc[(points_cross_section_repulsive['kappa']==2)], label='Non-Identical', s=7)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlabel(r'$\beta$')
#axs[1,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,1].set_title('$\kappa=2$', fontsize=10)
#axs[1,1].legend()
plt.tight_layout(pad=1.0)
pdf.savefig()

fig , axs = plt.subplots(2,2)
axs[0,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==5)],points_cross_section['difference'].loc[(points_cross_section['kappa']==5)], label='Non-Identical', s=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
#axs[0,0].set_xlabel(r'$\beta$')
axs[0,0].set_ylabel(r'$\frac{\sigma_v^\mathrm{non-ident} - \sigma_v^\mathrm{ident}}{\sigma_v^\mathrm{non-ident}}$')
axs[0,0].set_title('Attractive, $\kappa =5$', fontsize=10)
#axs[0,0].legend()




axs[0,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==5)],points_cross_section_repulsive['difference'].loc[(points_cross_section_repulsive['kappa']==5)], label='Non-Identical', s=7)

axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
#axs[0,1].set_xlabel(r'$\beta$')
#axs[0,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[0,1].set_title('Repuslive, $\kappa =5$', fontsize=10)
#axs[0,1].legend()

axs[1,0].scatter(points_cross_section['beta'].loc[(points_cross_section['kappa']==20)],points_cross_section['difference'].loc[(points_cross_section['kappa']==20)], label='Non-Identical', s=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r'$\beta$')
axs[1,0].set_ylabel(r'$\frac{\sigma_v^\mathrm{non-ident} - \sigma_v^\mathrm{ident}}{\sigma_v^\mathrm{non-ident}}$')
axs[1,0].set_title('$\kappa=20$', fontsize=10)
#axs[1,0].legend()

axs[1,1].scatter(points_cross_section_repulsive['beta'].loc[(points_cross_section_repulsive['kappa']==20)],points_cross_section_repulsive['difference'].loc[(points_cross_section_repulsive['kappa']==20)], label='Non-Identical', s=7)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlabel(r'$\beta$')
#axs[1,1].set_ylabel(r'$\sigma_v \kappa^2/4\pi$')
axs[1,1].set_title('$\kappa=20$', fontsize=10)
#axs[1,1].legend()

plt.tight_layout(pad=1.0)

pdf.savefig()

pdf.close()

        

