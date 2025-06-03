################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
###############################################
#imports
import katdal
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import astropy.coordinates as ac
import functools
import healpy as hp
import optparse
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import healpy as hp
from astropy import units as u
from matplotlib.offsetbox import AnchoredText
from matplotlib.colors import LogNorm
import time
import pickle
import sys
import os
Tcmb=2.725
import katcali
import katcali.visualizer as kv
import katcali.models as km
import katcali.rfi as kr
import katcali.solver as ks
import katcali.io as kio
import katcali.label_dump as kl
import katcali.diode as kd
import katcali.filter as kf
import katcali.beam_UHF as kb_u
from astropy.coordinates import SkyCoord
from astropy import units as u


fname = sys.argv[1]
output_file = sys.argv[2]

pols = ['h', 'v']
cali_output_file=f'/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/{output_file}/'
output_file = cali_output_file + 'cali_imgs/'
os.makedirs(output_file, exist_ok=True)

# UHF ch_plot: (272,2869)+(3133,3547)
data = kio.load_data(fname)
target_list, c0_list, bad_ants, flux_model_list = kio.check_ants(fname)
for i in range(64):
    ant = f'm{i:03d}'
    # if ant not in bad_ants:
    # if fname in ['1684087370', '1675643846'] and ant in ['m024', 'm037']:
    #     continue
    try:
        for pol in pols:
            recv = ant + pol
            with open(cali_output_file + str(fname) + '_' + str(recv) + f'/level2_data', 'rb') as f:
                d_r = pickle.load(f)
            # d_r = pickle.load(open(cali_output_file + str(fname) + '_' + str(recv) + f'/level2_data', 'rb'))
            print(fname, recv)
            T_map1=d_r['T_map']
            Tresi_map1=d_r['Tresi_map']
            gain_map1=d_r['gain_map']
            Tnd_ref_list1=d_r['Tnd_ref_list']
            Tnda_list1=d_r['Tnda_list']
            Tndb_list1=d_r['Tndb_list']
            Tnd_diff_ratio_list1=d_r['Tnd_diff_ratio_list']
            NRMSE1_list1=d_r['NRMSE1_list']
            NRMSE2_list1=d_r['NRMSE2_list']
            # d_r2=pickle.load(open(cali_output_file+str(fname)+'_'+str(recv)+f'_{ratios}_level2_Tnd_data','rb'))
            # Tnd_ref_list2=d_r2['Tnd_ref_list']
            # Tnda_list2=d_r2['Tnda_list']
            # Tndb_list2=d_r2['Tndb_list']

            plt.figure(figsize=(16,5))
            plt.subplots_adjust(hspace=0)
            gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
            plt.subplot(gs[0,0])
            plt.plot(Tnda_list1,'.',ms=4)
            plt.plot(Tndb_list1,'.',ms=4)
            plt.plot(Tnd_ref_list1, color='grey',zorder=0)
            if isinstance(target_list, list):
                plt.legend([f'Tnda ({target_list[0]})',f'Tndb ({target_list[1]})','Tnd_ref'],ncol=3)
            else:
                plt.legend([f'Tnda ({target_list})',f'Tndb ({target_list})','Tnd_ref'],ncol=3)
            #plt.xlabel('channel')
            plt.xticks([])
            plt.ylabel('Tnd (K)')
            plt.title('Tnd result '+fname+', '+recv)
            plt.subplot(gs[1,0])
            plt.plot(Tnd_diff_ratio_list1,'m.',ms=4)
            plt.xlabel('channel')
            plt.ylabel('Tnd_diff_ratio')
            plt.savefig(output_file+'Tnd_all_'+str(recv)+'.png', bbox_inches='tight')
            # plt.show()
            plt.close()

            del d_r, T_map1, Tresi_map1, gain_map1, Tnd_ref_list1, Tnda_list1, Tndb_list1
            del Tnd_diff_ratio_list1, NRMSE1_list1, NRMSE2_list1
    except Exception as e:
        print(f'{ant}: {e}')
        continue
