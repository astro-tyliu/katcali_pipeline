################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
###############################################
#imports
import katdal
import numpy as np
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
Tcmb=2.725
try:
    import katcali
except:
    import sys
    sys.path.insert(0, '/home/liutianyang/katcali')
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


fnames = ['1684087370', '1675643846']
pols = ['h', 'v']

py_output_file='/mnt/results/katcali_results/MeerKLASS_UHF/level2/py_results/cali_results/'
output_file='/mnt/results/katcali_results/MeerKLASS_UHF/level2/py_results/cali_imgs/'

ratio_t = 0.4
ratio_ch = 0.5
# ratio_ch = 0.171
ratios = "{:03d}".format(int(ratio_t * 100)) + "{:03d}".format(int(ratio_ch * 100))

ratios = "{:03d}".format(int(ratio_t * 100)) + "{:03d}".format(int(ratio_ch * 100))

# UHF ch_plot: (272,2869)+(3133,3547)
for fname in fnames:
    data = kio.load_data(fname)
    target, c0, bad_ants, flux_model = kio.check_ants(fname)
    for i in range(64):
        ant = f'm{i:03d}'
        if ant not in bad_ants:
            if fname in ['1684087370', '1675643846'] and ant in ['m024', 'm037']:
                continue
            for pol in pols:
                recv = ant + pol
                print(fname, recv)
                d_r = pickle.load(open(py_output_file + str(fname) + '_' + str(recv) + f'_{ratios}_level2_data', 'rb'))
                T_map1=d_r['T_map']
                Tresi_map1=d_r['Tresi_map']
                gain_map1=d_r['gain_map']
                Tnd_ref_list1=d_r['Tnd_ref_list']
                Tnda_list1=d_r['Tnda_list']
                Tndb_list1=d_r['Tndb_list']
                Tnd_diff_ratio_list1=d_r['Tnd_diff_ratio_list']
                NRMSE1_list1=d_r['NRMSE1_list']
                NRMSE2_list1=d_r['NRMSE2_list']
                # d_r2=pickle.load(open(py_output_file+str(fname)+'_'+str(recv)+f'_{ratios}_level2_Tnd_data','rb'))
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
                plt.legend(['Tnda','Tndb','Tnd_ref'],ncol=3)
                #plt.xlabel('channel')
                plt.xticks([])
                plt.ylabel('Tnd (K)')
                plt.title('Tnd result '+fname+', '+recv)
                plt.subplot(gs[1,0])
                plt.plot(Tnd_diff_ratio_list1,'m.',ms=4)
                plt.xlabel('channel')
                plt.ylabel('Tnd_diff_ratio')
                plt.savefig(output_file+'Tnd_all_'+str(fname)+'_'+str(recv)+'_'+ratios+'.png', bbox_inches='tight')
                # plt.show()
                plt.close()
