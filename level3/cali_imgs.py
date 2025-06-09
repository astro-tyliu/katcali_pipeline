################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
###############################################
#imports
import katdal
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
import gc
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

cali_output = sys.argv[1]
fname = sys.argv[2]

cali_output_file=f'/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/{cali_output}/'
output_file1 = cali_output_file + 'waterfall/'
output_file2 = cali_output_file + 'mean_std/'
os.makedirs(output_file1, exist_ok=True)
os.makedirs(output_file2, exist_ok=True)
num_ants = 10

pols = ['h', 'v']
# ch_plots = [800, 1300, 1800, 2300, 3300]
ch_plots = [300, 2800, 3200, 3500]
ch_ref = 3200

print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')
print (katcali.__version__)
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'], plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'], plt.rcParams['lines.linewidth'])

ants = [f'm{i:03d}' for i in range(num_ants)]
map_mean_std = np.ma.array(np.zeros((3, len(ch_plots), num_ants, len(pols))), mask=True)
resi_mean_std = np.ma.array(np.zeros((3, len(ch_plots), num_ants, len(pols))), mask=True)

for i, ant in enumerate(ants):
    
    for j, pol in enumerate(pols):
        
        recv = ant+pol
        try:
            with open(cali_output_file + str(fname) + '_' + str(recv) + '/level3_data', 'rb') as f:
                dd = pickle.load(f)
        except Exception as e:
            print(f'{ant}: {e}')
            continue
        
        data = kio.load_data(fname)
        #show the calibrator and bad ants information
        # target, c0, bad_ants, flux_model = kio.check_ants(fname)
        target_list,c0_list,bad_ants,flux_model_list=kio.check_ants(fname)
        data.select(ants=ant,pol=pol)

        vis,flags = kio.call_vis(fname,recv)
        vis_backup = vis.copy()
        ra, dec, az, el = kio.load_coordinates(data)
        timestamps, freqs = kio.load_tf(data)
        dump_period = data.dump_period
        # ang_deg = kio.load_ang_deg(ra, dec, c0)
        if isinstance(target_list, list):
            ang_deg=kio.load_ang_deg2(ra,dec,c0_list)
        else:
            ang_deg=kio.load_ang_deg(ra,dec,c0_list) #modeified for xcalib
        ang_deg=np.array(ang_deg)
        dp_tt, dp_ss, dp_f, dp_w, dp_t, dp_s, dp_slew, dp_stop = kl.cal_dp_label(data, flags, ant, pol, ch_ref, ang_deg)
        del data
        dp_sb, dp_se = dp_ss[0], dp_ss[-1]

        nd_on_time, nd_cycle, nd_set = kd.cal_nd_basic_para(fname)
        print (nd_on_time, nd_cycle, nd_set)
        nd_on_edge, nd_off_edge = kd.cal_nd_edges(timestamps, nd_set, nd_cycle, nd_on_time)
        print (len(nd_on_edge), len(nd_off_edge))
        nd_ratio, nd_0, nd_1x = kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)

        nd_t0, nd_t1x, nd_s0, nd_s1x, nd_t0_ca, nd_t0_cb, nd_t1x_ca, nd_t1x_cb = kl.cal_label_intersec(dp_tt, dp_ss, nd_0, nd_1x)
        labels_1x = kl.cal_label_intersec_complex(dp_tt, dp_ss, nd_0, nd_1x, nd_ratio)

        dd = pickle.load(open(cali_output_file+fname+'_'+str(recv)+'/level3_data','rb'))
        print (dd.keys())

        check_map1 = dd['T_map']-dd['Tsm_map']-dd['Tel_map']
        plt.figure(figsize=(12, 4))
        plt.subplot(121)
        plt.imshow(check_map1[nd_s0, 272:2869], aspect='auto')
        plt.colorbar()
        plt.xticks([0, 2868-272], [272, 2868])
        plt.xlabel("Channels")
        plt.ylabel("Timestamps")
        plt.subplot(122)
        plt.imshow(check_map1[nd_s0, 3133:3547], aspect='auto')
        plt.colorbar()
        plt.xticks([0, 3546-3133], [3133, 3546])
        plt.xlabel("Channels")
        plt.ylabel("Timestamps")
        plt.suptitle(f"Sky temperature of {fname} {recv}")
        plt.savefig(output_file1+'sky_temp_'+str(recv)+'.png', bbox_inches='tight')
        plt.close()

        check_map2 = dd['Tresi_map']
        plt.figure(figsize=(12, 4))
        plt.subplot(121)
        plt.imshow(check_map2[nd_s0, 272:2869], aspect='auto')
        plt.colorbar()
        plt.xticks([0, 2868-272], [272, 2868])
        plt.xlabel("Channels")
        plt.ylabel("Timestamps")
        plt.subplot(122)
        plt.imshow(check_map2[nd_s0, 3133:3547], aspect='auto')
        plt.colorbar()
        plt.xticks([0, 3546-3133], [3133, 3546])
        plt.xlabel("Channels")
        plt.ylabel("Timestamps")
        plt.suptitle(f"Residuals of {fname} {recv}")
        plt.savefig(output_file1+'residuals_'+str(recv)+'.png', bbox_inches='tight')
        plt.close()

#         map_plots = check_map1[nd_s0, ch_plots]
#         map_mean = np.mean(map_plots, axis=0)
#         map_mean_std[0, :, i, j] = map_mean
#         map_mean_std[1, :, i, j] = np.sqrt(np.mean((map_plots[map_plots > map_mean] - map_mean) ** 2, axis=0))
#         map_mean_std[2, :, i, j] = np.sqrt(np.mean((map_plots[map_plots < map_mean] - map_mean) ** 2, axis=0))

#         resi_plots = check_map2[nd_s0, ch_plots]
#         resi_mean = np.mean(resi_plots, axis=0)
#         resi_mean_std[0, :, i, j] = resi_mean
#         resi_mean_std[1, :, i, j] = np.sqrt(np.mean((resi_plots[resi_plots > resi_mean] - resi_mean) ** 2, axis=0))
#         resi_mean_std[2, :, i, j] = np.sqrt(np.mean((resi_plots[resi_plots < resi_mean] - resi_mean) ** 2, axis=0))

#         del target_list, c0_list, bad_ants, flux_model_list, vis, flags, vis_backup, ra, dec, az, el, timestamps, freqs, dump_period, ang_deg
#         del dp_tt, dp_ss, dp_f, dp_w, dp_t, dp_s, dp_slew, dp_stop, dp_sb, dp_se, nd_on_time, nd_cycle, nd_set, nd_on_edge, nd_off_edge, nd_ratio, nd_0, nd_1x
#         del nd_t0, nd_t1x, nd_s0, nd_s1x, nd_t0_ca, nd_t0_cb, nd_t1x_ca, nd_t1x_cb, labels_1x
#         del dd, check_map1, check_map2, map_plots, map_mean, resi_plots, resi_mean
#         gc.collect()

# for i, ch_plot in enumerate(ch_plots):
    
#     mean_values_1 = map_mean_std[0, i, :, :]
#     lower_std_1 = map_mean_std[1, i, :, :]
#     upper_std_1 = map_mean_std[2, i, :, :] 
    
#     mean_values_2 = resi_mean_std[0, i, :, :]
#     lower_std_2 = resi_mean_std[1, i, :, :]
#     upper_std_2 = resi_mean_std[2, i, :, :]
    
#     x = np.arange(mean_values_1.shape[0])
    
#     fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
#     axs[0].errorbar(x, mean_values_1[:, 0], yerr=[lower_std_1[:, 0], upper_std_1[:, 0]], 
#                 fmt='o', color='b', ecolor='b', capsize=5, capthick=1, elinewidth=1, label="h polarization")
#     axs[0].errorbar(x, mean_values_1[:, 1], yerr=[lower_std_1[:, 1], upper_std_1[:, 1]], 
#                 fmt='o', color='r', ecolor='r', capsize=5, capthick=1, elinewidth=1, label="v polarization")
#     axs[0].set_xlabel("No. ants")
#     axs[0].set_ylabel("Sky temperature (K)")
#     axs[0].legend()
#     axs[0].grid(True)
#     axs[1].errorbar(x, mean_values_2[:, 0], yerr=[lower_std_2[:, 0], upper_std_2[:, 0]], 
#                 fmt='o', color='b', ecolor='b', capsize=5, capthick=1, elinewidth=1, label="h polarization")
#     axs[1].errorbar(x, mean_values_2[:, 1], yerr=[lower_std_2[:, 1], upper_std_2[:, 1]], 
#                 fmt='o', color='r', ecolor='r', capsize=5, capthick=1, elinewidth=1, label="v polarization")
#     axs[1].set_xlabel("No. ants")
#     axs[1].set_ylabel("Residuals (K)")
#     axs[1].legend()
#     axs[1].grid(True)
#     fig.suptitle(f"Sky Temperature and Residuals of channel {ch_plot}", fontsize=16)
#     plt.tight_layout(rect=[0, 0, 1, 0.95])  # 调整布局，避免总标题和子图重叠
#     plt.savefig(output_file2+'mean_std_'+str(ch_plot)+'.png', bbox_inches='tight')
#     plt.close()
