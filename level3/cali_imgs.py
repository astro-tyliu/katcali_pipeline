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
import time
import pickle
import sys
import os
import gc

from matplotlib.offsetbox import AnchoredText
from matplotlib.colors import LogNorm
from astropy import units as u
from astropy.coordinates import SkyCoord

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

Tcmb = 2.725

# Input arguments
cali_output = sys.argv[1]
fname = sys.argv[2]
ant = sys.argv[3]
pol = sys.argv[4]

# Output paths
cali_output_file = f'/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/{cali_output}/'
output_file1 = os.path.join(cali_output_file, 'waterfall')
output_file2 = os.path.join(cali_output_file, 'mean_std')
os.makedirs(output_file1, exist_ok=True)
os.makedirs(output_file2, exist_ok=True)

# Constants
# num_ants = 64
# pols = ['h', 'v']
ch_plots = [300, 2800, 3200, 3500]
ch_ref = 3200

# Logging
print('Start @', time.asctime())
print('katcali version:', katcali.__version__)

# Set plotting parameters
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.linewidth'] = 1.5

# ants = [f'm{i:03d}' for i in range(num_ants)]
# ants = ['m005', 'm006', 'm007']

# Loop over receivers
data = kio.load_data(fname)
target_list, c0_list, bad_ants, flux_model_list = kio.check_ants(fname)

# import psutil, os
# print(f"Memory before loop: {psutil.Process(os.getpid()).memory_info().rss/1024/1024:.2f} MB")

# for ant in ants:
#     for pol in pols:
recv = ant + pol
data_path = os.path.join(cali_output_file, f'{fname}_{recv}', 'level3_data')

if not os.path.exists(data_path):
    print(f'{recv}: level3_data not found')

try:
    with open(data_path, 'rb') as f:
        dd = pickle.load(f)
        check_map1 = dd['T_map'] - dd['Tsm_map'] - dd['Tel_map']
        check_map2 = dd['Tresi_map']
        del dd
        gc.collect()
except Exception as e:
    print(f'{recv}: {e}')

try:
    data.select(ants=ant, pol=pol)
    vis, flags = kio.call_vis(fname, recv)
    vis_backup = vis.copy()
    ra, dec, az, el = kio.load_coordinates(data)
    timestamps, freqs = kio.load_tf(data)
    dump_period = data.dump_period

    if isinstance(target_list, list):
        ang_deg = kio.load_ang_deg2(ra, dec, c0_list)
    else:
        ang_deg = kio.load_ang_deg(ra, dec, c0_list)
    ang_deg = np.array(ang_deg)

    dp_tt, dp_ss, *_ = kl.cal_dp_label(data, flags, ant, pol, ch_ref, ang_deg)
    dp_sb, dp_se = dp_ss[0], dp_ss[-1]

    nd_on_time, nd_cycle, nd_set = kd.cal_nd_basic_para(fname)
    nd_on_edge, nd_off_edge = kd.cal_nd_edges(timestamps, nd_set, nd_cycle, nd_on_time)
    nd_ratio, nd_0, nd_1x = kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
    nd_t0, nd_t1x, nd_s0, *_ = kl.cal_label_intersec(dp_tt, dp_ss, nd_0, nd_1x)

    # Plot sky temperature
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    axs[0].imshow(check_map1[nd_s0, 272:2869], aspect='auto')
    axs[1].imshow(check_map1[nd_s0, 3133:3547], aspect='auto')
    for ax in axs:
        ax.set_xlabel("Channels")
        ax.set_ylabel("Timestamps")
        ax.figure.colorbar(ax.images[0], ax=ax)
    axs[0].set_xticks([0, 2868-272])
    axs[0].set_xticklabels([272, 2868])
    axs[1].set_xticks([0, 3546-3133])
    axs[1].set_xticklabels([3133, 3546])
    fig.suptitle(f"Sky temperature of {fname} {recv}")
    fig.tight_layout()
    fig.savefig(os.path.join(output_file1, f'sky_temp_{recv}.png'), bbox_inches='tight')
    plt.close(fig)

    # Plot residuals
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    axs[0].imshow(check_map2[nd_s0, 272:2869], aspect='auto')
    axs[1].imshow(check_map2[nd_s0, 3133:3547], aspect='auto')
    for ax in axs:
        ax.set_xlabel("Channels")
        ax.set_ylabel("Timestamps")
        ax.figure.colorbar(ax.images[0], ax=ax)
    axs[0].set_xticks([0, 2868-272])
    axs[0].set_xticklabels([272, 2868])
    axs[1].set_xticks([0, 3546-3133])
    axs[1].set_xticklabels([3133, 3546])
    fig.suptitle(f"Residuals of {fname} {recv}")
    fig.tight_layout()
    fig.savefig(os.path.join(output_file1, f'residuals_{recv}.png'), bbox_inches='tight')
    plt.close(fig)

except Exception as err:
    print(f'{recv}: Failed during processing - {err}')
finally:
    del vis, flags, vis_backup, ra, dec, az, el, timestamps, freqs, dump_period
    del ang_deg, dp_tt, dp_ss, dp_sb, dp_se, nd_on_time, nd_cycle, nd_set
    del nd_on_edge, nd_off_edge, nd_ratio, nd_0, nd_1x, nd_t0, nd_t1x, nd_s0
    del check_map1, check_map2
    gc.collect()

# print(f"Memory after loop: {psutil.Process(os.getpid()).memory_info().rss/1024/1024:.2f} MB")

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

        # del map_plots, map_mean, resi_plots, resi_mean

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
