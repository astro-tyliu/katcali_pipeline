#!/usr/bin/env python3
import faulthandler
faulthandler.enable()
import matplotlib
matplotlib.use('AGG')

print(matplotlib.__file__)

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
import healpy as hp
from astropy import units as u
from matplotlib.offsetbox import AnchoredText
import time
import pickle
import sys
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
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.interpolate import Rbf
from scipy.stats import mode
import os
import json

import sys
# sys.path.insert(0, '/home/liutianyang/katcali')
print(sys.path)
print(katcali.__file__)


# All scan modes (indices)
# dp_slew: slew
# dp_stop: stop
# dp_t: track
# dp_s: scan
# dp_w: slew and stop
# dp_tt: track and no flagging and close to
# dp_ss: scan and no flagging
# dp_f: flagging for scan or (flagging or far from) for track
# dp_sb: the first dp_ss
# dp_se: the last dp_ss
# dp_ptr_list: any points within 0.5 deg of any calibrators in dp_ss

# All labels (indices)
# nd_on_time: duration the diode is on each time (duration)
# nd_cycle: period each time on and off (duration)
# nd_set: startup (time)
# nd_on_edge: diode on begin (times)
# nd_off_edge: diode off end (times)
# nd_ratio: time ratios that diode is on each timestamps (ratios)

# nd_0: diode off
# nd_1x: diode on
# nd_t0: diode off and dp_tt (track and no flagging and close to)
# nd_t1x: diode on and dp_tt (track and no flagging and close to)
# nd_s0: diode off and dp_ss (scan and no flagging)
# nd_s1x: diode on and dp_ss (scan and no flagging)
# nd_t0_ca: nd_t0 before scan
# nd_t0_cb: nd_t0 after scan
# nd_t1x_ca: nd_t1x before scan
# nd_t1x_cb: nd_t1x after scan

print (katcali.__version__)
print (katdal.__version__)
print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')

plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5

p_radec=np.loadtxt('radio_source_fsky.txt')
# ra_a4059,dec_a4059=-0.74042,-34.76056

if len(sys.argv) == 1:
    # default: '1675632179', 'm060'
    sys.argv.extend(['1675632179', 'm060', '8.', '4.', '16.', '8.', '00000000_000000'])

fname=sys.argv[1]
ant=sys.argv[2]
Threshold_factor1, Threshold_factor2 = float(sys.argv[3]), float(sys.argv[4])  # diode off, on
Threshold_factor11,Threshold_factor22 = float(sys.argv[5]), float(sys.argv[6])  # diode off, on
base_dir = sys.argv[7] + '/'

# base_dir="/mnt/results/katcali_results/MeerKLASS_UHF/level1/py_results/"

ch_plot=3300
ch_ref=3300

data=kio.load_data(fname)
#print (data)

target_list,c0_list,bad_ants,flux_model_list=kio.check_ants(fname)

ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print (str(i) + ' is bad')
print (fname)
print (ants_good)

nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
print (nd_on_time,nd_cycle,nd_set)

# Select ant and polarization, then load data in

#select ant, polarization, and one channel to show data calibration
# Thresholds = []
# ratios = "{:03d}".format(int(ratio_t * 100)) + "{:03d}".format(int(ratio_ch * 100))

if ant in ants_good: 
    # output_dir = base_dir + fname + "_" + ant + "/"
    output_dir = base_dir
    # os.makedirs(output_dir, exist_ok=True)
    try:  
        for index, pol in enumerate(['h','v']):
            #load data, labels, and parameters
            data.select(ants=ant,pol=pol)
            recv=ant+pol
            # Thresholds.append(str(int(Threshold_factor1 * 100)) + str(int(Threshold_factor2 * 100)) + \
            #                   str(int(Threshold_factor11 * 100)) + str(int(Threshold_factor22 * 100)))

            ch_plot_list=[300,800,1300,1800,2300,2800,3200,3500]
            
            corr_id=kio.cal_corr_id(data,recv)
            assert(recv==data.corr_products[corr_id][0])
            assert(recv==data.corr_products[corr_id][1])
            print (corr_id,recv)

            ra,dec,az,el=kio.load_coordinates(data)
            timestamps,freqs=kio.load_tf(data)
            dump_period=data.dump_period

            vis,flags= kio.call_vis(fname,recv)
            vis_backup=vis.copy()

            vis=np.ma.array(vis_backup,mask=flags)
            print ('# origin mask ratio of '+recv+': '+str(np.shape(np.where(vis.mask==True))[1]/np.shape(vis)[0]/np.shape(vis)[1]))
            
            if isinstance(target_list, list):
                ang_deg=kio.load_ang_deg2(ra,dec,c0_list)
            else:
                ang_deg=kio.load_ang_deg(ra,dec,c0_list) #modeified for xcalib
            #p_radec=np.loadtxt('radio_source.txt')

            dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)

            ra1=ra.copy()
            if np.max(ra1[dp_ss]) - np.min(ra1[dp_ss]) > 300:
                ra1[ra>180]-=360
            # for i in range(len(ra)):
            #     if ra[i]>180:
            #         ra1[i]=ra[i]-360

            print (np.mean(ra),np.mean(ra1))
            ra=ra1
            print (np.mean(ra),np.mean(ra1))

            nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
            print (len(nd_on_edge),len(nd_off_edge))

            nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)

            nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)
            print (len(nd_t0),len(nd_t0_ca),len(nd_t0_cb))
            print (len(nd_t1x),len(nd_t1x_ca),len(nd_t1x_cb))
            print (len(nd_s0),len(nd_s1x))

            labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)

            dp_sb=dp_ss[0]
            dp_se=dp_ss[-1]

            p = SkyCoord(data.ra*u.deg,  data.dec*u.deg, frame='icrs')
            ang_lim=.5

            dp_ptr_list=kl.cal_ptr_mask(p,p_radec,nd_s0, dp_sb,dp_se,ang_lim)

            plt.figure(figsize=(24,20))
            plt.subplots_adjust(wspace=0.1,hspace=0.3)
            for i in range(len(ch_plot_list)):
                ch_plot1=ch_plot_list[i]
                p_data=np.ma.log10(vis[nd_s0,ch_plot1])
                plt.subplot(5,2,i+1)
                plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max()/1.05, s=14)
                if i>len(ch_plot_list)-3:
                    plt.xlabel('R.A.')
                if i%2==0:
                    plt.ylabel('Dec')
                plt.colorbar()
                plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.savefig(output_dir + fname + '_' + recv + '_scatter_map0.png')
            #plt.show()

            plt.figure(figsize=(24,20))
            plt.subplots_adjust(wspace=0,hspace=0.3)
            for i in range(len(ch_plot_list)):
                try:
                    ch_plot1=ch_plot_list[i]
                    p_data=np.ma.log10(vis[nd_s0,ch_plot1])
                    plt.subplot(5,2,i+1)
                    kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=60)
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    if i>len(ch_plot_list)-3:
                        plt.xlabel('R.A.')
                    if i%2==0:
                        plt.ylabel('Dec')
                    plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                except ValueError:
                    xi = np.linspace(ra[nd_s0].min(), ra[nd_s0].max(), 30)
                    yi = np.linspace(dec[nd_s0].min(), dec[nd_s0].max(), 30)
                    # zi = griddata((ra[nd_s0], dec[nd_s0]), p_data, (xi[None, :], yi[:, None]), method = 'linear')
                    zi = np.full((len(yi), len(xi)), np.nan)
                    CS = plt.imshow(zi[::-1, :], extent=(ra[nd_s0].min(),ra[nd_s0].max(),dec[nd_s0].min(),dec[nd_s0].max()), cmap=plt.get_cmap('jet',255), vmin=0, vmax=1, aspect='auto')
                    plt.colorbar(CS)
                    CS = plt.contour(xi, yi, zi, levels=15, linewidths=0.5, colors='k')
                    plt.xlim(ra[nd_s0].min(), ra[nd_s0].max())
                    plt.ylim(dec[nd_s0].min(), dec[nd_s0].max())
            plt.savefig(output_dir + fname + '_' + recv + '_map0.png')
            #plt.show()

            # RFI flagging
            ## Basic RFI flagging (all channels)
            dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
            nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
            nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)
            labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)

            flag_step=1
            vis_clean=kr.vis_flag_v2(data, flags, ch_ref, ant,pol,vis_backup,timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor1, Threshold_factor2, ratio_clean_key=0, plt_key=-1)

            ###second rfi
            rbf = Rbf(range(4096) ,np.ma.mean(vis[nd_0,:],axis=0), smooth=1) #vis has mask, flags has applied.
            cur=rbf(range(4096))  # cur has no mask

            vis_clean2=vis_clean.copy()
            for ch_i in range(4096):
                vis_clean2[:,ch_i]=vis_clean2[:,ch_i]/cur[ch_i] #rescaled!!!
            print (vis_clean2.mean())

            # print(11111111, vis)
            # print(22222222, np.sum(vis.mask))
            # plt.figure(figsize=(6, 5))
            # plt.subplot(211)
            # plt.plot(np.mean(vis[nd_0, :], axis=0), '.')
            # plt.plot(cur)
            # plt.subplot(212)
            # plt.plot(np.mean(vis_clean[nd_0, :], axis=0), 'g.')
            # plt.plot(cur)
            # plt.show()
            # sys.exit()

            vis_clean2_part1=vis_clean2

            vis_clean2_part1.mask[dp_ptr_list,:]=True

            flag_step=2
            #rfi flagging for raw vis data
            vis_clean2_part1=kr.vis_flag_v2(data, flags, ch_ref, ant,pol,vis_clean2_part1, timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor11, Threshold_factor22, ratio_clean_key=0, plt_key=-1)

            vis_clean2=vis_clean.copy() #reset vis_clean2
            vis_clean2.mask=vis_clean2_part1.mask
            vis_clean2.mask[dp_ptr_list,:]=vis_clean.mask[dp_ptr_list,:]

            vis_clean2=kr.clean_bad_ratio(vis_clean2)

            plt.figure(figsize=(24,20))
            plt.subplots_adjust(wspace=0.1,hspace=0.3)
            for i in range(len(ch_plot_list)):
                ch_plot1=ch_plot_list[i]
                p_data=vis_clean2[nd_s0,ch_plot1]
                plt.subplot(5,2,i+1)
                plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
                if i>len(ch_plot_list)-3:
                    plt.xlabel('R.A.')
                if i%2==0:
                    plt.ylabel('Dec')
                plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                #plt.xlim(-30,1)
                #plt.ylim(-37,-24)
                plt.colorbar()
            plt.savefig(output_dir + fname + '_' + recv + '_scatter_map2.png')
            #plt.show()

            plt.figure(figsize=(24,20))
            plt.subplots_adjust(wspace=0,hspace=0.3)
            for i in range(len(ch_plot_list)):
                try:
                    ch_plot1=ch_plot_list[i]
                    p_data=vis_clean2[nd_s0,ch_plot1]
                    plt.subplot(5,2,i+1)
                    kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=60)
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    if i>len(ch_plot_list)-3:
                        plt.xlabel('R.A.')
                    if i%2==0:
                        plt.ylabel('Dec')
                    plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    #plt.xlim(-30,1)
                    #plt.ylim(-37,-24)
                except ValueError:
                    xi = np.linspace(ra[nd_s0].min(), ra[nd_s0].max(), 30)
                    yi = np.linspace(dec[nd_s0].min(), dec[nd_s0].max(), 30)
                    # zi = griddata((ra[nd_s0], dec[nd_s0]), p_data, (xi[None, :], yi[:, None]), method = 'linear')
                    zi = np.full((len(yi), len(xi)), np.nan)
                    CS = plt.imshow(zi[::-1, :], extent=(ra[nd_s0].min(),ra[nd_s0].max(),dec[nd_s0].min(),dec[nd_s0].max()), cmap=plt.get_cmap('jet',255), vmin=0, vmax=1, aspect='auto')
                    plt.colorbar(CS)
                    CS = plt.contour(xi, yi, zi, levels=15, linewidths=0.5, colors='k')
                    plt.xlim(ra[nd_s0].min(), ra[nd_s0].max())
                    plt.ylim(dec[nd_s0].min(), dec[nd_s0].max())
            plt.savefig(output_dir + fname + '_' + recv + '_map2.png')
            #plt.show()

            vis_clean=vis_clean2.copy() #update vis_clean

            #compare vis before and after rfi flagging

            plt.figure(figsize=(14,5.4))
            plt.subplot(121)
            plt.imshow(np.ma.log10(vis_backup),aspect='auto',vmax=2.7)
            plt.ylabel('time dump')
            plt.xlabel('channel')
            plt.title('raw vis of '+str(fname)+', '+str(recv), y=1.12)
            plt.colorbar(label='Log10')
            plt.twiny()
            plt.imshow(np.ma.log10(vis_backup),aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps),0),vmax=2.7)
            plt.xlabel('Freq (MHz)')
            plt.subplot(122)
            plt.imshow(vis_clean,aspect='auto')
            #plt.ylabel('time dump')
            plt.xlabel('channel')
            plt.title('clean vis of '+str(fname)+', '+str(recv), y=1.12)
            plt.colorbar()
            plt.twiny()
            plt.imshow(vis_clean,aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps),0))
            plt.xlabel('Freq (MHz)')
            #plt.savefig(str(fname)+'raw_vis.pdf', bbox_inches='tight')
            plt.savefig(output_dir + fname + '_' + recv + '_raw_vis.png', bbox_inches='tight')
            #plt.show()

            plot_gsize=60

            #raw map in different resolution
            plt.figure(figsize=(24,9))
            plt.subplots_adjust(wspace=0.,hspace=.2)
            plt.subplot(221)
            p_data=np.ma.log10(vis_backup[nd_s0,ch_plot])
            plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, s=8)
            plt.title('Logged raw signal map of '+str(fname)+', '+str(recv)+ ' @'+ str(round(freqs[ch_plot]/1e6,0))+' MHz')
            #plt.xlabel('R.A. (J2000)')
            plt.ylabel('Dec (J2000)')
            clb = plt.colorbar()
            plt.subplot(222)
            p_data=vis_clean[nd_s0,ch_plot]
            plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, s=8)
            plt.title('RFI flagged raw signal map @'+ str(round(freqs[ch_plot]/1e6,0))+' MHz')
            #plt.xlabel('R.A. (J2000)')
            #plt.ylabel('Dec (J2000)')
            clb = plt.colorbar()
            plt.subplot(223)
            p_data=np.ma.array(np.ma.log10(vis_backup[nd_s0,ch_plot]),mask=False) #for plot_mdata mask
            kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
            #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.title('smoothed map of above')
            plt.xlabel('R.A. (J2000)')
            plt.ylabel('Dec (J2000)')
            plt.subplot(224)
            p_data=vis_clean[nd_s0,ch_plot]
            kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
            #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.title('smoothed map of above')
            plt.xlabel('R.A. (J2000)')
            #plt.ylabel('Dec (J2000)')
            plt.savefig(output_dir + 'F_' + fname + '_' + recv + '_rfi_comp_map.png', bbox_inches='tight')
            #plt.show()

            d={}
            d['mask']=vis_clean.mask
            fs=open(output_dir + fname + '_' + recv + '_mask','wb')
            pickle.dump(d,fs,protocol=2)
            fs.close()
            
            print ('# '+recv+' of '+fname+' finished successfully #')
            
            del vis
            del vis_backup
            del vis_clean
            del vis_clean2
            del vis_clean2_part1
            #end

        ## Check the final mask
        ## (mask_h and mask_v are both required)

        #get vis
        recv1=ant+'h'
        recv2=ant+'v'

        vis_h,flags_h= kio.call_vis(fname,recv1)
        vis_v,flags_v= kio.call_vis(fname,recv2)

        #get mask
        d1 = pickle.load(open(output_dir + fname + '_' + ant + 'h_mask','rb'))
        mask_h=d1['mask']
        d2 = pickle.load(open(output_dir + fname + '_' + ant + 'v_mask','rb'))
        mask_v=d2['mask']
        #vis_clean
        vis_clean_h=np.ma.array(vis_h,mask=mask_h)
        vis_clean_v=np.ma.array(vis_v,mask=mask_v)
        print (np.ma.mean(vis_clean_h),np.ma.mean(vis_clean_v))
        #get intersection mask
        vis_div=vis_clean_h/vis_clean_v #intersection mask
        vis_div=kr.clean_bad_ratio(vis_div) #clean bad ratio part
        #update vis_clean
        vis_clean_hh=np.ma.array(vis_clean_h,mask=vis_div.mask)
        vis_clean_vv=np.ma.array(vis_clean_v,mask=vis_div.mask)
        print (np.ma.mean(vis_clean_hh),np.ma.mean(vis_clean_vv))

        np.shape(np.where(vis_div.mask==True))
        mask_ratio_final=np.shape(np.where(vis_div.mask==True))[1]/np.shape(vis_div)[0]/np.shape(vis_div)[1]
        print ('# final mask ratio of '+ant+': '+ str(mask_ratio_final))

        d={}
        d['mask']=vis_div.mask
        fs=open(output_dir + fname + '_' + ant + '_mask','wb')
        pickle.dump(d,fs,protocol=2)
        fs.close()

        num_mask_t = np.mean(vis_div.mask, axis=0)
        mode_t = mode(num_mask_t)
        num_mask_ch = np.mean(vis_div.mask, axis=1)
        mode_ch = mode(num_mask_ch)
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        axes[0].plot(freqs / 1e6, num_mask_t)
        axes[0].set_xlabel('freqs')
        axes[0].set_ylabel('masked ratio')
        axes[0].set_ylim(top=1.13)
        axes[0].set_title('ratios')
        axes[0].text(580, 1.05, f'mode:{round(mode_t.mode[0], 3)}', fontsize=18, color='red')
        axes[0].text(830, 1.05, f'mode count:{mode_t.count[0]}', fontsize=18, color='red')
        # print(222222222, round(mode_t.mode[0], 3), mode_t.count[0])
        axes[1].plot(timestamps - timestamps[0], num_mask_ch)
        axes[1].set_xlabel('timestamps')
        axes[1].set_ylabel('masked ratio')
        axes[1].set_ylim(top=1.13)
        axes[1].set_title('ratios')
        axes[1].text(800, 1.05, f'mode:{round(mode_ch.mode[0], 3)}', fontsize=18, color='red')
        axes[1].text(4800, 1.05, f'mode count:{mode_ch.count[0]}', fontsize=18, color='red')
        # don't include the flags
        plt.savefig(output_dir + fname + '_' + ant + '_masked_ratio.png', bbox_inches='tight')

        plt.figure(figsize=(14,5.4))
        plt.subplot(121)
        plt.imshow(vis_clean_hh,aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps)*2,0))
        plt.xlabel('Freq (MHz)')
        plt.ylabel('time (s)')
        plt.title('clean vis of '+str(fname)+', '+str(ant)+'h', y=1.12)
        plt.colorbar()
        plt.twiny()
        plt.imshow(vis_clean_hh,aspect='auto',extent=(0,len(data.freqs),len(timestamps)*2,0))
        plt.xlabel('channel')
        plt.subplot(122)
        plt.imshow(vis_clean_vv,aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps)*2,0))
        plt.xlabel('Freq (MHz)')
        plt.ylabel('time (s)')
        plt.title('clean vis of '+str(fname)+', '+str(ant)+'v', y=1.12)
        plt.colorbar()
        plt.twiny()
        plt.imshow(vis_clean_vv,aspect='auto',extent=(0,len(data.freqs),len(timestamps)*2,0))
        plt.xlabel('channel')
        #plt.savefig(str(fname)+'raw_vis.pdf', bbox_inches='tight')
        plt.savefig(output_dir + fname + '_' + ant + '_clean_vis.png', bbox_inches='tight')
        #plt.show()

        #show channel changes

        for index, pol in enumerate(['h','v']):
            recv=ant+pol
            if pol=='h':
                vis_clean_pol=vis_clean_hh
            if pol=='v':
                vis_clean_pol=vis_clean_vv

            if pol == 'h':
                vis_div_x = np.ma.array(vis_div, mask=flags_h)
            if pol == 'v':
                vis_div_x = np.ma.array(vis_div, mask=flags_v)
            num_mask_tx = np.mean(vis_div_x.mask, axis=0)
            mode_tx = mode(num_mask_tx)
            num_mask_chx = np.mean(vis_div_x.mask, axis=1)
            mode_chx = mode(num_mask_chx)
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            axes[0].plot(freqs / 1e6, num_mask_tx)
            axes[0].set_xlabel('freqs')
            axes[0].set_ylabel('masked ratio')
            axes[0].set_ylim(top=1.13)
            axes[0].set_title('ratios')
            axes[0].text(580, 1.05, f'mode:{round(mode_tx.mode[0], 3)}', fontsize=18, color='red')
            axes[0].text(830, 1.05, f'mode count:{mode_tx.count[0]}', fontsize=18, color='red')
            # print(222222222, round(mode_t.mode[0], 3), mode_t.count[0])
            axes[1].plot(timestamps - timestamps[0], num_mask_chx)
            axes[1].set_xlabel('timestamps')
            axes[1].set_ylabel('masked ratio')
            axes[1].set_ylim(top=1.13)
            axes[1].set_title('ratios')
            axes[1].text(800, 1.05, f'mode:{round(mode_chx.mode[0], 3)}', fontsize=18, color='red')
            axes[1].text(4800, 1.05, f'mode count:{mode_chx.count[0]}', fontsize=18, color='red')
            # include the flags
            plt.savefig(output_dir + fname + '_' + recv + '_masked_ratio.png')

            plt.figure(figsize=(24,20))
            for i in range(len(ch_plot_list)):
                ch_plot1=ch_plot_list[i]
                p_data=vis_clean_pol[nd_s0,ch_plot1]
                plt.subplot(5,2,i+1)
                plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
                if i>len(ch_plot_list)-3:
                    plt.xlabel('R.A.')
                plt.ylabel('Dec')
                plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.savefig(output_dir + fname + '_' + recv + '_scatter_map.png')
            #plt.show()

            plt.figure(figsize=(24,20))
            plt.subplots_adjust(wspace =0)
            for i in range(len(ch_plot_list)):
                try:
                    ch_plot1=ch_plot_list[i]
                    p_data=vis_clean_pol[nd_s0,ch_plot1]
                    plt.subplot(5,2,i+1)
                    #kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=60)
                    kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    if i>len(ch_plot_list)-3:
                        plt.xlabel('R.A.')
                    plt.ylabel('Dec')
                    plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                    #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                except ValueError:
                    xi = np.linspace(ra[nd_s0].min(), ra[nd_s0].max(), 30)
                    yi = np.linspace(dec[nd_s0].min(), dec[nd_s0].max(), 30)
                    # zi = griddata((ra[nd_s0], dec[nd_s0]), p_data, (xi[None, :], yi[:, None]), method = 'linear')
                    zi = np.full((len(yi), len(xi)), np.nan)
                    CS = plt.imshow(zi[::-1, :], extent=(ra[nd_s0].min(),ra[nd_s0].max(),dec[nd_s0].min(),dec[nd_s0].max()), cmap=plt.get_cmap('jet',255), vmin=0, vmax=1, aspect='auto')
                    plt.colorbar(CS)
                    CS = plt.contour(xi, yi, zi, levels=15, linewidths=0.5, colors='k')
                    plt.xlim(ra[nd_s0].min(), ra[nd_s0].max())
                    plt.ylim(dec[nd_s0].min(), dec[nd_s0].max())
            plt.savefig(output_dir + fname + '_' + recv + '_map.png')
            #plt.show()
        print ('### '+ant+' of '+fname+' finished successfully ###')
        print ('###############################################')
        
        metadata = {
            "parameters": {
                "fname": fname,
                "ant": ant,
                "Threshold_factors": [Threshold_factor1, Threshold_factor2, Threshold_factor11, Threshold_factor22]
            },
            "ch_plot_list": ch_plot_list,
            "ch_plot": ch_plot,
            "ch_ref": ch_ref,
            "description": '### '+ant+' of '+fname+' finished successfully ###'
        }
        metadata_path = os.path.join(output_dir, f"{fname}_{ant}_metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=4)
            
    except(Exception):
        print ('*** '+ant+' of '+fname+' skipped becasue of some reason***')
        metadata = {
            "parameters": {
                "fname": fname,
                "ant": ant,
                "Threshold_factors": [Threshold_factor1, Threshold_factor2, Threshold_factor11, Threshold_factor22]
            },
            "description": '*** '+ant+' of '+fname+' skipped becasue of some reason***'
        }
        metadata_path = os.path.join(output_dir, f"{fname}_{ant}_metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=4)
        import traceback
        traceback.print_exc()
        raise
else:
    print ('*** '+ant+' of '+fname+' is a bad ant, skipped***')

del data
print ('# end @ ' + time.asctime(time.localtime(time.time())) +'#')

#END of Level 1 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
