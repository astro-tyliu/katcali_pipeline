#!/usr/bin/env python3
import matplotlib
matplotlib.use('AGG')

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
from matplotlib.colors import LogNorm
import time
import os
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
import katcali.beam_UHF as kb_u
from astropy.coordinates import SkyCoord
from astropy import units as u
import json
print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')
print (katcali.__version__)
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])


def cal_tag(target):
    #print (target)
    if target=='PKS1934-638':
        tag='1934-638'
    if target=='PKS2153-69':
        tag='2153-69'
    if target=='PictorA':
        tag='PictorA'
    if target=='3C273':
        tag='3C273'
    if target=='HydraA':
        tag='HydraA'
    if target=='3C353':
        tag='3C353'
    if target=='':
        tag=''
    return tag


# All scan modes (indices)
# dp_slew: slew
# dp_stop: stop
# dp_t: track
# dp_s: scan
# dp_w: slew and stop
# dp_tt: track and no flagging and close to the calibrator
# dp_ss: scan and no flagging
# dp_u: dp_tt + dp_ss
# dp_f: flagging for scan of (flagging or far from) for track
# dp_sb: the first dp_ss
# dp_se: the last dp_ss
# dp_ptr_list: any calibrators within ang_lim (0.5 deg) within dp_ss

# dp_c0a - dp_c4a: dp_tt when observing the target 0/1/2/3/4 before dp_sb (drop two edges), where dp_c1a equals to the union of dp_c1a1,dp_c1a2 and dp_c1a3
# dp_c0b - dp_c4b: dp_tt when observing the target 0/1/2/3/4 after dp_se (drop two edges), where dp_c1b equals to the union of dp_c1b1,dp_c1b2 and dp_c1b3
# dp_c1a1: before the earliest dp_c2a
# dp_c1a2: between dp_c1a1 and dp_c1a3
# dp_c1a3: after the latest dp_c3a
# dp_c1b1: before the earliest dp_c2b
# dp_c1b2: between dp_c1b1 and dp_c1b3
# dp_c1b3: after the latest dp_c3b
# dp_ca: the union of the dp_c0a - dp_c4a
# dp_cb: the union of the dp_c0b - dp_c4b
# dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b,dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3

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
# nd_t1x_ca: nd_t1 before scan
# nd_t1x_cb: nd_t1 after scan

# Select an observation block and load basic information in
fname=sys.argv[1]#'1684087370'#'1675632179'#
# Select ant and polarization, then load data in 
ant=sys.argv[2]#'m060'
pol=sys.argv[3]#'v'
input_file_name=sys.argv[4]
recv=ant+pol
output_file=sys.argv[5] + f'/{fname}_{recv}/'

input_file=f'/scratch3/users/liutianyang/katcali_pipeline/level1/py_results/{input_file_name}/'
# output_file=f'/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/cali_{file_timestamp}/{fname}_{recv}/'

data=kio.load_data(fname)
target_list,c0_list,bad_ants,flux_model_list=kio.check_ants(fname)

tag_list=[] #modeified for xcalib
if isinstance(target_list, list):
    for target in target_list:
        tag=cal_tag(target)
        tag_list.append(tag)
else:
    #print ('same calibrator for track-I/II')
    tag_list=cal_tag(target_list)
print (tag_list, flush=False)

target_start_list=[] #modeified for xcalib
if isinstance(target_list, list):
    for tag in tag_list:
        print (tag)
        if tag != '':
            data.select()
            data.select(targets=tag+'_u0.8')
            target_start=data.target_indices[0]
            target_start_list.append(target_start)
            data.select()
        else:
            target_start_list.append('')
else:
    data.select()
    data.select(targets=tag_list+'_u0.8')
    target_start_list=data.target_indices[0]
    data.select()
print (target_start_list)

ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print (str(i) + ' is bad')
        
print (fname)
print (ants_good, flush=False)

d={}
d2={}
count_ch=0

if ant in ants_good:

    # # ----------------------------------------------------------------------------------------- #
    # fs=open(output_file+'level2_data','rb')
    # d = pickle.load(fs)
    # fs.close()
    # T_map = d['T_map']
    # Tresi_map = d['Tresi_map']
    # gain_map = d['gain_map']
    # Tel_map = d['Tel_map']
    # Tnd_ref_list = d['Tnd_ref_list']
    # Tnda_list = d['Tnda_list']
    # Tndb_list = d['Tndb_list']
    # Tnd_diff_ratio_list = d['Tnd_diff_ratio_list']
    # NRMSE1_list = d['NRMSE1_list']
    # NRMSE2_list = d['NRMSE2_list']
    # gta_param = d['gta_param']
    # gtb_param = d['gtb_param']
    # sma_param = d['sma_param']
    # smb_param = d['smb_param']
    # timestamps = d['timestamps']
    # nd_0 = d['nd_0']
    # ra = d['ra']
    # dec = d['dec']

    # fs=open(output_file+'level2_Tnd_data','rb')
    # d2 = pickle.load(fs)
    # fs.close()
    # Tnd_ref_list = d2['Tnd_ref_list']
    # Tnda_list = d2['Tnda_list']
    # Tndb_list = d2['Tndb_list']
    # Tnd_diff_ratio_list = d2['Tnd_diff_ratio_list']
    # NRMSE1_list = d2['NRMSE1_list']
    # NRMSE2_list = d2['NRMSE2_list']

    # iii = np.where(np.array(Tnda_list)!=None)[0][-1] + 1
    # if iii < 2869:
    #     channels_cali=list(range(iii,2869))+list(range(3133,3547))
    # elif iii >= 2869 and iii < 3133:
    #     channels_cali=list(range(3133,3547))
    # elif iii >= 3133 and iii < 3547:
    #     channels_cali=list(range(iii,3547))
    # else:
    #     print(f"{recv} has been finished!")
    #     sys.exit()
    
    channels_cali=list(range(272,2869))+list(range(3133,3547))
    # channels_cali=list(range(3170,3230))
    # # ----------------------------------------------------------------------------------------- #

    nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
    print (nd_on_time,nd_cycle,nd_set, flush=False)

    #load data, labels, and parameters
    ch_ref=3200  #ch3608 is for 1023MHz, but Tnd and Trec models cut at 1015MHz
    ch_plot=3200
    check_ratio_limit=0.3
    data.select(ants=ant,pol=pol)
    corr_id=kio.cal_corr_id(data,recv)
    assert(recv==data.corr_products[corr_id][0])
    assert(recv==data.corr_products[corr_id][1])
    print (corr_id,recv)
    vis,flags= kio.call_vis(fname,recv)
    vis_backup=vis.copy()
    
    ra,dec,az,el=kio.load_coordinates(data)
    timestamps,freqs=kio.load_tf(data)
    dump_period=data.dump_period
    if isinstance(target_list, list):
        ang_deg=kio.load_ang_deg2(ra,dec,c0_list)
    else:
        ang_deg=kio.load_ang_deg(ra,dec,c0_list) #modeified for xcalib
    ang_deg=np.array(ang_deg)
    dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
    '''
    az_corr=az.copy()
    az_corr[az>180]-=360
    az=az_corr.copy()
    '''
    az[az>180]-=360
    
    nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
    print (len(nd_on_edge),len(nd_off_edge))
    nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
    
    # RFI flagging
    #check with .py result
    try:
        d3 = pickle.load(open(input_file+fname+'_'+ant+'/mask2', 'rb'))
        print ('mask2 loaded')
    except(Exception):
        d3 = pickle.load(open(input_file+fname+'_'+ant+'/mask', 'rb'))
        print ('mask loaded')

    os.makedirs(output_file, exist_ok=True)
        
    mask_inter=d3['mask']
    vis_clean=np.ma.array(vis,mask=flags)
    vis_clean=np.ma.array(vis_clean,mask=mask_inter)
    
    #load the scan and track labels
    dp_u=kl.cal_dp_u(dp_tt,dp_ss)
    
    ### full-band models prepare####
    
    #select beam pattern model
    beam_select='HiRes'
    #BM-III:pattern
    dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b,dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3=kl.cal_dp_c(fname,data,ant,pol,flags,ch_ref,dp_tt,dp_ss,ang_deg, target_start=target_start_list,n_src_off=4)
    Npix=256
    Ddeg=6
    pattern_fband=kb_u.load_pattern_fband(beam_select,pol)
    Aeff_max_fband=kb_u.load_Aeff_max_fband(beam_select,pol)
    # x_pix,y_pix=kb_u.cal_pix_params(data,c0,Npix,Ddeg)
    if isinstance(c0_list, list):
        x_pix,y_pix=[],[] #modified for xcalib
        for c0 in c0_list:
            if c0 is None:
                x_pix.append('')
                y_pix.append('')
            else:
                x_pix_i,y_pix_i=kb_u.cal_pix_params(data,c0,Npix,Ddeg)
                x_pix.append(x_pix_i)
                y_pix.append(y_pix_i)
            
        x_pix,y_pix=np.array(x_pix),np.array(y_pix)
    else:
        x_pix,y_pix=kb_u.cal_pix_params(data,c0_list,Npix,Ddeg)
    
    Trec_list=km.cal_Trec(data,ant,pol,freqs,band='UHF') #580-1015 MHz
    
    #Galactic model
    nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
    #gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, ch_plot, ch_plot+1, nside)
    gal_ori=km.cal_Gal_model_np2(vis, freqs, ra, dec, 0, len(freqs), nside, model_key=-1)
    print ('#Gal model is from Halsam!!!')
    gal_ori.flags.writeable=False #avoid change by mistake
    gal=gal_ori.copy()
    
    #select raw vis for track befor scan
    visa_ptr = vis_clean.copy()
    visb_ptr = vis_clean.copy()
    
    for i in range(len(timestamps)):
        if i not in dp_ca:
            visa_ptr.mask[i,:]=True
        if i not in dp_cb:
            visb_ptr.mask[i,:]=True
    
    #diode noise
    Tnd_std,Tnd_spl=km.Tnd_spl(data, ant,pol)
    
    #Spill
    Tspill_func=km.cal_Tspill_func(el,pol,freqs,band='UHF')
    
    ################################
    ####prepare for data storage#################
    T_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
    gain_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tel_map=np.ma.array(np.zeros_like(vis),mask=True)
    
    Tnd_ref_list=[None for i in range(len(freqs))]
    Tnda_list=[None for i in range(len(freqs))]
    Tndb_list=[None for i in range(len(freqs))]
    Tnd_diff_ratio_list=[None for i in range(len(freqs))]
    NRMSE1_list=[None for i in range(len(freqs))]
    NRMSE2_list=[None for i in range(len(freqs))]
    
    gta_param=[None for i in range(len(freqs))]
    gtb_param=[None for i in range(len(freqs))]
    
    sma_param=[None for i in range(len(freqs))]
    smb_param=[None for i in range(len(freqs))]
    ################################################
    
    for ch_cali in channels_cali:
    
        check_ratio=float(np.array(vis_clean.mask[:,ch_cali]==False).sum())/len(timestamps)
        print (ch_cali, check_ratio, flush=False)
        if check_ratio>check_ratio_limit: #0.6
            try:
       
                #spill model 
                ### Tspill=km.cal_Tspill(el,pol,freqs, ch_cali,'UHF') #works but not efficiency
                Tspill=Tspill_func((el,freqs[ch_cali]/1e6))[:,0]
            
                #atmosphere emission model
                Tatmo=km.calc_atmosphere_model_1ch(data,ch_cali)
                
                #atmosphere emission model
                Tatmo_abs_1ch=km.calc_atmosphere_trans_factor_1ch(data,ch_cali) #old name:calc_atmosphere_abs_factor_1ch
                
                #elevation related emission model
                Tel=Tspill+Tatmo 
                
                #load the diode injection model and get a reference value
                #note: diode version are different dish by dish!
                #Tnd_std,Tnd_ref,noise,Tnd_spl= km.call_Tnd(data, ant, pol,freqs,ch_cali,1) #works but not efficiency
                Tnd_ref=Tnd_spl(freqs[ch_cali]/1e9)
                print (Tnd_std,Tnd_ref, flush=False)
                
                # T_ptr2,pattern,x_pix_max,y_pix_max=kb_u.cal_BMIII_1ch(data,ch_cali,flux_model, dp_ca,dp_cb,pattern_fband,x_pix,y_pix,Aeff_max_fband)
                if isinstance(target_list, list):
                    #modified for xcalib
                    T_ptr2a=np.full(len(timestamps),np.NaN)
                    T_ptr2b=np.full(len(timestamps),np.NaN)
                    for i in range(np.shape(x_pix)[0]):
                        if x_pix[i]=='':
                            print ('# track-'+str(i+1)+' is empty')
                        else:
                            if i==0:
                                T_ptr2a,pattern,x_pix_max,y_pix_max=kb_u.cal_BMIII_1ch(data,ch_cali,flux_model_list[0], dp_ca,[],pattern_fband,x_pix[0],y_pix[0],Aeff_max_fband)
                            if i==1:
                                T_ptr2b,pattern,x_pix_max,y_pix_max=kb_u.cal_BMIII_1ch(data,ch_cali,flux_model_list[1], [],dp_cb,pattern_fband,x_pix[1],y_pix[1],Aeff_max_fband)
                    
                    T_ptr2=np.nanmax([T_ptr2a,T_ptr2b],axis=0)
                else:
                    T_ptr2,pattern,x_pix_max,y_pix_max=kb_u.cal_BMIII_1ch(data,ch_cali,flux_model_list, dp_ca,dp_cb,pattern_fband,x_pix,y_pix,Aeff_max_fband)
                
                Tgal=gal[:,ch_cali]
                
                #####choose beam model
                T_ptr=T_ptr2 #BM-III 
                
                ####set input parameters
                # ga0,gb0=ks.cal_gain0(fname,data,ant,pol,flags,ch_cali,dp_tt,dp_ss,ang_deg,T_ptr,vis_clean,n_src_off=4,target_start=target_start_list) #gain level ##need to check the label still works for 2021 data!!!!!!
                ga0,gb0=ks.cal_gain0(fname,data,ant,pol,flags,ch_cali,dp_tt,dp_ss,ang_deg,T_ptr,vis_clean,n_src_off=4,target_start=target_start_list)
                print (ga0,gb0, flush=False)
                assert(isinstance(ga0,np.float))
                assert(isinstance(gb0,np.float))
                
                Trec0=Trec_list[ch_cali]
                
                if target_list[0]!='':
                    print (Trec0)
                    eta_p0=1.0
                    func_sm_param0=[Trec0]
                    func_gt_param0=[ga0,0,0,0,0]
                    print(np.mean(timestamps), np.std(timestamps))
                    print(np.ma.mean(visa_ptr), np.ma.std(visa_ptr))
                    print(ch_cali, eta_p0, Tnd_ref, Tnd_std)
                    print(np.ma.mean(nd_ratio), np.ma.std(nd_ratio))
                    print(np.ma.mean(T_ptr), np.ma.std(T_ptr))
                    print(np.ma.mean(Tel), np.ma.std(Tel))
                    print(np.ma.mean(Tgal), np.ma.std(Tgal))
                    print(func_gt_param0)
                    print(func_sm_param0)
                    print(np.ma.mean(nd_0), np.ma.std(nd_0))
                    print(np.mean(nd_1x), np.std(nd_1x))
                    print(np.mean(dp_ca), np.std(dp_ca))
                    print(np.ma.mean(vis_clean), np.ma.std(vis_clean))
                    print (Tnd_ref, 11111, flush=False)
                    ####fitting
                    instru_pa=ks.solve_params0_v3(timestamps, visa_ptr, ch_cali, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1x, band='UHF')
                    ####get fitting result
                    print(22222, flush=False)
                    Tnda=instru_pa[0]
                    eta_pa=instru_pa[1]
                    sma=instru_pa[2]
                    gta=instru_pa[3:]
                    print (Tnda, eta_pa, sma, gta, 33333, flush=False)

                    m1=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pa, Tnda, Tel, Tgal, gta, sma, nd_0, nd_1x)
                    ma=np.ma.array(m1,mask=visa_ptr[:,ch_cali].mask)    
                    ga=ks.func_gt(timestamps,gta)
                    resia=(visa_ptr[:,ch_cali]-ma)/ga

                if target_list[-1]!='':
                    ####set input parameters
                    #Trec0=Trec_list[ch_plot]
                    print (Trec0, flush=False)
                    eta_p0=1.0
                    func_sm_param0=[Trec0]
                    func_gt_param0=[gb0,0,0,0,0]
                    
                    print (Tnd_ref, flush=False)
                    
                    ####fitting######
                    instru_pb=ks.solve_params0_v3(timestamps, visb_ptr, ch_cali, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal,
                                          func_gt_param0, func_sm_param0, nd_0, nd_1x, band='UHF')
                    
                    ######get fitting result#####
                    Tndb=instru_pb[0]
                    eta_pb=instru_pb[1]
                    smb=instru_pb[2]
                    gtb=instru_pb[3:]
                    print (Tndb, eta_pb, smb, gtb, flush=False)

                    m2=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pb, Tndb, Tel, Tgal, gtb, smb, nd_0, nd_1x)
                    mb=np.ma.array(m2,mask=visb_ptr[:,ch_cali].mask)
                    gb=ks.func_gt(timestamps,gtb)
                    resib=(visb_ptr[:,ch_cali]-mb)/gb
                
                if ch_cali == ch_plot:
                    ##show model and raw vis
                    plt.figure(figsize=(20,8))
                    plt.subplots_adjust(wspace=0.1,hspace=0.25)
                    plt.subplot(221)
                    if target_list[0]!='':
                        m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pa, Tnda, Tel, Tgal, gta, sma, nd_0, nd_1x)
                        ma=np.ma.array(m,mask=visa_ptr[:,ch_cali].mask)
                        plt.plot(timestamps[dp_c0a]-timestamps[0],ma[dp_c0a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c0a]-timestamps[0],visa_ptr[dp_c0a,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],ma[dp_c1a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],visa_ptr[dp_c1a,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2a]-timestamps[0],ma[dp_c2a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],ma[dp_c3a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],ma[dp_c4a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2a]-timestamps[0],visa_ptr[dp_c2a,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],visa_ptr[dp_c3a,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],visa_ptr[dp_c4a,ch_cali],'g--',drawstyle='steps-mid')
                        #plt.xlabel('time (s)')
                        plt.ylabel('raw signal')
                        plt.legend(['model','raw signal'],ncol=2)
                        plt.ylim(np.nanmin(visa_ptr[dp_c0a,ch_cali])-20,np.nanmax(visa_ptr[dp_c1a,ch_cali])+30)
                        plt.title('Fitting curve of track-I at '+ str(int(freqs[ch_cali]/1e6)) +' MHz, '+str(fname)+', '+str(recv))
                    plt.subplot(223)
                    if target_list[0]!='':
                        ga=ks.func_gt(timestamps,gta)
                        resia=(visa_ptr[:,ch_cali]-ma)/ga
                        plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c0a][-1]-timestamps[0]],[0,0],'k-.')
                        plt.plot(timestamps[dp_c0a]-timestamps[0],resia[dp_c0a],'.-')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],resia[dp_c1a],'.-')
                        plt.plot(timestamps[dp_c2a]-timestamps[0],resia[dp_c2a],'.-')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],resia[dp_c3a],'.-')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],resia[dp_c4a],'.-')
                        plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c4a][-1]-timestamps[0]],[0,0],'k-.')
                        plt.xlabel('time (s)')
                        plt.ylabel('Tres (K)')
                        plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'], fontsize=12, ncol=3)
                        plt.ylim(-0.19,0.19)
                        plt.title('$T_{res}$ of the fitted curve above')
                    ##show model and raw vis
                    plt.subplot(222)
                    if target_list[-1]!='':
                        m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pb, Tndb, Tel, Tgal, gtb, smb, nd_0, nd_1x)
                        mb=np.ma.array(m,mask=visb_ptr[:,ch_cali].mask)
                        plt.plot(timestamps[dp_c0b]-timestamps[0],mb[dp_c0b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c0b]-timestamps[0],visb_ptr[dp_c0b,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],mb[dp_c1b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],visb_ptr[dp_c1b,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2b]-timestamps[0],mb[dp_c2b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],mb[dp_c3b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],mb[dp_c4b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2b]-timestamps[0],visb_ptr[dp_c2b,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],visb_ptr[dp_c3b,ch_cali],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],visb_ptr[dp_c4b,ch_cali],'g--',drawstyle='steps-mid')
                        plt.ylim(np.nanmin(visb_ptr[dp_c0b,ch_cali])-20,np.nanmax(visb_ptr[dp_c1b,ch_cali])+30)
                        #plt.xlabel('time (s)')
                        #plt.ylabel('raw signal')
                        plt.legend(['model','raw signal'],ncol=2)
                        plt.title('Fitting curve of track-II at '+ str(int(freqs[ch_cali]/1e6)) +' MHz')
                    plt.subplot(224)
                    if target_list[-1]!='':
                        gb=ks.func_gt(timestamps,gtb)
                        resib=(visb_ptr[:,ch_cali]-mb)/gb
                        plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c0b][-1]-timestamps[0]],[0,0],'k-.')
                        plt.plot(timestamps[dp_c0b]-timestamps[0],resib[dp_c0b],'.-')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],resib[dp_c1b],'.-')
                        plt.plot(timestamps[dp_c2b]-timestamps[0],resib[dp_c2b],'.-')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],resib[dp_c3b],'.-')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],resib[dp_c4b],'.-')
                        plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c4b][-1]-timestamps[0]],[0,0],'k-.')
                        plt.xlabel('time (s)')
                        #plt.ylabel('residual (K)')
                        plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'],fontsize=12, ncol=3)
                        #plt.title('calibrator Part II of '+str(fname)+', '+str(recv)+ ' @'+ str(round(freqs[ch_cali]/1e9,3)) +' GHz')
                        plt.ylim(-0.19,0.19)
                        plt.title('$T_{res}$ of the fitted curve above')
                    plt.savefig(output_file + 'F_cali_ch'+str(ch_cali) + '.pdf', bbox_inches='tight')
                    #plt.show()
                
                
                print ('Tnd_ref: ', Tnd_ref)
                ####data need to storage######
                Tnd_ref_list[ch_cali]=Tnd_ref
                if target_list[0]!='':
                    print ('Tnda: ', Tnda)

                    NRMSE1=ks.cal_NRMSE(ma/ga,resia)
                    print ('NRMSE1 =', NRMSE1)

                    if Tnda>0 and Tnda<100:
                        Tnda_list[ch_cali]=Tnda
                        NRMSE1_list[ch_cali]=NRMSE1
                    
                if target_list[-1]!='':
                    print ('Tndb: ', Tndb)

                    NRMSE2=ks.cal_NRMSE(mb/gb,resib)
                    print ('NRMSE2 =', NRMSE2)

                    if Tndb>0 and Tndb<100:
                        Tndb_list[ch_cali]=Tndb
                        NRMSE2_list[ch_cali]=NRMSE2
                    
                if target_list[0]!='' and target_list[-1]!='':
                    print ('average: ', (Tnda+Tndb)/2.)
                    Tnd_diff_ratio=abs(Tnda-Tndb)/np.max([Tnda,Tndb])
                    print ('Tnd_diff_ratio:', str(round(Tnd_diff_ratio*100,2)),'%')

                    if Tnda >0 and Tndb >0 and Tnda<100 and Tndb<100:
                        Tnd_diff_ratio_list[ch_cali]=Tnd_diff_ratio
                
                gta_param[ch_cali]=gta
                gtb_param[ch_cali]=gtb
                
                gain_map[dp_ca,ch_cali]=ga[dp_ca]
                gain_map[dp_cb,ch_cali]=gb[dp_cb]
                
                calT_tra=visa_ptr[:,ch_cali]/ga
                calT_trb=visb_ptr[:,ch_cali]/gb
                if Tnda>0:
                    assert((abs(calT_tra[dp_ca]-ma[dp_ca]/ga[dp_ca]-resia[dp_ca])<1e-10).all()==True)
                if Tndb>0:
                    assert((abs(calT_trb[dp_cb]-mb[dp_cb]/gb[dp_cb]-resib[dp_cb])<1e-10).all()==True)
                T_map[dp_ca,ch_cali]=calT_tra[dp_ca]
                T_map[dp_cb,ch_cali]=calT_trb[dp_cb]
                
                Tresi_map[dp_ca,ch_cali]=resia[dp_ca]
                Tresi_map[dp_cb,ch_cali]=resib[dp_cb]
                
                sma_param[ch_cali]=sma
                smb_param[ch_cali]=smb
                
                Tel_map[:,ch_cali]=Tel
                
                print (ch_cali, Tnd_ref, Tnda, Tndb)

                if ch_cali%50==0:
                    #d={}
                    d['T_map']=T_map
                    d['Tresi_map']=Tresi_map
                    d['gain_map']=gain_map
                    d['Tel_map']=Tel_map
                    d['Tnd_ref_list']=Tnd_ref_list
                    d['Tnda_list']=Tnda_list
                    d['Tndb_list']=Tndb_list
                    d['Tnd_diff_ratio_list']=Tnd_diff_ratio_list
                    d['NRMSE1_list']=NRMSE1_list
                    d['NRMSE2_list']=NRMSE2_list
                    d['gta_param']=gta_param
                    d['gtb_param']=gtb_param
                    d['sma_param']=sma_param
                    d['smb_param']=smb_param
                    d['timestamps']=timestamps
                    d['nd_0']=nd_0
                    d['ra']=ra
                    d['dec']=dec
                    fs=open(output_file+'level2_data','wb')
                    pickle.dump(d,fs,protocol=2)
                    fs.close()
                    
                    #d2={}
                    d2['Tnd_ref_list']=Tnd_ref_list
                    d2['Tnda_list']=Tnda_list
                    d2['Tndb_list']=Tndb_list
                    d2['Tnd_diff_ratio_list']=Tnd_diff_ratio_list
                    d2['NRMSE1_list']=NRMSE1_list
                    d2['NRMSE2_list']=NRMSE2_list
                    fs=open(output_file+'level2_Tnd_data','wb')
                    pickle.dump(d2,fs,protocol=2)
                    fs.close()

                print ('***channel '+ str(ch_cali) +' finished')    
                count_ch+=1

            except Exception as error:
                print("An error occurred:", error) 
                print ('***channel '+ str(ch_cali) +' failed...')
                print (np.shape(data))
                data.select()
                data.select(ants=ant,pol=pol)
                print (np.shape(data))
                print ('data reset applied.')
        else:
            print ('***channel '+ str(ch_cali) +' skipped')   
        print(11111, flush=False)
            
    ####save data####
    #d={}
    d['T_map']=T_map
    d['Tresi_map']=Tresi_map
    d['gain_map']=gain_map
    d['Tel_map']=Tel_map
    d['Tnd_ref_list']=Tnd_ref_list
    d['Tnda_list']=Tnda_list
    d['Tndb_list']=Tndb_list
    d['Tnd_diff_ratio_list']=Tnd_diff_ratio_list
    d['NRMSE1_list']=NRMSE1_list
    d['NRMSE2_list']=NRMSE2_list
    d['gta_param']=gta_param
    d['gtb_param']=gtb_param
    d['sma_param']=sma_param
    d['smb_param']=smb_param
    d['timestamps']=timestamps
    d['nd_0']=nd_0
    d['ra']=ra
    d['dec']=dec
    fs=open(output_file+'level2_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()
    
    #d2={}
    d2['Tnd_ref_list']=Tnd_ref_list
    d2['Tnda_list']=Tnda_list
    d2['Tndb_list']=Tndb_list
    d2['Tnd_diff_ratio_list']=Tnd_diff_ratio_list
    d2['NRMSE1_list']=NRMSE1_list
    d2['NRMSE2_list']=NRMSE2_list
    fs=open(output_file+'level2_Tnd_data','wb')
    pickle.dump(d2,fs,protocol=2)
    fs.close()

    ###begin: reference plot only###

    Tnda_list=kio.mask_None(Tnda_list)
    Tndb_list=kio.mask_None(Tndb_list)

    import matplotlib.gridspec as gridspec
    plt.figure(figsize=(16,5))
    plt.subplots_adjust(hspace=0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0,0])
    plt.plot(Tnda_list,'.',ms=4)
    plt.plot(Tndb_list,'.',ms=4)
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    plt.legend(['Tnda','Tndb','Tnd_ref'],ncol=3)
    #plt.xlabel('channel')
    plt.xticks([])
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.subplot(gs[1,0])
    plt.plot(Tnd_diff_ratio_list,'m.',ms=4)
    plt.xlabel('channel')
    plt.ylabel('Tnd_diff_ratio')
    plt.savefig(output_file+'Tnd_all.pdf', bbox_inches='tight')

    metadata = {
        "parameters": {
            "fname": fname,
            "ant": ant,
            "pol": pol,
            "input_file": input_file_name
        },
        "ch_ref": ch_ref,
        "ch_plot": ch_plot,
        "check_ratio_limit": check_ratio_limit,
        "channels_cali": channels_cali,
        "description": '### '+recv+' of '+fname+' finished successfully ###'
    }
    metadata_path = os.path.join(output_file, f"metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

else:
    metadata = {
        "parameters": {
            "fname": fname,
            "ant": ant,
        },
        "description": '***BAD ANT, Level2 calibartion for '+fname+' '+ant+pol+' skipped'
    }
    metadata_path = os.path.join(output_file, f"metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)
    print ('***BAD ANT, Level2 calibartion for '+fname+' '+ant+pol+' skipped')
    print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata') 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
