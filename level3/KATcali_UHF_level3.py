#!/usr/bin/env python3

#imports
import os
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
import json
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
import katcali.beam as kb

from astropy.coordinates import SkyCoord
from astropy import units as u

print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')

print (katcali.__version__)

print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])

# Select an observation block and load basic information in
fname=sys.argv[1] #'1675632179'
#select ant, polarization, and one channel to show data calibration
ant=sys.argv[2] #'m059'
pol=sys.argv[3] #'v'
input_file1_name=sys.argv[4]
input_file2_name=sys.argv[5]
file_timestamp=sys.argv[6]
recv=ant+pol
ch_ref=3200

input_file1=f'/scratch3/users/liutianyang/katcali_pipeline/level1/py_results/{input_file1_name}/'
input_file2=f'/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/{input_file2_name}/'
output_file=f'/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/{file_timestamp}/'

d2=pickle.load(open(input_file2+fname+'_'+str(recv)+'_level2_Tnd_data','rb'),encoding='latin-1')
print(d2.keys())
Tnd_ref_list=d2['Tnd_ref_list']
Tnda_list=d2['Tnda_list']
Tndb_list=d2['Tndb_list']
Tnd_diff_ratio_list=d2['Tnd_diff_ratio_list']
NRMSE1_list = d2['NRMSE1_list']
NRMSE2_list = d2['NRMSE2_list']

Tnda_list=np.array(Tnda_list,dtype=np.float64)
Tnda_list=np.ma.masked_invalid(Tnda_list)
Tndb_list=np.array(Tndb_list,dtype=np.float64)
Tndb_list=np.ma.masked_invalid(Tndb_list)
Tnd_diff_ratio_list=np.array(Tnd_diff_ratio_list,dtype=np.float64)
Tnd_diff_ratio_list=np.ma.masked_invalid(Tnd_diff_ratio_list)
NRMSE1_list=np.array(NRMSE1_list,dtype=np.float64)
NRMSE1_list=np.ma.masked_invalid(NRMSE1_list)
NRMSE2_list=np.array(NRMSE2_list,dtype=np.float64)
NRMSE2_list=np.ma.masked_invalid(NRMSE2_list)

Tnd_list=np.ma.mean([Tnda_list,Tndb_list],axis=0)

data=kio.load_data(fname)

#show the calibrator and bad ants information
# target,c0,bad_ants,flux_model=kio.check_ants(fname)
target_list,c0_list,bad_ants,flux_model_list=kio.check_ants(fname)

ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print (str(i) + ' is bad')
        
print (fname)
print (ants_good)

# desi1_list = ['1678743988', '1682448988', '1678295187', '1676313206', '1677020482', '1675643846', '1678122565', '1676657789', '1677777992', '1675623808', '1677002481', '1677195529', '1684087370', '1675210948', '1678726283', '1678467685', '1679605292', '1677174749', '1678381591', '1675021905', '1677795989', '1678899080', '1675816512', '1675106912']
# desi2_list = ['1679592842', '1678734987', '1689003684', '1689176790', '1688399183', '1677183387', '1679247986', '1677011008', '1680626188', '1681143685', '1680798562', '1683492604', '1681920680', '1685641589', '1675632179', '1679419886', '1689090392', '1684781618', '1679615321', '1679333668', '1680644082', '1681229848']
# if fname in desi1_list:
#     classify_dir = '/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/imgs_desi1/'
# elif fname in desi2_list:
#     classify_dir = '/scratch3/users/liutianyang/katcali_pipeline/level2/py_results/imgs_desi2/'
# tmp = os.listdir(classify_dir + 'good/' + input_file2_name + '/')
# level2_good = [f[8:13] for f in tmp if f.endswith('.png')]
# print('good receivers in level 2: ', level2_good)
# del tmp

if ant in ants_good:
# if recv in level2_good:
    
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
    # ang_deg=kio.load_ang_deg(ra,dec,c0)
    if isinstance(target_list, list):
        ang_deg=kio.load_ang_deg2(ra,dec,c0_list)
    else:
        ang_deg=kio.load_ang_deg(ra,dec,c0_list) #modeified for xcalib
    ang_deg=np.array(ang_deg)
    dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
    dp_sb,dp_se=dp_ss[0],dp_ss[-1]
    
    nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
    print (nd_on_time,nd_cycle,nd_set)
    nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
    print (len(nd_on_edge),len(nd_off_edge))
    nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
    
    nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)
    labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)
    
    # RFI flagging
    try:
        d3 = pickle.load(open(input_file1+fname+'_'+ant+'_mask2','rb'))
        print ('# mask2 loaded')
    except(Exception):
        d3 = pickle.load(open(input_file1+fname+'_'+ant+'_mask','rb'))
        print ('# mask loaded')
                              
    mask_inter=d3['mask']
    vis_clean=np.ma.array(vis,mask=mask_inter)
    
    ###raw vis preparsion
    vis_clean_tt=vis_clean.copy()
    vis_clean_tt.mask[:dp_sb,:]=True
    vis_clean_tt.mask[dp_se+1:,:]=True
    
    ## load the foreground models
    dp_u=kl.cal_dp_u(dp_tt,dp_ss)
    
    ### full-band models prepare####
    
    #Trec_list=km.cal_Trec(data,ant,pol,freqs)
    Trec_list=km.cal_Trec(data,ant,pol,freqs,band='UHF') #580-1015 MHz
    
    #Galactic model
    nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
    #gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, ch_plot, ch_plot+1, nside)
    gal_ori=km.cal_Gal_model_np2(vis, freqs, ra, dec, 0, len(freqs), nside, model_key=-1)
    print ('#Gal model is from Halsam!!!')
    gal_ori.flags.writeable=False #avoid change by mistake
    gal=gal_ori.copy()
    
    #spill model 
    Tspill_func=km.cal_Tspill_func(el,pol,freqs,band='UHF')

    ####prepare for data storage#################
    gt_param=[None for i in range(len(freqs))]
    sm_param=[None for i in range(len(freqs))]
    T_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
    gain_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tsm_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tel_map=np.ma.array(np.zeros_like(vis),mask=True)

    channels_cali = list(range(272,2869))+list(range(3133,3547))
    
    if isinstance(target_list, str):
        target_list = [target_list, target_list]

    Tnd_diff_ratio_limit = None
    NRMSE_limit = None
    if target_list[0] == target_list[-1]:
        Tnd_diff_ratio_limit = 1.0  # 0.05
    elif target_list[0] == '' or target_list[-1] == '':
        NRMSE_limit = 1.0  # 0.004
    else:
        Tnd_diff_ratio_limit = 1.0  # 0.05
    
    ch_count=0
    
    d={}
    #####################################
    
    os.makedirs(output_file, exist_ok=True)
    for ch_plot in channels_cali:

        if target_list[0] == '':
            NRMSE2 = NRMSE2_list[ch_plot]
            print('target_list: ', target_list)
            print('NRMSE_limit: ', NRMSE_limit)
            print('NRMSE2: ', NRMSE2)
            judge = NRMSE2 < NRMSE_limit
        elif target_list[-1] == '':
            NRMSE1 = NRMSE1_list[ch_plot]
            print('target_list: ', target_list)
            print('NRMSE_limit: ', NRMSE_limit)
            print('NRMSE1: ', NRMSE1)
            judge = NRMSE1 < NRMSE_limit
        else:
            Tnd_diff_ratio=Tnd_diff_ratio_list[ch_plot]
            print('target_list: ', target_list)
            print('Tnd_diff_ratio_limit: ', Tnd_diff_ratio_limit)
            print('Tnd_diff_ratio: ', Tnd_diff_ratio)
            judge = Tnd_diff_ratio < Tnd_diff_ratio_limit
        
        if judge:
    
            try:
    
                Tspill=Tspill_func((el,freqs[ch_plot]/1e6))[:,0]
                Tatmo=km.calc_atmosphere_model_1ch(data,ch_plot)
                Tel=Tspill+Tatmo 
                
                Tgal=gal[:,ch_plot]
                
                ###set input params
                Tnd=Tnd_list[ch_plot]
                print (Tnd)
                assert(isinstance(Tnd,np.float))
                
                # calibration for scan part
                ####param0
                g0=10.
                Tptr=0 #no point source
                eta_p0=1.0
                Trec0=Trec_list[ch_plot]
                func_sm_param0=[Trec0,0,0,0]
                func_gt_param0=[g0,0,0,0,0]#must be [-5:] from func_obj_sm
                
                ##fitting
                instru_p=ks.solve_params_sm_v3(timestamps, vis_clean_tt, ch_plot, nd_ratio, Tptr, eta_p0, Tnd, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1x, band='UHF')
                
                ###output
                eta_p=instru_p[0]
                sm=instru_p[1:-5]
                gt=instru_p[-5:] 
                print (eta_p, sm, gt)
                
                m=ks.calc_total_model_sm_v3(timestamps, nd_ratio, Tptr, eta_p, Tnd, Tel, Tgal,  gt, sm, nd_0, nd_1x)
                residual=(vis_clean[:,ch_plot]-m)/ks.func_gt(timestamps,gt)
                T=vis_clean[:,ch_plot]/ks.func_gt(timestamps,gt)
                
                gain=ks.func_gt(timestamps,gt)
                Tsm=ks.func_sm(timestamps,sm)
                            
                #print (residual[nd_s0].mean(),residual[nd_s0].std())
                print ('residual in ch'+str(ch_plot)+': '+str(residual[nd_s0].mean())+' +/- '+str(residual[nd_s0].std()))
                ####data need to storage######
                sm_param[ch_plot]=sm
                gt_param[ch_plot]=gt
                T_map[nd_s0,ch_plot]=T[nd_s0]
                Tresi_map[nd_s0,ch_plot]=residual[nd_s0]
                gain_map[dp_sb:dp_se+1,ch_plot]=gain[dp_sb:dp_se+1]
                Tsm_map[dp_sb:dp_se+1,ch_plot]=Tsm[dp_sb:dp_se+1]
                Tel_map[:,ch_plot]=Tel

                print ('***channel'+ str(ch_plot) +' finished') 
                ch_count+=1

                if ch_plot%50==0:
                    ####save data####
                    d['gt_param']=gt_param
                    d['sm_param']=sm_param
                    d['T_map']=T_map
                    d['Tresi_map']=Tresi_map
                    d['gain_map']=gain_map
                    d['Tel_map']=Tel_map
                    d['Tgal_map']=gal_ori
                    d['Tsm_map']=Tsm_map
                    d['nd_s0']=nd_s0
                    d['timestamps']=timestamps
                    d['ra']=ra
                    d['dec']=dec
                    fs=open(output_file+f'{fname}_{recv}_level3_data','wb')
                    pickle.dump(d,fs,protocol=2)
                    fs.close()
                
            except Exception as error:
                print("An error occurred:", error) 
                print ('***channel'+ str(ch_plot) +' failed...')
                print (np.shape(data))
                data.select()
                data.select(ants=ant,pol=pol)
                print (np.shape(data))
                print ('data reset applied.')
        else:
            print ('***channel'+ str(ch_plot) +' skipped')
            
    metadata = {
        "parameters": {
            "fname": fname,
            "ant": ant,
            "pol": pol,
            "level1_input_file": input_file1_name,
            "level2_input_file": input_file2_name,
        },
        "ch_ref": ch_ref,
        "Tnd_diff_ratio_limit": Tnd_diff_ratio_limit,
        "NRMSE_limit": NRMSE_limit,
        "channels_cali": channels_cali,
        "description": '### '+recv+' of '+fname+' finished successfully ###'
    }
    metadata_path = os.path.join(output_file, f"{fname}_{recv}_metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

####save data####
d['gt_param']=gt_param
d['sm_param']=sm_param
d['T_map']=T_map
d['Tresi_map']=Tresi_map
d['gain_map']=gain_map
d['Tel_map']=Tel_map
d['Tgal_map']=gal_ori
d['Tsm_map']=Tsm_map
d['nd_s0']=nd_s0
d['timestamps']=timestamps
d['ra']=ra
d['dec']=dec
fs=open(output_file+f'{fname}_{recv}_level3_data','wb')
pickle.dump(d,fs,protocol=2)
fs.close()

print ('# total fitted channel number: '+str(ch_count))
print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata')
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
