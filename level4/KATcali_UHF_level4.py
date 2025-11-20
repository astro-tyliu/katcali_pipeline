#!/usr/bin/env python3
import matplotlib
matplotlib.use('AGG')
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
from datetime import datetime
import pickle
import sys
import json
import os
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
from matplotlib.gridspec import GridSpec


current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
Tcmb=2.725

print('start @ ' + time.asctime(time.localtime(time.time())) +'#')

print(katcali.__version__)

print(plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print(plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])

fname=sys.argv[1]
ant=sys.argv[2]
input_file3_name = sys.argv[3]
file_timestamp = sys.argv[4]
print(fname,ant)

p_radec=np.loadtxt('radio_source_fsky.txt')

# ch_ref = 3200

sigma = 4.
input_file=f'/scratch3/users/liutianyang/katcali_pipeline/level3/py_results/{input_file3_name}/'
output_file=f'/scratch3/users/liutianyang/katcali_pipeline/level4/py_results/{file_timestamp}/'

def cal_map_I(map_h, map_v):    
    assert(np.shape(map_h)==np.shape(map_v))
    map_h.mask[map_h<0]=True
    map_v.mask[map_v<0]=True
    map_m=(map_h+map_v)/2.   
    
    diff_ratio=np.ma.exp(abs(np.ma.log(map_h)-np.ma.log(map_v))) #difference betwwen map h,v; 
                                                        #function to make sure maps are exchangable
    return map_m, diff_ratio

try:
    dh=pickle.load(open(input_file+fname+'_'+ant+'h_level3_data','rb'), encoding='latin-1')
    dv=pickle.load(open(input_file+fname+'_'+ant+'v_level3_data','rb'), encoding='latin-1')
    os.makedirs(output_file, exist_ok=True)

    assert((dh['ra']==dv['ra']).all()==True)
    assert((dh['dec']==dv['dec']).all()==True)
    assert((dh['timestamps']==dv['timestamps']).all()==True)
    assert((dh['Tgal_map']==dv['Tgal_map']).all()==True)
    #assert((dh['nd_s0']==dv['nd_s0']).all()==True)

    ra=dh['ra']
    dec=dh['dec']
    timestamps=dh['timestamps']
    Tgal_map=dh['Tgal_map']
    nd_s0=dh['nd_s0']

    gain_map_h=dh['gain_map']
    T_map_h=dh['T_map']
    Tresi_map_h=dh['Tresi_map']
    Tsm_map_h=dh['Tsm_map']
    nd_s0_h=dh['nd_s0'] 
    Tel_map_h=dh['Tel_map']
    Tsky_map_h=T_map_h-Tel_map_h-Tsm_map_h

    gain_map_v=dv['gain_map']
    T_map_v=dv['T_map']
    Tresi_map_v=dv['Tresi_map']
    Tsm_map_v=dv['Tsm_map']
    nd_s0_v=dv['nd_s0'] 
    Tel_map_v=dv['Tel_map']
    Tsky_map_v=T_map_v-Tel_map_v-Tsm_map_v

    d_checker={}
    d_checker['Tsm_map_h']=Tsm_map_h
    d_checker['Tsm_map_v']=Tsm_map_v
    d_checker['Tgal_map']=Tgal_map
    d_checker['Tcmb']=Tcmb

    freqs=kio.cal_freqs_UHF(range(4096))
    freqs = np.array(freqs)

    #freq ranges are from Marta[Martha]
    ch_Vod1 = kio.ch_finder(768,freqs)
    ch_Vod2 = kio.ch_finder(778,freqs)
    ch_MTN1 = kio.ch_finder(801,freqs)
    ch_MTN2 = kio.ch_finder(811,freqs)
    ch_Telk1 = kio.ch_finder(811,freqs)
    ch_Telk2 = kio.ch_finder(821,freqs)
    ch_tower=list(range(ch_Vod1,ch_Vod2))+list(range(ch_MTN1,ch_MTN2))+list(range(ch_Telk1,ch_Telk2))
    ch_tower.sort()

    ch_e=3547  #UHF_band
    
    assert((T_map_h.mask==Tresi_map_h.mask).all()==True)
    assert((T_map_v.mask==Tresi_map_v.mask).all()==True)
    div=T_map_h/T_map_v
    mask_clean=div.mask.copy()
    mask_clean_backup=mask_clean.copy()

    l=[T_map_h,T_map_v,Tresi_map_h,Tresi_map_v,Tsm_map_h,Tsm_map_v,gain_map_h,gain_map_v]
    l_copy=[T_map_h.copy(),T_map_v.copy(),Tresi_map_h.copy(),Tresi_map_v.copy(),Tsm_map_h.copy(),Tsm_map_v.copy(),gain_map_h.copy(),gain_map_v.copy()]

    l_str=[r'$T_{cal}(\nu)$'+', pol=HH',r'$T_{cal}(\nu)$'+', pol=VV',r'$T_{res}(\nu)$'+', pol=HH',r'$T_{res}(\nu)$'+', pol=VV',
           r'$T_{rec}(\nu)$'+', pol=HH',r'$T_{rec}(\nu)$'+', pol=VV','gain'+r'$(\nu)$'+', pol=HH','gain'+r'$(\nu)$'+', pol=VV']

    del_point0,del_point1,del_point2,del_point3,del_point4,del_point5,del_point6,del_point7=[],[],[],[],[],[],[],[]

    #####filter######
    niter=6
    
    plt.figure(figsize=(14,len(l)*niter*1.6))
    plt.subplots_adjust (wspace=0.2, hspace=0.35) 
    for i in range(niter):
        print ('> The iteration '+str(i+1)+' is in progress...')
        mask_clean2=mask_clean.copy()
        
        del_point=[]
        
        for l_i in range(len(l)):
            a=l[l_i]
            plt.subplot(len(l)/2*niter,2,i*len(l)+l_i+1)
            a_ch=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0)
            plt.plot(freqs/1e6,a_ch,'g.', ms=4)
            if l_i !=2 and l_i!=3:
                plt.ylim(a_ch.min()*0.95,a_ch.max()*1.05)
            
            try:
                ax1=kf.curve_Rbf_filter_ma_array(range(272,2869),a_ch[272:2869],sigma=sigma)
                #ax1=kf.curve_filter_ma_array(range(272,2869),a_ch[272:2869],sigma=sigma)
                plt.plot(freqs[ax1]/1e6,a_ch[ax1],'.',color='r',ms=10) #deleted data
                mask_clean2[:,ax1]=True ###apply to the mask
                del_point+=list(ax1)
            except:
                ax1=[]
                print ('* 272-2869 filter not applied')
            #'''
            try:
                ax2=kf.curve_Rbf_filter_ma_array(range(3133,3547),a_ch[3133:3547],sigma=sigma)
                #ax2=kf.curve_filter_ma_array(range(3133,3547),a_ch[3133:3547],sigma=sigma)
                plt.plot(freqs[ax2]/1e6,a_ch[ax2],'.',color='r',ms=10) #deleted data
                mask_clean2[:,ax2]=True ###apply to the mask
                del_point+=list(ax2)
            except:
                ax2=[]
                print ('* 3133-3547 filter not applied')
            
            try:    
                ax3=kf.curve_Rbf_filter_ma_array(range(272,3547),a_ch[272:3547],sigma=sigma)
                #ax3=kf.curve_filter_ma_array(range(272,3547),a_ch[272:3547],sigma=sigma,k=1)
                plt.plot(freqs[ax3]/1e6,a_ch[ax3],'.',color='r',ms=10) #deleted data
                mask_clean2[:,ax3]=True ###apply to the mask
                del_point+=list(ax3)
            except:
                ax3=[]
                print ('* 272-3547 filter not applied')
            #'''
            if l_i==0:
                del_point0+=list(ax1)+list(ax2)+list(ax3)
            if l_i==1:
                del_point1+=list(ax1)+list(ax2)+list(ax3)
            if l_i==2:
                del_point2+=list(ax1)+list(ax2)+list(ax3)
            if l_i==3:
                del_point3+=list(ax1)+list(ax2)+list(ax3)
            if l_i==4:
                del_point4+=list(ax1)+list(ax2)+list(ax3)
            if l_i==5:
                del_point5+=list(ax1)+list(ax2)+list(ax3)
            if l_i==6:
                del_point6+=list(ax1)+list(ax2)+list(ax3)
            if l_i==7:
                del_point7+=list(ax1)+list(ax2)+list(ax3)
        
            
            plt.title('time mean '+ l_str[l_i]+', iter= '+str(i+1))
            
            if (l_i+2) % 2 == 0 and (l_i+1) % 7 !=0:
                plt.ylabel('Temperature (K)')
                
            if (l_i+1) % 7 ==0:
                plt.ylabel('gain')    
                
                
            if (l_i+2) % 8 == 0 or (l_i+1) % 8 == 0:
                plt.xlabel('Frequency (MHz)')
                
        print ('masks before and after filter are same:',)
        print ((mask_clean==mask_clean2).all())
        mask_clean=mask_clean2.copy()
        
        print (str(len(set(del_point))) + ' channels masked by '+str(len(del_point))+' times')
        if len(del_point)==0:
            print ('***filter done after iter time '+str(i+1)) 
            break
    true_iter = i + 1
    len_del = len(del_point)

    print (del_point0)
    print (del_point1)
    print (del_point2)
    print (del_point3)
    print (del_point4)
    print (del_point5)
    print (del_point6)
    print (del_point7)

    del_point_all=del_point0+del_point1+del_point2+del_point3+del_point4+del_point5+del_point6+del_point7
    del_point_all=list(set(del_point_all))
    print (len(del_point_all))
    
    plt.figure(figsize=(14,80))
    plt.subplots_adjust (wspace=0.2, hspace=0.3) 
    
    for l_i in range(len(l)):
        a=l[l_i]
        plt.subplot(len(l)/2*niter,2,i*len(l)+l_i+1)
        
        if l_i==0:
            del_point_plot=del_point0
        if l_i==1:
            del_point_plot=del_point1
        if l_i==2:
            del_point_plot=del_point2
        if l_i==3:
            del_point_plot=del_point3
        if l_i==4:
            del_point_plot=del_point4
        if l_i==5:
            del_point_plot=del_point5
        if l_i==6:
            del_point_plot=del_point6
        if l_i==7:
            del_point_plot=del_point7
            
        a_ch=np.ma.mean(a,axis=0) #raw
        plt.plot(freqs[del_point_all]/1e6,a_ch[del_point_all],'.',color='gold',ms=6) #all deleted data
        plt.plot(freqs[del_point_plot]/1e6,a_ch[del_point_plot],'.',color='r',ms=6) #deleted data
        
        a_ch1=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0) #clean 
        plt.plot(freqs/1e6,a_ch1,'.', ms=4) #clean data
        
        if l_i !=2 and l_i!=3:
            plt.ylim(a_ch.min()*0.95,a_ch.max()*1.05)
     
        plt.title('time mean '+ l_str[l_i])
    
        if (l_i+2) % 2 == 0 and (l_i+1) % 7 !=0:
            plt.ylabel('Temperature (K)')
    
        if (l_i+1) % 7 ==0:
            plt.ylabel('gain')    
    
        if (l_i+2) % 8 == 0 or (l_i+1) % 8 == 0:
            plt.xlabel('Frequency (MHz)')
            
        if l_i==1:
            plt.legend(['deleted data due to other parameters', 'deleted data due to this parameter','clean data'], fontsize=11.5)
    plt.savefig(output_file+'F_'+fname+'_'+ant+'_ch_filter_final.png', bbox_inches='tight')
    #plt.show()

    mask_clean[:,ch_tower]=True

    plt.figure(figsize=(14,80))
    plt.subplots_adjust (wspace=0.2, hspace=0.3) 
    
    for l_i in range(len(l)):
        a=l[l_i]
        plt.subplot(len(l)/2*niter,2,i*len(l)+l_i+1)
              
        a_ch1=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0) #clean 
        plt.plot(freqs/1e6,a_ch1,'.', ms=4) #clean data
    
       
        if l_i !=2 and l_i!=3:
            plt.ylim(a_ch1.min()*0.95,a_ch1.max()*1.05)
     
        plt.title('time mean '+ l_str[l_i])
    
        if (l_i+2) % 2 == 0 and (l_i+1) % 7 !=0:
            plt.ylabel('Temperature (K)')
    
        if (l_i+1) % 7 ==0:
            plt.ylabel('gain')    
    
        if (l_i+2) % 8 == 0 or (l_i+1) % 8 == 0:
            plt.xlabel('Frequency (MHz)')
            
        if l_i==1:
            plt.legend(['clean data'], fontsize=11.5)
    plt.savefig(output_file+'F_'+fname+'_'+ant+'_ch_filter_final_clean.png', bbox_inches='tight')
    # plt.show()

    def mask2int(mask):
        int_mask=np.zeros(np.shape(mask))
        #print (int_mask)
        int_mask[np.where(mask==True)]=1
        int_mask[np.where(mask==False)]=0
        return int_mask

    Tsky_map_h=np.ma.array(Tsky_map_h,mask=mask_clean)
    Tsky_map_v=np.ma.array(Tsky_map_v,mask=mask_clean)

    Tresi_map_h=np.ma.array(Tresi_map_h,mask=mask_clean)
    Tresi_map_v=np.ma.array(Tresi_map_v,mask=mask_clean)

    Tsky_map,Tsky_ratio=cal_map_I(Tsky_map_h, Tsky_map_v)
    Tresi_map=(Tresi_map_h+Tresi_map_v)/2.
    assert((Tsky_map.mask==Tresi_map.mask).all()==True)

    plt.figure(figsize=(14,5.4))
    #plt.subplots_adjust (wspace=0.2, hspace=0) 
    plt.subplot(121)
    plt.imshow(Tsky_map[:,272:3547],aspect='auto',extent=(freqs[272]/1e6,freqs[3547]/1e6,len(timestamps)*2,0))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('time (s)')
    #plt.xlim(550,ch_e)
    clb = plt.colorbar()
    clb.set_label('Kelvin', labelpad=-35, y=1.06, rotation=0, fontsize=12)
    plt.twiny()
    plt.imshow(Tsky_map[:,272:3547], aspect='auto',extent=(272,3547,len(timestamps)*2,0))
    plt.xlabel('channel')
    plt.subplot(122)
    plt.imshow(Tresi_map[:,272:3547],aspect='auto',vmax=0.3, extent=(freqs[272]/1e6,freqs[3547]/1e6,len(timestamps)*2,0))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('time (s)')
    #plt.xlim(550,ch_e)
    clb = plt.colorbar()
    clb.set_label('Kelvin', labelpad=-45, y=1.06, rotation=0, fontsize=12)
    plt.twiny()
    plt.imshow(Tresi_map[:,272:3547], aspect='auto',vmax=0.3, extent=(272,3547,len(timestamps)*2,0))
    plt.xlabel('channel')
    plt.savefig(output_file+'F_'+fname+'_'+ant+'_Inten.png', bbox_inches='tight')
    # plt.show()

    assert((Tsky_ratio.mask==mask_clean).all()==True)
    print (np.ma.max(Tsky_ratio))

    assert((Tsky_map.mask==mask_clean).all()==True)
    assert((Tresi_map.mask==mask_clean).all()==True)
    assert((Tsky_ratio.mask==mask_clean).all()==True)

    # data=kio.load_data(fname)
    # j = 0
    # for i in range(64):
    #     ant = 'm' + '%03d' %i
    #     for pol in ['h', 'v']:
    #         recv = ant + pol
    #         try:
    #             vis,flags= kio.call_vis(fname,recv)
    #             target_list,c0_list,bad_ants,flux_model_list=kio.check_ants(fname)
    #             if isinstance(target_list, list):
    #                 ang_deg=kio.load_ang_deg2(ra,dec,c0_list)
    #             else:
    #                 ang_deg=kio.load_ang_deg(ra,dec,c0_list)
    #             dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
    #             j = 1
    #             break
    #         except Exception as e:
    #             continue
    #     if j:
    #         break
    # del data, vis

    # ra1=ra.copy()
    # if np.max(ra1[dp_ss]) - np.min(ra1[dp_ss]) > 300:
    #     ra1[ra>180]-=360

    d={}
    d['Tsky_map']=Tsky_map
    d['Tresi_map']=Tresi_map
    d['Tsky_ratio']=Tsky_ratio

    d['timestamps']=timestamps
    d['ra']=ra
    d['dec']=dec
    #assert((nd_s0_h==nd_s0_v).all())
    d['nd_s0_h']=nd_s0_h
    d['nd_s0_v']=nd_s0_v
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()

    d1={}

    d1['Inten_mask']=mask_clean
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_mask','wb')
    pickle.dump(d1,fs,protocol=2)
    fs.close()

    ch_plot=3300
    plot_gsize=90
    ra1=ra.copy()
    if np.max(ra1[nd_s0]) - np.min(ra1[nd_s0]) > 300:
        ra1[ra>180]-=360

    plt.figure(figsize=(14,5))
    ax=plt.subplot()
    p_data=Tresi_map[nd_s0,ch_plot]
    if (p_data.mask==True).all()==False:
        kv.plot_mdata(ra1[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize, grid_method='linear', levels=6, x_mask=1, y_mask=2, cmap=kv.cmap2())
    ax.invert_xaxis()
    plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    #plt.text(146.3, 8.5, 'Kelvin', rotation=0)
    plt.xlabel('R.A. (J2000)')
    plt.ylabel('Dec (J2000)')
    #plt.savefig('F_calibrated_map.pdf', bbox_inches='tight')
    plt.title(fname+', '+ant+', Tres-map, ch'+str(ch_plot))
    plt.savefig(output_file+'F_'+fname+'_'+ant+'_Tres.png', bbox_inches='tight')
    #plt.show()

    metadata = {
        "parameters": {
            "fname": fname,
            "ant": ant,
            "level3_input_file": input_file3_name,
            "file_timestamp": file_timestamp
        },
        "ch_e": ch_e,
        "ch_plot": ch_plot,
        "sigma": sigma,
        "max_iterations": niter,
        "true_iterations": true_iter,
        "final_num_del": len_del,
        "run_time": current_time,
        "description": '### '+ant+' of '+fname+' finished successfully ###'
    }
    metadata_path = os.path.join(output_file, f"{fname}_{ant}_metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

except IOError:
    if os.path.exists(output_file):
        metadata = {
            "parameters": {
                "fname": fname,
                "ant": ant,
                "level3_input_file": input_file3_name,
                "file_timestamp": file_timestamp
            },
            "run_time": current_time,
            "description": '### failed for '+fname+', '+ant
        }
        metadata_path = os.path.join(output_file, f"{fname}_{ant}_metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=4)
    print ('### failed for '+fname+', '+ant)

print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata')
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
