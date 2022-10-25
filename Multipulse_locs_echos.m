function [SW_click_inds, SW_click_pks,IPI]=Multipulse_locs_echos(Y_zoom,ey_norm,locs,pks,Fs,F_ds,W_seg,MP_thresh,Plot_flag)

% locs=Locs; pks=Pks;
IPI_auto=[]; IPI_cep=[];

    d_P=[];f_P=[]; MP_thresh_vals=zeros(1,length(locs));
    SW_click_inds=[]; SW_click_pks=[]; d_P1=0; d_P0=[];f_P0=[]; f_P1=[];
    Feature_ind=0; E_ratio=[]; IPI=[];
    
    seg_ds=round(W_seg*F_ds); seg_ds2=round(W_seg*Fs);
    locs_samples=locs*F_ds; locs_samples2=locs*Fs;
    time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];
    
    for ind=1:length(locs)

        Second_pulse=[];
        if  (locs_samples(ind)+seg_ds)<length(ey_norm)
           Clicks_bank= ey_norm(int32(locs_samples(ind)-round(seg_ds)):int32(locs_samples(ind)+seg_ds));   % pick region of analysis
           Y_bank= Y_zoom(int32(locs_samples(ind)-round(seg_ds)):int32(locs_samples(ind)+seg_ds)); 
        else
           Clicks_bank= ey_norm(int32(locs_samples(ind)-round(seg_ds)):int32(length(ey_norm)));   % pick region of analysis
           Y_bank= Y_zoom(int32(locs_samples(ind)-round(seg_ds)):int32(length(ey_norm)));   % pick region of analysis
        end  
        

        
        Clicks_bank=Clicks_bank/max(Clicks_bank);
        Clicks_bank(1:15)=0; Clicks_bank(end-15:end)=0;
        
        exit=0; limit=0.95; go=1;
        while exit==0 && go<10
            go=go+0.1;
            [pks_zoom,locs_zoom] =findpeaks(Clicks_bank,F_ds,'MinPeakHeight',1e-3,'MinPeakDistance',0.2e-3*go*(1+rand(1))+1e-3);  % Apply instantaneous energy detector (find peaks)
             Main_peak=find(pks_zoom==1);
             if length(Main_peak)==1
                 exit=1;
             end
             tmp = rand;
             if tmp > limit
                break
             end
        end

%         figure; plot(time_C,Clicks_bank); hold on; plot(locs_zoom,pks_zoom,'x');
     IPI_max=18e-3;      


     if ~isempty(Main_peak)       
            Lim_left=locs_zoom(Main_peak)-IPI_max;
            Lim_right=locs_zoom(Main_peak)+IPI_max;              
            pks_zoom(locs_zoom<Lim_left)=0;
            pks_zoom(locs_zoom>Lim_right)=0;
            Lz=length(nonzeros(pks_zoom));
            if Lz>1
                [pk1,I1]=sort(pks_zoom);
                tk1=locs_zoom(I1);
                if Lz>2
                    RP=(pk1(end-1)-pk1(end-2))/pk1(end-1);
                    if abs(tk1(end-1)-tk1(end))>abs(tk1(end-2)-tk1(end)) 
                        Second_pulse=find(pks_zoom==pk1(end-1));                   
                    elseif RP<0.3 
                        Second_pulse=find(pks_zoom==pk1(end-2));
                    else
                        Second_pulse=find(pks_zoom==pk1(end-1));
                    end
                else                   
                   Second_pulse=find(pks_zoom==pk1(end-1));
                end
                W=abs(locs_zoom(Main_peak)-locs_zoom(Second_pulse));
                if W>1.8e-3 && W<9e-3
                    W_size=W;
                else
                W_size=1.8e-3;        %set min window size
                end
            else
                locs_zoom=locs_zoom(find(pks_zoom~=0));
                pks_zoom=nonzeros(pks_zoom);
                W_size=1.8e-3;  
            end

            if Lz>1
                Make_order=sort([Main_peak Second_pulse]);
                First_peak=Make_order(1);
                Second_peak=Make_order(2);

                Mask1=zeros(1,length(Clicks_bank));  Mask2=zeros(1,length(Clicks_bank));
                Mask1(int32(locs_zoom(First_peak)*F_ds-10):int32(locs_zoom(First_peak)*F_ds+10))=1;
                Mask2(int32(locs_zoom(Second_peak)*F_ds-10):int32(locs_zoom(Second_peak)*F_ds+10))=1;
                Mask_signal=Mask1+Mask2;
                Mask_noise=ones(1,length(Mask_signal))-Mask_signal;
                UnMask_signal=Mask_signal.*Clicks_bank';            
                UnMask_noise=Mask_noise.*Clicks_bank';
                UnMask_noise_filt=medfilt1(UnMask_noise,25);

                Clicks_bank=UnMask_noise_filt+UnMask_signal;            
                window=round(F_ds*(W_size-0.3e-3));
                
                if locs_zoom(First_peak)*F_ds-0.2*window>0 
                    a1=int32(ceil(locs_zoom(First_peak)*F_ds-0.2*window));
                else
                    a1=int32(1);
                end
                if locs_zoom(Second_peak)*F_ds+1.2*window<length(Clicks_bank)
                    a2=int32(locs_zoom(Second_peak)*F_ds+1.2*window);
                else
                    a2=int32(length(Clicks_bank));
                end
                xxx=Clicks_bank(a1:a2); 
                                              
            else
                window=round(F_ds*(W_size-0.3e-3));
                if locs_zoom*F_ds-0.5*window>0 
                    a1=int32(locs_zoom*F_ds-0.5*window);
                else
                    a1=int32(1);
                end
                if locs_zoom*F_ds+2*window<length(Clicks_bank)
                    a2=int32(locs_zoom*F_ds+2*window);
                else
                    a2=int32(length(Clicks_bank));
                end
                                    
                xxx=Clicks_bank(a1:a2);                          
            end
            
            [ey_zoom,~]=energyop(xxx,0);
             check=length(xxx)-window-2;
             if check>100
                 for n=1:(length(xxx)-window-2)
                    xn=ey_zoom(int32(n):int32(window+n));
                    [Tau_w, ~]=grpdelay(xn,F_ds); 
                    Tau_filt= medfilt1(Tau_w,3);
                    PSF(n)=-mean(Tau_filt);    
                 end

                 PSF_filt=lowpass(PSF,20,F_ds);
                 M=mean(PSF_filt);
                 PSF_filt2=PSF_filt+abs(M)*ones(1,length(M));
                 PSF_filt_norm=PSF_filt2/max(PSF_filt2);

        %% Detection stage

                 ROI=PSF_filt_norm;
                 [pks_F,locs_F] =findpeaks(ROI,F_ds,'MinPeakHeight',MP_thresh,'MinPeakDistance',1e-3);  % Apply instantaneous energy detector (find peaks)
                 [PK_order,I_PK]=sort(pks_F);
                 locs_F=locs_F(I_PK);
                 if length(PK_order)>1
                    MP_thresh_vals(ind)=PK_order(end-1);
                 end
                 
                 
%              time_x=[0:1/F_ds:(1/F_ds)*(length(xxx)-1)];
%              time_psf=[0:1/F_ds:(1/F_ds)*(length(PSF_filt_norm)-1)];
% %              figure; plot(time_x,xxx); hold on; 
%              figure;
%              Y_seg=Y_bank(a1:a2)/max(Y_bank(a1:a2)); t_seg=[0:1/F_ds:(1/F_ds)*(length(Y_seg)-1)];
%              plot(1e3*t_seg,Y_seg);hold on;
%              plot(1e3*time_psf,PSF_filt_norm); hold on;
% %              plot(locs_F,pks_F,'x');
% %                  


                 Detector=length(pks_F);
                 if Detector>1                
                    Feature_ind=Feature_ind+1;
                    if ~isempty(Second_pulse)
                         IPI(Feature_ind)=abs(locs_zoom(Main_peak)-locs_zoom(Second_pulse));
                    else
                         IPI(Feature_ind)=abs(locs_F(end)-locs_F(end-1));
                    end

                    SW_click_inds(Feature_ind)=locs(ind);
                    SW_click_pks(Feature_ind)=pks(ind);
                 end
             end
       end
    end
    
    
        if Plot_flag==1
            subplot(4,1,3);
            if ~isempty(SW_click_inds)
                plot(time,ey_norm); hold on; plot(SW_click_inds,SW_click_pks,'co','Linewidth',2);
                xlabel('time [sec]'); ylabel('TKEO'); title('Multipulse Detection'); ylim([0 1]);
            end
        end
               
    
end
