function [SW_click_inds, SW_click_pks,d_P,f_P,E_ratio]=Multipulse_locs3(Y_zoom,ey_norm,locs,pks,F_ds,W_seg,MP_thresh,Plot_flag,fois,wind)

% locs=Locs; pks=Pks;

 d_P=[];f_P=[]; MP_thresh_vals=zeros(1,length(locs));
    SW_click_inds=[]; SW_click_pks=[]; d_P1=0; d_P0=[];f_P0=[]; f_P1=[];
    Feature_ind=0; E_ratio=[]; IPI=[];
    
    
    seg_ds=round(W_seg*F_ds);
    locs_samples=locs*F_ds;
    time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];
  
  
    for ind=1:length(locs)
        Save_C=[];
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
     IPI_max=2*10e-3;      


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
                Y_roi=Y_bank(a1:a2);                                             
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
                Y_roi=Y_bank(a1:a2);
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

    if length(pk1)>3
                                
        IP_candidates=sort(tk1(end-3:end));
        S=0; C_ind=0;
         for n=1:length(IP_candidates)-1
           for m=n+1:length(IP_candidates)
             for l=m+1:length(IP_candidates)
                S=S+1;
                C(S,:)=[log((IP_candidates(l)-IP_candidates(m))/(IP_candidates(m)-IP_candidates(n))) n m l];               
                if abs(C(S,1))<0.17
                    C_ind=C_ind+1;
                    Triple_ind(C_ind,:)= [n m l];
                    Save_C(C_ind,:)=[abs(C(S,1)) C(S,2:end)];
                end
             end
           end
         end
        
        
        if ~isempty(Save_C) 
            Feature_ind=Feature_ind+1;
           
            
    %% Detection stage
             
             ROI=PSF_filt_norm;
             [pks_F,locs_F] =findpeaks(ROI,F_ds,'MinPeakHeight',MP_thresh,'MinPeakDistance',1e-3);  % Apply instantaneous energy detector (find peaks)
             if length(pks_F)>1
                 Order_pks=sort(pks_F);
                 locs_F_m=[locs_F(pks_F==Order_pks(end)) locs_F(pks_F==Order_pks(end-1))];
                 locs_F_m=sort(locs_F_m);
             end
             E_analysis=Y_bank;
             ED = movmean(E_analysis.^2,75);
             ED_norm=ED/max(ED);
             [pks_ED,~] =findpeaks(ED_norm,F_ds,'MinPeakHeight',0.02,'MinPeakDistance',2e-3);  % Apply instantaneous energy detector (find peaks)
             Search=find(pks_ED==1);   
             Relevant_pks=pks_ED(Search:end);
             Order_pks=sort(Relevant_pks);

             Detector=length(pks_F);
             if Detector>1
%                 Feature_ind=Feature_ind+1;
                if length(Relevant_pks)>1 
                    E_ratio(Feature_ind)=Order_pks(end)/Order_pks(end-1);
                else
                    E_ratio(Feature_ind)=0;
                end
                if F_ds*locs_F_m(1)-wind<1 && F_ds*locs_F_m(2)+wind<length(Y_roi)
                      P0=Y_roi(1:int32(F_ds*locs_F_m(1)+wind));
                      P1=Y_roi(int32(F_ds*locs_F_m(2)-wind):int32(F_ds*locs_F_m(2)+wind));
                elseif F_ds*locs_F_m(1)-wind<1 && F_ds*locs_F_m(2)+wind>length(Y_roi)
                      P0=Y_roi; P1=Y_roi;
                elseif F_ds*locs_F_m(1)+wind>length(Y_roi)
                      P0=Y_roi; P1=Y_roi;
                elseif F_ds*locs_F_m(2)+wind>length(Y_roi)
                      P0=Y_roi(int32(F_ds*locs_F_m(1)-wind):int32(F_ds*locs_F_m(1)+wind));
                      P1=Y_roi(int32(F_ds*locs_F_m(2)-wind):end);              
                else
                    P0=Y_roi(int32(F_ds*locs_F_m(1)-wind):int32(F_ds*locs_F_m(1)+wind));
                    P1=Y_roi(int32(F_ds*locs_F_m(2)-wind):int32(F_ds*locs_F_m(2)+wind));                   
                end
                                                           
                
                [d_P0(Feature_ind),f_P0(Feature_ind)]=superlets2(P0,F_ds,fois);
                if ~isempty(P1)
                   [d_P1(Feature_ind),f_P1(Feature_ind)]=superlets(P1,F_ds,fois);
                else
                    d_P1(Feature_ind)=0;
                    f_P1(Feature_ind)=0;
                end
                
                if d_P1(Feature_ind)==0 || f_P1(Feature_ind)==0
                    f_P(Feature_ind)=f_P0(Feature_ind);
                    d_P(Feature_ind)=d_P0(Feature_ind);               
                elseif f_P0(Feature_ind)>f_P1(Feature_ind)
                    f_P(Feature_ind)=f_P0(Feature_ind);
                    d_P(Feature_ind)=d_P0(Feature_ind);
                else
                    f_P(Feature_ind)=f_P1(Feature_ind);
                    d_P(Feature_ind)=d_P1(Feature_ind);
                end

                SW_click_inds(Feature_ind)=locs(ind);
                SW_click_pks(Feature_ind)=pks(ind);
            end
        end
    end
             end
     end
    
                    
    end
    
            if Plot_flag==1
            subplot(4,1,3);
            if ~isempty(SW_click_inds)
                plot(time',ey_norm); hold on; plot(SW_click_inds,SW_click_pks,'co','Linewidth',2);
                xlabel('time [sec]'); ylabel('TKEO'); title('Multipulse Detection'); ylim([0 1]);
            end
        end

end
