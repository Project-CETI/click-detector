function  [Coda_save,Detected_pattern_t]=Coda_click_Detector(SNR_window,SNR_thresh,F_ds,Y_filtered,Plot_flag,MP_thresh,W_seg,Dt_coda,fois,wind,ICI_max_coda,ICI_min_coda,Th_coda,E_th,consistency_T_coda)

% fois=fois_coda; 
% SNR_window=SNR_window_coda;
% SNR_thresh=10;
% consistency_T=consistency_T_coda;

            Coda_save=[]; Pks=[]; Locs=[]; Detected_pattern_t=[];
            c=0;                   % counter for the transient detection stage

            Y_zoom=Y_filtered;                   

            t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
            [ey,ex]=energyop(Y_zoom,0);             % Apply TKEO (ey is the output signal with enhanced SNR)
            ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
            ED_thresh=10*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
            Th_min=0.03;
            time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];   
            [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',Th_min);  % Apply instantaneous energy detector (find peaks)          
            locs_samples=locs*F_ds;                 % store indices of time readings of the detected impulses

            %% Eliminate transients bellow the minimum allowed SNR

            for i=1:length(locs)
                if locs_samples(i)>SNR_window && (locs_samples(i)+SNR_window)<length(ey_norm)
                    tmp=ey_norm(int32(locs_samples(i)-SNR_window):int32(locs_samples(i)+SNR_window));
                    SNR(i)=pks(i)/median(tmp); 
                    if SNR(i)>SNR_thresh
                        c=c+1;
                        Locs(c)=locs(i);
                        Pks(c)=pks(i);
                    end
                end
            end
            
%             if length(Locs)>2 && length(Locs)<100
            if length(Locs)>2 

               %% Show plots for the transient detection stage

                if Plot_flag==1
                    figure;set(gcf, 'Position', get(0,'Screensize'));
                    subplot(4,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal');
                    subplot(4,1,2); plot(time,ey_norm); xlabel('time [sec]'); ylabel('TKEO'); title('Transient');
                    hold on; plot(Locs,Pks,'x')    
                end
            
                %% Multipulse Detection Stage
                [MP_t,MP_p,d_P_coda,f_P_coda,E_ratio]=Multipulse_locs3(Y_zoom,ey_norm,Locs,Pks,F_ds,W_seg,MP_thresh,Plot_flag,fois,wind);  % Run multipulse detector

               %% Coda Detection Stage

               [Detected_pattern_t,Detected_pattern_p]=ICI4(time,ey_norm,MP_t,MP_p,Y_zoom,F_ds,f_P_coda,d_P_coda,E_ratio,Dt_coda,ICI_max_coda,ICI_min_coda,Th_coda,E_th,consistency_T_coda); % Run Fuzzy-Logic based coda detector
               if ~isempty(Detected_pattern_t)
                   Eliminate_inds=Eliminate_MultiPath2(Detected_pattern_t,Detected_pattern_p);
                    Eliminate_inds=Eliminate_inds-ones(1,length(Eliminate_inds));
                    Detected_pattern_t(Eliminate_inds)=[];
                    Detected_pattern_p(Eliminate_inds)=[];
               end
               if ~isempty(Detected_pattern_t)
                   [Coda_times,Coda_Type]=Classify_Coda_type(Detected_pattern_t,Detected_pattern_p,time,ey_norm,Plot_flag);  % Apply Coda type classifier       
                   Coda_save=[Coda_times ; Coda_Type];  
               end
            end
           
end