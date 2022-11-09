function [TOA_tag,TOA_other]=EL_click_Detector_tags(F_ds,Y_filtered,Plot_flag,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo)

% SNR_window=SNR_window_echo;
% SNR_thresh=SNR_thresh_echo;
% consistency_T=consistency_T_echo;

TOA_tag=[]; TOA_other=[];
locs=[]; pks=[];
L_other=[]; L_tag=[];
P_other=[]; P_tag=[];
        
        Y_zoom=Y_filtered;     % Aply band pass filter               
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
        [ey,~]=energyop(Y_zoom,0);              % Apply TKEO (ey is the output signal with enhanced SNR)
        ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
        ED_thresh=0.1*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
        time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)]; 
        [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',15e-3);  % Apply instantaneous energy detector (find peaks)          

        %% Eliminate transients bellow the minimum allowed SNR

        [Array,I]=sort(pks);
       [SDCM_All, GF] = get_jenks_interface(Array);
%        figure; 
%        subplot(2,1,1); plot(Array,'rx','Linewidth',2); grid on; ylabel('Normalized peak');
%        subplot(2,1,2); plot(GF,'gx','Linewidth',2);  grid on; ylabel('GF');
        Seg_T=find(GF==max(GF))-0;
       
                 
        P_tag=pks(I(Seg_T:length(pks)));
        L_tag=locs(I(Seg_T:length(pks)));
        
 [Sorted,I2]=sort(L_tag);  
 ICIs=diff(Sorted);
 Stab_T=mean(ICIs)-2*std(ICIs);
%  figure; 
%  subplot(2,1,1); plot(Sorted,'rx','Linewidth',2); grid on; ylabel('Normalized peak');
%  subplot(2,1,2); plot(ICIs,'gx','Linewidth',2);  grid on; ylabel('ICI'); hold on;
%                  plot([1:length(ICIs)],Stab_T*ones(1,length(ICIs)),'k-','Linewidth',2);  grid on; ylabel('ICI'); hold on;
 
 Ref_std=std(ICIs);

 
 while (1)     
     trial=find(ICIs<Stab_T);
     if isempty(trial)
         break
     end
     Sorted(trial(1)+1)=[];
     Pass=I2(trial(1)+1);
     P_other=[pks(I(1:Seg_T-1))' P_tag(Pass)];
     L_other=[locs(I(1:Seg_T-1))' L_tag(Pass)];
     L_tag(Pass)=[];
     P_tag(Pass)=[];
     ICIs=[];
     ICIs=diff(Sorted);
     Stab_T=mean(ICIs)-2*std(ICIs);
%      figure; 
%      subplot(2,1,1); plot(Sorted,'rx','Linewidth',2); grid on; ylabel('Normalized peak');
%      subplot(2,1,2); plot(ICIs,'gx','Linewidth',2);  grid on; ylabel('ICI'); hold on;
%                      plot([1:length(ICIs)],Stab_T*ones(1,length(ICIs)),'k-','Linewidth',2);  grid on; ylabel('ICI'); hold on;
 end
                                     
        TOA_tag=L_tag'; TOA_other=L_other;
               
        if Plot_flag==1
            figure;set(gcf, 'Position', get(0,'Screensize'));
            subplot(3,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal');
            subplot(3,1,2); plot(time,ey_norm); hold on; 
            xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
            hold on; plot(L_tag,P_tag,'ro','Linewidth',2);           
            legend('','Tagged whale','Fontsize',10);
            subplot(3,1,3); plot(time,ey_norm); hold on;            
            hold on; plot(L_other, P_other,'go','Linewidth',2);
            legend('','Other whale','Fontsize',10);
        end
end

