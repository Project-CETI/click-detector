function [TOA,TOA_tag,TOA_other]=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag,Tag_flag)
%function TOA=Run_Detector(Detector_flag,F_ds,Enhanced_sampled,Detection_flag)
%
%Description:
%This function gets a signal buffer of M secons and outputs a vector
%contains the time of arrivals of detected sperm whale Codas/Echolocation clicks, 
%
%Input:
%Y_filtered - Buffer (of Tsec [seconds]) of the filtered signal
%F_ds - Sample frequency (after downsampling)
%Detector_flag - a flag {0,1}. If 1 - apply coda detector, If 0 - apply echolocation clicks detector
%Plot_flag - a flag {0,1}. If 1 - Show figure of detections, If 0 - disable figures
%
%Output:
% TOA - Vector of MX1 with time of arrival in seconds for M detections


%% Run detector over each analysis window
    
      %% Determine operational parameters (Local parameters(i.e. parameters for each detector))
        
        Dt_coda=0.4;                                  % set threshold of fuzzy-logic coda detector             
        MP_thresh=0.5;                                % set threshold for the multipulse detector (for peaks of the slope phase function)
        W_seg=28e-3;                                  % [sec] -set window length in for multipulse analysis

        %% Coda detector parameters:
        ICI_max_coda=1;                             % set maximal allowed ICI
        ICI_min_coda=0.05;                          % set minimum allowed ICI 
        consistency_T_coda=0.2;                     % set maximal allowd consistency parameter        
        SNR_window_coda=F_ds*50e-3;                 % Define time window for SNR calculation of each detected transient
        SNR_thresh_coda=10;                         % Define minimum allowed SNR
        fois_coda= linspace(2e3,8e3,100);           % set spectogram bounds
        Th_coda=0.8;                                % set threshold for maximal allowed diversity in click's amplitude
        E_th=0.7;                                   % set threshold for maximal allowed diversity in click's energy
        wind=F_ds*0.3e-3;                           % set window for analysis of the multipulsed components of a click
        

        %% Echolocation detector parameters:
        
        ICI_max_echo=1.8;                           % set maximal allowed ICI
        ICI_min_echo=0.42;                          % set minimum allowed ICI
        consistency_T_echo=0.22;                    % set maximal allowd consistency parameter        
        SNR_window_echo=F_ds*50e-3;                 % Define time window for SNR calculation of each detected transient
        SNR_thresh_echo=100;                        % Define minimum allowed SNR
        fois_echo= linspace(2e3,18e3,100);          % set spectogram bounds
        Th_echo=0.7;                                % set threshold for maximal allowed diversity in click's amplitude
        Gather_TOA=[];
        TOA=[];
        TOA_tag=[];
        TOA_other=[];
        
       %% Select and run detector
       
        if Detector_flag
            [Coda_save,TOA]=Coda_click_Detector(SNR_window_coda,SNR_thresh_coda,F_ds,Y_filtered,Plot_flag,MP_thresh,W_seg,Dt_coda,fois_coda,wind,ICI_max_coda,ICI_min_coda,Th_coda,E_th,consistency_T_coda);
            Codas_info={Coda_save};       
        else
            if Tag_flag
            [TOA_tag,TOA_other]=EL_click_Detector_tags(F_ds,Y_filtered,Plot_flag,consistency_T_echo,ICI_max_echo,ICI_min_echo,Th_echo,MP_thresh,W_seg);
            else
            detections=EL_click_Detector_server(SNR_window_echo,SNR_thresh_echo,Fs,F_ds,Y_filtered,Plot_flag,MP_thresh,W_seg,consistency_T_echo,ICI_max_echo,ICI_min_echo,Th_echo); 
            TOA=cell2mat(detections);
            end                           
        end

    
end

