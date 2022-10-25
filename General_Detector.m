function Run_Detector(SNR_window_coda,SNR_thresh_coda,T,ii,F_ds,Enhanced_sampled,Detection_flag,MP_thresh,W_seg,Dt_coda,fois_coda,wind,ICI_max_coda,ICI_min_coda,Th_coda,E_th,consistency_T_coda)

    %% Run detector over each analysis window


for ii=1:1%NOI

    if Detector_flag
        [Coda_save,Detected_codas]=Coda_click_Detector(SNR_window_coda,SNR_thresh_coda,T,ii,F_ds,Enhanced_sampled,Detection_flag,MP_thresh,W_seg,Dt_coda,fois_coda,wind,ICI_max_coda,ICI_min_coda,Th_coda,E_th,consistency_T_coda);
        Codas_info(Timer)={Coda_save};
        Gather_Codas=[Gather_Codas (Timer-1)*T_sec+Detected_codas]                
    else
        [Detected_Echos,Train_percentage,NOD,Features,Raw_Features]=EL_click_Detector_ROC2(SNR_window_echo,SNR_thresh_echo,T,ii,F_ds,Enhanced_sampled,Detection_flag,MP_thresh,W_seg,consistency_T_echo,ICI_max_echo,ICI_min_echo,Th_echo,Thresholds);                
        TMP=cell2mat(Detected_Echos);                   
        Gather_Echos={[cell2mat(Gather_Echos) (Timer-1)*T_sec+TMP]};
    end

end
        






