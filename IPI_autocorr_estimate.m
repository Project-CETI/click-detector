function IPI_auto=IPI_autocorr_estimate(Sig,Fs)
    [autocor,lags] = xcorr(Sig,Fs,'coeff');
    tt=1e3*lags/Fs;
%     figure; plot(tt,autocor);
%     xlabel('Lag (ms)')
%     ylabel('Autocorrelation')
    [Ps,Ls] =findpeaks(autocor,Fs,'MinPeakHeight',1e-3,'MinPeakDistance',2e-3);  
    [~,I] = sort(Ps,'descend');
    Ls_sorted=Ls(I);
    IPI_auto=Ls_sorted(1)-Ls_sorted(2);
end






