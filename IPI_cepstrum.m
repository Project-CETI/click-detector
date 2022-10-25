function IPI_cep=IPI_cepstrum(x,fs)

% x=Y_filtered(5.3e4:5.6e4);
% fs=F_ds;

% cepstral analysis
[C, q] = cepstrum(x, fs);
                      
Ceps=abs(C);
[pks,locs] =findpeaks(Ceps,fs,'MinPeakDistance',2.2e-3);  
inds=find(locs>2e-3);
Locs=locs(inds); Pks=pks(inds);

% figure;
%     plot(1e3*q, abs(C))
%     hold on; plot(1e3*Locs,Pks,'rx') 
      
    [P,I] = sort(Pks,'descend');
    IPI_cep=Locs(I(1));
end
