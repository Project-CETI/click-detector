clear all; clc;

                                                                            
Program_folder=pwd;                                     % current folder
Data_folder=[Program_folder '\Recordings'];             % folder that contains the recordings
% Data_folder=Program_folder;

cd(Data_folder);
% files = dir('10M_ch2_0-5*.wav');
% files = dir('10M_ch2*.wav');
files = dir('AMAR*.flac');
% files = dir('SW1*.flac');
% files = dir('*filtered*.wav');

filename=files.name;
[y,Fs] = audioread(filename);                 % load recordings
cd(Program_folder)  
T_sec=10;
S_factor=floor(Fs/48e3);                      % Define factor for resampling to 48khz
Y = decimate(y,S_factor);      % Resample recording to 48khz
Fs=Fs/S_factor;
File_duration=(1/Fs)*(length(Y)-1);           % Calculate duration of the loaded recording
NOI=floor(File_duration/T_sec);               % Calculate the number of windows in the current recording
% Thresholds=linspace(0.6,1,1e3);
Thresholds=0.748;
FA=zeros(NOI*size(files,1),length(Thresholds));
c=0;

for jj=1:size(files,1)
    jj
for ii=1:NOI
    c=c+1;
    T=T_sec*Fs;
    Y_zoom=bandpass(Y(int32((ii-1)*T+1):int32((ii-1)*T+T)),[2e3 22e3],Fs);     % Aply band pass filter               
    N=45;
%     figure; plot(Y_zoom)


    FWHM = 0.4e-3;
    Ts=1/Fs;
    Sigma_G= FWHM/(4*sqrt(log(2)));
    if Fs==48e3
       N=41;
    elseif Fs==96e3
       N=82;
    elseif Fs==156250
        N=132;
    elseif Fs==192e3
        N=164;
    elseif Fs==250e3
        N=213;
    end

    i=[-N:N];
    Gaussian_Template=(Ts/(Sigma_G*sqrt(2*pi)))*exp(-((i*Ts).^2)/(2*Sigma_G^2));

            [ey,ex]=energyop(Y_zoom,0);

    MAF1=conv(Gaussian_Template,ex); MAF1(1:N-1)=[]; MAF1(end-N:end)=[];
    MAF2=movmean(ex,2*N+1);
    FDR=(MAF1-MAF2)./(MAF1);
    FDR(isnan(FDR))=0;
    FDR(FDR<0)=0;

    [Pks,Locs] =findpeaks(FDR,'MinPeakHeight',0.6,'MinPeakDistance',Fs*25e-3);  % Apply instantaneous energy detector (find peaks)          
    Decision=max(Pks);
    for q=1:length(Thresholds)
        Th=Thresholds(q);
        if Decision>Th
            FA(c,q)=FA(c,q)+1;
        end
    end
%     figure;
%     subplot(4,1,1); plot(ex);
%     subplot(4,1,2); plot(MAF1);
%     subplot(4,1,3); plot(MAF2);
%     subplot(4,1,4); plot(FDR);
%     ylim([0.5 1]);
end
end



Pd=sum(FA,1)/(NOI*size(files,1));
figure; plot(Thresholds,Pd,'-x','Linewidth',2)

FA_Gabor=FA;
Pd_Gabor=Pd;
save('Pd_Gabor.mat','Pd_Gabor')

load('Pd_Gabor.mat')
figure; plot(FA_Gabor,Pd_Gabor,'-x','Linewidth',2)




[Pks,Locs] =findpeaks(FDR,'MinPeakHeight',0.6,'MinPeakDistance',Fs*25e-3);  % Apply instantaneous energy detector (find peaks)          

%     for q=1:length(Thresholds)
%         Th=Thresholds(q);
%         if Decision>Th
%             FA(q)=FA(q)+1;
%         end
%     end

id=find(Pks>1);
Pks(id)=1;

[a,b]=xlsread('10M_ch2_0-5.xls');
  
  
  t_GT=a(:,2);
  tolerance=25e-3;
  for q=1:20
      TP=0;
      p_inds=find(Pks>Thresholds(q));
      t_capture=Locs(p_inds)/Fs;

      for i=1:length(t_GT)
         n=t_GT(i);
         [val,idx]=min(t_capture-n);
         if val<tolerance
             TP=TP+1;
         end        
      end
      Fa(q)= length(t_capture)-length(t_GT);
      Pd(q)=TP/length(t_GT);
  end

Fa(Fa<0)=0;  
figure; plot(Fa,Pd,'-.','Linewidth',2)


