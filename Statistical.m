

clear all; clc;

                                                                            
Program_folder=pwd;                                     % current folder
%Data_folder=[Program_folder '\Recordings'];             % folder that contains the recordings
Data_folder=[Program_folder];
cd(Data_folder);
% files={'10M_ch2_0-5*.wav','SW1*.flac','AMAR2*.flac'};
% files = dir('10M_ch2_0-5*.wav');
files = dir('AMAR*.flac');
% files = dir('*filtered*.wav');


DR=[]; SIGMA=[];
for q=1:length(files)
    q
cd(Data_folder);
% SF=dir([cell2mat(files(q))]);
filename=files(q).name;
[y,Fs] = audioread(filename);                 % load recordings
cd(Program_folder)  

 S_factor=floor(Fs/48e3);                 % Define factor for resampling to 48khz  
 Y_decimated = decimate(y,S_factor);
 y=[]; Fs=Fs/S_factor;
 y=Y_decimated;
Tsec=2;
T=Tsec*Fs;
     File_duration=(1/Fs)*(length(y)-1);           % Calculate duration of the loaded recording
     NOI=floor(File_duration/Tsec);               % Calculate the number of windows in the current recording
dr=[]; Sigma=[];

for ii=1:NOI
    
    Y_zoom=y(int32((ii-1)*T+1):int32((ii-1)*T+T));     % Aply band pass filter               
%     figure; plot(Y_zoom)

    [s,f,t]=spectrogram(Y_zoom,hann(128),64,128,Fs);
    % figure;
    % imagesc(t,f,abs(s).^2);

    F_i=find(f<3.375e3 | f>16.875e3);
    s(F_i,:)=[];

    E_spec=abs(s).^2;
    Win=24;
    Yf = movmedian(E_spec',2*Win+1);
    %  figure;
    % imagesc(Yf')
    Ymu=Yf';

    Yn=(E_spec-Ymu)./Ymu;
    % figure;
    % imagesc(Yn)

    E_f=sum(Yn);
    E_f_sorted=sort(E_f,'descend');
    % figure; subplot(2,1,1); plot(E_f,'.'); subplot(2,1,2); plot(E_f_sorted,'.');

    xmax=mean(E_f_sorted(1:30));
    xref=median(E_f_sorted);
    dr(ii)=20*log10(xmax/xref);
    Sigma(ii)=std(E_f);

end

DR=[DR dr];
SIGMA=[SIGMA Sigma];
% hold on; plot(Sigma,dr,'x','Linewidth',2); 
% xlabel('\sigma','Fontsize',14); ylabel('dr','Fontsize',14); grid on;
% ylim([0 100]);

end

dr_noise=DR;
sigma_noise=SIGMA;
% dr_SW=DR;
% sigma_SW=SIGMA;


figure;
plot(sigma_noise,dr_noise,'bx','Linewidth',2); hold on;
plot(sigma_SW,dr_SW,'rx','Linewidth',2);  grid on;
xlabel('\sigma','Fontsize',18,'FontWeight','bold'); ylabel('dr','Fontsize',18,'FontWeight','bold');


 figure;
 subplot(1,2,1); 
 plot(F_phase_noise2(:,1),F_phase_noise2(:,2),'bx','Linewidth',1); hold on; 
%  plot(Gather_Features_phase(:,1),Gather_Features_phase(:,2),'rx','Linewidth',1);
 plot(F_phase_SW1(:,1),F_phase_SW1(:,2),'rx','Linewidth',1);
 axis([0 2 0 10]);grid on;
 xlabel('m','Fontsize',14);
 ylabel('\sigma_{IPI}','Fontsize',14);
 legend('Noise transients','Sperm whales clicks','Fontsize',12);
subplot(1,2,2); 
plot(sigma_noise,dr_noise,'bx'); hold on;
plot(sigma_SW,dr_SW,'rx');
xlabel('\sigma','Fontsize',14); ylabel('dr','Fontsize',14); grid on;
legend('Noise recordings','Sperm whales recordings','FontSize', 12)

C1=[sigma_SW'  dr_SW']; L1=size(C1,1);
C2=[sigma_noise'  dr_noise']; L2=size(C2,1);

waveform=[C1 ; C2];
Group3=[ones(1,L1) zeros(1,L2)];

Criterions={'ttest','entropy','bhattacharyya','roc','wilcoxon'};
for i=1:5
    [~, Zp] = rankfeatures(phase', Group,'Criterion',cell2mat(Criterions(i)),'CrossNorm', 'meanvar');
%     [~, Zs] = rankfeatures(waveform', Group2,'Criterion',cell2mat(Criterions(i)),'CrossNorm', 'meanvar');
    [~, Zpr] = rankfeatures(phase_raw', Group3,'Criterion',cell2mat(Criterions(i)),'CrossNorm', 'meanvar');        
    SI_p(i)=mean(Zp);
%     SI_s(i)=mean(Zs);
    SI_pr(i)=mean(Zpr);
end

inds_n=find(w1_ceps(:,1)<0.2 & w1_ceps(:,2)<0.9);
inds_sw=find(w2_ceps(:,1)<0.2 & w2_ceps(:,2)<0.9);

W1_phase=w1_phase(inds_n,:);
W2_phase=w2_phase(inds_sw,:);

W1_autocor=w1_autocor(inds_n,:);
W2_autocor=w2_autocor(inds_sw,:);

W1_ceps=w1_ceps(inds_n,:);
W2_ceps=w2_ceps(inds_sw,:);


w1_phase=F_phase_noise2;
w2_phase=F_phase_SW1;

w2_raw=F_raw_SW;
w1_raw=F_raw_n;

w1_autocor=F_auto_noise2;
w2_autocor=F_auto_SW1;

w1_ceps=F_cepstrum_noise2;
w2_ceps=F_cepstrum_SW1;

w1_waveform=[sigma_noise' dr_noise'];
w2_waveform=[sigma_SW' dr_SW'];

X=w1_waveform(:,1);

W1_waveform(:,1)=w1_waveform(:,1)/max(w1_waveform(:,1));
W1_waveform(:,2)=w1_waveform(:,2)/max(w1_waveform(:,2));
W2_waveform(:,1)=w2_waveform(:,1)/max(w2_waveform(:,1));
W2_waveform(:,2)=w2_waveform(:,2)/max(w2_waveform(:,2));

figure; 
subplot(1,2,1);
plot(w1_waveform(:,1),w1_waveform(:,2),'x'); hold on;
plot(w2_waveform(:,1),w2_waveform(:,2),'x'); 
subplot(1,2,2);
plot(W1_waveform(:,1),W1_waveform(:,2),'x'); hold on;
plot(W2_waveform(:,1),W2_waveform(:,2),'x'); 

% Classes={w1_phase w2_phase ; W1_phase W2_phase ; w1_raw w2_raw ; w1_autocor w2_autocor ; w1_ceps w2_ceps ; W1_waveform W2_waveform}
Classes={W1_phase W2_phase ;  W1_autocor W2_autocor ; w1_ceps w2_ceps}

for i=1:3
    w1=cell2mat(Classes(i,1));
    w2=cell2mat(Classes(i,2));    
   [J1(i),J2(i),J3(i),FDR(i)]=SI(w1,w2);
end
    

[FDR ; J1 ; J2 ; J3 ]





 plot(F_phase_noise2(:,1),F_phase_noise2(:,2),'bx','Linewidth',1); hold on; 
 plot(F_phase_SW1(:,1),F_phase_SW1(:,2),'rx','Linewidth',1);

plot(sigma_noise,dr_noise,'bx'); hold on;
plot(sigma_SW,dr_SW,'rx');

