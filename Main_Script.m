clear all; clc;

Program_folder=pwd;                                     % current folder
% files = dir('SW*.flac');
% files = dir('AMAR*.flac');
files = dir('10M*.wav');

%% Pick Detector

Detector_flag=0; % 1- apply coda detector | 0- apply echolocation clicks detector
Plot_flag=1;     % 1- show detection figures | 0- dont show

Other_whale_Detector_flag=0; %1- apply other whale detector | 0- disable other whale detector

%% Determine operational parameters (Global parameters)

F_low_coda = 3e3;
F_high_coda = 7e3;
F_low_echo = 2e3;
F_high_echo = 22e3;
FsAnalyze = 48e3;
T_sec=10; %[sec]                         % Define duration of window for analysis in seconds
%%

if Detector_flag
    F_low=F_low_coda;
    F_high=F_high_coda;
else
    F_low=F_low_echo;
    F_high=F_high_echo;
end

Gather_IPI=[]; Gather_IPI_auto=[]; Gather_Features=[];

for file_ind=1:length(files)
    
    filename=files(file_ind).name;
    filesave=[Program_folder '\' filename '.xls'];
    
    Audio_name=[strfind(filename, 'flac') strfind(filename, 'wav')];
    if ~isempty(Audio_name)
        Timer=0;
        [y,Fs] = audioread(filename);                 % load recordings
        Y=y(:,1);                                     % Choose chanel one from the WRU
        if Fs < FsAnalyze
            S_factor = 1;
        elseif Detector_flag
            S_factor=floor(Fs/FsAnalyze);                 % Define factor for resampling to 48khz
        else
            S_factor=floor(Fs/FsAnalyze);
        end        
        File_duration=(1/Fs)*(length(Y)-1);           % Calculate duration of the loaded recording
        Y_decimated = decimate(Y,S_factor);           % Resample recording to 48khz
        F_ds=Fs/S_factor;                             % Sample frequency of the decimated recording (48khz by default)       
        T=F_ds*T_sec;                                 % Define duration of window for analysis in samples
        T_raw=Fs*T_sec;
        NOI=floor(File_duration/T_sec);               % Calculate the number of windows in the current recording
        buffer_index=0; Gather_TOA=[]; 
        
        for Buffer_ind=1:NOI
            [file_ind Buffer_ind]
            Y_filtered_raw=bandpass(Y(int32((Buffer_ind-1)*T_raw+1):int32((Buffer_ind-1)*T_raw+T_raw)),[F_low, F_high],Fs);
            Y_filtered=bandpass(Y_decimated(int32((Buffer_ind-1)*T+1):int32((Buffer_ind-1)*T+T)),[F_low, F_high],F_ds);     % Aply band pass filter and extract buffer             
            [TOA,IPI,IPI_auto,Features]=Run_Detector(Y_filtered_raw,Y_filtered,Fs,F_ds,Detector_flag,Plot_flag);  % Run detection
            Gather_TOA=[Gather_TOA (Buffer_ind-1)*T_sec+TOA];  % Acumulate time of arrivals from each buffer  
            Gather_IPI=[Gather_IPI IPI];
            Gather_IPI_auto=[Gather_IPI_auto IPI_auto];
            Gather_Features=[Gather_Features ; Features];
        end
        
         if Other_whale_Detector_flag
            DetectOtherWhale=find_whale_neighbours(Y,Gather_TOA,Fs); % Run 'other whale' detector
            if ~isempty(DetectOtherWhale)
                writematrix(DetectOtherWhale',filesave); % Save time of arrivals of the detected 'other whales clicks'
            end
         elseif ~isempty(Gather_TOA)               
                 writematrix(Gather_TOA',filesave);  % Save time of arrivals of the detected clicks                                                              
         end 

    end
end




 
% IPI_SW_m=1e3*Gather_IPI;
% % 
% % IPI_SW_m=Gather_IPI*1e3;
% save('IPI_SW_m.mat','IPI_SW_m')
% Model=18*rand(1,1e8);
% % 


% 
% figure;
% h1=histogram(Model,50,'Normalization','probability'); hold on;
% h2=histogram(IPI_noise,50,'Normalization','probability'); hold on; 
% h3=histogram(IPI_SW,50,'Normalization','probability'); hold on;
% h4=histogram(IPI_SW_m,50,'Normalization','probability'); 
% 
% dist_SW=round(KLDiv(h1.Values,h3.Values),2)
% dist_SW_m=round(KLDiv(h1.Values,h4.Values),2)
% dist_noise=round(KLDiv(h1.Values,h2.Values),2)
% legend('Noise model',['Noise transients: KLD=' num2str(dist_noise)],['Sperm whales: KLD=' num2str(dist_SW)],'Fontsize',12);
% ylabel('PDF','Fontsize',14); xlabel('IPI [ms]','Fontsize',14); 
% % figure;
% % plot(h1.Values,'-.','Linewidth',2); hold on;
% % plot(h2.Values,'-.','Linewidth',2);
% % plot(h3.Values,'-.','Linewidth',2);
% legend('Noise model',['Noise transients: KLD=' num2str(dist_noise)],['Single sperm whales: KLD=' num2str(dist_SW)],['Multiple sperm whales: KLD=' num2str(dist_SW_m)],'Fontsize',12);
% ylabel('PDF','Fontsize',14); xlabel('IPI [ms]','Fontsize',14);
% 

