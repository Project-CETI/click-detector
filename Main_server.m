clear all; clc;

    Program_folder=pwd;                                     % current folder
    Data_folder=Program_folder;             % folder that contains the recordings
    cd(Data_folder);
    files = dir('SW*.flac');

    Results_folder=Data_folder;             % folder that contains the recordings


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
    
        
    for file_ind=1:length(files)

        cd(Data_folder);
        filename=files(file_ind).name;
        Audio_name=[strfind(filename, 'flac') strfind(filename, 'wav')];
        if ~isempty(Audio_name)
            Timer=0;
            [y,Fs] = audioread(filename);                 % load recordings
            cd(Program_folder);
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
                Y_filtered=bandpass(Y_decimated(int32((Buffer_ind-1)*T+1):int32((Buffer_ind-1)*T+T)),[F_low, F_high],F_ds);     % Aply band pass filter and extract buffer             
                tic
                TOA=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag);  % Run detection
                toc
            end
        end
        
    end
    



