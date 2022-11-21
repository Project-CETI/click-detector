clear all; clc;

Program_folder=pwd;
% files = dir('twowhales*.flac');
files = dir('twowhales*.flac');

 
    %% Pick Detector

    Detector_flag=0; % customized for buoy recievers: 1- apply coda detector | 0- apply echolocation clicks detector
    Plot_flag=1;     % 1- show detection figures | 0- dont show
    Tag_flag=1;    % 1- customized for Dtags: set to 1 for recordings from tags.
    Error_flag=length(files);

    %% Check filename:
    if Error_flag==0
        error('File not found. Please check the file name (line 4)');        
    end
        
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
            buffer_index=0; 
            Gather_TOA=[]; 
            TOA_other_whale=[];
            TOA_other=[];

            for Buffer_ind=1:NOI
                Y_filtered=bandpass(Y_decimated(int32((Buffer_ind-1)*T+1):int32((Buffer_ind-1)*T+T)),[F_low, F_high],F_ds);     % Aply band pass filter and extract buffer                             
                if Tag_flag
                    [~,TOA,TOA_other]=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag,Tag_flag);  % Run detection
                else
                    [TOA,~,~]=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag,Tag_flag);  % Run detection                    
                end
                Gather_TOA=[Gather_TOA (Buffer_ind-1)*T_sec+sort(TOA)];
                TOA_other_whale=[TOA_other_whale (Buffer_ind-1)*T_sec+sort(TOA_other)];
            end
            
            if Tag_flag
                  filesave=[Program_folder '\' 'Tag_detected_' filename '.xls'];
                  writematrix(Gather_TOA',filesave);
                  filesave=[Program_folder '\' 'Other_whale_detected_' filename '.xls'];
                  writematrix(TOA_other_whale',filesave);                                                    
            else                
                 writematrix(Gather_TOA',filesave);  % Save time of arrivals of the detected clicks                                                              
            end 
        end
        
    end
    



