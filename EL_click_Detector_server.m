function TOA=EL_click_Detector_server(SNR_window,SNR_thresh,Fs,F_ds,Y_filtered,Plot_flag,MP_thresh,W_seg,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo)

% SNR_window=SNR_window_echo;
% SNR_thresh=SNR_thresh_echo;
% consistency_T=consistency_T_echo;


        FA=zeros(9,20);
        Valid_flag=0;
        TOA={};
        Groups={}; 
        Groups_pks={};
        
        c=0;                   % counter for the transient detection stage
        Y_zoom=Y_filtered;     % Aply band pass filter               
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
        [ey,~]=energyop(Y_zoom,0);              % Apply TKEO (ey is the output signal with enhanced SNR)
        ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
        ED_thresh=1*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
        time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];   
        [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',25e-3);  % Apply instantaneous energy detector (find peaks)          
        locs_samples=locs*F_ds;                 % store indices of time readings of the detected impulses

        %% Eliminate transients bellow the minimum allowed SNR
        
        Locs=[]; Pks=[];
        for i=1:length(locs)
            if locs_samples(i)>SNR_window && (locs_samples(i)+SNR_window)<length(ey_norm)
                tmp=ey_norm(int32(locs_samples(i)-SNR_window):int32(locs_samples(i)+SNR_window));
                SNR(i)=pks(i)/median(tmp); 
                if SNR(i)>SNR_thresh
                    c=c+1;
                    Locs(c)=locs(i);
                    Pks(c)=pks(i);
                end
            end
        end

        if length(Pks)>60
           [Pks2,I] = maxk(Pks,60);
           Locs2=Locs(I);
           Pks=Pks2;
           Locs=Locs2;
        end

        if length(Locs)>2 && length(Locs)<100
           %% Show plots for the transient detection stage

        if Plot_flag==1
            figure;set(gcf, 'Position', get(0,'Screensize'));
            subplot(4,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal');
            subplot(4,1,2); plot(time,ey_norm); hold on; 
            xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
            hold on; plot(Locs,Pks,'x')    
        end
        

            %% Multipulse Detection Stage

            [MP_t,MP_p,IPI]=Multipulse_locs_echos(Y_zoom,ey_norm,Locs,Pks,Fs,F_ds,W_seg,MP_thresh,Plot_flag);  % Run multipulse detector
        

          if length(IPI)>2
            [C_inds,Features_phase]=segment_IPIs_revised(IPI);   

            %% Detection of Click Trains 
                Params={'T0','T3','T4','T5','T6','lim_u','lim_d','L_shape','l5','l6'};

                for lp=1:length(Params)
                    load(cell2mat(Params(lp)));    
                end
                
                l1=lim_d(1);
                l2=lim_d(1);                
                Features=Features_phase;
                
                for Li=1:9                   
                    l1=lim_d(1);
                    l2=lim_d(Li);
                    l3=lim_u(Li);
                    l4=lim_u(1);                
                    y0=polyval(T0,Features(1));
                    y1=L_shape(1)*Features(1);
                    y2=L_shape(Li)*Features(1);
                    y3=polyval(T3,Features(1));
                    y4=polyval(T4,Features(1));
                    y5=T5*Features(1);
                    y6=T6*Features(1);
                    D_thresh_vals=linspace(0.1,0.95,20);

                    if length(MP_t)>2 
                         Reject=0;
                         if Features(1)>l1 && Features(1)<l2 && Features(2)>y1  && Features(2)<y0  
                             Reject=1;                            
                         elseif Features(1)>l2 && Features(1)<l3 && Features(2)>y1 && Features(2)<y2
                             Reject=1;
                         elseif Features(1)>l3 && Features(1)<l4 && Features(2)>y1 && Features(2)<y3
                             Reject=1;
                         elseif Features(1)>l5 && Features(1)<l6 && Features(2)>y4 && Features(2)<y6
                             Reject=1;  
                         elseif Features(1)>l6  && Features(2)>y5 && Features(2)<y6
                             Reject=1;                                
                         end

                       if Reject==0
                           [Times,I_Times]=sort(MP_t(cell2mat(C_inds)));
                           Po=MP_p(cell2mat(C_inds));
                           Powers=Po(I_Times);
                           [~,Final_seq]=ICI_extract_Sequence2(time,ey_norm,Times,Powers,Y_zoom,F_ds,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo);   % Run click trains detector 
                           Groups={Times(Final_seq)};
                           Groups_pks={Powers(Final_seq)};                                                                                
                           Detected_Echos=sort(cell2mat(Groups));
                           NOD=length(cell2mat(Groups));            
                           Train_percentage2=NOD/length(Times);
                          Groups={}; Groups_pks={};
                       else 
                           Train_percentage2=0;
                       end
                     
                       for qq=1:length(D_thresh_vals)
                           D_thresh=D_thresh_vals(qq);   
                          if Train_percentage2>D_thresh
                              Valid_flag=1;
                              FA(Li,qq)=1;
                          end                                     
                       end
                     end
                    
                end
            
          end
        end
        
                            if Valid_flag==1
%                      if Train_percentage>-1
                          [MP_t,mem]=sort(MP_t);
                           MP_p=MP_p(mem);
                           Th_selected=Adapt_threshold(MP_t,MP_p);
                            if Th_selected<1
                                eliminate_pks=find(MP_p<Th_selected);
                                MP_t(eliminate_pks)=[];
                                MP_p(eliminate_pks)=[];
                            end
                           [~,Final_seq]=ICI_extract_Sequence2(time,ey_norm,MP_t,MP_p,Y_zoom,F_ds,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo);   % Run click trains detector 
                           Detected_Echos=MP_t(Final_seq);
                           Detected_Echos_pks=MP_p(Final_seq);
                           TOA={Detected_Echos};
%                       end
%                     end
%                                         
%               
              if Plot_flag
                  subplot(4,1,4);
                  Time_of_arrival=cell2mat(Groups);
                  Peak_of_arrival=cell2mat(Groups_pks);
                  subplot(4,1,4); plot(time,ey_norm); hold on;
                  if ~isempty(Detected_Echos)
                     plot(Detected_Echos,Detected_Echos_pks,'go','Linewidth',2);
                  end
%                   plot(Time_of_arrival,Peak_of_arrival,'go','Linewidth',2)
                  xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
                  legend('','Echolocation clicks');
              end
              
            end
        
        
end
               
           

