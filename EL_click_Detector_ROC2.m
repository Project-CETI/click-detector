function [TOA,IPI,IPI_auto,Features]=EL_click_Detector_ROC2(SNR_window,SNR_thresh,Fs,F_ds,Y_filtered_raw,Y_filtered,Plot_flag,MP_thresh,W_seg,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo)

% SNR_window=SNR_window_echo;
% SNR_thresh=SNR_thresh_echo;
% consistency_T=consistency_T_echo;

IPI=[]; IPI_auto=[];
FA=zeros(20,20);
Valid_flag=0;
      TOA={};
      Detected_Echos_pks=[]; Detected_Echos=[];
      Train_percentage=0; Train_percentage2=0;
       Features=[10 10]; Raw_Features=[10 10];

        Groups={}; 
        Groups_pks={};
        
        c=0;                   % counter for the transient detection stage
        Y_zoom=Y_filtered;     % Aply band pass filter               
        Y_zoom_raw=Y_filtered_raw;
        
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
        [ey,~]=energyop(Y_zoom,0);             % Apply TKEO (ey is the output signal with enhanced SNR)
        ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
        ED_thresh=1*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
%         ED_thresh=0.05;
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

            [MP_t,MP_p,MP_thresh_vals,IPI,IPI_auto]=Multipulse_locs_echos(Y_zoom_raw,Y_zoom,ey_norm,Locs,Pks,Fs,F_ds,W_seg,MP_thresh,Plot_flag);  % Run multipulse detector
            IPI_auto=abs(IPI_auto);
        
%        figure; 
%        subplot(4,1,4); plot(time,ey_norm); hold on; plot(MP_t(cell2mat(C_inds(1))),MP_p(cell2mat(C_inds(1))),'og','Linewidth',2); 
%        plot(MP_t(cell2mat(C_inds(1))),MP_p(cell2mat(C_inds(1))),'og','Linewidth',2); 
         
            if isempty(MP_t)
                MP_t=0;
            elseif length(MP_t)>4
                 [in,~] = sort(1e3*IPI,'descend');
                  xin=[1:length(in)];
                  Pin = polyfit(xin,in,1);
                  yfit = polyval(Pin,xin);
                  mi=yfit(1)-yfit(2);
                  Raw_Features=[mi std(1e3*IPI)];
            end            
 
          S_test=std(IPI)-std(IPI_auto);
          if S_test>0
              IPI=IPI_auto;
          end
          if length(IPI)>2
%             [m,C_std,C_inds,Features,Cut]=segment_IPIs(IPI);
            [C_inds,Features]=segment_IPIs_revised(IPI);    

                %% Detection of Click Trains 
                Line_parms={'ms.mat','r.mat','p.mat','lims.mat','lim2.mat','lim3.mat'};
                for q=1:size(Line_parms,2)
                   load(cell2mat(Line_parms(q)));
                end
                
                for Li=9:9
                    y1=ms(Li)*Features(1);
                    y2=polyval(p,Features(1));
                    y3=ms(52)*Features(1)+r(52);              

                    D_thresh_vals=linspace(0.1,0.95,20);


                    if Features(1)<0.04
                        Valid_flag=1;
                        FA(Li,:)=ones(1,20);
                    else
                         if length(MP_t)>2 
                             Reject=0;
                             if Features(1)>lims(Li) && Features(1)<lim2 && Features(2)>y1  && Features(2)<y2  
                                 Reject=1;                            
                             elseif Features(1)>lim3 && Features(1)<lim2 && Features(2)>3.8 && Features(2)<y3
                                 Reject=1;
                             elseif Features(1)>lim2 && Features(2)<y3 && Features(2)>y1
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
                         end

                       for qq=12:12%length(D_thresh_vals)
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


                      
                Valid_flag=1;     
                    if Valid_flag==1
%                      if Train_percentage>-1
                          [MP_t,mem]=sort(MP_t);
                           MP_p=MP_p(mem);
%                            Th_selected=Adapt_threshold(MP_t,MP_p);
%                             if Th_selected<1
%                                 eliminate_pks=find(MP_p<Th_selected);
%                                 MP_t(eliminate_pks)=[];
%                                 MP_p(eliminate_pks)=[];
%                             end
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
%               
           
 end
