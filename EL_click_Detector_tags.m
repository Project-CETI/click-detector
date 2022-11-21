function [TOA_tag,TOA_other]=EL_click_Detector_tags(F_ds,Y_filtered,Plot_flag,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo,MP_thresh,W_seg)


TOA_tag=[]; TOA_other=[];
locs=[]; pks=[];
L_other=[]; L_tag=[];
P_other=[]; P_tag=[];
Valid_flag=0;
        
        Y_zoom=Y_filtered;     % Aply band pass filter               
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
        [ey,~]=energyop(Y_zoom,0);              % Apply TKEO (ey is the output signal with enhanced SNR)
        ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
        ED_thresh=0.1*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
        time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)]; 
        [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',15e-3);  % Apply instantaneous energy detector (find peaks)          

        %% Eliminate transients bellow the minimum allowed SNR

        [Array,I]=sort(pks,'descend');
       [SDCM_All, GF] = get_jenks_interface(Array);
%        figure; 
%        subplot(2,1,1); plot(Array,'rx','Linewidth',2); grid on; ylabel('Normalized peak');
%        subplot(2,1,2); plot(SDCM_All,'gx','Linewidth',2);  grid on; ylabel('GF');
       
              
    T_locs=locs; M_bank=[]; BANK=[];
    locs_inds=T_locs*F_ds;
    slice=0.004*F_ds;
    for ij=1:length(locs_inds)
          BANK(ij,:)=Y_zoom(int32(locs_inds(ij)-slice):int32(locs_inds(ij)+slice));
          M_bank(ij)=max(BANK(ij,:));
    end
    T_pks=M_bank/max(M_bank);
       
       [Array2,I2]=sort(T_pks,'descend');
       [SDCM_All2, GF2] = get_jenks_interface(Array2);
%        figure; 
%        subplot(2,1,1); plot(Array2,'rx','Linewidth',2); grid on; ylabel('Normalized peak');
%        subplot(2,1,2); plot(SDCM_All2,'gx','Linewidth',2);  grid on; ylabel('GF');
           
                
       Seg_T=find(GF2==max(GF2));
                        
        P_tag=pks(I(1:Seg_T));
        L_tag=locs(I(1:Seg_T));

[l_cands,Id]=sort(L_tag);
p_cands=P_tag(Id);
    [~,Final_seq]=ICI_extract_Sequence2(time,ey_norm,l_cands,p_cands,Y_zoom,F_ds,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo);   % Run click trains detector 

    if ~isempty(Final_seq)
        L_tag=l_cands(Final_seq);          
        P_tag=p_cands(Final_seq);      
        P_other=pks(I(Seg_T+1:length(pks)));
        L_other=locs(I(Seg_T+1:length(pks)));

        Memberhood=ismember([1:Seg_T],Final_seq);
        Pass_ind=find(Memberhood==0);
        L_other=[L_other ; l_cands(Pass_ind)]; 
        P_other=[P_other ; p_cands(Pass_ind)]; 
        
        %% detection and verification of other whale clicks
       [MP_t,MP_p,IPI]=Multipulse_locs_echos(Y_zoom,ey_norm,L_other,P_other,F_ds,W_seg,MP_thresh,0);  % Run multipulse detector

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
                
                for Li=9:9                   
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
                    
              
               if Valid_flag
                          [MP_t,mem]=sort(MP_t);
                           MP_p=MP_p(mem);
                           [~,Other_seq]=ICI_extract_Sequence2(time,ey_norm,MP_t,MP_p,Y_zoom,F_ds,consistency_T,ICI_max_echo,ICI_min_echo,Th_echo);   % Run click trains detector 
                           L_other=MP_t(Other_seq);
                           P_other=MP_p(Other_seq);
               else
                  L_other=[];  
                  P_other=[];             
               end
               
        else
           L_other=[];  
           P_other=[];             
        end
               
        
    else
        L_tag=[]; P_tag=[];
        L_other=[];  P_other=[];
    end
        
                
        TOA_tag=L_tag'; TOA_other=L_other;
               
        if Plot_flag
            figure;set(gcf, 'Position', get(0,'Screensize'));
            subplot(3,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal'); 
            subplot(3,1,2); plot(time,ey_norm); hold on; 
            xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
            hold on; plot(L_tag,P_tag,'ro','Linewidth',2);           
            legend('','Tagged whale','Fontsize',10);
            subplot(3,1,3); plot(time,ey_norm); hold on;            
            hold on; plot(L_other, P_other,'go','Linewidth',2);
            legend('','Other whale','Fontsize',10);
        end
end

