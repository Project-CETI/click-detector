function [Coda_times,Coda_Type]=Classify_Coda_type(Detected_pattern_t,Detected_pattern_p,time,ey_norm,Detection_flag)

           Sequences=diff(Detected_pattern_t);
           S_ind=find(Sequences>1);
           Coda_Length=0;
           Coda_Type={};  Coda_times={}; Coda_pks={}; 
           if Detection_flag
             subplot(5,1,5); plot(time,ey_norm); ylim([0 1]); xlabel('time [sec]'); ylabel('TKEO'); title('Detection of Codas');   
             legendInfo{1}=['']; 
           end

           
           for seqs=1:length(S_ind)+1
               Coda_ICIs=[];

                   if seqs==length(S_ind)+1
                       Coda_ICIs=Sequences(seqs+Coda_Length(seqs):end);
                       Decide= Coda_ICIs>0.04 & Coda_ICIs<1;
                       if sum(Decide)>=1
                           Coda_times(seqs)={Detected_pattern_t(seqs+Coda_Length(seqs):end)};
                           Coda_pks(seqs)={Detected_pattern_p(seqs+Coda_Length(seqs):end)};
                       else
                           Coda_times(seqs)={[]};
                           Coda_pks(seqs)={[]};
                       end
                   else
                       Coda_ICIs=Sequences(seqs+Coda_Length(seqs):S_ind(seqs)-1);
                       Decide= Coda_ICIs>0.04 & Coda_ICIs<1;
                       if sum(Decide)>=1
                           Coda_times(seqs)={Detected_pattern_t(seqs+Coda_Length(seqs):S_ind(seqs))};
                           Coda_pks(seqs)={Detected_pattern_p(seqs+Coda_Length(seqs):S_ind(seqs))};
                       else
                           Coda_times(seqs)={[]};
                           Coda_pks(seqs)={[]};
                       end
                       Coda_Length(seqs+1)=Coda_Length(seqs)+length(Coda_ICIs);
                   end
               if ~isempty(Coda_times)
                   Pattern=zeros(1,9);                                 
                   Pattern(1:length(Coda_ICIs))=Coda_ICIs;            
                   Coda_Type(seqs)={Coda_Type_clas(Pattern)};

                   if Detection_flag
                       hold on; plot(cell2mat(Coda_times(seqs)),cell2mat(Coda_pks(seqs)),'o','Linewidth',2); 
                       legendInfo{seqs+1} = ['Coda ' num2str(seqs) ' = ' cell2mat(Coda_Type(seqs))];
                   end
               end
           end   
           
           if Detection_flag
               legend(legendInfo)
           end
end



