function [Send_seq,Send_seq2]=ICI_extract_Sequence2(time,ey_norm,locs,pks,Y_zoom,F_ds,consistency_T,ICI_max,ICI_min,Th)
% %  v=1;
% % locs=sort(MP_t(cell2mat(C_inds(v))));
% % pks=MP_p(cell2mat(C_inds(v)));
% locs=MP_t;
% pks=MP_p;

% [locs,Id]=sort(L_tag);
% pks=P_tag(Id);
% Th=Th_echo;
% ICI_max=ICI_max_echo;
% ICI_min=ICI_min_echo;
% ICI_max=1.8;


    Final_seq={}; Final_ICI={}; Final_L={};  Compare=[];  Candidate_Trains={}; Send_seq=[];
    Send_seq2=[]; Compare2=[];  Candidate_Trains2={};
    T_locs=locs;
    locs_inds=T_locs*F_ds;

    slice=0.015*F_ds;
    for ij=1:length(locs_inds)
          BANK(ij,:)=Y_zoom(int32(locs_inds(ij)-slice):int32(locs_inds(ij)+slice));
          M_bank(ij)=max(BANK(ij,:));
    end
    T_pks=M_bank;

    S=0; C_ind=0; C=[]; Triple=[];  Save_C=[]; Triple_ind=[];
    for n=1:length(T_locs)-1
        for m=n+1:length(T_locs)
             for l=m+1:length(T_locs)
                S=S+1;
                C(S,:)=[log((T_locs(l)-T_locs(m))/(T_locs(m)-T_locs(n))) n m l];
                AD(S)=(T_pks(l)+T_pks(m)+T_pks(n))^2/(3*(T_pks(l)^2+T_pks(m)^2+T_pks(n)^2)) ;
                if abs(C(S,1))<consistency_T && AD(S)>Th && T_locs(l)-T_locs(m)<ICI_max && T_locs(l)-T_locs(m)>ICI_min
                    C_ind=C_ind+1;
                    Triple(C_ind,:)= [T_locs(n) T_locs(m) T_locs(l)];
                    Triple_ind(C_ind,:)= [n m l];
                    Save_C(C_ind,:)=[abs(C(S,1)) C(S,2:end) T_locs(l)-T_locs(m) AD(S)];
                end
             end
        end
    end

    Tmp=Save_C;
    Flag=size(Save_C,1);

    if  Flag>0 && Flag<3
        Send_seq=unique(Tmp(:,2:4));
        Send_seq=Send_seq';
        Send_seq2=Send_seq;    
    elseif Flag>0
        if Flag<5
            En=1;
        else
            En=Flag-4;
        end
       for Q=1:En 
            Save_C=Tmp(Q:end,:);
            SEQ=Save_C(:,2:4);
            seq_save=[];
            i=1;
            Final_seq={};
            Final_flag=0;
            q=1; init=0;
            %%
            while(Final_flag==0)
                Candidates=[];

                for T=seq_save
                    [eliminate_inds , ~]=find(Save_C(:,2:4)==T);
                    Save_C(eliminate_inds,:)=[];
                    eliminate_inds=[]; a1=[];
                end

                if isempty(seq_save) && init>0
                    Save_C(1,:)=[];
                end

                if size(Save_C,1)<1
                    Final_flag=1;
                end

                seq_save=[];
                SEQ=Save_C(:,2:4);
                i=1;
                j=i+1; k=1;
                Exit_flag=0;
                while(i<(size(SEQ,1)) && Exit_flag==0)

                  s=1;
                    while(j<size(SEQ,1)+1)
                        Inter=intersect(SEQ(i,:),SEQ(j,:));
                        if length(Inter)==2
                            DIF=abs(Save_C(i,5)-Save_C(j,5));
                            if Inter(2)==SEQ(i,3) && Inter(2)~=SEQ(j,3) && DIF<0.2
                                Candidates(s,:)=[SEQ(j,:) Save_C(j,1) Save_C(j,5) j Save_C(j,6)];
                                s=s+1;
                            end
                        end
                        j=j+1;        
                    end

                    if length(Candidates)>1
                        Candidates_diversity=Candidates(:,7);
                        Choose=find(Candidates_diversity==max(Candidates_diversity));
                        seq_save(k,:)=Candidates(Choose(1),:);
                        i= seq_save(k,6);
                        j=i+1;
                        k=k+1;
                        Candidates=[];
                    elseif isempty(Candidates)
                        Exit_flag=1;
                    else
                        seq_save(k,:)=Candidates;
                        i= seq_save(k,6);
                        j=i+1;
                        k=k+1;
                        Candidates=[];
                    end
                end

                contents=seq_save;

                if ~isempty(seq_save)
                 seq_save=[SEQ(1,1) seq_save(:,1)' seq_save(end,2:3)];
                 seq_save=unique(seq_save);
                end

                if length(seq_save)<4
                    Final_seq(q)={[]};
                    Final_ICI(q)={[]};
                    Final_L(q)={[]};
                    L1(q)=0;
                    q=q+1;
                else
                    Final_seq(q)={seq_save};
                    Final_ICI(q)={ones(1,length(seq_save))*mean(contents(:,5))};
                    Final_L(q)={ones(1,length(seq_save))*length(seq_save)};
                    L1(q)=length(seq_save);
                    q=q+1;
                end
                        contents=[];
                        SEQ=[];
                        init=init+1;
            end
            
            S_index=find(L1==max(L1));
            Candidate_Trains(Q)={cell2mat(Final_seq(S_index(1)))};
            Compare(Q)=max(L1);
            
            Candidate_Trains2(Q)={cell2mat(Final_seq)};
            Compare2(Q)=sum(L1);
            L1=[];
       end
      
       Final_selection=find(Compare==max(Compare));
       Final_selection2=find(Compare2==max(Compare2));
       if ~isempty(Final_selection)
          Send_seq=Candidate_Trains(Final_selection(1));
       end
       if ~isempty(Final_selection2)
          Send_seq2=Candidate_Trains2(Final_selection2(1));
       end       
    
    Send_seq2=cell2mat(Send_seq2);
    end

    %% Extract the taged clicks
% %         peaks_tag=find(Compare==max(Compare));
%         for i=1:length(Compare)
%             Score(i)=mean(pks(cell2mat(Candidate_Trains(i))));           
%         end
%         
%         Eli=find(Compare<10); Score(Eli)=0;
%         Tag_choice=find(Score==max(Score));
%         if Score(Tag_choice(1))>0.5
%             Send_seq2=[];
%             Send_seq2=cell2mat(Candidate_Trains(Tag_choice(1)));
%         end
    
end
    
