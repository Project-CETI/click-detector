function [T_plot,P]=ICI4(time,ey_norm,locs,pks,Y_zoom,F_ds,f_P_coda,d_P_coda,E_ratio,D_th,ICI_max,ICI_min,Th,E_th,consistency_T)
                    
% D_th=0.45;
% ICI_max=1;
% ICI_min=0.05;
% Th=0.8;
% % consistency_T=consistency_T_coda;
% consistency_T=0.2;
% locs=MP_t;
% pks=MP_p;


    Final_seq={}; Final_ICI={}; Final_L={}; Validation=[]; M_bank=[]; E=[];
    T_locs=locs;
    locs_inds=T_locs*F_ds;

    slice=0.015*F_ds;
    Acomulate=2e-3*F_ds;
%     figure;
    for ij=1:length(locs_inds)
        if locs_inds(ij)-slice>0
          BANK(ij,:)=Y_zoom(int32(locs_inds(ij)-slice):int32(locs_inds(ij)+slice));
          M_bank(ij)=max(abs(BANK(ij,:)));
          E(ij)=sum(Y_zoom(int32(locs_inds(ij)-Acomulate):int32(locs_inds(ij)+Acomulate)).^2);
%           subplot(6,6,ij);plot(BANK(ij,:)); 
%           title(['E = ' num2str(round(1e4*E(ij)),3) 'K']); 
%           ylim([-0.05 0.05]);
        end
    end
    
    T_pks=M_bank;
    T_E=1e3*E;
    
    if length(T_E)>10
        Keep1=GetCodas(T_E); 
        Eliminate_inds=Eliminate_MultiPath(locs,pks);
        for i=1:length(Eliminate_inds)
           Keep1(Keep1==Eliminate_inds(i))=[];
        end
    else
        Keep1=[1:length(T_E)];
    end
    
    S=0; C_ind=0; C=[]; Triple=[];  Save_C=[]; Triple_ind=[];

    for n=1:length(T_locs)-1
        for m=n+1:length(T_locs)
             for l=m+1:length(T_locs)
                S=S+1;
                C(S,:)=[log((T_locs(l)-T_locs(m))/(T_locs(m)-T_locs(n))) n m l];
                AD(S)=(T_pks(l)+T_pks(m)+T_pks(n))^2/(3*(T_pks(l)^2+T_pks(m)^2+T_pks(n)^2)) ;
                E_test=sum([sum(n==Keep1) sum(m==Keep1) sum(l==Keep1)]);
                V_test=JF(T_E(n),T_E(m),T_E(l));
                if abs(C(S,1))<consistency_T && AD(S)>Th && T_locs(l)-T_locs(m)<ICI_max && T_locs(l)-T_locs(m)>ICI_min && E_test==3 && V_test>E_th
                    C_ind=C_ind+1;
                    Triple(C_ind,:)= [T_locs(n) T_locs(m) T_locs(l)];
                    Triple_ind(C_ind,:)= [n m l];
                    Save_C(C_ind,:)=[abs(C(S,1)) C(S,2:end) T_locs(l)-T_locs(m) AD(S) mean([T_E(n) T_E(m) T_E(l)]) JF(T_E(n),T_E(m),T_E(l))];

                end
             end
        end
    end
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Coda_detector_fileName='Coda_detector_revised2.fis'; 
%     Coda_detector = readfis(Coda_detector_fileName);

    Coda_detector_fileName='Coda_validation_via_duration.fis'; 
    Coda_detector = readfis(Coda_detector_fileName);
    
   if ~isempty(Save_C)
        Selected=Save_C;
        for i=1:size(Selected,1)
             trip= Selected(i,2:4);
             TIMES(i,:)=locs(trip);
             PEAKS(i,:)=pks(trip);
             Fcoda=f_P_coda(trip);
             Dcoda=d_P_coda(trip);
             Er=E_ratio(trip);
             for ind=1:3                             
%                    output(ind) = evalfis(Coda_detector,[Fcoda(ind) Dcoda(ind) Er(ind)]) ;  
                   output(ind) = evalfis(Coda_detector,[Dcoda(ind)]) ;                                

             end

             Validation(i)=mean(output)>D_th;
             M(i)=mean(output);
        end
   end
    if ~isempty(Validation)
            valid_inds=find(Validation==1);           
            T_valid=TIMES(valid_inds,:);          
            P_valid=PEAKS(valid_inds,:);          

            T_plot=unique(T_valid(:));
            P_inds=zeros(length(T_plot),1);
            T_plot=T_plot';
          for i=1:length(T_plot)
             n=T_plot(i);
             [~,P_inds(i)]=min(abs(locs-n));       
          end
    else
        T_plot=[];
        P_inds=[];
    end
    P=pks(P_inds);
%     M

end