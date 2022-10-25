function [C_inds,Reduced]=segment_IPIs_revised(IPI)

          id_save={}; Group={};
          Reduced=[0 0]; mf=[]; ms=[]; C_inds={}; 
          [input,I] = sort(1e3*IPI,'descend');

          q=0; DIF=1;
          Original_input=input;     

           x=[1:length(input)];
           P = polyfit(x,input,1);
           yfit = polyval(P,x);
           if length(yfit)>1
               m=yfit(1)-yfit(2);
           else
                m=yfit(1);
           end

           Reduced=[m std(input)];
           C_inds={I};
           
%           figure; 
%           subplot(3,1,1); plot(input,'bo','Linewidth',2); hold on;  plot(x,yfit,'r-.','Linewidth',2);
%           xlabel('Transient index','Fontsize',14); ylabel('IPI','Fontsize',14); ylim([0 20]); xlim([0 60]);
% %           legend('',['m=' num2str(round(m,2)), ', \sigma_{IPI}=' num2str(round(std(input)),2) '[ms]'],'Fontsize',14);
%           legend('','$\hat{m}=0.22,~\sigma_{z}=4[ms]$','interpreter','latex','Fontsize',16);
% 
%           subplot(3,1,3); plot([12:27],input(SW_candidate_inds),'bo','Linewidth',2); hold on;  plot([12:27],yfit,'r-.','Linewidth',2);
%           xlabel('Transient index','Fontsize',14); ylabel('IPI','Fontsize',14); ylim([0 16]);
%           xlim([0 54]);         
% %           legend('',['m=' num2str(round(mf,2)), ', \sigma_{IPI}=' num2str(round(std(input_f),2)) '[ms]'],'Fontsize',14);
%           legend('','$\hat{m}^s=0.02,~\sigma_{z}^s=0.15[ms]$','interpreter','latex','Fontsize',16);


%         figure; subplot(2,1,1);
        samples=3;
        input=[ input zeros(1,2*samples)];
        for q=1:length(input)-samples
          x=[q:q+samples];
          P = polyfit(x,input(q:q+samples),1);
          yfit = polyval(P,x);
          ms(q)=yfit(1)-yfit(2);
%           plot(input,'x','Linewidth',2); grid on; xlabel('Click #'); ylabel('IPI [ms]');
%           hold on; plot(x,yfit,'r-.','Linewidth',2);
        end


%        subplot(3,1,2);
%        plot(ms,'-x','Linewidth',2); grid on; xlabel('q','Fontsize',14); ylabel('Shape','Fontsize',14);
% xlabel('$q$','interpreter','latex','Fontsize',16);
% ylabel('$\hat{\textbf{m}_q}$','interpreter','latex','Fontsize',16);

       id=[]; g=0;
    for q=1:length(Original_input)-samples
        if abs(ms(q))<0.15
            id=[id q];
        elseif length(id)>1
            g=g+1;
            Group(g)={[id:id(end)+samples]};
            id=[];    
        else
            id=[];
        end    
    end

    if ~isempty(Group)
        Le=[];
        for i=1:size(Group,2)
            GL=length(cell2mat(Group(i)));
            Le=[Le GL];
        end

        SW_candidate=Group(find(Le==max(Le)));
        if size(SW_candidate,2)>1
            for q=1:size(SW_candidate,2)
                Test(q)=std(Original_input(cell2mat(SW_candidate(q))));
            end
        else
            Test=1;
        end
        SW_candidate_inds=cell2mat(SW_candidate(find(Test==min(Test))));
       input_f=Original_input(SW_candidate_inds);
       xf=[1:length(input_f)];
       Pf = polyfit(xf,input_f,1);
       yfit = polyval(Pf,xf);           
       mf=yfit(1)-yfit(2);
       Reduced=[mf std(input_f)];
       C_inds={};
       C_inds={I(SW_candidate_inds)};       
    end
 
end


        