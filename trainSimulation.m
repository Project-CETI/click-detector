

consistency_T=0.22;
NOS_all=round(linspace(5,60,52));
F3=[];
NOS_all=5:60;
% figure;
for j=1:52     
    j    
    NOS=NOS_all(j);
    for i=1:1e3      
%           TEST=rand(1,NOS_all(j))*18e-3;
% %           [~,~,C_inds,Reduced,Cut]=segment_IPIs(TEST);
%           [C_inds,Reduced]=segment_IPIs_revised(TEST);
% %           Reduced=[Reduced length(cell2mat(Cut)) NOS];
%           Reduced=[Reduced length(cell2mat(C_inds))];
%           F3=[F3 ; Reduced];
          T_pks=rand(1,NOS);   TEST=rand(1,NOS)*10;
         [Send_seq,Send_seq2]=ICI_sim(TEST,T_pks,consistency_T);
         N(j,i)=length(Send_seq2);
    end    
end



Gt=0.6;

for j=1:29
    Pfa(j)=sum((N(j,:)/(NOS_all(j)))>Gt)/size(N,2);
end

figure; plot(5:33,Pfa,'x--','Linewidth',2)









figure; 
subplot(1,2,1); plot(F3(:,1),F3(:,2),'.'); axis([0 10 0 10]);
subplot(1,2,2); plot(FS(:,1),FS(:,2),'.');  axis([0 10 0 10]);


[~,I]=sort(F3(:,3));

X=F3(I,:);
A=unique(X(:,3));

    inds=find(X(:,3)==A(1));
    mu_m=mean(X(inds,1)); s_m=std(X(inds,1));
    mu_sigma=mean(X(inds,2)); s_sigma=std(X(inds,2));
    Shapes(1)=(mu_sigma-s_sigma)/(mu_m-s_m);


q=[1:52];
% we=5.5*(1-exp(-q*0.01-1));
we=5.5*(1-exp(-q*0.02-0.5));
gamma=we;
figure;
colorstflex = {'#FF0000','#00FF00','#0000FF','#00FFFF','#FF00FF','#FFFF00','#000000','#D95319','#7E2F8E','#77AC30','#A2142F'};

for i=5:8%length(A)
    inds=find(X(:,3)==A(i));
    mu_m=mean(X(inds,1)); s_m=std(X(inds,1));
    mu_sigma=mean(X(inds,2)); s_sigma=std(X(inds,2));
    Shapes(i)=(mu_sigma-s_sigma)/(mu_m-s_m);
    L_shape(i)=0.5*(Shapes(i)+Shapes(i-1));
    T(i,:)=[mu_m-gamma(i)*s_m mu_sigma-gamma(i)*s_sigma];
%     plot(X(inds,1),X(inds,2),'x','color',cell2mat(colorstflex(i)),'Linewidth',2); hold on; axis([0 3.5 0 8]);
    plot(X(inds,1),X(inds,2),'.','Linewidth',2); hold on; axis([0 3.5 0 8]);
    xxx=linspace(0.8,3,1e6);
    yl=L_shape(i)*xxx;
    plot(xxx,yl,'k-','Linewidth',3); 
    inds=[];
     legendInfo{i}=[num2str(i+A(1)-1), '-residual transients: '];
     legend(legendInfo,'Fontsize',12);
end

x=T(:,1);
P = polyfit(x,T(:,2),2);
m_T=linspace(0.21,1.2,1e6);
yfit = polyval(P,m_T);
hold on;  plot(m_T,yfit,'k-','Linewidth',3);  axis([0 3.5 0 8]);
xu=linspace(0.2,0.4,1e6); xd=linspace(1.1,3,1e6);
y0=L_shape(2)*xd; y1=(2*L_shape(end)-0.98*L_shape(end-1))*xu;
m_f=2*L_shape(end)-0.98*L_shape(end-1);


hold on;  plot(xd,y0,'k-','Linewidth',3); 
hold on;  plot(xu,y1,'k-','Linewidth',3); 

yfit0=polyval(P,xd);
yfit1=polyval(P,xu);

[~,d2]=min(abs(y0-yfit0));
[~,u2]=min(abs(y1-yfit1));
lim2_d=xd(d2);
lim2_u=xu(u2);

[lim2_d lim2_u]


T4=P;

T0=polyval(P,linspace(lim2_u,lim2_d,1e2));
T1=L_shape(2)*linspace(lim2_d,3.5,1e2);
T2=m_f*linspace(lim2_u,1.3,1e2);
hold on;  plot(linspace(lim2_u,lim2_d,1e2),T0,'k-','Linewidth',3); 
hold on;  plot(linspace(lim2_d,3.5,1e2),T1,'k-','Linewidth',3); 
hold on;  plot(linspace(lim2_u,1.3,1e2),T2,'k-','Linewidth',3); 







%%
elim=[3 14 24 32 40 41 52 53 65 72 78 79 88 94 108 109 117 118 131 138 139];
Gather_Features_phase(elim,:)=[];
hold on;  plot(Gather_Features_phase(:,1),Gather_Features_phase(:,2),'rx','Linewidth',2);  


%%


X=F3(I,:);
X2_inds=find(X(:,1)<0.3 & X(:,2)<0.7);
X2=X(X2_inds,:);
A2=unique(X2(:,3));
q2=[10:-1:1];
we=3.5*(1-exp(-q2*0.25));
% we=5.5*(1-exp(-q*0.05-1));
gamma=we;
gamma_up=2.5;
figure;

    inds=find(X2(:,3)==A2(1));
    mu_m=mean(X2(inds,1)); s_m=std(X2(inds,1));
    mu_sigma=mean(X2(inds,2)); s_sigma=std(X2(inds,2));
    Shapes(1)=(mu_sigma-s_sigma)/(mu_m-s_m);


for i=2:10%length(A2)
    inds=find(X2(:,3)==A2(i));
    mu_m=mean(X2(inds,1)); s_m=std(X2(inds,1));
    mu_sigma=mean(X2(inds,2)); s_sigma=std(X2(inds,2));
    T(i,:)=[mu_m-gamma(i)*s_m mu_sigma-gamma(i)*s_sigma];
    T_up(i,:)=[mu_m+gamma_up*s_m mu_sigma+gamma_up*s_sigma];
    Shapes(i)=(mu_sigma-s_sigma)/(mu_m-s_m);
    L_shape(i)=0.5*(Shapes(i)+Shapes(i-1));
    Lines(i,:)=linspace(lim_d(i-1),lim_u(i-1),100)*0.5*(Shapes(i)+Shapes(i-1));
    plot(X2(inds,1),X2(inds,2),'.','Linewidth',2); hold on; 
    inds=[];
%     plot(linspace(lim_d(i-1),lim_u(i-1),100),Lines(i,:),'k-','Linewidth',3); hold on;
end

 x=T(:,1);
P = polyfit(x,T(:,2),2);
m_T=linspace(0.02,0.08,1e2);
yfit = polyval(P,m_T);

T_up(1,:)=[];
 x_up=T_up(:,1);
P_up = polyfit(x_up,T_up(:,2),2);
m_T_up=linspace(0.125,0.165,1e2);
yfit_up = polyval(P_up,m_T_up);


hold on;  plot(m_T,yfit,'k-','Linewidth',3);
% hold on;  plot(T(:,1),T(:,2),'kx','Linewidth',3);
hold on;  plot(m_T_up,yfit_up,'k-','Linewidth',3);



N_tau=6;
yd=polyval(P,linspace(lim_d(1),lim_d(N_tau),1e6));
yu=polyval(P_up,linspace(lim_u(1),lim_u(N_tau),1e6));

hold on;  plot(linspace(lim_d(1),lim_d(N_tau),1e6),yd,'k-','Linewidth',3); 
hold on;  plot(linspace(lim_u(1),lim_u(N_tau),1e6),yu,'k-','Linewidth',3); 
hold on;  plot(linspace(lim_d(1),lim_u(1),1e6),yl(1,:),'k-','Linewidth',3); 
hold on;  plot(linspace(lim_d(N_tau),lim_u(N_tau),1e6),yl(N_tau,:),'k-','Linewidth',3); 


hold on;  plot(linspace(0,0.2,1e6),L_shape(i)*linspace(0,0.2,1e6),'k-','Linewidth',3); 

%%
xl=linspace(0.02,0.2,1e6);
L_shape_rel=L_shape(1:10);

for q=1:9
yl(q,:)=L_shape_rel(q+1)*linspace(lim_d(q),lim_u(q),1e6);
% yl(q,:)=L_shape_rel(q+1)*xl;

end

yd=polyval(P,linspace(lim_d(1),lim_d(end),1e6));
yu=polyval(P_up,linspace(lim_u(1),lim_u(end),1e6));

for i=1:9
[~,d]=min(abs(yl(i,:)-yd));
[~,u]=min(abs(yl(i,:)-yu));
lim_d(i)=xl(d);
lim_u(i)=xl(u);
end

[lim_d' lim_u']

%%



 hold on;   plot(Gather_Features(:,1),Gather_Features(:,2),'ro','Linewidth',2); 

 Fratures1_4=Gather_Features;
 
 
 



figure;
for q=1:52
    plot(F3(1e3*(q-1)+1:1e3*q,1),F3(1e3*(q-1)+1:1e3*q,2),'.'); hold on;
end


% save('Ntrain.mat','N')
All_res=unique(F3(:,3));
figure;
for q=1:length(All_res)
    id=find(F3(:,3)==All_res(q));
    plot(F3(id,1),F3(id,2),'.'); hold on;
    id=[];
end


for q=1:length(All_res)
    id=find(F3(:,3)==All_res(q));
    Tmp=F3(id,1:2);
    Nan_inds=find(isnan(Tmp(:,1)) | isnan(Tmp(:,2)));
    Tmp(Nan_inds,:)=[];
    [M(q),I]=min(Tmp(:,1));
    P(q)=Tmp(I,2);
    Mu(q,:)=[mean(Tmp(:,1)) mean(Tmp(:,2))];
    Sig(q,:)=[std(Tmp(:,1)) std(Tmp(:,2))];
    id=[]; Tmp=[];
end

InitP=[M' P'];

% Lambda=3*ones(19,1);
Lambda=linspace(2,5,52);
    Shapes=(Mu(:,2)-Lambda'.*(Mu(:,2)-Sig(:,2)))./(Mu(:,1)-Lambda'.*(Mu(:,1)-Sig(:,1)));
    InitP=[Mu(:,1)-Lambda'.*Sig(:,1) Mu(:,2)-Lambda'.*Sig(:,2)];

    
    
    lim1=0.05; lim2=0.23; lim3=0.21;
    
        p = polyfit(InitP(:,1),InitP(:,2),2);
        x2 = linspace(0.05,0.23,1e3);
        y2 = polyval(p,x2);
   figure; plot(InitP(:,1),InitP(:,2),'kx','Linewidth',2); hold on;
            plot(x2,y2,'r-','Linewidth',2); grid on;
        xlabel('m','FontSize', 14); grid on; ylabel('\sigma_{IPI}','FontSize', 14);       
        legend('Marginal points','Quadratic fit','FontSize', 14)
        
        
        ind=1
        x = linspace(0,0.25,1e5); 
        Par = @(x) 3.8;
        Quad = @(x) p(1)*x.^2 + p(2)*x + p(3);
        for ind=1:52
            Line = @(x) ms(ind)*x + r(ind);        
            f = @(x)Quad(x)-Line(x);
            [Intersections,I]=min(abs(f(x)));
            lims(ind)=x(I);
        end
        figure; plot(lims,'o--')
        
         
         xq=find(lims<0.02);
         vx=find(lims>0.02);
         vy=lims(vx);
         vq = interp1(vx,vy,xq)
        lims(xq)=vq;
         
        figure; plot(lims,'o--')

        
        
        figure;
        plot(F3(:,1),F3(:,2),'g.'); hold on;
        plot(SW(:,1),SW(:,2),'r.'); hold on;
        plot(x2,y2,'k-','Linewidth',3); ylabel('\sigma_{IPI}','Fontsize',12); xlabel('m'); hold on; grid on;
        plot(x3,y3,'k-','Linewidth',3); ylabel('\sigma_{IPI}','Fontsize',12); xlabel('m'); hold on; grid on;
        plot(x4,y4,'k-','Linewidth',3); ylabel('\sigma_{IPI}','Fontsize',12); xlabel('m'); hold on; grid on;
         xlabel('m','Fontsize',14); 
        for i=1:5:20
        ind=52;
        x1 = linspace(lims(i),0.6,1e3);
        y1=ms(i)*x1;
        x2 = linspace(lims(i),lim2,1e3);
        y2 = Quad(x2);
        x3 = linspace(lim3,0.6,1e3);
        y3 = Line(x3);        
        x4 = linspace(lim3,lim2,1e3);
        y4 = Par(x4)*ones(1,length(x4)); 
        plot(x1,y1,'-','Linewidth',3); ylabel('\sigma_{IPI}','Fontsize',14); xlabel('m','Fontsize',14); hold on; grid on;                
        end



legend('Noise','Sperm whales','','','','N_{\tau}=7','N_{\tau}=12','N_{\tau}=17','N_{\tau}=22','FontSize', 14)


figure; subplot(3,2,1); plot(All_res,Mu(:,1),'-.','Linewidth',2); ylabel('\mu_{m}','Fontsize',12); xlabel('N_{res}');
 subplot(3,2,2); plot(All_res,Mu(:,2),'-.','Linewidth',2); ylabel('\mu_{\sigma_{IPI}}','Fontsize',12); xlabel('N_{res}');grid on;
 subplot(3,2,3); plot(All_res,Sig(:,1),'-.','Linewidth',2); ylabel('\sigma_{m}','Fontsize',12); xlabel('N_{res}');grid on;
 subplot(3,2,4); plot(All_res,Sig(:,2),'-.','Linewidth',2); ylabel('\sigma_{\sigma_{IPI}}','Fontsize',12); xlabel('N_{res}');grid on;
 subplot(3,2,5); plot(All_res,Shapes,'-.','Linewidth',2); ylabel('Shape','Fontsize',12); xlabel('N_{res}');grid on;
subplot(3,2,6); plot(F3(:,1),F3(:,2),'g.'); hold on;
plot(InitP(5:end,1),InitP(5:end,2),'-.','Linewidth',2); ylabel('\sigma_{IPI}','Fontsize',14); xlabel('m'); hold on; grid on;
  

%%
figure;
% plot(F(:,1),F(:,2),'.'); hold on; 
plot(F2(:,1),F2(:,2),'.');
ylabel('\sigma_{IPI}','FontSize', 14);
xlabel('m','FontSize', 14); grid on;

in=find(F2(:,3)<58);
in2=find(F2(:,3)>58);

figure; h=histogram2(F2(in,1),F2(in,2),[100 100],'FaceColor','flat');
colorbar; hold on;

figure; h2=histogram2(F2(in2,1),F2(in2,2),[100 100],'Normalization','probability');
colorbar

figure; h=histogram2(F2(in,1),F2(in,2),[100 100],'Normalization','probability','DisplayStyle','tile','ShowEmptyBins','off');
colorbar; hold on;

All_res=round(linspace(7,60,52));
for i=1:52    
    res=All_res(i);
    res_inds=find(F3(:,3)==res);
    P_res(i)=length(res_inds)/size(F3,1);
    res_inds=[];
end

figure;
plot(All_res,P_res,'x-.','Linewidth',2);
 ylabel('P(N_{res})','FontSize', 14); xlabel('N_{res}','FontSize', 14); ylim([0 1]);

figure; semilogy(All_res,P_res,'x-.','Linewidth',2);
 ylabel('P(N_{res})','FontSize', 14); xlabel('N_{res}','FontSize', 14); ylim([0 1]);


save('P_res.mat','P_res')


% N=Nn;
figure;
% subplot(1,2,1);
Thresholds=linspace(0.1,0.95,20);
% Thresholds=0.6;
for j=1:length(Thresholds)
    Th=Thresholds(j);
    for i=1:20
        R=N(i,:)/All_res(i);
       FAs(i,j)=sum(R>Th)/length(R);
    end
    plot(All_res(1:20),FAs(:,j),'x-.','Linewidth',2); grid on; hold on;
    xlabel('N_{res}','FontSize', 14); ylabel('P( P_{fa} | N_{res} ) ','FontSize', 14); ylim([0 1]);
end
legend('\delta=0.6','FontSize', 14)


% subplot(1,2,2);
% for j=1:4
%    plot(Thresholds,FAs(j,:),'x-.','Linewidth',2); grid on; hold on;
%     xlabel('\delta'); ylabel('P_{fa}'); 
% end

for i=1:20
    Pfa(:,i)=FAs(1:i,:)'*P_res(1:i)';
end

figure;
semilogy(Thresholds,Pfa(:,1:20)*6*60*60,'-.','Linewidth',2); xlabel('\delta','Fontsize',14); ylabel('FA / hour','Fontsize',12); grid on;
legend('','N_{\tau}=8','N_{\tau}=9','N_{\tau}=10','N_{\tau}=11','N_{\tau}=12','N_{\tau}=13','FontSize', 14)



figure; yyaxis left
plot(All_res(1:1),Pfa,'o'); xlabel('N_{\tau}'); ylabel('P_{fa}'); grid on; yaxis([0 Pfa(end)]);
yyaxis right
plot(All_res(1:20),Pfa*6*60*60,'-.','Linewidth',2); xlabel('N_{\tau}'); ylabel('FAR [FA/hour]'); grid on;
yaxis([0 Pfa(end)*6*60*60]);


figure; plot(Thresholds,Pfa,'x-.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=[];
F=F2;

figure;
Reses=[10 20 30 40 50 60];
for i=1:4
    res_m=find(F(:,3)<Reses(i));
    subplot(2,2,i); plot(F(res_m,1),F(res_m,2),'.'); grid on; axis([0 0.5 0 4]); 
end

figure;
for j=1:4
res_m=find(F(:,3)==Reses(j));
%   figure;
  plot(F(res_m,1),F(res_m,2),'.'); grid on; axis([0 0.45 0 6]); hold on;
  ylabel('\sigma_{IPI}','FontSize', 14);
xlabel('m','FontSize', 14); grid on; 
end
legend('10','20','30','40','FontSize', 14)

m40=F(res_m,1);
Sig40= F(res_m,2);

figure;
subplot(2,2,1);
h1=histogram(m40); h1.Normalization = 'probability'; h1.BinWidth = 0.01;
grid on;  ylabel('Pdf','FontSize', 14); xlabel('m','FontSize', 14); xlim([0.1 0.4]);

subplot(2,2,2);
h1=histogram(Sig40); h1.Normalization = 'probability'; h1.BinWidth = 0.08;
grid on;  ylabel('Pdf','FontSize', 14); xlabel('\sigma_{IPI}','FontSize', 14); xlim([1 4]);

subplot(2,2,3);
xm=linspace(0.1,0.4,1e3);
pdm = fitdist(m40,'Normal'); Pdm_est = normpdf(xm,pdm.mu,pdm.sigma);
plot(xm,Pdm_est/max(Pdm_est)); xlim([0.1 0.4]);
grid on;  ylabel('Pdf','FontSize', 14); xlabel('m','FontSize', 14);
subplot(2,2,4);
xs=linspace(1,4,1e3);
pds = fitdist(Sig40,'Normal'); Pds_est = normpdf(xs,pds.mu,pds.sigma);
plot(xs,Pds_est/max(Pds_est)); xlim([1 4]);
grid on;  ylabel('Pdf','FontSize', 14); xlabel('\sigma_{IPI}','FontSize', 14);



figure;
Lines=[10 20 30 40 50 60];
for i=1:6
%     res_m=find(F(:,3)==Lines(i));
    res_m=find(F(:,3)==21);

    m40=F(res_m,1);
    Sig40= F(res_m,2);
    xm=linspace(0.05,0.55,1e3);
    pdm = fitdist(m40,'Normal'); Pdm_est = normpdf(xm,pdm.mu,pdm.sigma);
    xs=linspace(0,7,1e3);
    pds = fitdist(Sig40,'Normal'); Pds_est = normpdf(xs,pds.mu,pds.sigma);
    M(2)=(pds.mu-pds.sigma)/(pdm.mu-pdm.sigma);
    Xm=linspace(0,0.5,1e3);
    C = cov(F(res_m,1:2));
    Ys=0.5*(M(1)+M(2))*Xm;
     plot(F(res_m,1),F(res_m,2),'.'); grid on; axis([0 0.45 0 7]); hold on;
     hold on; plot(Xm,Ys,'k-','Linewidth',3); grid on; axis([0 0.5 0 4]); 
end
% legend('','15','','20','','25','','30','','35','','60')
legend('10','20','30','40','50','60','FontSize', 14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; h = histogram2(F(:,1),F(:,2),[100 48],'Normalization','probability')
xlabel('m')
ylabel('\sigma_{IPI}')
h.FaceColor = 'flat';
colorbar



F=F2;
%%
DELTAS=[10 13 16 19];
THR=linspace(0.1,0.9,10);
for qq=1:4
    DEL=DELTAS(qq);
    for q=1:10
    SW=SW_All; SW_train3=SW_train_All;
    SW_All=SW; SW_train_All=SW_train3;
%     EL2=find(SW(:,1)>0.32);
%     SW(EL2,:)=[]; 
%     SW_train(EL2)=[];
    SW_All2=SW;

    TH=THR(q);


    Eliminate=[];
    Eliminate=find(SW_train3<TH); 
    SW(Eliminate,:)=[]; 

    gamma=2.5;

        res_m=find(F(:,3)==DEL);
        m40=F(res_m,1);
        Sig40= F(res_m,2);
        xm=linspace(0,2,1e3);
        pdm = fitdist(m40,'Normal'); Pdm_est = normpdf(xm,pdm.mu,pdm.sigma);
        xs=linspace(0,18,1e3);
        pds = fitdist(Sig40,'Normal'); Pds_est = normpdf(xs,pds.mu,pds.sigma);
        M(1)=(pds.mu-pds.sigma)/(pdm.mu-pdm.sigma);
        Mem(1,:)=[pdm.mu pds.mu];

    figure;
    plot(F(:,1),F(:,2),'.');

    Lines=[11:60];
    for i=1:50
        res_m=find(F(:,3)==Lines(i));
        if ~isempty(res_m)
            m40=F(res_m,1);
            Sig40= F(res_m,2);
            xm=linspace(0,2,1e3);
            pdm = fitdist(m40,'Normal'); Pdm_est = normpdf(xm,pdm.mu,pdm.sigma);
            xs=linspace(0,18,1e3);
            pds = fitdist(Sig40,'Normal'); Pds_est = normpdf(xs,pds.mu,pds.sigma);
            M(i+1)=(pds.mu-pds.sigma)/(pdm.mu-pdm.sigma);
            Mem(i+1,:)=[pdm.mu pds.mu];
            Boundry(i)=norm(Mem(i+1,:) - Mem(i,:));
%             M2(i,:)=[pds.mu-(gamma+0.03*i)*pds.sigma pdm.mu-(gamma+0.03*i)*pdm.sigma];
            M2(i,:)=[pds.mu-(gamma)*pds.sigma pdm.mu-(gamma)*pdm.sigma];

            Xm=linspace(0,2,1e3);
        %     C = cov(F(res_m,1:2));
            Ys=0.5*(M(i)+M(i+1))*Xm;
    %          plot(F(res_m,1),F(res_m,2),'.'); grid on; axis([0 0.5 0 6]); 
             hold on; 
             if i==1 
             Ys1=Ys;
             end
        end
    end

            Xm=linspace(0,2,1e3);
            Xm_pol=Xm;
            plot(Xm,Ys1,'k-','Linewidth',3); grid on; axis([0 2 0 18]);


    El=find(M2(:,1)==0); 
    M2(El,:)=[];  
        p = polyfit(M2(:,2),M2(:,1),3);
        x1 = linspace(0,2,1e3);
        y1 = polyval(p,x1);
%         figure; plot(M2(:,2),M2(:,1),'kx','Linewidth',2)
%         hold on; plot(x1,y1,'r-','Linewidth',2)
        xlim([0.11 0.2]);
    [val,lim1]= min(abs(y1-Ys1(1)));

        res_m=find(F(:,3)==60);
        m40=F(res_m,1);
        Sig40= F(res_m,2);
        xm=linspace(0,2,1e3);
        pdm = fitdist(m40,'Normal'); Pdm_est = normpdf(xm,pdm.mu,pdm.sigma);
        xs=linspace(0,18,1e3);
        pds = fitdist(Sig40,'Normal'); Pds_est = normpdf(xs,pds.mu,pds.sigma);
        M(50)=(pds.mu-(gamma+0.02*i)*pds.sigma)/(pdm.mu-(gamma+0.02*i)*pdm.sigma);
    %     M(50)=(pds.mu-(gamma)*pds.sigma)/(pdm.mu-(gamma)*pdm.sigma);

        Xm=linspace(0,2,1e3);
        Ys=M(i)*Xm+0.035*Boundry(i);
        plot(Xm,Ys,'k-','Linewidth',3); grid on; axis([0 2 0 18]); 
        Ys=4.13;
        [val2,lim2]= min(abs(y1-Ys(1)));

    %     hold on; plot(M2(2:end-1,2),M2(2:end-1,1),'k-','Linewidth',3)
        hold on; plot(x1(lim1:lim2),y1(lim1:lim2),'k-','Linewidth',3)
        hold on; plot([0.227 x1(lim2)] ,[y1(lim2) 4.13],'k-','Linewidth',3)

            hold on; plot(SW(:,1),SW(:,2),'r.','Linewidth',7)



        id1=find(SW(:,1)>0.07 & SW(:,1)<0.25);
        Y1_est1=M(1)*SW(id1,1);
        Y2_est=polyval(p,SW(id1,1));

        id2=find(SW(:,1)>0.25);
        Y1_est2=M(1)*SW(id2,1);
        Y3_est=M(i)*SW(id2,1)+0.035*Boundry(i);


        term11=find(SW(id1,2)>Y1_est1);
        term2=find(SW(id1,2)<Y2_est);
        term12=find(SW(id2,2)>Y1_est2);
        term3=find(SW(id2,2)<Y3_est);

        IT1=intersect(term11,term2);
        IT2=intersect(term12,term3);

        Pd(qq,q)=(size(SW,1)-length([IT1' IT2']))/size(SW_All2,1);
    end   
end
    figure; plot(THR,Pd,'x')
 %% 
    
    
    figure;
% subplot(1,2,1);
Thresholds=linspace(0.1,0.9,10);
% Thresholds=0.6;
for j=1:length(Thresholds)
    Th=Thresholds(j);
    for i=1:30
        R=N(i,:)/All_res(i);
       FAs(i,j)=sum(R>Th)/length(R);
    end
    plot(All_res(1:30),FAs(:,j),'x-.','Linewidth',2); grid on; hold on;
    xlabel('N_{res}','FontSize', 14); ylabel('P( P_{fa} | N_{res} ) ','FontSize', 14); ylim([0 1]);
end
legend(['\delta=' num2str(round(Thresholds(1),2))],['\delta=' num2str(round(Thresholds(2),2))],['\delta=' num2str(round(Thresholds(3),2))],['\delta=' num2str(round(Thresholds(4),2))],['\delta=' num2str(round(Thresholds(5),2))],['\delta=' num2str(round(Thresholds(6),2))],['\delta=' num2str(round(Thresholds(7),2))],['\delta=' num2str(round(Thresholds(8),2))],['\delta=' num2str(round(Thresholds(9),2))],['\delta=' num2str(round(Thresholds(10),2))],'FontSize', 14)



for i=1:4
    Pfa(i,:)=FAs(1:(DELTAS(i)-9),:)'*P_res(1:(DELTAS(i)-9))';
end

    figure; 
    subplot(1,2,1); plot(Pfa'*6*60*60,Pd','-o','Linewidth',2); grid on;
    ylabel('Pd'); xlabel('FAR [FA/hour]'); grid on;
    legend('N_{\tau}=10','N_{\tau}=13', 'N_{\tau}=16', 'N_{\tau}=19','FontSize', 14);
    subplot(1,2,2); semilogx(Pfa'*6*60*60,Pd','-o','Linewidth',2); grid on;
    ylabel('Pd'); xlabel('FAR [FA/hour]'); grid on;
    legend('N_{\tau}=10','N_{\tau}=13', 'N_{\tau}=16', 'N_{\tau}=19','FontSize', 14);


figure; yyaxis left
plot(All_res(1:1),Pfa,'o'); xlabel('N_{\tau}'); ylabel('P_{fa}'); grid on; yaxis([0 Pfa(end)]);
yyaxis right
plot(All_res(1:1),Pfa*6*60*60,'o'); xlabel('N_{\tau}'); ylabel('FAR [FA/hour]'); grid on;
yaxis([0 Pfa(end)*6*60*60]);



%%

Sel=[1 2 3 4 5 8 15];

figure; 
subplot(1,2,1);
    plot(Pfa(:,Sel)*6*60*60,100*Pd(Sel,:)','-o','Linewidth',2); grid on;
    xlabel('FA / hour','FontSize', 14); ylabel('Pd [%]','FontSize', 14);
    legend('','N_{\tau}=8','N_{\tau}=9','N_{\tau}=10','N_{\tau}=11','N_{\tau}=15','N_{\tau}=22','FontSize', 14)

subplot(1,2,2);
    semilogx(Pfa(:,Sel)*6*60*60,100*Pd(Sel,:)','-o','Linewidth',2); grid on;
    xlabel('FA / hour','FontSize', 14); ylabel('Pd [%]','FontSize', 14);
legend('','N_{\tau}=8','N_{\tau}=9','N_{\tau}=10','N_{\tau}=11','N_{\tau}=15','N_{\tau}=22','FontSize', 14)


