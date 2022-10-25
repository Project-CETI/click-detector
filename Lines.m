
%% Up

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
% figure;
colorstflex = {'#FF0000','#00FF00','#0000FF','#00FFFF','#FF00FF','#FFFF00','#000000','#D95319','#7E2F8E','#77AC30','#A2142F'};

for i=2:length(A)
    inds=find(X(:,3)==A(i));
    mu_m=mean(X(inds,1)); s_m=std(X(inds,1));
    mu_sigma=mean(X(inds,2)); s_sigma=std(X(inds,2));
    Shapes(i)=(mu_sigma-s_sigma)/(mu_m-s_m);
    L_shape(i)=0.5*(Shapes(i)+Shapes(i-1));
%     T(i,:)=[mu_m-gamma(i)*s_m mu_sigma-gamma(i)*s_sigma];
%     plot(X(inds,1),X(inds,2),'x','color',cell2mat(colorstflex(i)),'Linewidth',2); hold on; axis([0 3.5 0 8]);
    plot(X(inds,1),X(inds,2),'.','Linewidth',2); hold on; axis([0 3.5 0 8]);
    inds=[];
%      legendInfo{i}=[num2str(i+A(1)-1), '-residual transients: '];
%      legend(legendInfo,'Fontsize',12);
end

y4=polyval(T4,linspace(lim2_u,lim2_d,1e2));
y5=T5*linspace(lim2_d,3.5,1e2);
y6=T6*linspace(lim2_u,1.3,1e2);
hold on;  plot(linspace(lim2_u,lim2_d,1e2),y4,'k-','Linewidth',3); 
hold on;  plot(linspace(lim2_d,3.5,1e2),y5,'k-','Linewidth',3); 
hold on;  plot(linspace(lim2_u,1.3,1e2),y6,'k-','Linewidth',3); 



%% Down

X=F3(I,:);
X2_inds=find(X(:,1)<0.3 & X(:,2)<0.7);
X2=X(X2_inds,:);
A2=unique(X2(:,3));
q2=[10:-1:1];
we=3.5*(1-exp(-q2*0.25));
% we=5.5*(1-exp(-q*0.05-1));
gamma=we;
gamma_up=2.5;
% figure;

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

N_tau=4;
yd=polyval(T0,linspace(lim_d(1),lim_d(N_tau),1e6));
yu=polyval(T3,linspace(lim_u(1),lim_u(N_tau),1e6));

hold on;  plot(linspace(lim_d(1),lim_d(N_tau),1e6),yd,'k-','Linewidth',3); 
hold on;  plot(linspace(lim_u(1),lim_u(N_tau),1e6),yu,'k-','Linewidth',3); 
hold on;  plot(linspace(lim_d(1),lim_u(1),1e6),yl(1,:),'k-','Linewidth',3); 
hold on;  plot(linspace(lim_d(N_tau),lim_u(N_tau),1e6),yl(N_tau,:),'k-','Linewidth',3); 


