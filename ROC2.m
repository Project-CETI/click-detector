

% f=dir('Features_Phase*.mat');
% FP=[];
% f=dir('*FA*.mat');
% FA=[];
f=dir('Pd*.mat');
Pd=[];

for i=1:length(f)
   L=load([f(i).name]);
%    FA=[FA ; L.Gather_FA];
%    FP=[FP ; L.Gather_Features_phase];
   Pd=[Pd ; L.Gather_FA];

end

figure; plot(FP(:,1),FP(:,2),'.')

%% SW
Gd=zeros(9,20);
Ld=size(Pd,1)/9;
for i=1:Ld
   Gd=Gd+Pd(9*(i-1)+1:9*i,:); 
end
pd=Gd/Ld;
% pd=Gd/124;

%% Noise:
G=zeros(9,20);
L=size(FA,1)/9
for i=1:L
   G=G+FA(9*(i-1)+1:9*i,:); 
end
Gn=G/(L/360);

%% ROC
Th=linspace(0.1,0.95,20);

figure; 
for i=[1 2 3 9]
  plot(Gn(i,:),pd(i,:),'-*','Linewidth',2); hold on; grid on;
end
xlabel('FA/hour','Fontsize',14);
ylabel('Pd','Fontsize',14);
legend('N_{\tau}=13','N_{\tau}=12','N_{\tau}=11','N_{\tau}=5','Fontsize',12);

Pd_10=pd(9,:);
Pd_20=pd(9,:);
Pd_30=pd(9,:);
Pd_60=pd(9,:);

PDs=[Pd_10 ; Pd_20 ; Pd_30 ; Pd_60];

figure; 
for i=1:4
  plot(Gn(9,:),PDs(i,:),'-*','Linewidth',2); hold on; grid on;
end
xlabel('FA/hour','Fontsize',14);
ylabel('Pd','Fontsize',14);
legend('T_{buf}=10sec','T_{buf}=20sec','T_{buf}=30sec','T_{buf}=60sec','Fontsize',12);
xlim([0 20])

