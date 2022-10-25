
GT_all=xlsread('10M_ch2_0-5.xls');
GT=GT_all(:,2);


Det=TOA;
% Det=xlsread('10M_ch2_0-5.wav.xls');
tolerence=linspace(10e-3,100e-3,100);
for j=1:100
    tol=tolerence(j);
    TP=0; 
    for i=1:length(Det)
        [val,id]=min(abs(GT-Det(i)));
        if val<tol
            TP=TP+1;
        end
    end
    Pd(j)=TP/length(Det);
%     Pd(j)=TP/length(GT);

end

figure;
plot(1e3*tolerence,Pd,'rx-','Linewidth',2); grid on; hold on; xlabel('Tolaerance [msec]'); ylabel('Pd');
load Pd_phase
plot(1e3*tolerence,Pd_phase,'bx-','Linewidth',2); grid on; hold on; xlabel('Tolaerance [msec]'); ylabel('Pd');
load Pd_ATKEO
plot(1e3*tolerence,Pd_ATKEO,'gx-','Linewidth',2); xlabel('Tolaerance [msec]'); ylabel('Pd');
% axis([0 1 0 2]);


Pd_phase=Pd;
save('Pd_phase.mat','Pd_phase')

Pd_ATKEO=Pd;
save('Pd_ATKEO.mat','Pd_ATKEO')