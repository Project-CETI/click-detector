
function [duration,f_centroid]=superlets2(Section,F_ds,fois)

Test=0; stop=0;
init_thresh=0.65;
srord   = [1, 3]; 
Clicks_bank= Section;   % pick region of analysis
time_pulse= [0:1/F_ds:(1/F_ds)*(length(Clicks_bank)-1)]'; 

xS_slt=aslt(Clicks_bank, F_ds, fois, 3, srord, 0);                % calculate spectogram 

M_val_sel=max(max(xS_slt));                                         % calculate the maximum value within the spectogram
slt_NORM=xS_slt./M_val_sel; 

while Test<0.05 && stop<100
    stop=stop+1;
    Binary=double(slt_NORM>init_thresh);                                    
    init_thresh=init_thresh-0.01;
    Binary(end,1:end)=zeros(1,size(Binary,2));     
    Binary(1,1:end)=zeros(1,size(Binary,2));
    Binary(1:end,end)=zeros(1,size(Binary,1));
    Binary(1:end,1)=zeros(1,size(Binary,1));
    Test=mean(mean(Binary));
end 

if Test>0.05
    m=iblobs(Binary);
    m_inds=m.area>90 & m.area<7000;
    M=m(m_inds);

    roi=find(M.vc==max(M.vc));

    duration=1e3*(M(roi).umax-M(roi).umin)/F_ds;         % calculate the pulse duration [ms]
    f_centroid=1e-3*fois(round(M(roi).vc))   ;                 % calculate the frequency centroid [khz]

    if length(duration)>1
        duration=0.005;
        f_centroid=2;
    end
else
        duration=0;
        f_centroid=2;
end
% figure;
% subplot(2,1,1);
%         imagesc(1e3*time_pulse, fois, xS_slt ); xlabel('t[ms]');
%         set(gca, 'ydir', 'normal');
%         colormap jet;
% subplot(2,1,2);
%         imagesc(1e3*time_pulse, fois, Binary ); xlabel('t[ms]');
%         set(gca, 'ydir', 'normal');
%         colormap jet;
end

