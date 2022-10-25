function Th_selected=Adapt_threshold(Locs,Pks)

% Locs=MP_t; Pks=MP_p;
    L=100; Th_selected=[];
    steps=linspace(min(Pks),max(Pks),L);
    edges=linspace(0,2,20);
    for i=1:L
        Th=steps(i);
        Select_inds=Pks>Th;
        Select=Locs(Select_inds);
        D=diff(Select);
        D_relevant=D(D>0.4);
        [Counts,edges] = histcounts(D_relevant,edges);
        N(i)=max(Counts);
    end

    Th_max=steps(N==max(N));
    Th_selected=Th_max(1);
    if max(N)<4
        Th_selected=max(Pks)+eps;
    end
    
end

