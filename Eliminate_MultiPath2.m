function Eliminate_inds=Eliminate_MultiPath2(locs,pks)
   
    Eliminate_inds=[];
    
    Multi_Path=diff(locs);
    Multi_Path2=diff(pks);
    
     MPa=[Multi_Path ; Multi_Path2]; 
%     figure; plot(Multi_Path,Multi_Path2,'x')
%     MP_inds=find(Multi_Path<70e-3 && Multi_Path2<0);
     MP_inds=find(MPa(1,:)<70e-3 & MPa(2,:)>0);

    Seqs=diff(MP_inds);
    Decision=find(Seqs==1, 1);
    if length(Decision)<3
        Eliminate_inds=MP_inds+ones(1,length(MP_inds));
    end
    
    
end
    
    




