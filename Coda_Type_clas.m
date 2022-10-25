
function Coda_Type=Coda_Type_clas(Detected_pattern)

    data = readcell('Codas_Reference.xlsx');
    D=data(2:end,1:12);
    Analysis=data(2:end,3:11);
    Coda_type=data(2:end,12);

    ICI=zeros(size(data,1)-1,9);
    for i=1:size(Analysis,1) 
        NOC(i)=cell2mat(D(i,1));
        for j=1:size(Analysis,2)      
            ICI(i,j)=cell2mat(Analysis(i,j));       
        end    
    end


    tokeep = Detected_pattern ~= 0;   %logical array indicating which elements to keep
    workingP = Detected_pattern(tokeep);
    D_NOC = length(workingP);
    A_inds=find(NOC==D_NOC+1);

    Correspondance=[];
    for i=A_inds    
            Correspondance(i)=norm(ICI(i,:)-Detected_pattern);               
    end

    toset = Correspondance == 0;
    Correspondance(toset)=100;

    T=0.23;
    Type_ind=Correspondance==min(Correspondance);

    if min(Correspondance)<T
        Coda_Type=cell2mat(D(Type_ind,12));
    else
        Coda_Type=['Unseen'];
    end

end









