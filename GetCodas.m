function Keep=GetCodas(Intensities)

    input = sort(Intensities,'descend');
    
    for loop=1:4

        classes = 2;
        for i = 1 : classes-1
        if i == 1
        data = input;
        elseif i > 1
        data = remaining_elements;
        end
        total = length (data);
        [~, GF] = get_jenks_interface(data);
        [~, I1] = max(GF);
        sub_array{i} = data(I1+1:total);
        remaining_elements = data (1:I1);
        end
        output = vertcat({data(1:I1)}, flipud(sub_array'));
        output{:};

        Break(loop)=output(1);    
        input=cell2mat(output(2));

    end


    for i=1:3
        G=cell2mat(Break(i));
        for j=1:length(G)
           [~,in(j)]= min(abs(Intensities-G(j)));
        end
        Break_inds(i)={in};
        in=[];
    end
    
    Keep=sort([cell2mat(Break_inds(1)) cell2mat(Break_inds(2)) cell2mat(Break_inds(3))]);
end