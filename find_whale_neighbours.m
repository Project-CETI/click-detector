function Detect_low=find_whale_neighbours(Y,Gather_TOA,Fs)

    %% tag decoding
    Detect_low=[];
    Detect = Gather_TOA;
    Power = zeros(1, length(Detect));
    for ind = 1: length(Power)
        Current = Y(round(Detect(ind)*Fs-0.1*Fs): round(Detect(ind)*Fs+0.1*Fs));
        Power(ind) = max(abs(Current));
    end
    Th = mean(Power) - 2*std(Power);
    if Th < mean(Power)/2
        Detect_low = Detect(find(Power < Th));
    %     for ind = 1: length(Detect_low)
    %         figure;
    %         plot(Y(round(Detect_low(ind)*Fs-0.1*Fs): round(Detect_low(ind)*Fs+0.1*Fs)));
    %     end
    else
        disp('No other whales');
    end

end