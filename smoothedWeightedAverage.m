function whaleOut = smoothedWeightedAverage(whaleIn, twin)
whaleOut = whaleIn;
spd = 60*60*24;

for wn = 1:numel(whaleIn)

    %initialize wloc estimate and velocity estimate:
    whaleOut{wn}.wlocSmooth = nan(size(whaleOut{wn}.wloc));
    whaleOut{wn}.velocityEst = whaleOut{wn}.wlocSmooth;

    Iuse = find(~isnan(whaleIn{wn}.wloc(:,1)));
    if ~isempty(Iuse)
        weights = sum(~isnan(whaleIn{wn}.TDOA(Iuse, :)), 2);

        v = zeros(length(Iuse), 3);
        wlocEst = v;

        % iterate through localized detections and calculate weighted
        % average for new whale position estimate:
        for i = 1:length(Iuse)
            %             Ind = max([i-numPoints, 1]):min([i+numPoints, length(Iuse)]); % indices to use in weighted line fitting

            % find indices of points within time window twin:
            Ind = find(abs(whaleIn{wn}.TDet(Iuse) - whaleIn{wn}.TDet(Iuse(i)) )<=twin./spd);

            for idim = 1:3 % iterate through each dimension and calculate weighted average whale location

                wi = whaleIn{wn}.wloc(Iuse(Ind), idim); % whale locations in this dimension and window
                wi0 = sum(wi.*weights(Ind))/sum(weights(Ind)); % weighted average
                wlocEst(i, idim) = wi0;
            end

        end
        
        v = diff(wlocEst)./(diff(whaleOut{wn}.TDet(Iuse)).*spd);

        whaleOut{wn}.wlocSmooth(Iuse, :) = wlocEst;
        whaleOut{wn}.velocityEst(Iuse(2:end), :) = v;
    end
end

end