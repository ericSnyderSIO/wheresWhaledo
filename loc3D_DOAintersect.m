function whaleLoc = loc3D_DOAintersect(DET, hydLoc, paramFile)
% produces a rough 3-D track estimate from the period of time with
% overlapping detections on both four-channel arrays.
% whaleLoc is a struct with the 3-D tracks
% DET is a struct containing detection tables
% hydLoc is a struct with the hydrophone locations in lat, lon, z

global brushing
loadParams(paramFile)
spd = 60*60*24;

colorNums = unique([DET{1}.color; DET{2}.color]); % find all unique labels
colorNums(colorNums==2) = []; % remove unlabeled points

h0 = mean([hydLoc{1}; hydLoc{2}]);

% convert hydrophone locations to meters:
[h1(1), h1(2)] = latlon2xy_wgs84(hydLoc{1}(1), hydLoc{1}(2), h0(1), h0(2));
h1(3) = abs(h0(3))-abs(hydLoc{1}(3));

[h2(1), h2(2)] = latlon2xy_wgs84(hydLoc{2}(1), hydLoc{2}(2), h0(1), h0(2));
h2(3) = abs(h0(3))-abs(hydLoc{2}(3));

figure(31)
scatter3(h1(1), h1(2), h1(3), 24, 'k^', 'filled')
hold on
scatter3(h2(1), h2(2), h2(3), 24, 'k^', 'filled')

for wn = 1:length(colorNums) % iterate through each whale number
    I1 = find(DET{1}.color==colorNums(wn)); % indices on array 1 labeled as whale wn
    I2 = find(DET{2}.color==colorNums(wn)); % indices on array 2 labeled as whale wn

    if ~isempty(I1) && ~isempty(I2) % make sure there are detections on both arrays with this label

        t1 = DET{1}.TDet(I1); % times of detections on array 1
        t2 = DET{2}.TDet(I2); % times of detections on array 2

        % period of overlapping detections
        tstart = max([min(DET{1}.TDet(I1)), min(DET{2}.TDet(I2))]);
        tend = min([max(DET{1}.TDet(I1)), max(DET{2}.TDet(I2))]);

        ti = tstart:1/spd:tend;

        if ~isempty(ti) % will be empty if there aren't overlapping detections
            doa1 = DET{1}.DOA(I1, :);
            doa2 = DET{2}.DOA(I2, :);

            doa1i = interp1(t1, doa1, ti);
            doa2i = interp1(t2, doa2, ti);

            for i = 1:length(ti)
                D = [doa1i(i, :); -doa2i(i, :)];
                R = D.'\(h2-h1).';

                w1 = R(1).*doa1i(i,:) + h1;
                w2 = R(2).*doa2i(i,:) + h2;

                w(i, :) = mean([w1; w2]);
                werr(i) = sqrt(sum((w1-w2).^2));

            end
            whaleLoc{wn}.xyz = w;

            [lat, lon] = xy2latlon_wgs84(w(1, :), w(2, :), h0(1), h0(2));
            z = w(3,:) - abs(h0(3));
            whaleLoc{wn}.LatLonDepth = [lat; lon; z];

            scatter3(whaleLoc{wn}.xyz(:, 1), whaleLoc{wn}.xyz(:, 2), whaleLoc{wn}.xyz(:, 3), ...
                24, brushing.params.colorMat(wn+2, :), 'filled')
        end
    end
end
hold off

