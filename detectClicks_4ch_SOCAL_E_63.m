function [detTable] = detectClicks_4ch_SOCAL_E_63(tstart, tend, XH, H, c, paramFile)
% the typical detector, but with the ADCP pings removed
global detParam
loadParams(paramFile)

spd = 60*60*24;

t1 = tstart;
t2 = t1 + detParam.twin/spd;

if isfield(XH, 'deploymentName')
    wbtext = ['Running detector on ', XH.deploymentName];
else
    wbtext = 'Running detector';
end
wb = waitbar(0, wbtext);

% design filter:

if length(detParam.fc)==1
    [b, a] = ellip(4,0.1,40,detParam.fc*2/detParam.fs,'high');
else
    [b, a] = ellip(4,0.1,40,detParam.fc.*2/detParam.fs);
end

% Add echosounder detector
[bping, aping] = ellip(4,0.1,40,[20e3 30e3].*2/detParam.fs)
pingThresh = 70;
detTable = table;

idet = 1; % counter for number of detections

while t2<=tend
    [x, t] = readxwavSegment(t1, t2, XH);

    % detect clicks
    xf = filtfilt(b, a, x);

    if max(xf)>=detParam.th

        [pks, ind] = findpeaks(xf(:,1), 'minPeakHeight', detParam.th, 'minPeakDistance', detParam.minPkDist);

        tdet = t1 + ind/detParam.fs/spd; % times of detections

        % detect and remove ADCP pings:
        xfping = filtfilt(bping, aping, x(:,1));
        if max(xfping)<pingThresh
            [pingpks, pingind] = findpeaks(xfping, 'minPeakHeight', pingThresh);

            iping = find(diff(pingind)<100);
            trem = median(t(pingind(iping)));

            % since there shouldn't be more than 1 ADCP ping in a 30 sec window, I
            % can just say anything in between +/- .02 sec of mean point should be removed

            Irem = find(abs(tdet-trem)<(.02/spd));

            tdet(Irem) = [];
            ind(Irem) = [];
            pks(Irem) = [];
        end
        % calculate TDOA for each detection
        for i = 1:length(ind)
            i1 = max([1, ind(i) - detParam.maxdn]);
            i2 = min([length(xf), ind(i) + detParam.maxdn]);

            xclk = xf(i1:i2, :);
            [xc, lags] = xcov(xclk);

            tdoa = [0,0,0,0,0,0];
            xamp = tdoa;
            for pn = 1:length(detParam.xcRow) % iterate through each hydrophone pair

                % find 3 biggest peaks in xcov
                [xcpks, Nloc] = findpeaks(xc(:, detParam.xcRow(pn)), 'MinPeakDistance', 12, 'NPeaks', 3', 'SortStr', 'descend');

                if (0.8*xcpks(1))>xcpks(2) % largest peak is significantly bigger than 2nd largest
                    tdoa(pn) = lags(Nloc(1))/detParam.fs;
                    xamp(pn) = xcpks(1);
                else % largest peak is not much bigger than 2nd largest - reflection likely causing ambiguity
                    % sort peaks chronologically
                    [NlocSort, IND] = sort(Nloc, 'ascend');
                    tdoa(pn) = lags(NlocSort(2))/detParam.fs;
                    xamp(pn) = xcpks(IND(2));
                end

            end
            doa = H\(tdoa.'.*c);
            doa = doa./sqrt(sum(doa.^2));

            el = 180 - acosd(doa(3));
            az = atan2d(doa(2), doa(1));

            %         detTable.('XAmp')(idet, :) = xamp;
            %         detTable.('TDOA')(idet, :) = tdoa;
            %         detTable.('DOA')(idet, :) = doa.';
            %         detTable.('Ang')(idet, :) = [az, el];
            %         detTable.('Label')(idet) = '0';
            %         detTable.('DAmp')(idet:idet+length(pks)-1) = pks;
            %         detTable.('TDet')(idet:idet+length(pks)-1) = tdet;

            tempTable = table(tdet(i), pks(i), [az, el], doa.', tdoa, xamp, ...
                'VariableNames', {'TDet', 'DAmp', 'Ang', 'DOA', 'TDOA', 'XAmp'});

            detTable = [detTable; tempTable];

            idet = idet+1;
        end
    end
    calcPerc = (t2-tstart)/(tend-tstart);
    waitbar(calcPerc, wb, wbtext);

    t1 = t2;
    t2 = t1 + detParam.twin/spd;
end

close(wb)
detTable.('color') = 2.*ones(size(detTable.('TDet')));
detTable.('Label') = num2str(zeros(size(detTable.('TDet'))));