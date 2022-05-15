function [detTable] = detectClicks_1ch(tstart, tend, XH, paramFile)

global detParam
loadParams(paramFile)

spd = 60*60*24;
detParam.twintwin = 30; % 30 second window;

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

detTable = table;

idet = 1; % counter for number of detections



while t2<=tend
    [x, t] = readxwavSegment(t1, t2, XH);
    
    xf = filtfilt(b, a, x);

    if max(xf)>=detParam.th

        [pks, ind] = findpeaks(xf, 'minPeakHeight', detParam.th, 'minPeakDistance', detParam.minPkDist);

        tdet = t1 + ind/detParam.fs/spd; % times of detections
        
        tempTable = table(tdet, pks, 'VariableNames', {'TDet', 'DAmp'});
        
%         detTable.('DAmp')(idet:idet+length(pks)-1) = pks;
%         detTable.('TDet')(idet:idet+length(pks)-1) = tdet;
        detTable = [detTable; tempTable];

        idet = idet + length(pks);
        
    end
    
    calcPerc = (t2-tstart)/(tend-tstart);
    waitbar(calcPerc, wb, wbtext);

    t1 = t2;
    t2 = t1 + detParam.twin/spd;
end

close(wb)

detTable.('Label') = num2str(zeros(size(detTable.('TDet'))));
detTable.('color') = 2.*ones(size(detTable.('TDet')));