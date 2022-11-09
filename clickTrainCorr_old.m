function [CTC, DETout] = clickTrainCorr(DETin, whaleNum, kerInst, labeledInst, unlabeledInst, varargin)
% whale = clickTrainCorr(DETin, whaleNum, kerInst, labeledInst)
% whale = clickTrainCorr(DETin, whaleNum, kerInst, labeledInst, paramFile)
% Performs click train correlation (CTC) on DET tables.
%
% INPUTS:
% -DETin: struct of det tables, where DETin{1} is instrument 1's DET table,
% etc
% -whaleNum: whale number being used in CTC
% -kerInst: instrument(s) being used as kernel of CTC. If it's a vector,
% code will iterate over all instruments listed and perform CTC then
% associate clicks. To use both instruments 1 and 2 as kernels, use
% kerInst=[1,2]
% -labeledInst: the instruments with detections already assigned whale
% labels. For instance, if instrument 1 and 2 have already been labeled
% using brushDOA, use labeledInst=[1,2]
% -unlabeledInst: the unlabled instruments to be used in CTC
% -paramFile (OPTIONAL): the parameter file with variables to be used in CTC
%
% Outputs:
% -'CTC': a table of only the detections/TDOAs associated with each labeled whale
% - DETout: same as DETin but with the labels updated to reflect output of CTC

global CTCparam

% load in params
if nargin == 5 % param file specified
    loadParams(varargin{1})
else % no param file specified, load default file
    loadParams('CTCparams.txt')
end

% set variables:
spd = 60*60*24; % seconds per day, for converting between seconds and datenum

DETout = DETin; % assign output DET table to input one, labels will be changed as code performs CTC

% determine start/end time to perform CTC:
TSTART = nan;       % initialize start time of CTC
TEND = nan;         % initialize end time of CTC

ind1 = 1;           % index used to keep track of table concatenation
for ninst = 1:length(kerInst)
    
    if ~isempty(DETout{kerInst(ninst)})
        wnInd{ninst} = find(DETout{kerInst(ninst)}.color==whaleNum);
        NumDet = length(wnInd{ninst});
        
        % make a temporary table to concatenate with CTC
        tempCTC = table(nan(NumDet, 1), nan(NumDet,CTCparam.numInst), nan(NumDet,CTCparam.numTDOA), nan(NumDet,CTCparam.numInst),...
            'VariableNames', {'TDet', 'DAmp', 'TDOA', 'DETind'});
        tempCTC.TDet = DETout{kerInst(ninst)}.TDet(wnInd{ninst});
        tempCTC.DAmp(:, ninst) = DETout{kerInst(ninst)}.DAmp(wnInd{ninst});
        tempCTC.TDOA(:, CTCparam.indTDOA{kerInst(ninst)}) = DETout{kerInst(ninst)}.TDOA(wnInd{ninst}, :);
        tempCTC.DETind(:, ninst) = wnInd{ninst};

        ind2 = ind1 + NumDet - 1; 
        CTC(ind1:(ind1+NumDet-1), :) = tempCTC;
        ind1 = ind2+1;

        if ~isempty(wnInd)
            TSTART = min([TSTART, min(DETout{kerInst(ninst)}.TDet(wnInd{ninst}))]);
            TEND = min([TEND, max(DETout{kerInst(ninst)}.TDet(wnInd{ninst}))]);
        end
    end
end

% make CTC.detInd, which indicates the index of DET where this detection
% came from


% sort CTC by detection time
[~, Isort] = sort(CTC.TDet);
CTC = CTC(Isort, :);
% NOTE: CTC likely contains two separate rows for the same detection
% received on different instruments. If the code determines that they are
% the same detection, they will be combined into the same line later

for ndet = 1:length(CTC.TDet) % iterate through each detection and perform CTC
    thisArr = find(~isnan(CTC.DETind(ndet, :))); % which 4ch instrument was this detection found on
    otherArr = 3-thisArr; % other 4ch inst. 

    tstart = CTC.TDet(ndet) - CTCparam.twin/(2*spd); % beginning of window 
    tend = tstart + CTCparam.twin/spd; % end of window
    
    t = tstart:1/(CTCparam.fsct*spd):tend; % time vector of click trains
    X = zeros(length(t), CTCparam.numInst); % initialize click train of deltas
    Xh = zeros(length(t) + CTCparam.Nhann-1, CTCparam.numInst); % initialize click train of hanning windows
    for ninst = 1:length(kerInst)
        % find detections within time window on kerInst(ninst)
        Idet = find(~isnan(CTC.DETind(:, kerInst(ninst))) & ... % detections on inst. kerInst(ninst)
            CTC.TDet>=tstart & CTC.TDet<=tend); % detections within time window
        CTCInd{ninst} = Idet;
        % iterate through each detection and find the time in click train
        % which most closely matches detection time:
        for idet = 1:length(Idet)
            [~, Imatch] = min(abs(t-CTC.TDet(Idet(idet)))); % closest value in t to current detection
            X(Imatch, kerInst(ninst))=1; % set this index to 1
        end
        Xh(:, kerInst(ninst)) = conv(X(:, kerInst(ninst)), CTCparam.Wk); % convolve w/ hanning window
    end

    for ninst = 1:length(unlabeledInst)
        % find detections within time window:
        Idet = find(DETout{unlabeledInst(ninst)}.TDet>=tstart & DETout{unlabeledInst(ninst)}.TDet<=tend); 
           
        for idet = 1:length(Idet)
            [~, Imatch] = min(abs(t-DETout{unlabeledInst(ninst)}.TDet(Idet(idet)))); % closest value in t to current detection
            X(Imatch, unlabeledInst(ninst))=1; % set this index to 1
        end
        Xh(:, unlabeledInst(ninst)) = conv(X(:, unlabeledInst(ninst)), CTCparam.Wk); % convolve w/ hanning window
    end

    

    if length(labeledInst)==2 % both arrays are labeled
        [XC, lags] = xcorr(Xh(:,1), Xh(:,2), CTCparam.maxLag);
        [pks, Ipk] = findpeaks(XC, 'Npeaks', 2, 'SortStr', 'descend');
        if length(pks)==2
            if CTCparam.peakRatio*pks(1)>pks(2) % sufficiently unique peak in xcorr
                bestLag = lags(Ipk(1)); % delay which best aligns click trains

                % Find the closest detection from the other array to this
                % one
                [tdif, I] = min(abs(CTC.TDet(CTCInd{otherArr})-CTC.TDet(ndet) + bestLag/(CTCparam.fsct*spd))); 
                
                % if detections on two arrays are close enough to be
                % assumed the same click, combine them in CTC
                if tdif<=(CTCparam.tCloseEnough/spd)
                    CTC.DAmp(ndet, otherArr) = CTC.DAmp(CTCInd{otherArr}(I), otherArr);
                    CTC.DETind(ndet, otherArr) = CTC.DETind(CTCInd{otherArr}(I), otherArr);
                    CTC.TDOA(ndet, CTCparam.indTDOA{otherArr}) = CTC.TDOA(CTCInd{otherArr}(I), CTCparam.indTDOA{otherArr});
                    CTC.TDOA(ndet, CTCparam.indTDOA{3}(1)) = bestLag/CTCparam.fsct;
                    ok = 1;
                end

            end
        end

        for ninst = 1:length(unlabeledInst)
            if thisInst==1
                [XC, lags] = xcorr(Xh(:,1), Xh(:,unlabeledInst(ninst)), CTCparam.maxLag);
                [pks, Ipk] = findpeaks(XC, 'Npeaks', 2, 'SortStr', 'descend');
                if CTCparam.peakRatio*pks(1)>pks(2) % sufficiently unique peak in xcorr
                
                end
            else
                [XC, lags] = xcorr(Xh(:,2), Xh(:,unlabeledInst(ninst)), CTCparam.maxLag);
            end
        end
    elseif length(labeledInst)==1 && labeledInst(1) == 1 % only inst 1 is labeled
    
    elseif length(labeledInst)==1 && labeledInst(1) == 2 % only inst 2 is labeled

    end

    ok = 1;

end


% FOR TUESDAY:
% I think there's a faster way to manage the data. If I just keep all the
% DET files separated (don't combine 4ch's yet), I can iterate through each
% detection on Inst 1 and do whale association. I then have a flag that
% says whether a detection was already associated.

% like DET{3}.associated(idet) = 1 if this detection has been copied into
% CTC.

% Steps for a detection would be:
% 1. Do CTC for arr 1
%   a. Make click train (could be a function: feed in TDet and t and output
%      a click train? [Xct] = makeCT(TDet, t, CTCparam.Nhann))
%   b. iterate through other instruments, make click trains (using same
%      function)
%   c. xcorr click trains from other arrays with currentArray
%   d. peak ratio comparison to determine how well it found a unique
%      correlation
% 2. Associate clicks
%   a. if a click train corr passed test in 1.d., apply expected delay. 
%   b. if a click on other inst is close enough to click on current array
%      (after delay is applied), then assume it is the same click
%   c. Transfer relevant info from DET to CTC.
%   d. DET{instNo}.moved(I) = 1; 
%   e. CTC.TDOA(??) = bestLag/fsct;
% 3. Repeat for array 2
% 4. Remove redundant entries in CTC. If CTC.DETind(:,instNo) has repeat values
% (likely), determine which is more likely correct (based on adjacent
% TDOAs? or ratio of peak(2) to peak(1)? or number of clicks used in click
% train?)

% above algorithm requires fewer if statements and triple-layer for-loops, I
% think
