function whaleOut = localize(whaleIn, hloc, H1, H2, drift, c, M)
% whaleOut = localize(whaleIn, hloc, H1, H2)
% takes the TDOAs from a "whale" struct (output of calcTDOAfromCTC) and
% performs iterative model generation to determine the most likely whale location.
% INPUTS:
% -whaleIn = the "whale" struct
% -hloc = hydrophone locations (m).
% -H1 and H2 = H matrices for 4ch instruments (m).
% -drift is a 1xNinst array of the drift for each instrument in hloc (s).
% -c = speed of sound (m/s).
% -M = base model (coarse grid). This will be used for the first iteration
% in localization, then each subsequent localization will be calculated

% Overall process:
% 1) Iterate through each detection
%    a) Use LMSE to determine most likely whale location based on TDOA in
%    courseModel 

whaleOut = whaleIn;

for wn = 1:numel(whaleOut)
    if ~isempty(whaleOut{wn})
        for ndet = 1:length(whaleOut{wn}.TDet)
            
        end
    end
end




end

function wloc = localize_1Det(TDOA, sig2, x, y, z, hloc, H1, H2)
LMSE = (M.TDOA - TDOA(ndet, :)).^2    
end

function [TDOA, wloc] = makeModel(xv, yv, zv, h, H1, H2, c)

% make wloc (matrix of whale positions)
[cx, cy, cz] = ndgrid(xv, yv, zv);
wloc = [cx(:), cy(:), cz(:)];

s1 = wloc-h(1, :);
r1 = sqrt(sum(s1.^2, 2)); % range to instrument 1
s1 = s1./r1; % direction vector to instrument 1

s2 = wloc-h(2, :);
r2 = sqrt(sum(s2.^2, 2)); % range to instrument 2
s2 = s2./r2; % direction vector to instrument 2

s3 = wloc-h(3, :);
r3 = sqrt(sum(s3.^2, 2)); % range to instrument 3

s4 = wloc-h(4, :);
r4 = sqrt(sum(s4.^2, 2)); % range to instrument 4

% small aperture TDOAs
TDOA(:, 1:6) = (s1*H1.')./c;
TDOA(:, 7:12) = (s2*H2.')./c;

% large aperture TDOAs
TDOA(:, 13) = (r1-r2)./c;
TDOA(:, 14) = (r1-r3)./c;
TDOA(:, 15) = (r1-r4)./c;
TDOA(:, 16) = (r2-r3)./c;
TDOA(:, 17) = (r2-r4)./c;
TDOA(:, 18) = (r3-r4)./c;

end