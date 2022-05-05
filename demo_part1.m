% close all force
% clear all
clear inDet DATA

%% DEMO, PART 1: functionality and variable receiver numbers

Nrec = 8;

N = 10000;
fs = 10000;

for i = 1:Nrec
    DATA{i}(1,:) = (1:N)./fs;
    DATA{i}(2,:) = randn(N,1).*10;

    I = find(DATA{i}(2,:)>30);
    inDet{i}(1,:) = DATA{i}(1, I);
    inDet{i}(2,:) = DATA{i}(2, I);
end

[outDet, labels] = brushDet(DATA, inDet, 'brushDet.params');