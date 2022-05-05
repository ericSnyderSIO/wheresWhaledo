% Detector experiments

[detA] = detectClicks_4ch_SOCAL_E_63_reflectRem(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 65, 1, 'detClicks_4ch.params');

[detB] = detectClicks_4ch_SOCAL_E_63_reflectRem(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 65, .9, 'detClicks_4ch.params');

[detC] = detectClicks_4ch_SOCAL_E_63_reflectRem(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 65, .8, 'detClicks_4ch.params');

[detD] = detectClicks_4ch_SOCAL_E_63_reflectRem(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 65, .7, 'detClicks_4ch.params');


%%

figure(1)
subplot(4,1,1)
plot(detA.('TDet'), detA.('Ang')(:,1), '.')
subplot(4,1,2)
plot(detB.('TDet'), detB.('Ang')(:,1), '.')
subplot(4,1,3)
plot(detC.('TDet'), detC.('Ang')(:,1), '.')
subplot(4,1,4)
plot(detD.('TDet'), detD.('Ang')(:,1), '.')
