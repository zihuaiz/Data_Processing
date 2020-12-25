%Brendon-This function adjusts for the nonlinearity in the piezo for each y
%trace individually

clear all; clc; close all;
fold = 'C:\Users\deleonlab\Documents\data_mat\';
fname = 'PLtlb1piezoitr_43399.mat'; 
% The scans and regions for average. Usually the number of scan and scan
% points.
subDomAve = 1:4;
peakDomain = 1:400; 
FreqAbs = 36151;
dat = load([fold fname]);
%%
PL = dat.pl';
ref = dat.ref';
max_yval = length(dat.yvals);
%dat.yvals = 1:10;
dat.yvals = 1:max_yval;
periodAdjust = 1.5*[0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0];
cleanSubIndex = 1;  %starting index (frequency) for subarray of scan where the number of peaks is constant
peakThresh = 0.003;
peakHeight = 0.005;
delta = zeros(length(dat.yvals),1);
for ii=1:length(dat.yvals)
    [pks,locs] = findpeaks(-ref(ii,cleanSubIndex:end),dat.xvals(cleanSubIndex:end),'MinPeakHeight',peakHeight);
    delta(ii,1) = locs(1); 
end


figure

delta = delta-delta(1);
deltaS = smooth(delta,20);
plot(dat.yvals,deltaS,'r-',dat.yvals,delta,'k-')
xlabel('Iteration');
ylabel('delta');

%%%%%%%%now correct the plots for the erorrs
PLEAdjust = zeros(length(dat.yvals),length(dat.xvals));
refAdjust = zeros(length(dat.yvals),length(dat.xvals));
for ii=1:length(dat.yvals)
    xTemp = dat.xvals-delta(ii);
    PLEAdjust(ii,:) = interp1(xTemp,PL(ii,:),dat.xvals);
    refAdjust(ii,:) = interp1(xTemp,ref(ii,:),dat.xvals);
end

for ii=1:length(dat.yvals)  
    if ii==1  %fit nonlinearity of piezo in this range with first sweep  
        [pks,locs] = findpeaks(-refAdjust(ii,cleanSubIndex:end),'MinPeakHeight',peakHeight);
        pksI = 1.5*(1:length(locs));
        pksI = pksI + 1.5/(locs(2)-locs(1))*locs(1) - 1.5;
      freq = interp1(locs,pksI,1:length(dat.xvals),'linear','extrap'); %interpolation including endpoints

    else
      [pks,locs] = findpeaks(-refAdjust(ii,cleanSubIndex:end),'MinPeakHeight',peakHeight);  
      pksI = 1.5*(1:length(locs));
      pksI = pksI + 1.5/(locs(2)-locs(1))*locs(1) - 1.5;
      freqTemp = interp1(locs,pksI,1:length(dat.xvals),'linear','extrap'); %interpolation including endpoints
      PLEAdjust(ii,:) = interp1(freqTemp+periodAdjust(ii),PLEAdjust(ii,:),freq,'linear');
      refAdjust(ii,:) = interp1(freqTemp+periodAdjust(ii),refAdjust(ii,:),freq,'linear');
    end
end

% figure
% subplot(2,1,1)
% pcolor(dat.xvals,dat.yvals,PL);
% shading flat
% shading interp
% %view(0,90)
% xlabel('Frequency (V)')
% ylabel('Iteration')
% %fprintf('Max Count Rate = %.3f',max
% subplot(2,1,2)
% pcolor(dat.xvals,dat.yvals,ref);
% shading flat
% shading interp
% xlabel('Frequency (V)')
% ylabel('Iteration')
%dat.yvals = 1:max_yval;
figure
subplot(3,1,1)
pcolor(dat.xvals(cleanSubIndex:end),dat.yvals,ref(dat.yvals,cleanSubIndex:end));
shading flat
%shading interp
title('Raw data with clean SubIndex')  
xlabel('Frequency (V)')
ylabel('Iteration')
subplot(3,1,2)
subIndex = 1; %index (iteration) to plot 1D trace of 
[pks,locs] = findpeaks(-ref(subIndex,cleanSubIndex:end),dat.xvals(cleanSubIndex:end),'MinPeakHeight',peakHeight);
plot(dat.xvals(cleanSubIndex:end),ref(subIndex,cleanSubIndex:end),locs,-pks,'r+')
title('Ref signal peak finding check')
xlabel('Frequency (V)')
ylabel('Voltage')
subplot(3,1,3)
title('Corrected axis')
plot(freq,ref(subIndex,:),'k-')
xlabel('Frequency (GHz)');


% %%%%%%%%now correct the plots for the erorrs
% PLEAdjust = zeros(length(dat.yvals),length(dat.xvals));
% refAdjust = zeros(length(dat.yvals),length(dat.xvals));
% for(ii=1:length(dat.yvals))
%     xTemp = dat.xvals-delta(ii);
%     PLEAdjust(ii,:) = interp1(xTemp,PL(ii,:),dat.xvals);
%     refAdjust(ii,:) = interp1(xTemp,ref(ii,:),dat.xvals);
% end

figure(89)
subplot(3,1,1)
pcolor(freq,dat.yvals,PLEAdjust);
shading flat
%shading interp
%view(0,90)
title('PLE scans with correction')
xlabel('Frequency (GHz)')
ylabel('Iteration')

subplot(3,1,2)
pcolor(freq,dat.yvals,refAdjust);
shading flat
shading interp
title('Corrected cavity signal')
xlabel('Frequency (GHz)')
ylabel('Iteration')
subplot(3,1,3)
indexP = 1;
title('First scan')
plot(freq,PLEAdjust(indexP,:))
xlabel('Frequency (GHz)')
ylabel('Intensity')


%%
PLEAdjust = PLEAdjust'; %transpose again to be consistent with other array conventions (in python)
PLEAdjust(isnan(PLEAdjust)) = 0;%basically we lose some of the domain from the adjustment and setting nan to 0

refAdjust = refAdjust';
refAdjust(isnan(refAdjust)) = 0; 


yvals = dat.yvals;


PLEAverage = zeros(length(freq),1);
for ii=subDomAve
   PLEAverage = PLEAdjust(:,ii)+PLEAverage;
end
PLEAverage = PLEAverage/length(subDomAve);
size(PLEAverage)
%%%%%%%%%%%%%%%%%%%No Fit%$$$$$$$$$$$$$$$$$$$
figure(87)
plot(freq(peakDomain),PLEAverage(peakDomain))
title('Average PLE scan')
xlabel('Frequency (GHz)');
ylabel('cps')

 
data.xvals = freq';
data.pl = PLEAdjust;
data.yvals = yvals';
data.ave = PLEAverage;
save(strcat(fname(1:end-4),'corrected','.mat'),'-struct','data');


