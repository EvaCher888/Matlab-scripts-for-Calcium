function [pow]=Analyse_FFT(x,y,pow,autoLevel)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - Variables that can be changed - - - - - -  %

% Set the accepted time error in seconds for each sample
% This variable checks that the time steps are evenly spaced
timeError = 0.25;

% Set the lowest accepted frequency in periods to be analyzed
% This variable excludes frequencies lower than lowPeriod periods
lowPeriod = 1.0;

% Set the level for auto-storage of the strongest peaks
% This variable stores peaks with powers greater than autoLevel 
% of the highest/strongest peak.
%autoLevel = 0.75;

% Set increment of the power spectrum 
% This variable determines how many values that will be saved in 
% the power spectrum file. If set to 1 all values are saved
pow.increment = 1;



        
%         pow.counter = i_emd; 
%   fprintf('%5i (%i)\r',pow.counter-1,pow.noCells);
  
  % Extract single cell calcium recording to analyze
  
  %y =detrend(y);%.0*hann(256);
  
  % -*- Key -*-
  % Trend elimination and centering of graph
%   if trendDegree > 0
%     pow.trend = polyfit(xdata,y,trendDegree);
%     y         = y - polyval(pow.trend,pow.timeVec);
%   end
  % -*- Key -*-
  % Fast Fourier Transform
  Y = fft(y,pow.N);
  
  % -*- Key -*-
  % Power spectral density (the Power Spectrum)
  % The factor 2 compensates for missing negative values
  P            = Y.*conj(Y)/pow.N;
  P            = P(1:pow.N/2+1);
  P(2:pow.N/2) = 2*P(2:pow.N/2);
  pow.B        = [pow.B P];
  
  % -*- Key -*_
  % Determine all frequency peaks
  dP        = gradient(P);
  allPeaksP = find(dP(1:end-1)>0 & dP(2:end)<0);
  for peakN = 1:length(allPeaksP)
    if P(allPeaksP(peakN))<P(allPeaksP(peakN)+1)
      allPeaksP(peakN) = allPeaksP(peakN)+1;
    end
  end
  allValleysP = find(dP(1:end-1)<0 & dP(2:end)>0);
  for valleyN = 1:length(allValleysP)
    if P(allValleysP(valleyN))>P(allValleysP(valleyN)+1)
      allValleysP(valleyN) = allValleysP(valleyN)+1;
    end
  end
  allValleysP = [1; allValleysP];
  pow.peaks   = [pow.freqVec(allPeaksP)' P(allPeaksP)];

  % Exclude frequencies lower than lowPeriod periods  
  
  
  timeSpan  = x(end)-x(1);
  pow.peaks = pow.peaks(find(pow.peaks(:,1)>lowPeriod*1e3/timeSpan),:);
  
  % -*- Key -*_
  % Store analyzed data in pow.peaks with three columns
  % Frequency PSDheight RelativePower(in %)
  
  areaP = trapz(pow.freqVec,P);
  for j = 1:size(pow.peaks,1)
    [tmp,I] = sort(abs(pow.freqVec(allValleysP)-pow.peaks(j,1)));
    if length(allValleysP) > 1
      lowLim  = min(allValleysP(I(1:2)));
      highLim = max(allValleysP(I(1:2)));
    else
      lowLim  = 1;
      highLim = length(pow.freqVec);
    end
    pow.peaks(j,3) = 100*trapz(pow.freqVec(lowLim:highLim),P(lowLim:highLim))/areaP;
  end

  % Sort peaks depending on relative power or PSD height

    pow.peaks = flipud(sortrows(pow.peaks,2));

  
  % Automatically stores all the highest peaks determined by autoLevel
 
      %pow.noPeaks = size(find(pow.peaks(:,2)>autoLevel*pow.peaks(1,2)),1);
           pow.noPeaks = size(find(pow.peaks(:,2)>autoLevel),1);
  
  % Reduce the pow.peaks matrix to pow.noPeaks
  pow.peaks = pow.peaks(1:pow.noPeaks,:);

  % Store the most dominant frequency for mean calculation
  %pow.freq.vec = [pow.freq.vec pow.peaks(1,1)];
  
  