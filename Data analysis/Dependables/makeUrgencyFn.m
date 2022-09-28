%% Create urgency signal
% Oscillations in the higher Beta band (16:30) are taken as indicators of motor preperations
% and are use  to constrain our neurally-informed models, the grand-average Beta amplitude 
% waveform in the 8 sec ITI was taken (Calculated with dataAnalysis.m) and normalised according
% to the motor execution threshold indexed by average pre-response beta in a given condition,
% and a 2nd-order polynomial was then fitted to capture the smooth urgency trend relative to threshold
%
% NOTE for future data analysis or Neurally-informed modelling this is intergrated and adapted 
% to be used flexibiliy in dataAnalysis.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------    Threshold   -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first get grandAverageBeta response-locked data to determine the threshold.

load('grandAverageBetaERD.mat') % find in /Context-Dependent-Detection/Data analysis/Dependables/
timeSeries = grandAverageResponse.timeSerie; 
response = [grandAverageResponse.WeakConstant grandAverageResponse.StrongConstant grandAverageResponse.WeakMixed grandAverageResponse.StrongMixed];

figure; 
plot(timeSeries, response) 

for indCond = 1:4
    threshold(indCond) = response(find(timeSeries>=-67,1), indCond);
end

threshold(3) = mean(threshold(3:4)); threshold(4)=[]; % collapse across the two threshold levels for Mixed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------    ITI dynamics   --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First get grandAverageBeta target-locked data to created the urgency signal.
 load('grandAverageBetaERD_8secsITI.mat')

timeSeries = grandAverageERD.timeSeries;
betaITI = [grandAverageERD.WeakContext-threshold(1) grandAverageERD.StrongContext-threshold(2) grandAverageERD.MixedContext-threshold(3)];

figure; hold on;
plot(timeSeries, betaITI)

% remove NaNs and everything after the target-onset.
del = find(isnan(betaITI(:,1)) | timeSeries > 0); 
timeSeries(del) = [];
betaITI(del,:)  = [];

% Fit 2-order polynomial.
b = [];
for indCond = 1:3
    b(:,indCond) = polyfit(ti,betaITI(:,indCond),2);
    plot(timeSeries, b(1,indCond)*timeSeries.^2 + b(2,indCond)*timeSeries + b(3,indCond))
end


% Make a normalised urgency function then for the model:
f  = 60;
dt = 1/f;
timeSeries = [-7500:dt*1000:1500];
for indCond = 1:3
    urg(:,c) = b(1,indCond)*timeSeries.^2 + b(2,indCond)*timeSeries + b(3,indCond);
    urg(find(timeSeries>0), indCond) = urg(find(timeSeries>=0, 1),indCond);
end

urg = urg./max(max(urg));
urg = 1-urg;

figure; hold on
plot(timeSeries,urg)

% and finally just to play around with pink noise:
% urg = urg+pinknoise(size(urg))
% plot(timeSeries,urg)
