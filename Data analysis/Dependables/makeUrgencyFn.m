% Make urgency function
load('grandAverageBetaERD.mat')
tr = grandAverageResponse.timeSerie; 
betaResp = [grandAverageResponse.WeakConstant grandAverageResponse.StrongConstant grandAverageResponse.WeakMixed grandAverageResponse.StrongMixed];
figure; plot(tr,betaResp) 

for c=1:4
    th(c) = betaResp(find(tr>=-67,1),c);
end
th(3)=mean(th(3:4)); th(4)=[]; % collapse across the two threshold levels for Mixed.

load('grandAverageBetaERD_8secsITI.mat')

ti = grandAverageERD.timeSeries;
betaITI = [grandAverageERD.WeakContext-th(1) grandAverageERD.StrongContext-th(2) grandAverageERD.MixedContext-th(3)];
% betaITI = [grandAverageResponse.WeakConstant-th(1) grandAverageResponse.StrongConstant-th(2) nanmean([grandAverageResponse.WeakMixed grandAverageResponse.StrongMixed],2)-th(3)];

figure; plot(ti,betaITI)
legend
grid
del = find(isnan(betaITI(:,1)) | ti>0);
ti(del)=[];
betaITI(del,:)=[];
b = [];
for c=1:3
    b(:,c) = polyfit(ti,betaITI(:,c),2);
end

hold on
for c=1:3
    plot(ti,b(1,c)*ti.^2+b(2,c)*ti+b(3,c))
end

% Make a normalised urgency function then for the model:
f=60;
dt=1/f;
tm = [-7500:dt*1000:1500];
for c=1:3
    urg(:,c) = b(1,c)*tm.^2+b(2,c)*tm+b(3,c);
    urg(find(tm>0),c) = urg(find(tm>=0,1),c);
end

urg = urg./max(max(urg));
urg = 1-urg;

figure; hold on
plot(tm,urg)

% save urg urg tm

% and finally just to play around with pink noise:
urg = urg+pinknoise(size(urg))
plot(tm,urg)