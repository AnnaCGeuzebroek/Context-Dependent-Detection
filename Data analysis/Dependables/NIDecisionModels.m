%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------  supporting functions     ----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Leaky accumulation  ---------------------------
function [err, simdat, pred, DecValue] = NIDecisionModels(obj, pm, parNames, urgency, noise, getCPP, reflectingBound)
% Simulates beh matrix for continuous leaky accumulator with certain given parameters
% It has 3 bound parameters (one for W one for S one for M) but single leak parameter and single non-decision time (tnd)
% assumed to be before decision process (think of it as a delay from changes on the screen to being counted in the accumulator)
% and 1 free drift rate parameter for lo coh scaled up for high coh
% datsum = obj.modelBehaviour.datsum;
% Simulates beh matrix for continuous leaky accumulator with certain given parameters

% It has 3 bound parameters (one for W one for S one for M) but single leak parameter and single non-decision time (tnd)
% assumed to be before decision process (think of it as a delay from changes on the screen to being counted in the accumulator)
% and 1 free drift rate parameter for lo coh scaled up for high coh

if ~exist('getCPP','var'); getCPP = 0; end
if ~exist('reflectingBound','var'); reflectingBound = 0; end

            
% pm are parameters
TOW    = obj.modelBehaviour.datsum.TOW;
bt     = obj.modelBehaviour.datsum.bt;
trialMatrix = obj.modelBehaviour.datsum.trialMatrix;

bound = []; drift = []; leak = []; tnd = []; boost = []; criteria = []; trialnoise = [];
for indParName = 1:size(parNames,1)
    eval(sprintf('%s(end+1) = pm(:, %i);', lower(parNames(indParName,:)), indParName))
end

if isempty(bound)
    bound = 1;
end


dt = 1/obj.stim.refreshRate; % sample period
maxRTms = obj.stim.RTdeadLine(end); % the upper limit for RTs to be classed as hits
minRTms = obj.stim.RTCutOff; % minium allowable RT to call a response a 'hit', in msec. There was only one RT in the S blocks across all subjects that was shorter than this, and next shortest was nearly 100 ms later

postRespPause = 1*obj.stim.refreshRate; % don't start accumulating until this number of samples after response. This is simply a reasonable guess

if size(drift,2) == 1
    % if only one drift rate is fitted we are going to assume that the
    % drift rate scale depending on the coherence levels.
    
    % however, this is how it is implemented now, making it more
    % transfable (for example with more coherence levels):
    drift = drift*(obj.modelBehaviour.cohs./min(obj.modelBehaviour.cohs)); % scale up the drift rate parameters for other conditions
elseif size(drift,2) > length(obj.modelBehaviour.cohs)
    drift = reshape(drift, length(obj.modelBehaviour.cohs), []);
end

% Without loss of generality, and to facilitate DV simulation and comparison with NI models later, will assign 67 ms of
% the nondecision time to the motor end, and the free part will be for PRE decisional NDT:
MT = 4*dt; % 4 refreshes = one ssvep cycle, 66.7 ms
simdat = [];
acc = 1;
hit = 1; miss = 2; false_alarm = 3; % codes for response types
waitAfterTarget = obj.stim.refreshRate/2;  % start the urgency signal half a second into the ITI
continAccDur = obj.stim.refreshRate/10; % for how many samples should the CPP keep accumulating after reaching commitment? From R-locked CPPs it looks like around 100 ms

for n = 1:size(noise,3) % for robustness, it can help to simulate more trials than there are in the real data...
    for indBlock = 1:length(TOW)
        
        if length(bound) > 1; currBound = bound(bt(indBlock)); else; currBound = bound; end
        if isempty(leak); currLeak = 0; elseif length(leak) > 1; currLeak = leak(bt(indBlock)); else; currLeak = leak; end
        if isempty(trialnoise); currTrialnoise = 1; elseif length(trialnoise) > 1; currTrialnoise = trialnoise(bt(indBlock)); else; currTrialnoise = trialnoise; end
        if isempty(criteria); currCriteria = 0; elseif length(criteria) > 1; currCriteria = criteria(bt(indBlock)); else; currCriteria = criteria; end % TODO add the criteria
        
        if isempty(boost)
            currBoost = 1;
        elseif ~isempty(boost)
            if obj.modelBehaviour.learnBoost == bt(indBlock)
                currBoost = boost;
            elseif obj.modelBehaviour.learnBoost ~= bt(indBlock)
                currBoost = 1;
            end
        end
        
        if size(drift,1) > 1
            if obj.modelBehaviour.learnBoost == bt(indBlock,1)
                currDrift = drift(2,:);
            else
                currDrift = drift(1,:);
            end
        else
            currDrift = drift;
        end
        
        % make sensory evidence waveform - the target-on pulse train plus Gaussian noise
        % changed this to make it more transferable, e.g. get the timeline
        % which should include the
        
        timeLine = TOW{indBlock}';
        
        timeLine(timeLine ~= 0) = currDrift(timeLine(timeLine ~= 0)).*currBoost;
        sensEv = timeLine + currTrialnoise.*noise(indBlock,1:length(TOW{indBlock}), n);
        
        DV(1:round(tnd/dt)) = 0; % initialize DV to a single 0 for the first time points of the block, up until the pre-decision nondecision time
        
        lastresp = -postRespPause; % we only accumulate if it's a certain time since last response. Initialise to this value so accumulation begins right away at start of block
        
        RespT   = []; % keep track of all responses
        targT   = []; % keep track of all targets
        
        DVendval = nan; % **************

        % now simulate the block by looping through all sample points:
        for indTime = round(tnd/dt)+1:length(sensEv)
            
            if  indTime <= lastresp+continAccDur | indTime > lastresp + postRespPause %  indTime <= lastresp+continAccDur || don't accumulate unless it has been a sufficient time since last response
                % this gives us the 'reflecting bound'
                if reflectingBound  
                    DV(indTime) = max((1-currLeak)*DV(indTime-1) + sensEv(indTime-round(tnd/dt)) - currCriteria, 0);% There is the main model equation! (see e.g. Ossmy et al 2013)
                else
                    DV(indTime) = (1-currLeak)*DV(indTime-1) + sensEv(indTime-round(tnd/dt)) - currCriteria; % There is the main model equation! (see e.g. Ossmy et al 2013)
                end
            else
                
                if isnan(DVendval) % I'm using this as a way to know when it is time to linearly ramp the CPP down to zero
                    DVendval = DV(indTime-1);
                end
                DV(indTime) = max(0,DVendval*(1-(indTime-(lastresp+continAccDur))/25)); % linear decrease to a floor of zero
                if DV(indTime) == 0, DVendval = nan; end % when the CPP has reached back down to zero, turn off the linear ramp-down
                % DV(indTime) = 0;
            end
            
            % detect target transitions
            if TOW{indBlock}(indTime) > TOW{indBlock}(indTime-1)
                targT = [targT indTime*dt];
            end
            
            % detect responses
            if indTime > waitAfterTarget
                if indTime-waitAfterTarget <= size(urgency,1)
                    urgencyDV(indTime) = DV(indTime) + urgency(indTime-waitAfterTarget,bt(indBlock)); % total decision signal including urgency
                else  % when there's a miss, hard to know what to do with urgency - most neutral thing might be to keep it steady where it is
                    urgencyDV(indTime) = DV(indTime) + urgency(end, bt(indBlock));
                end
            else
                urgencyDV(indTime) = DV(indTime);
            end
            
              % Now check for bound crossing. Only do this if it is not within postRespPause since last bound crossing
              if indTime > lastresp+postRespPause  % **************
                  
                  if urgencyDV(indTime) > currBound % TODO for disc task < -currBound beside RespT check if correct or mistake.
                      RespT = [RespT indTime*dt+MT]; % log the response after adding the non-decision time
                      lastresp = indTime; % and now this is the last response that happened, at sample n
                      
                      if ~isempty(targT)
                          waitAfterTarget = round(targT(end)/dt)+maxRTms*obj.stim.refreshRate; % start counting the urg again at 0.5 sec after target offset
                      end
                  end
              end
        end
        
        % now make the output matrix that has all the RTs w.r.t. target onset
        % and false alarms. This has to be equivalent to how the real data were
        % analysed! Deals in seconds
        
        for indTarget = 1:length(targT)
            nextrespind = find(RespT > targT(indTarget)+minRTms & RespT < targT(indTarget)+maxRTms, 1); % find the index of the next response which is within the allowable @hit@ window
            if ~isempty(nextrespind) % TODO give additional Hit -- correct or wrong.
                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) hit RespT(nextrespind)-targT(indTarget)]; % 1 = hit
            else % if there WAS no next response, set the response parameters for this trial as 'not a number'
                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) miss nan]; % 2 = miss
            end
            
%             figure, plot(DV), hold on, plot(targT/dt, repmat(drift(1),size(targT)), '*')
%             plot(timeLine)
%             plot(RespT/dt, repmat(drift(1),size(RespT)), '*')
            
%             floor((obj.stim.targetEpoch(1)./1000)/dt)
%             ceil((obj.stim.targetEpoch(2)./1000)/dt)
%              timeFrame = [-6:60];
            timeFrame = [floor((obj.stim.targetEpoch(1)./1000)/dt):ceil((obj.stim.targetEpoch(end)./1000)/dt)];
              if getCPP == 1
                 if (targT(indTarget)/dt) + ceil((obj.stim.targetEpoch(end)./1000)/dt) < length(DV)
                     DecValue.Target(:,acc) = DV(round(targT(indTarget)/dt) + timeFrame);
                   if ~isempty(nextrespind)
                       DecValue.Response(1:length((floor((obj.stim.responseEpoch(1)./1000)/dt):ceil((obj.stim.responseEpoch(end)./1000)/dt))),acc) = DV(round(RespT(nextrespind)/dt) + [floor((obj.stim.responseEpoch(1)./1000)/dt):ceil((obj.stim.responseEpoch(end)./1000)/dt)]);
                   else
                       DecValue.Response(1:length((floor((obj.stim.responseEpoch(1)./1000)/dt):ceil((obj.stim.responseEpoch(end)./1000)/dt))),acc) = nan;
                   end
                 else
                     DecValue.Target(:,acc) = nan;
                     if ~isempty(nextrespind)
                         DecValue.Response(:,acc) = nan;
                     end
                 end
                 acc = acc + 1;
             end
        end
        % TODO check what to do with FA in disc task
        ITIstartT = [0 targT+1];
        for indITI = 1:length(ITIstartT)
            nexttargind = find(targT > ITIstartT(indITI),1); % index of next target
            % from this establish the end of the window starting from the current ITI start where we will check for FAs
            if ~isempty(nexttargind)
                endtime = targT(nexttargind) + minRTms;
            else
                endtime = length(DV)*dt-.125; % like in real data, if we didn't find a next target then this must be the end of the block, so check up as far as we ould possibly extract an erpr
            end
            % now find responses in this ITI window:
            nextrespind = find(RespT > ITIstartT(indITI)+maxRTms-(obj.stim.duration/1000) & RespT < endtime); % find indices of responses the ITI window, ruling out any at the very start that are within the hit window from the previous target. Target duration is 1 sec, so fs in sample points
            
            for m = 1:length(nextrespind)
                simdat  = [simdat; trialMatrix{indBlock}(indITI,:) false_alarm RespT(nextrespind(m))-ITIstartT(indITI)];
                
                DecValue.Target(:,acc) = nan;
                % TODO here we can start seeing psychophysical
                % kernel if we have the actual traces
                DecValue.Response(:,acc) = nan;
                
                acc = acc + 1;
            end
        end
        
    end
end
% It's possible that with certain parameters the DV NEVER crosses the bound, so need a quick fix so no error happens:
% this happened once by total fluke - there was no 'E' below, so there must have been a false alarm before the first target (cond=nan) and nothing else!

if isempty(simdat)
    tmpCond = unique(reshape([trialMatrix{:}]', size(trialMatrix{1},2), [])', 'rows');
    simdat = [tmpCond repmat(3, size(tmpCond,1),1) repmat(10, size(tmpCond,1),1)];
end

if obj.modelBehaviour.ChiOrG == 1
    [err,pred] = Chisquared(obj, simdat);
elseif obj.modelBehaviour.ChiOrG == 2
    [err,pred] = Gsquared(obj,simdat);
end
end