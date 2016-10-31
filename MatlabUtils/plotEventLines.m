% Plot vertical event lines
% -------------------------
% all = 1: plot all events (including general events)
% all = 0: plot only foot events (all except general events)
% 
% stanceShade = 'R': shade right leg stance
% stanceShade = 'L': shade left leg stance
% stanceShade = 'B': shade both leg stance
% ------------------------------------------------------------------

function plotEventLines(C3Dkey, timeVectorToUse, all, stanceShade)

if nargin == 2
    all = 1;
    stanceShade = '';
elseif nargin == 3
    stanceShade = '';
elseif nargin ~= 4
    disp('Error plotting event lines: USAGE: plotEventLines(C3Dkey, timeVectorToUse, all*\n')
    return
end

col = 'm:';
ratioTxtAbove = 1.05;
timeVec = C3Dkey.event.(timeVectorToUse);


% Only plot events within the time frame of the x axes plot window
xl = xlim;
events2use = [];
for i = 1:length(timeVec)
    if timeVec(i) >= xl(1) && timeVec(i) <= xl(end)
       events2use = [events2use i]; 
    end
end
% disp(events2use)


% Only shade stance events within the time frame of the x axes plot window
if strcmpi(stanceShade, 'R')
%     tmp = C3Dkey.stanceFrames.right;
    tmp = getStance(C3Dkey, stanceShade, events2use);

elseif strcmpi(stanceShade, 'L')
%     tmp = C3Dkey.stanceFrames.left;
    tmp = getStance(C3Dkey, stanceShade, events2use);
    
elseif strcmpi(stanceShade, 'B')
%     tmp = [C3Dkey.stanceFrames.right ; C3Dkey.stanceFrames.left];
    tmp = [getStance(C3Dkey, 'r', events2use) ; getStance(C3Dkey, 'l', events2use)];

elseif strcmpi(stanceShade, 'I')
    tmp = getStance(C3Dkey, stanceShade, events2use);
    
elseif strcmpi(stanceShade, 'C')
    tmp = getStance(C3Dkey, stanceShade, events2use);    
    
else 
    tmp = [];
end
stanceInd = [];
if ~isempty(tmp) 
    for i = 1:length(tmp(:,1))
        if ~isempty(find(events2use, tmp(1))) && ~isempty(find(events2use, tmp(2)))
            stanceInd = [stanceInd ; tmp(i,:)];
        end
    end
end


% Plot event lines
try
    for i = events2use
        % all = 0
        if ~all
            if ~strcmp(C3Dkey.event.txt{i}, 'rGEN') && ...
               ~strcmp(C3Dkey.event.txt{i}, 'lGEN') && ...
               ~strcmp(C3Dkey.event.txt{i}, 'GEN')
                h = vline(timeVec(i), col, C3Dkey.event.txt{i}, ratioTxtAbove);
                hasbehavior(h, 'legend', false);
            end
            
        else
            h = vline(timeVec(i), col, C3Dkey.event.txt{i}, ratioTxtAbove);
            hasbehavior(h, 'legend', false);
        end
        
    end
    
    % shade stance leg
    for i = 1:length(stanceInd(:,1))
        ShadePlotForEmpahsisVert([timeVec(stanceInd(i,1)) timeVec(stanceInd(i,2))], 'm', 0.1);
    end
catch
    return
end
