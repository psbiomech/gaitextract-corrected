% Extract, filter, and plot data from a saved OpenSim data file (*.mot or *.sto)
% into a Matlab structure
% Tim Dorn
% June 2010
% 
% --------------------------------------------------------------------
% Usage: [out, plots2make, labels, dataFile, hp, xl] = extractMotFile(ProcessTask1, Value1, ProcessTask2, Value2, ...)
% --------------------------------------------------------------------
% 
% Inputs:   ProcessTaskX = Processing task X
% 
%       Task: 'FILE'        Value = filename string (default = choose file manually)
%       Task: 'FILEPRELOADED'        Value = preloaded out structure from extractMotFile
%       Task: 'XAXIS'       Value = field for x axis (default = 'time')
%       Task: 'XAXISLAB'    Value = x axis label (if different to XAXIS)
%       Task: 'YAXIS'       Value = y axis label
%       Task: 'FILT'        Value = low pass filter frequency (default = 0)
%       Task: 'PLOT'        Value = plot variables
%                                -1 = off
%                                 0 = user selectable (default)
%                                 [n1 n2] = plot the columns specified by the array
%       Task: 'FIG'         Value = figure to plot into (default = 123)
%       Task: 'VERTLINES'   Value = C3D event key (from getEvents.m)
%       Task: 'EMG'         Value = EMG mot file (from batchEMGprocess.m)
%       Task: 'GRAPHMIN'    Value = Minimum value for vertical axis of graphs
%       Task: 'GRAPHMAX'    Value = Maximum value for vertical axis of graphs
%       Task: 'HLINE'       Value = Vector of horizontal lines for each index
%       Task: 'FIBRELEN'    Value = OpenSim model structure (from Load_OSIM)
%       Task: 'ALLON1FIG'   Value = 0 = seperate figures, 
%                                   1 = one figure (default)
%                                   2 = overlay data on same figure 
%       Task: 'SUPERTITLE'  Value = label for big title of plot
%       Task: 'MAX'         Value = 1(maximize), 0(don't maximize)
%       Task: 'FIGPREP'     Value = 1(default), 0(don't prepare figure)
%       Task: 'WRITETITLE'  Value = 1(write title), 0(don't write title)
%       Task: 'LINEW'       Value = line width as a number
%       Task: 'LINEC'       Value = line colour as a 1x3 vector
%       Task: 'LINES'       Value = line style as a string, i.e. '-'
%       Task: 'HLINE0'      Value = 0(noline), 1(zero horizontal line)
%       Task: 'MULTIPLY'    Value = multiply factor
%       Task: 'ADD'         Value = addition number
%       Task: 'DIFF'        Value = 0(no differentiate), 1(differentiate)
%       Task: 'STDDEV'      Value = indices from which to add stddev data to plot
% 
% 
% 
% 
% Outputs:  out.name:     trial name
%           out.labels:   extracted data labels
%           out.data:     extracted data (after filtering)
%           out.dataFile: filename that the data cames from
%           out.filtFreq: low pass filter frequency
%           performPlots: indices of plots extracted
%           labels:       all labels from dataFile
%           dataFile:     filename used
%           hp:           plot handle
%           
% 
% Notes:    If no input arguments are given i.e. data = extractMotFile, 
%           the function allows you to select a file for plotting and will
%           superimpose over existing plots if a common variable is being
%           plotted.
% 
% --------------------------------------------------------------------
% 
% Copyright (c)  2008 Tim Dorn
% Use of the GaitExtract Toolbox is permitted provided that the following
% conditions are met:
% 	1. The software is not distributed or redistributed.  Software distribution is allowed 
%     only through https://simtk.org/home/c3dtoolbox.
% 	2. Use of the GaitExtract Toolbox software must be acknowledged in all publications,
%      presentations, or documents describing work in which the GaitExtract Toolbox was used.
% 	3. Credits to developers may not be removed from source files
% 	4. Modifications of source code must retain the above copyright notice, this list of
%     conditions and the following disclaimer. 
% 
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
%  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
%  SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
%  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
%  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
%  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR BUSINESS INTERRUPTION) OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
%  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% --------------------------------------------------------------------

function [out, plots2make, labels, dataFile, hp, xl] = extractMotFile(varargin)

larg = length(varargin);
usage = 'Usage: [out, plots2make, labels, dataFile] = extractMotFile(ProcessTask1, Value1, ProcessTask2, Value2, ...)';

if mod(larg, 2) == 1,
    disp(usage)
    return
end
labels = [];

% Collate tasks
% -------------
numTasks = length(varargin)/2;
for i=1:numTasks
   tasks{i,1} = varargin{2*(i-1)+1};    % task name
   tasks{i,2} = varargin{2*(i-1)+2};    % task value
end


% Set defaults
% ------------
file = [];
xaxisLab = [];
dataFile = [];
dataFilePreLoaded = [];
xaxis = 'time';
yaxis = '';
performPlots = 0;
freq = 0;
figNum = 123;
C3Dkey = [];
emgFile = [];
EMGdataset = [];
ALLON1FIG = 1;
hl = [];
model_addfibrelen = [];
penn = 0;
supT = [];
maxPlot = 1;
figPrep = 1;
dataInput = [];
saveFile = 0;
excitation = [];
activation = [];
addkinematics = [];
model_addtendonstrain = [];
writeTitle = 1;
LineW = '';
LineC = '';
LineS = '';
hline0 = 1;
model_maxfom = [];
calcmean = 0;
mult = 1;
add = 0;
differentiate = 0;
addcontrolconstraints = [];
normModel = [];
addwork = [];
netAreaTOTAL = 0;
timecrop = [];
plotArea = 0;
labelMinMax = 0;
normSampFreq = 0;
generateMotionFile = [];
otherstring = '';
stddevInd = [];
includeStddev = 0;


% Set Variables
% -------------
for i=1:numTasks
    
    switch upper(tasks{i,1})
        case 'FILE'
            dataFile = tasks{i,2};
            out.dataFile = dataFile;

        case 'FILEPRELOADED'
            dataFilePreLoaded = tasks{i,2};            
            
        case 'XAXIS'
            xaxis = tasks{i,2};
            
        case 'XAXISLAB'
            xaxisLab = tasks{i,2};    
            
        case 'YAXIS'
            yaxis = tasks{i,2};            
            
        case 'FILT'
            freq = tasks{i,2};

        case 'PLOT'
            performPlots = tasks{i,2};

        case 'FIG'
            figNum = tasks{i,2};
            
        case 'VERTLINES'
            C3Dkey = tasks{i,2}{1};
            stanceShade = tasks{i,2}{2};
            
        case 'EMG'
            emgFile = tasks{i,2};
            if isstruct(emgFile)
                EMGdataset = emgFile;
                timeEMG = getData('time', EMGdataset);
            elseif ischar(emgFile)
                EMGdataset = extractMotFile('file', emgFile, 'plot', -1);
                timeEMG = getData('time', EMGdataset);
            end
            
        case 'GRAPHMIN'
            graphMin = tasks{i,2};

        case 'GRAPHMAX'
            graphMax = tasks{i,2};
            
        case 'HLINE'
            hl = tasks{i,2};
            
        case 'FIBRELEN'
            model_addfibrelen = tasks{i,2};
            
        case 'PENNATION'
            penn = tasks{i,2};     
            
        case 'ALLON1FIG'
            ALLON1FIG = tasks{i,2};

        case 'SUPERTITLE'
            supT = tasks{i,2};
            
        case 'MAX'
            maxPlot = tasks{i,2};    
            
        case 'FIGPREP'
            figPrep = tasks{i,2};    
            
        case 'DATAINPUT'
            dataInput = tasks{i,2};      
            
        case 'SAVEFILE'
            saveFile = tasks{i,2}; 
            
        case 'ADDEXCITATION'
            excitation = tasks{i,2}; 
            
        case 'ADDACTIVATION'
            activation = tasks{i,2}; 

        case 'ADDTENDONSTRAIN'
            model_addtendonstrain = tasks{i,2};
            
        case 'ADDCONTROLCONSTRAINTS'
            addcontrolconstraints = tasks{i,2};
            if ~isempty(addcontrolconstraints)
                cc = xml_read(addcontrolconstraints);
            end
            
        case 'ADDKINEMATICS'
            addkinematics = tasks{i,2}{1};
            Qind = tasks{i,2}{2};
              
        case 'WRITETITLE'
            writeTitle = tasks{i,2}; 

        case 'LINEW'
            LineW = tasks{i,2}; 

        case 'LINEC'
            LineC = tasks{i,2}; 

        case 'LINES'
            LineS = tasks{i,2}; 
            
        case 'HLINE0'
            hline0 = tasks{i,2}; 
            
        case 'MAXFOM'
            model_maxfom = tasks{i,2};
            
        case 'CALCMEAN'
            calcmean = tasks{i,2};
 
        case 'MULTIPLY'
            mult = tasks{i,2};
            
        case 'ADD'
            add = tasks{i,2};    
            
        case 'DIFF'
            differentiate = tasks{i,2};
            
        case 'NORMALIZE'
            normModel = tasks{i,2}{1};
            normField = tasks{i,2}{2};

        case 'ADDWORK'
            addwork = tasks{i,2};
            bigWorkData = [];
            bigWorkCols = {};
            if ~isfield(addwork, 'workDir')
                addwork.workDir = 1;
            end
             
        case 'TIMECROP'
            timecrop = tasks{i,2};
            
        case 'PLOTAREA'
            plotArea = tasks{i,2};
            
        case 'LABELMINMAX'
            labelMinMax = tasks{i,2};
            
        case 'OTHERSTRING'
            otherstring = tasks{i,2};

        case 'NORMALIZETIME'
            normTimeKey = tasks{i,2}{1};
            verboseOutput = tasks{i,2}{5};
            if isempty(normTimeKey)
                [file, pathStr] = uigetfile({'*.mat', 'C3DKey MAT File (*.mat)'}, ...
                    'Select C3Dkey *.mat file');
                C3DfileLoc = [pathStr, file];
                load(C3DfileLoc);
                normTimeKey = C3Dkey;
            else
                C3Dkey = normTimeKey;
                pathStr = fileparts(tasks{i,2}{4}); % {4} is path of c3d mat file in GaitExtract directory
            end

            try % to get body mass
                bodyMass = load(sprintf('%s\\..\\..\\Static\\tmass.data', pathStr));
            catch
                fprintf('Could not locate body mass file\n');
                bodyMass = 9999;
            end
            
            try % to get body height
                height = load(sprintf('%s\\..\\..\\Static\\height.data', pathStr));
            catch
                fprintf('Could not locate body height file\n');
                bodyHeight = 9999;
            end
            
            try % to get test leg
                fid = fopen(sprintf('%s\\..\\testLeg.data', pathStr), 'rb');
                testLeg = fread(fid, 1, 'uint8=>char');
                fclose(fid);
            catch
                fprintf('Unknown test leg\n');
            end

            if verboseOutput
                try
                    fprintf('Trial: %s\t(mass = %.2f kg, height = %.1f cm)\n', ...
                    C3Dkey.c3dFile, bodyMass, height);
                catch
                    height = 0;
                end
                
                if strcmpi(testLeg, 'r')
                    fprintf('Test Leg: RIGHT\n')
                else
                    fprintf('Test Leg: LEFT\n')
                end
            end
            
            normSampFreq = tasks{i,2}{2};       % normalizing sample frequency
            eventIndex = tasks{i,2}{3};         % event indices to normalize between
            
        case 'GENERATEMOTFILE'
            generateMotionFile = tasks{i,2};
            
        case 'STDDEV'
            stddevInd = tasks{i,2};   % indices of stddev data
            
        case 'INCLUDESTDDEV'
            includeStddev = tasks{i,2};   % include plot of stddev (with the stddev_index = mean_index+1
            
        otherwise
            error('Error: UNKNOWN TASK (%s)...\n', upper(tasks{i,1}));
    end
end

if freq < 0, 
    error('Error: filtering frequency must be > 0');
end

if isempty(dataInput) && isempty(dataFilePreLoaded)
    
    if isempty(dataFile) || ~exist(dataFile, 'file')
        fprintf('Warning: dataFile [%s] is not defined or does not exist...\n', dataFile');
        [file, pathStr] = uigetfile({'*.mot;*.sto', 'OpenSim Motion Files (*.mot, *.sto)'; ...
           '*.*',  'All Files (*.*)'}, 'Select Motion File');
        dataFile = [pathStr, file];
    else
        file = dataFile;
    end



    % Extract motion data from external file
    % --------------------------------------
    nlrows = 1;
    nhead = 1;

    
% I decided to remove this part because ssr.exe is an additional file 
% required in the system path. All it does is replace the OpenSim NaN
% description with 'NaN' (TD Sept2010)
%     try
%         % rename nans to be detectable. REQUIRES ssr.exe in system32 directory
%         eval(sprintf('!ssr 0 "1.#QNAN000000000000000" "NaN" "%s"', dataFile));
%         eval(sprintf('!ssr 0 "1.#IND0000000000000000" "NaN" "%s"', dataFile));
%     catch
%         fprintf('NaNs in the file will cause an error...\n')
%     end


    fid = fopen(dataFile, 'r');
    if fid < 0
        fprintf('\nERROR: %s could not be opened for reading...\n\n', dataFile);
        return
    end

    buffer = fgetl(fid);
    while strcmpi(strtrim(buffer), 'endheader') == 0,
        if nhead == 1,
            out.name = buffer;
        end
        nhead = nhead + 1;
        if ~isempty(findstr(buffer, 'nColumns')),
            % get number of columns from the header
            ncols = strread(buffer, '%*s%d', 'delimiter', '=');
        end
        if ~isempty(findstr(buffer, 'nRows')),
            % get number of rows from the header
            nrows = strread(buffer, '%*s%d', 'delimiter', '=');
        end
        buffer = fgetl(fid);
    end
    fclose(fid);

    [labels, x, y] = readColData(dataFile, ncols, nhead, nlrows);
    out.labels = labels;
    out.data = [x, y];
else
    out = dataInput;
    if ~isempty(out)
        [nrows, ncols] = size(out.data);
    end
end

if ~isempty(dataFilePreLoaded)
   out = dataFilePreLoaded;
   [nrows, ncols] = size(out.data);
   labels = out.labels;
end


% perform multiplication/addition operations here
out.data = [out.data(:,1) out.data(:,2:end) * mult + add];
xData = getData(xaxis, out);
sampFreq = round(1 / (xData(end) - xData(end-1)));
out.ts = xData(2)-xData(1);



% Filter Data (optional)
% ----------------------

if freq > 0,
%     fprintf('Filtering Data @ %d Hz\n', freq);
    for i = 1:ncols
        if ~strcmpi(out.labels{i}, 'time'),
            out.data(:,i) = smooth(out.data(:,i), freq, sampFreq);
        end
    end
end
out.freq = freq;


% Differentiate (optional)
% ------------------------
if differentiate
    for j = 2:ncols
        pp = spline(xData, out.data(:,j));
        ppd = diffpp(pp);
        out.data(:,j) = ppval(ppd, xData);
    end
end


% Normalize Time
% ---------------

if normSampFreq > 0
    PERCENT = 1;
    out.C3Dkey = C3Dkey;
    out.bodyMass = bodyMass;
    out.height = height;
    
    timeIndex = getIndex(out.labels, 'time');
    timeVec = out.data(:, timeIndex);
    startTime = timeVec(1);
    endTime = timeVec(end);
    
    % look at the time column and determine which events are in the selection
    allEvents = normTimeKey.event.txt;
    events2use = [];
    for i = 1:length(allEvents)
        if startTime <= normTimeKey.event.times0(i) && endTime+0.01 >= normTimeKey.event.times0(i)
           events2use = [events2use i]; 
        end
    end

    event.txt = {'START', normTimeKey.event.txt{events2use}, 'END'};
    event.time = [startTime, normTimeKey.event.times0(events2use), endTime];
    
    for j = 1:length(event.time)
        listLabels{j} = sprintf('%s   (t=%.2f)', event.txt{j}, event.time(j));
    end
    
    % load event details from C3Dkey and get user to select two events to normalize between  
    if isempty(eventIndex)
        [eventIndex, v] = listdlg('PromptString', 'Select 2 events to normalize in between:',...
                          'SelectionMode', 'multiple',...
                          'ListString', listLabels);
        if v == 0, return; end      % cancel button pressed
    end
    
    if length(eventIndex) ~= 2
        error('Select only 2 events (start and end event) to normalize in between!\n')
    end
    
    % We provide an automatic detection of events here so you don't need to
    % always manually work out the event numbers
    if eventIndex(1) == -1              % Right foot full gait cycle (rFS to rFS)
        eventIndex = getEventIndex(event.txt, 'rFS', 'rFS');
        
    elseif eventIndex(1) == -2          % Left foot full gait cycle (lFS to lFS)
        eventIndex = getEventIndex(event.txt, 'lFS', 'lFS');
        
    elseif eventIndex(1) == -3          % Right foot stance (rFS to rFO)
        eventIndex = getEventIndex(event.txt, 'rFS', 'rFO');
            
    elseif eventIndex(1) == -4          % Left foot stance (lFS to lFO)
        eventIndex = getEventIndex(event.txt, 'lFS', 'lFO');
        
    end
    
    
    eventIndex = sort(eventIndex);
    if verboseOutput
        fprintf('List Of Labels:');
        disp(listLabels)
        l = length(event.txt);
        if eventIndex(1) > l || eventIndex(2) > l
            fprintf('EventIndex:');
            disp(eventIndex)
            fprintf('EventTxt:');
            disp(event.txt)
            error('EventIndex is out of range!');
        end
        fprintf('Normalizing from:  %s(t=%.3f)[%d]  to  %s(t=%.3f)[%d]\n', ...
            event.txt{eventIndex(1)}, event.time(eventIndex(1)), eventIndex(1), ...
            event.txt{eventIndex(2)}, event.time(eventIndex(2)), eventIndex(2));
        
        % out.normalizedTime1 = [beginning of stance, end of stance, midswing, beginning of next stance]
        out.normalizedTime1 = [event.time(eventIndex(1)) event.time(eventIndex(1)+1) ...
            event.time(eventIndex(1)+1) + 0.5*(event.time(eventIndex(2))-event.time(eventIndex(1)+1)) event.time(eventIndex(2))];
    end
    
    % Store normalized event times
    event.eventIndex = eventIndex;
    event.normalizedTime = event.time(eventIndex(1):eventIndex(2))-event.time(eventIndex(1));
    event.normalizedText = event.txt(eventIndex(1):eventIndex(2));
    out.event = event;
    out.testLeg = testLeg;
    out.absoluteTime = event.time(eventIndex(2))-event.time(eventIndex(1));
    
    z=abs(timeVec-event.time(eventIndex(1)));
    tmp = find(min(z)==z);
    frames(1) = tmp(1);
    
    z=abs(timeVec-event.time(eventIndex(2)));
    tmp = find(min(z)==z);
    frames(2) = tmp(end);

    % normalize into a fixed number of points and replace
    newTime = 1:normSampFreq;
    if PERCENT
       newTime = 100 * ((newTime-1) / newTime(end-1));
       out.event.normalizedTime =  100 * (out.event.normalizedTime/out.event.normalizedTime(end));
    end
    data2resamp = out.data(frames(1):frames(end), :);
    out.data = resamp2(data2resamp', normSampFreq)';    % do the normalization
    out.data(:,timeIndex) = newTime;
end



% Generate new MOT file
% ----------------------
if ~isempty(generateMotionFile)
    try
        newFileName = sprintf('%s_%s-%s.mot', strtok(generateMotionFile, '.'), ...
            event.txt{eventIndex(1)}, event.txt{eventIndex(end)});
    catch
        newFileName = generateMotionFile;
    end
    generateMotFile(out.data, out.labels, newFileName);
end








% Plot extracted data
% -------------------
if ~isempty(performPlots) && performPlots(1) >=0, 
    
%     ALLON1FIG = 1;          % 0 = seperate figures, 1 = one figure
    
    % Range of plot colours
    col = [0 0 1;
           1 0 0;
           0 0 0;
           0.2 0.5 0.1;
           0.5 0.4 0.68;
           1 0.7 0.2;
           0.7 0.7 0.3
           0.31 0.50 0.67
           0.23 0.62 0.62
           0.75 0.50 0.41
           0.39 0.31 0.10];
       
    lw = [2 2 3];                             % Range of line widths
    ls = [{'-'}, {'--'}, {'-'}];     % ['.', 'x', '*', '+', 'o',];      % Range of line styles
    
    for i = 1:length(out.labels)
        listLabels{i} = sprintf('[%d]  %s', i, out.labels{i});
        
        if isstruct(normModel)
            try
                if strcmpi(normField, 'Vmax')
                    out.data(:,i) = out.data(:,i) / ...
                        (normModel.Muscles.(strtok(out.labels{i}, '.')).(normField) * ...
                        normModel.Muscles.(strtok(out.labels{i}, '.')).optimal_fiber_length);
                else
                    out.data(:,i) = out.data(:,i) / ...
                        normModel.Muscles.(strtok(out.labels{i}, '.')).(normField);
                end
            catch
            end
        end

    end
    
    if performPlots(1) == 0,
        [plots2make,v] = listdlg('PromptString', 'Variables to plot:',...
                      'SelectionMode', 'multiple',...
                      'ListString', listLabels);
        if v == 0, return; end      % cancel button pressed
    else
        plots2make = performPlots;
    end
    
    l = length(plots2make);
    pl = ceil(sqrt(l));
                    
    % Main plotting loop
    for i = 1:l
        j = plots2make(i);

        
        if figPrep == 1
            if ALLON1FIG == 0
                pl = 1;
                figure(100+j)
                
            elseif ALLON1FIG == 1
                h = figure(figNum);
                set(h,'PaperType','a4', 'PaperPositionMode','manual',...
                    'PaperOrientation','landscape', 'PaperUnits','centimeters',...
                    'PaperPosition', [0,0,29,20]);
                
                hs = subplot(pl, pl, i);
                out.subplotHandles(i) = hs;
%                 maximize
                hold on
                
            elseif ALLON1FIG == 2   % overlay on same figure
                h = figure(figNum);
%                 set(h,'PaperType','a4', 'PaperPositionMode','manual',...
%                     'PaperOrientation','landscape', 'PaperUnits','centimeters',...
%                     'PaperPosition', [0,0,29,20]);
%                 maximize
                hold on
                
            elseif length(ALLON1FIG) == 1       % num of columns in plot
                h = figure(figNum);
                numR = ceil(l/ALLON1FIG);
                hs = subplot(numR, ALLON1FIG, i);
                out.subplotHandles(i) = hs;
                
            elseif length(ALLON1FIG) == 2
                h = figure(figNum);
                hs = subplot(ALLON1FIG(1), ALLON1FIG(1), ALLON1FIG(1)*(ALLON1FIG(2)-1)+i);
                out.subplotHandles(i) = hs;
                hold on
                
            elseif length(ALLON1FIG) == 3
                h = figure(figNum);
                hs = subplot(ALLON1FIG(1), ALLON1FIG(2), ALLON1FIG(1)*(ALLON1FIG(2)-1)+i);
                out.subplotHandles(i) = hs;
                hold on
            end
        end
        
        if maxPlot ~= 0
            maximizecustom(maxPlot);
        end
        
        hold on
%         , 'Name', sprintf('%s', out.labels{j}), 'NumberTitle', 'off');
        [a,b,c,d] = legend;
        
        if ~isempty(LineW)
            lw = LineW;
        end
        if ~isempty(LineC)
            col = LineC;
        end
        if ~isempty(LineS)
            ls = LineS;
        end
        
        if plotArea
           hp = area(xData, out.data(:,j), ...
                'LineWidth', lw(1+mod(length(d), length(lw))), ...
                'LineStyle', ls{1+mod(length(d), length(ls))}, ...
                'FaceColor',[0.9 0.9 0.9]); 
                     
            ts = xData(2) - xData(1);
            areaAboveYaxis = trapz(max(out.data(:,j),0))*ts;
            areaBelowYaxis = trapz(min(out.data(:,j),0))*ts;
            netArea1 = areaAboveYaxis+areaBelowYaxis;
            netArea2 = areaAboveYaxis-areaBelowYaxis;  
            
%             uistack(hp, 'bottom');

        elseif ~isempty(stddevInd)
            cc = 0.9-0.05*i;
            if cc < 0.5
                cc = 0.5;
            end
                
            stddev = out.data(:,stddevInd(i));
            hand = errorshade3(xData, out.data(:,j), stddev, '', cc*ones(1,3));  % used to be confplot
            hp = hand(1);
            hasbehavior(hand(2), 'legend', false);     % dont show stddev in legend entries!
            uistack(hp, 'top');
            set(hp, 'Color', col(1+mod(length(d), size(col,1)), :), ...
                'LineWidth', lw(1+mod(length(d), length(lw))), ...
                'LineStyle', ls{1+mod(length(d), length(ls))})
            
            
        else
            hp = plot(xData, out.data(:,j), ...
                'Color', col(1+mod(length(d), size(col,1)), :), ...
                'LineWidth', lw(1+mod(length(d), length(lw))), ...
                'LineStyle', ls{1+mod(length(d), length(ls))});
            uistack(hp, 'top');
        end
        
%     col{1+mod(length(d), length(col))})
        
        avg = mean(out.data(:,j));
        avgabs = mean(abs(out.data(:,j)));
        maxval = max(out.data(:,j));
        meanrms = rms(out.data(:,j));
        
%             leg = sprintf('%s (avg = %.2f, avgabs = %.2f, max = %.2f)', ...
%                 file, avg, avgabs, maxval);

        if calcmean
            hline(meanrms, 'r--', sprintf('RMS = %.3f\n', meanrms));
            out.meanrms(j) = meanrms;
        end

        % The legend value denotes the legend in the matlab plotbrowser
        leg = out.labels{j};
        lh = legend([d';{leg}], 'Location', 'Best', 'Interpreter', 'none');
        legend hide
%         set(lh, 'Visible', 'off')
%         legend boxoff
        
        if(plotArea)
            xlabel(sprintf('%s, PosArea=%f, NegArea=%f, TotArea=%f', xaxis, areaAboveYaxis, areaBelowYaxis, netArea2), ...
                'FontSize', 8, 'Interpreter', 'none')
        else
            xlabel(xaxis, 'FontSize', 8, 'Interpreter', 'none')
        end
        
        if ~isempty(xaxisLab)
            xlabel(xaxisLab, 'FontSize', 8, 'Interpreter', 'none')
        end

        % plot yaxis only on left hand subplots
        if mod(i-1, pl) == 0
            ylabel(yaxis, 'FontSize', 8, 'Interpreter', 'none')
        end
        
        yltmp = ylim;
        axis tight
        xl = xlim;
        ylim(yltmp);
        
        if hline0
            h_h = hline(0, 'k:');
            hasbehavior(h_h, 'legend', false);     % dont show stddev in legend entries!
        end
        
        if ~isempty(hl)
            if length(hl) == 1
                hline(hl, 'k:');
            else
                hline(hl(j), 'k:');
            end
        end
        

        if isstruct(model_maxfom)
            try
                fom = model_maxfom.Muscles.(strtok(labels{j}, '.')).max_isometric_force * mult;   % OSim <2.0
            catch
                fom = model_maxfom.Forces.(strtok(labels{j}, '.')).max_isometric_force * mult;    % OSim >2.0
            end
            hline(fom, 'k:');
        end
        
        
        axis tight
        croptime(timecrop, xl)
        yl = ylim;
        diff = 0.1*(yl(2) - yl(1));
        newyl = [yl(1)-diff, yl(2)+diff];
        ylim(newyl);

        if exist('graphMin', 'var')
            ylim([graphMin, yl(2)]);
            yl = ylim;  % refresh yl
        end
        
        if exist('graphMax', 'var')
            ylim([yl(1), graphMax]);
        end        
        

        % add stddev
        if includeStddev == 1
            stddev = out.data(:,j+1);
            hand = errorshade3(xData, out.data(:,j), stddev, '', 0.8*ones(1,3));  % used to be confplot
            hp = hand(1);
            hasbehavior(hand(2), 'legend', false);     % dont show stddev in legend entries!
            uistack(hp, 'top');
            set(hp, 'Color', col(1+mod(length(d), size(col,1)), :), ...
                'LineWidth', lw(1+mod(length(d), length(lw))), ...
                'LineStyle', ls{1+mod(length(d), length(ls))})
        end
        
        
        % add kinematics
        if ~isempty(addkinematics)
            hold on
            xData2 = addkinematics.data(:, getIndex(addkinematics.labels, 'time'));
            yData2 = addkinematics.data(:, Qind);
            Ucol = [.3 .8 .3];

            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData2, yData2, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 2, 'Color', Ucol)

            hold on
            out.haxes = haxes;

            axis tight
            croptime(timecrop, xl)
            ylabel('MomArm / JntAng', 'FontSize', 8, 'Interpreter', 'none')
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            ylim(newyl);
            
            axes(haxes(2));
            axis tight
            croptime(timecrop, xl)
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl2 = ylim;
            diff2 = 0.1*(yl2(2) - yl2(1));
            newyl2 = [yl2(1)-diff2, yl2(2)+diff2];
            ylim(newyl2);
        end
        
        
        
        
        
        % add excitation / activation values
        if ~isempty(excitation)
            hold on
            xData2 = excitation.data(:, getIndex(excitation.labels, 'time'));
            try
                yData2 = excitation.data(:, getIndex(excitation.labels, [strtok(out.labels{j}, '.') '.excitation']));
                Ucol = [.3 .8 .3];
            catch
                try
                    yData2 = excitation.data(:, getIndex(excitation.labels, strtok(out.labels{j}, '.')));
%                     Ucol = [1 .64 .57];   % red 
                    Ucol = [.3 .8 .3];      % green
                catch
                    yData2 = excitation.data(:, getIndex(excitation.labels, [strtok(out.labels{j}, '.') '.activation']));
%                     Ucol = [.4 .6 .7];    % red
                    Ucol = [.3 .8 .3];      % green
                end
            end
            
            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData2, yData2, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 1, 'Color', Ucol)

            hold on
            out.haxes = haxes;

%                 axes(haxes(1));  
            axis tight
            croptime(timecrop, xl)
            if mod(i-1, pl) == 0
                ylabel('Force (N) / Excitation', 'FontSize', 8, 'Interpreter', 'none')
            end
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl = ylim;
            ylim([0, newyl(2)]);
%             ylim(newyl);
            
            axes(haxes(2));
            axis tight
            croptime(timecrop, xl)
            ylim([0, 1.1]);
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
        end
        
        
        if ~isempty (addcontrolconstraints)
            % get muscle control constraints
            for k = 1:length(cc.objects.ControlLinear)
                t = []; minP = []; maxP = [];
                if strcmpi(cc.objects.ControlLinear(k).ATTRIBUTE.name, ...
                        sprintf('%s.excitation', strtok(out.labels{j}, '.')))
                    try 
                        for m = 1:length(cc.objects.ControlLinear(k).min_nodes.ControlLinearNode)
                            t(m) = cc.objects.ControlLinear(k).min_nodes.ControlLinearNode(m).t;
                            minP(m) = cc.objects.ControlLinear(k).min_nodes.ControlLinearNode(m).value;
                            maxP(m) = cc.objects.ControlLinear(k).max_nodes.ControlLinearNode(m).value;
                        end
                    catch

                    end
                    break
                end
            end

            t = [xl(1) t xl(2)];
            if isempty(minP)
                minP = [0.01 minP 0.01];
            else
                minP = [minP(1) minP minP(end)];
            end
            if isempty(maxP)
                maxP = [1 maxP 1];
            else
                maxP = [maxP(1) maxP maxP(end)];
            end               

            % plot shaded control constraints
            hold on
            [ha hb hc] = shadedplot(t, minP, maxP, [1 1 0.8]);
            set(hb, 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.93 0.8])
            set(hc, 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.93 0.8])
            uistack(ha(2), 'bottom')
            alpha(0.2)
        end
            
                    
        if ~isempty(activation)
            xData3 = activation.data(:, getIndex(activation.labels, 'time'));
            try
                yData3 = activation.data(:, getIndex(activation.labels, [strtok(out.labels{j}, '.') '.activation']));
                Acol = [.71 .87 .5];
            catch
                yData3 = activation.data(:, getIndex(activation.labels, strtok(out.labels{j}, '.')));
                Acol = [1 .64 .57];
            end
            hold on
            plot(xData3, yData3, 'LineStyle', '-', 'LineWidth', 1, 'Color', Acol)
        end       
        
        
        % add tendon strain values
        % (assumes the plotting of tendon length with model as second variable input)
        if ~isempty(model_addtendonstrain)
%             hold on
%             lst = model_addtendonstrain.Muscles.(strtok(out.labels{j}, '.')).tendon_slack_length;
%             hline(lst, 'k:');
%             % NOTE THAT Thelen defines stretch not from lom, or lst, but rather from
%             % the length of the muscle/tendon in a standing posture (all angles = 0)
%             tendonStrain = (out.data(:,j) - lst)./lst;
%             if any(tendonStrain<0)
%                 fprintf('Warning: Tendon strain < 0 (slack tendon) at some stage for %s\n', out.labels{j})  
%             end
%             tendonStrain(tendonStrain<0) = 0;
%             
%             [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData, tendonStrain, ...
%                 'plot', 'plot');
%             set(hline1, 'LineStyle', 'none')
%             set(hline1, 'LineWidth', 2)
%             set(hline2, 'LineWidth', 1)
% %             set(hline1, 'LineStyle', 'none')       % remove tendon strain line
%             
%             hold on
%             out.haxes = haxes;
% 
% %                 axes(haxes(1));  
%             axis tight
%             croptime(timecrop, xl)
%             if mod(i-1, pl) == 0
%                 ylabel('Lengths (m) / TendonStrain', 'FontSize', 8, 'Interpreter', 'none')
%             end
%             set(gca, 'YTickLabelMode', 'auto')
%             set(gca, 'YTickMode', 'auto')
%             yl = ylim;
%             ylim([0, 1.1*yl(2)]);
            
            % assume here we are also plotting muscle fiblen (plotVarArms(9) && plotVarArms(20))
            if isfield(model_addtendonstrain, 'replaceForceName')
                newMuscName = regexprep(regexprep(strtok(out.labels{j}, '.'), '(mean))', '_r'), '(', '');
                ofl = model_addtendonstrain.Forces.(newMuscName).optimal_fiber_length;
            else
                ofl = model_addtendonstrain.Forces.(strtok(out.labels{j}, '.')).optimal_fiber_length;
            end
            ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
            hline(ofl, 'k:');
            
%             axes(haxes(2));
%             axis tight
            ylim([0 yltmp(2)])

            croptime(timecrop, xl)
%             ylim([0, 0.1]);
%             set(gca, 'YTickLabelMode', 'auto')
%             set(gca, 'YTickMode', 'auto')
        end

        
        % add work = integral of power
        % (assumes the plotting of power - this will overlay the integral)
        if ~exist('hs', 'var')
            hs = gca;
        end
        
        if ~isempty(addwork)
            hold on
            power = out.data(:,j);
            ts = xData(3)-xData(2);

            % addwork.workDir takes into account the 'direction' of work
            % for our positive Y axis. e.g. for muscles, addwork.workDir = -1
            % because legnthening (a positive movement) represents a negative work.
            netWork = cumtrapz(xData, power);  % work = trapezoidal integration of power
            areaAboveYaxis = addwork.workDir*trapz(max(power,0))*ts;   % eccentric contractions
            areaBelowYaxis = addwork.workDir*trapz(min(power,0))*ts;   % concentric contractions
            netArea = areaAboveYaxis+areaBelowYaxis;
            netAreaTOTAL = netAreaTOTAL + netArea;

%             h1 = hline(rms(power(find(power>0))), 'k:', sprintf('NEGWORK: %.3fJ', negWork));
%             set(h1, 'Visible', 'off');
%             h2 = hline(-rms(power(find(power<0))), 'k:', sprintf('POSWORK: %.3fJ', posWork));
%             set(h2, 'Visible', 'off');
%             h3 = hline(netWorkMuscle, 'k:', sprintf('NET: %.3fJ', netWorkMuscle));
%             set(h3, 'Visible', 'off');

            if ~isfield(addwork, 'shading')
                addwork.shading = 1;
            end
            
            if addwork.workDir < 0          % MUSCLES
                if addwork.shading == 1
                    h_area(1) = areashade(xData, power, 0, 'r');            % negative power (red) y>0
                    set(h_area(1), 'FaceAlpha', 0.3);
                    h_area(2) = areashade(xData, power, 0, 'g','h');        % positive power (green) y<0
                    set(h_area(2), 'FaceAlpha', 0.3);
                    hasbehavior(h_area(1), 'legend', false);     % dont show shading in legend entries!
                    hasbehavior(h_area(2), 'legend', false);     % dont show shading in legend entries!
                end
                aboveBelowLabels = {'NEGWORK(absorption)', 'POSWORK(generation)'};
            else                            % JOINTS
                if addwork.shading == 1
                    h_area(1) = areashade(xData, power, 0, 'g');            % positive power (green) y>0
                    set(h_area(1), 'FaceAlpha', 0.3);
                    h_area(2) = areashade(xData, power, 0, 'r','h');        % negative power (red) y<0
                    set(h_area(2), 'FaceAlpha', 0.3);
                    hasbehavior(h_area(1), 'legend', false);     % dont show shading in legend entries!
                    hasbehavior(h_area(2), 'legend', false);     % dont show shading in legend entries!
                end
                aboveBelowLabels = {'POSWORK(generation)', 'NEGWORK(absorption)'};
            end


%             set(h_area(2), 'FaceColor', [0 0 1]);
%             plot(xData, netWork, 'LineWidth', 2', 'Color', [.4 .6 .3]);
            
            % Detect stance and swing times in the plotted area
            % --------------------------------------------------
            xl = xlim;   % time axis of plotted area
            xl(2) = xl(2) + 0.02;
            if ~isfield(addwork, 'event')
                addwork.event.times0 = out.data(:,1);
            end
            timeVec = addwork.event.times0;
                
            events2use = [];
            for i = 1:length(timeVec)
                if timeVec(i) >= xl(1) && timeVec(i) <= xl(end)
                   events2use = [events2use i]; 
                end
            end
            % event vector indices (one gait cycle)
            if ~isfield(addwork, 'testLeg')
                addwork.testLeg = 'r';
            end
            eventIndStanceTMP = getStance(addwork, addwork.testLeg, events2use);
            eventIndSwingTMP = getSwing(addwork, addwork.testLeg, events2use);
            
            % use only one row (as there may be multiple stances and swings detected)
            eventIndStance = eventIndStanceTMP(1,:);  % first instance of stance
            % now we need to find a swing phase that follows immediately after stance has ended
            tmpInd = [];
            for k = 1:size(eventIndSwingTMP,1)
                if eventIndSwingTMP(k,1) >= eventIndStance(end)
                    tmpInd = [tmpInd ; eventIndSwingTMP(k,:)];
                end
            end
            eventIndSwing = tmpInd(1,:);
            % check swing follows stance
            if eventIndStance(end) ~= eventIndSwing(1)
                error('SWING DOES NOT IMMEDIATELY FOLLOW STANCE!\n')
            end
            
            out.gaitCycleEvents = addwork.event.txt(eventIndStance(1):eventIndSwing(end));
            if ~strcmpi(addwork.event.txt{eventIndStance(1)}, addwork.event.txt{eventIndSwing(end)})
                error('CAN''T DETECT A FULL GAIT CYCLE BEGINNING WITH STANCE ON LEG: %s\n\n', upper(addwork.testLeg));
            end
            
            % time boundaries (one gait cycle)
            timeStance = timeVec(eventIndStance);
            timeSwing = timeVec(eventIndSwing);
            timemidSwing = timeSwing(1) + ((timeSwing(2) - timeSwing(1)) / 2);
            timeSwing1 = [timeSwing(1) timemidSwing];
            timeSwing2 = [timemidSwing timeSwing(2)];
            
            % time vector indices
            indStance = findClosest(xData,timeStance(1)) : findClosest(xData,timeStance(2));
            indSwing1 = findClosest(xData,timeSwing1(1)) : findClosest(xData,timeSwing1(2));
            indSwing2 = findClosest(xData,timeSwing2(1)) : findClosest(xData,timeSwing2(2));
            
            % work calculations
            stanceAboveYaxis = addwork.workDir*trapz(max(power(indStance),0))*ts;   % eccentric contractions (stance)
            stanceBelowYaxis = addwork.workDir*trapz(min(power(indStance),0))*ts;   % concentric contractions (stance)
            swing1AboveYaxis = addwork.workDir*trapz(max(power(indSwing1),0))*ts;   % eccentric contractions (swing1)
            swing1BelowYaxis = addwork.workDir*trapz(min(power(indSwing1),0))*ts;   % concentric contractions (swing1)
            swing2AboveYaxis = addwork.workDir*trapz(max(power(indSwing2),0))*ts;   % eccentric contractions (swing2)
            swing2BelowYaxis = addwork.workDir*trapz(min(power(indSwing2),0))*ts;   % concentric contractions (swing2)
            
            % normalize with time (if expressed as a percentage)
            if isfield(addwork, 'meanAbsoluteTime')
                netWork = netWork*(addwork.meanAbsoluteTime/xData(end));
                areaAboveYaxis = areaAboveYaxis*(addwork.meanAbsoluteTime/xData(end));
                areaBelowYaxis = areaBelowYaxis*(addwork.meanAbsoluteTime/xData(end));
                netArea = netArea*(addwork.meanAbsoluteTime/xData(end));
                netAreaTOTAL = netAreaTOTAL*(addwork.meanAbsoluteTime/xData(end));
                stanceAboveYaxis = stanceAboveYaxis*(addwork.meanAbsoluteTime/xData(end));
                stanceBelowYaxis = stanceBelowYaxis*(addwork.meanAbsoluteTime/xData(end));
                swing1AboveYaxis = swing1AboveYaxis*(addwork.meanAbsoluteTime/xData(end));
                swing1BelowYaxis = swing1BelowYaxis*(addwork.meanAbsoluteTime/xData(end));
                swing2AboveYaxis = swing2AboveYaxis*(addwork.meanAbsoluteTime/xData(end));
                swing2BelowYaxis = swing2BelowYaxis*(addwork.meanAbsoluteTime/xData(end));
%                 supT = sprintf('%s: Absolute Time = %.3f s\n', addwork.name, addwork.meanAbsoluteTime);
            end
            
            if ~isempty(findstr(pwd, 'NORMALIZED'))
                units = 'J/kg';
            else
                units = 'J';
            end
            
            oneGaitCycleAboveYaxis = stanceAboveYaxis+swing1AboveYaxis+swing2AboveYaxis;
            oneGaitCycleBelowYaxis = stanceBelowYaxis+swing1BelowYaxis+swing2BelowYaxis;
            oneGaitCycleNet = stanceAboveYaxis+swing1AboveYaxis+swing2AboveYaxis+stanceBelowYaxis+swing1BelowYaxis+swing2BelowYaxis;
            
            % set titles of subplots
%             title_desc = sprintf('%s: %.3f + %.3f == %.3f (%s)', otherstring, areaAboveYaxis, areaBelowYaxis, netArea, units); 
            title_desc = sprintf('%s: %.3f + %.3f == %.3f (%s)', otherstring, oneGaitCycleAboveYaxis, oneGaitCycleBelowYaxis, oneGaitCycleNet, units);
%             if ~exist('hs', 'var')
%                 hs = gca;
%             end
            if writeTitle == 1
                title(hs, sprintf('%s\n%s\n', leg, title_desc), 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'Bold')
%                 set(get(gca, 'Title'), 'String', str);
            elseif writeTitle == 2
                oldtitle = get(get(hs, 'Title'), 'String');
                if size(oldtitle,1) == 3
                    newtitle = sprintf('%s\n%s\n', strtrim(oldtitle(1,:)), strtrim(oldtitle(2,:)), title_desc);
                elseif size(oldtitle,1) == 4
                    newtitle = sprintf('%s\n%s\n%s\n', strtrim(oldtitle(1,:)), strtrim(oldtitle(2,:)), strtrim(oldtitle(3,:)), title_desc);
                end
                set(get(hs, 'Title'), 'String', newtitle);
                set(get(hs, 'Title'), 'FontSize', 8);
            end
            
            % save to array to eventually output to file
            bigWorkCols = [bigWorkCols leg];
            
            bigWorkData = [bigWorkData ...
                [stanceAboveYaxis ; stanceBelowYaxis ; stanceAboveYaxis+stanceBelowYaxis ; ...
                swing1AboveYaxis ; swing1BelowYaxis ; swing1AboveYaxis+swing1BelowYaxis ; ...
                swing2AboveYaxis ; swing2BelowYaxis ; swing2AboveYaxis+swing2BelowYaxis ; ...
                   stanceAboveYaxis+swing1AboveYaxis+swing2AboveYaxis ; stanceBelowYaxis+swing1BelowYaxis+swing2BelowYaxis ; ...
                   stanceAboveYaxis+swing1AboveYaxis+swing2AboveYaxis+stanceBelowYaxis+swing1BelowYaxis+swing2BelowYaxis]];
            
            % print to screen
            if addwork.shading == 1
                fprintf('STANCE PHASE: t = %.2f -> %.2f\n', timeStance(1), timeStance(2));
                fprintf('%s (STANCE) --> %s: %.5f %s\n', leg, aboveBelowLabels{1}, stanceAboveYaxis, units)
                fprintf('%s (STANCE) --> %s: %.5f %s\n', leg, aboveBelowLabels{2}, stanceBelowYaxis, units)
                fprintf('%s (TOTALSTANCE): %.5f %s\n', leg, stanceAboveYaxis+stanceBelowYaxis, units)

                fprintf('FIRST HALF OF SWING PHASE: t = %.2f -> %.2f\n', timeSwing1(1), timeSwing1(2));
                fprintf('%s (SWING1) --> %s: %.5f %s\n', leg, aboveBelowLabels{1}, swing1AboveYaxis, units)
                fprintf('%s (SWING1) --> %s: %.5f %s\n', leg, aboveBelowLabels{2}, swing1BelowYaxis, units)
                fprintf('%s (TOTALSWING1): %.5f %s\n', leg, swing1AboveYaxis+swing1BelowYaxis, units)

                fprintf('SECOND HALF OF SWING PHASE: t = %.2f -> %.2f\n', timeSwing2(1), timeSwing2(2));
                fprintf('%s (SWING2) --> %s: %.5f %s\n', leg, aboveBelowLabels{1}, swing2AboveYaxis, units)
                fprintf('%s (SWING2) --> %s: %.5f %s\n', leg, aboveBelowLabels{2}, swing2BelowYaxis, units)
                fprintf('%s (TOTALSWING2): %.5f %s\n', leg, swing2AboveYaxis+swing2BelowYaxis, units)

                fprintf('FULL CYCLE OF GAIT: t = %.2f -> %.2f\n', timeStance(1), timeSwing2(2));
                fprintf('%s (TOTALGAITCYCLE) --> %s: %.5f %s\n', leg, aboveBelowLabels{1}, oneGaitCycleAboveYaxis, units)
                fprintf('%s (TOTALGAITCYCLE) --> %s: %.5f %s\n', leg, aboveBelowLabels{2}, oneGaitCycleBelowYaxis, units)
                fprintf('%s (TOTALGAITCYCLE): %.5f %s\n\n', leg, oneGaitCycleNet, units)
            end
        end        
        
        
        % plot fibre length / pennation angles
        % (assumes the plotting of muscle fibre length with model as second variable input)
        if ~isempty(model_addfibrelen)
            ofl = model_addfibrelen.Muscles.(strtok(out.labels{j}, '.')).optimal_fiber_length;
            ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
            hline(ofl, 'k:', sprintf('%f\n', ofl));
            graphMin = 0;
            ylim([graphMin, 1.5*ofl]);
            
            if penn
                hold on
                muscWidth = model_addfibrelen.Muscles.(out.labels{j}).optimal_fiber_length * ...
                    sin(model_addfibrelen.Muscles.(out.labels{j}).pennation_angle);
                pennAng = asin(muscWidth./out.data(:,j)) * 180/pi;
                
                if abs(sum(imag(pennAng))) > 0
                   fprintf('Warning: 90deg pennation angle for %s\n', out.labels{j})  
                end
                
                [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData, ...
                    real(pennAng), 'plot', 'plot');
                set(hline1, 'LineStyle', 'none')
                set(hline2, 'LineStyle', '--')
                set(hline1, 'LineWidth', 2)
                set(hline2, 'LineWidth', 1)
                
                hold on
                ofl = model_addfibrelen.Muscles.(out.labels{j}).optimal_fiber_length;
%                 ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
%                 hline(ofl, 'b--', sprintf('%f\n', ofl));
                out.haxes = haxes;

%                 axes(haxes(1));  
                axis tight
                xlim(xl)
                if mod(i-1, pl) == 0
                    ylabel('FibLen (m) / PenAng (deg)', 'FontSize', 8, 'Interpreter', 'none')
                end
                set(gca, 'YTickLabelMode', 'auto')
                set(gca, 'YTickMode', 'auto')
                yl = ylim;
                ylim([0, 1.1*yl(2)]);
                axes(haxes(2));
                axis tight
                xlim(xl)                
                set(gca, 'YTickLabelMode', 'auto')
                set(gca, 'YTickMode', 'auto')
            end
        end


        

%             plotbrowser
        warning off
        if writeTitle && isempty(addwork)
            if ~isempty(C3Dkey)
                title(gca, sprintf('%s\n', strtok(out.labels{j}, '.')), ...
                    'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'Bold')
            else
                title(gca, sprintf('%s\n', strtok(out.labels{j}, '.')), ...
                    'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'Bold')
            end
        end
        warning on

        % Check for EMG to superimpose...
        if ~isempty(EMGdataset)
            EMG_heightRatio = 0.3;
            EMGsingleMusc = getData(sprintf('%s__Processed', strtok(out.labels{j}, '.')), EMGdataset);
            if length(EMGsingleMusc) ~= 1
%                 fprintf('Found %s in EMG dataset!\n', out.labels{j})
                maxEMG = max(EMGsingleMusc);
                GraphLimit = max(ylim) - min(ylim);
                maxEMGline = (GraphLimit*EMG_heightRatio) + min(ylim);
                scaledEMG = EMGsingleMusc .* (GraphLimit*EMG_heightRatio/maxEMG) + min(ylim);
                hold on
                hline(maxEMGline, 'c:');
                plot(timeEMG, scaledEMG, 'Color', [0.678 0.922 1.0])
            end
        end

        if saveFile && ~ALLON1FIG 
            mysave(gcf, sprintf('Figure%s', num2str(gcf)));
            close(100+j)
        end
    
        croptime(timecrop, xl)
        
        
        if ~isempty(C3Dkey)
            plotEventLines(C3Dkey, 'times0', 0, stanceShade);
        end
        
        if xData(1)==0 && xData(end)==100
            set(hs, 'XTick', [0 25 50 75 100]);
        end
        
    end
end   % END MAIN PLOTTING LOOP


if labelMinMax
    yData = out.data(:,j);
    Low = find(yData == min(yData));
    Low = Low(1);

    High = find(yData == max(yData));
    High = High(1);

    hold on
    plot(xData(Low), yData(Low), 'rs', xData(High), yData(High), 'g^');
    text(xData(Low), yData(Low), sprintf('Min: %.2f', yData(Low)), 'FontSize', 8);
    text(xData(High), yData(High), sprintf('Max: %.2f', yData(High)), 'FontSize', 8);
end


if ~isempty(addwork)
%    suptitle(sprintf('Total net work: %.3fJ', netWorkTOTAL));
   bigWorkRows = {'STANCE(+)', 'STANCE(-)', 'STANCE(net)', ...
       'SWING1(+)', 'SWING1(-)', 'SWING1(net)', ...
       'SWING2(+)', 'SWING2(-)', 'SWING2(net)', ...
       'GAITCYCLE(+)', 'GAITCYCLE(-)', 'GAITCYCLE(net)'};

    out.bigWorkRows = bigWorkRows;
    out.bigWorkCols = bigWorkCols;
    out.bigWorkData = bigWorkData;
    out.netAreaTOTAL = netAreaTOTAL;
    out.timeStance = timeStance;
    out.timeSwing1 = timeSwing1;
    out.timeSwing2 = timeSwing2;
    
end


if ~isempty(supT)
    suptitle(supT); 
end


if saveFile && ALLON1FIG 
    mysave(gcf, sprintf('Figure%s', num2str(gcf)));
end


% turn off legends for subplots (all except 1)
% if ALLON1FIG == 1
%     for i = 2:l
%         subplot(pl, pl, i)
%         legend off
%     end
% end




function eventIndex = getEventIndex(labs, e1, e2)

eventIndex(1,1) = getIndex(labs, e1);
eventIndex(1,2) = getIndex(labs(eventIndex(1,1)+1:end), e2) + eventIndex(1,1);

if eventIndex(1,1) == eventIndex(1,2)
    % check that END event label is not the one we are after...
    if eventIndex(1,1) + 4 == length(labs)
        eventIndex(1,2) = length(labs);
    % check that START event label is not the one we are after...
    elseif eventIndex(1,1) - 4 == 1
        eventIndex(1,1) = 1;    
    else
        fprintf('ERROR... can not detect the events from %s to %s\n', e1, e2);
        fprintf('List Of Events:');
        disp(labs)
        [eventIndex, v] = listdlg('PromptString', 'Select 2 events to normalize in between:',...
                              'SelectionMode', 'multiple',...
                              'ListString', labs);
         if v == 0, return; end      % cancel button pressed
    end
end



function croptime(timecrop, xl)

if ~isempty(timecrop)
    for k = 1:2
        if timecrop(k) < 0
            timecrop(k) = xl(k);
        end
    end
    xlim(timecrop);
else
    xlim(xl);
end


