    % Extract and process kinetic data (GRF, CoP, GRMo, GRMx) from a C3D file.
% Tim Dorn
% Last Updated: September 2010
% 
% ----------------------------------------------------------------------
% Usage: [GRF, CoP, GRMo, GRMx] = getKinetics(C3Dkey, plottog*, filterFreq*, markersDyn*)
% ----------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
% 
%           plottog* = toggles the display of various plots
%               0 = no plots (default)
%               1 = plot kinetic graphs only (all in MODEL coordinates)
%               2 = plot kinetic graphs & kinetic verification (all in MODEL coords)
%               4 = plot kinetic verification only (all in MODEL coords)
% 
%           filterFreq* = filter frequency for GRF, GRMo (optional)
%                   < 0 means no filtering is done
%                   Should be used with caution.
% 
%           markersDyn* = optional dynamic markers data structure (from
%                   getMarkers.m) to aid with the verification of extracted
%                   kinetics (used in verifyKinetics.m)
% 
% 
% Notes:    1) Events MUST be labeled in VICON (and hence the c3d file)
%
%           2) Corners must be defined in VICON properly.
%                  Looking down onto the plate from above:
%                       corner 1: where X, Y are both most positive in FP coordinates
%                       corners 2, 3, 4 going CLOCKWISE
%                   (refer to Vicon manual for more information)
% 
%           3) Force plate origins must be defined properly
%                  from force plate true origin to the center of the
%                  plate surface... in FP coordinate system (given in FP manual)
% 
%           4) Output is saved as a *.mot file used by OpenSim
% 
%           5) smoothAll.on is an option (within this function - not an input)
%              which defines the region of filtering the kinetics (if filtering
%              is specified). 
% 
%              smoothAll.on = 1: smooth the whole force/moment trajectory over all time
%              smoothAll.on = 0: smooth the force/moment trajectory over the time
%                             that the foot is on the plate (helps prevent some
%                             under/overshoot prior to foot contact
%
% ----------------------------------------------------------------------
% 
% The output is set up as follows:
% 
% GRF(1,:)  Right Foot --> GRF X             CoP(1,:)  Right Foot --> CoP X
% GRF(2,:)  Right Foot --> GRF Y             CoP(2,:)  Right Foot --> CoP Y
% GRF(3,:)  Right Foot --> GRF Z             CoP(3,:)  Right Foot --> CoP Z
% GRF(4,:)  Left Foot  --> GRF X             CoP(4,:)  Left Foot  --> CoP X
% GRF(5,:)  Left Foot  --> GRF Y             CoP(5,:)  Left Foot  --> CoP Y
% GRF(6,:)  Left Foot  --> GRF Z             CoP(6,:)  Left Foot  --> CoP Z
% 
% GRMo(1,:)  RF -> GRM X about FP origin     GRMx(1,:)  RF -> GRM X about CoP
% GRMo(2,:)  RF -> GRM Y about FP origin     GRMx(2,:)  RF -> GRM Y about CoP
% GRMo(3,:)  RF -> GRM Z about FP origin     GRMx(3,:)  RF -> GRM Z about CoP
% GRMo(4,:)  LF -> GRM X about FP origin     GRMx(4,:)  LF -> GRM X about CoP
% GRMo(5,:)  LF -> GRM Y about FP origin     GRMx(5,:)  LF -> GRM Y about CoP
% GRMo(6,:)  LF -> GRM Z about FP origin     GRMx(6,:)  LF -> GRM Z about CoP
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

function [GRF, CoP, GRMo, GRMx] = getKinetics(C3Dkey, plottog, filterFreq, markersDyn)

usage = 'Usage: [timeVec, GRF, CoP, GRMo, GRMx] = getKinetics(C3Dkey, plottog*, filterFreq*, markersDyn*)';

smoothAll.on = 0;      % smoothAll.on = 1: smooth the whole force/moment trajectory over
                       %                all time
                       % smoothAll.on = 0: smooth the force/moment trajectory over the time
                       %                that the foot is on the plate (helps with
                       %                under/overshoot prior to foot contact

fixCoP = 1;            % Attempt to render CoP of spikes


                    
if nargin == 1,
    plottog = 0;
    filterFreq = 0;
    markersDyn = [];
elseif nargin == 2,
    filterFreq = 0;
    markersDyn = [];
    if plottog > 1,
        plottog = 1;
        disp('No markers given for verification: only plotting graphs (plottog = 1)')
    end
elseif nargin == 3,
    markersDyn = [];
    if plottog > 1,
        plottog = 1;
        disp('No markers given for verification: only plotting graphs (plottog = 1)')
    end
elseif nargin ~= 4,
    disp(usage);
    return
end
if plottog < 0 || plottog > 5 || plottog == 3 || plottog == 5
    plottog = 0;
end

if C3Dkey.allowed.kinetics == 0
   error('C3Dkey was generated with no force plate order input. Can not extract kinetics. Regenerate the C3Dkey with the force plate order variable and try again...\n'); 
end


    


% -------------------------------------------
% Load important parameters from the c3d file
% -------------------------------------------

x = 1; y = 2; z = 3;
numVer = 15;         % num of plots to use when visually verifying
loadLabels;

itf = c3dserver();
openc3d(itf, 0, C3Dkey.c3dFileFULL);
numSeq = length(C3Dkey.sequence.frames(:,1));

if glab.storeInfo == 1 && exist(glab.infoDirectory, 'dir') ~= 7,
    mkdir(glab.infoDirectory);
end



% ------------------------
% Display trial statistics
% ------------------------

fprintf('\n**********************************************\n')
fprintf('DATA OUTPUT AT %d HZ (ANALOG FREQUENCY)\n', C3Dkey.aFreq) 
fprintf('**********************************************\n')

fprintf('Trial Type: %s\n', C3Dkey.trialType)
fprintf('Processing Dynamic Task File: %s\n', C3Dkey.c3dFile)
fprintf('Subject Name: %s\n\n', C3Dkey.name)

fprintf('File Stats (uncropped): %d analog frames (%d video frames), ratio = %.1f\n\n', ...
    C3Dkey.numFrames.uncroppedA, C3Dkey.numFrames.uncroppedV, C3Dkey.r);

fprintf('Motion Capture Frequency: %d Hz\n', C3Dkey.vFreq)
fprintf('Motion Capture Period:    %.5f sec\n\n', 1/C3Dkey.vFreq)

fprintf('Analog Capture Frequency: %d Hz\n', C3Dkey.aFreq)
fprintf('Analog Capture Period:    %.5f sec\n\n', 1/C3Dkey.aFreq)

fprintf('File Stats (cropped to event cycle): %d analog frames (%d video frames), ratio = %.1f\n\n', ...
    C3Dkey.numFrames.croppedA, C3Dkey.numFrames.croppedV, C3Dkey.r);

fprintf('------------------\n')
fprintf('EVENT DESCRIPTIONS\n')
fprintf('------------------\n')
numEvents = length(C3Dkey.event.txt);

for i = 1:numEvents
    timeCycle = C3Dkey.event.times(i)-C3Dkey.event.times(1);
    
    eventTxt = sprintf('%s: %.4f%% of labeled event cycle  \t(%.4f sec, #%d ANALOG (#%d VIDEO)\n', ...
        C3Dkey.event.txt{i}, C3Dkey.event.percent(i), timeCycle, C3Dkey.event.Aframes(i), C3Dkey.event.Vframes(i));
    fprintf('%s', eventTxt)
end
fprintf('\n');



% ------------------------------------------------------------
% Prepare Analog Channels for Extraction (GAIT LAB DEPENDANT!)
% ------------------------------------------------------------
lp = length(glab.FP.suffix);
for i = 1:C3Dkey.numPlatesUsed
   for j = 1:lp
      plateLabel{(i-1)*lp+j} = sprintf(glab.FP.string, ...
          glab.FP.prefix{j}, abs(C3Dkey.FP_order(i)), glab.FP.suffix{j});
   end
   plateSign(6*(i-1)+1:6*i) = plateDirVec(C3Dkey.FP_order(i));
end

if isempty(C3Dkey.FP_order),
    error(sprintf('No force plate order defined in the event key.\nRecreate the key using a force plate order input: C3Dkey = getEvents(c3dFile, direction, FP_order)'));
end

fprintf('Number of Force Plates Available in Trial: %d\n', C3Dkey.numPlatesTotal)
fprintf('Number of Force Plates Used For Extraction: %d\n', C3Dkey.numPlatesUsed)
disp(['Force Plate Order: ', num2str(C3Dkey.FP_order)])
fprintf('\n');



% ----------------------------------------------------
% Retrieve physical force plate data from the c3d file
% and optionally display to the screen
% ----------------------------------------------------
FP = calcFPData(itf, C3Dkey);

header = '      X         Y         Z';
for i = 1:C3Dkey.numPlatesUsed
    disp('-----------------------------------------------')
    disp(['-------- Stepping Order PLATE #', num2str(i), ' (FP', ...
        num2str(C3Dkey.FP_order(i)), ') --------'])
    disp('-----------------------------------------------')
    disp(['XYZ corner coordinates (GLOBAL VICON COORDS) for plate #', num2str(i)])
    disp(header)
    disp(FP.corners(:,:,i))
    disp(['XYZ coordinates (GLOBAL VICON COORDS) for CENTER of plate #', num2str(i)])
    disp(header)
    disp(FP.viconOrig2FPcenterSurfaceGLOB(i,:))
    disp(['XYZ coordinates (GLOBAL VICON COORDS) for ORIGIN of plate #', num2str(i)])
    disp(header)
    disp(FP.viconOrig2FPorigGLOB(i,:))
end
   


% --------------------------------------------------------
% Determine force plate sequences from the labelled events
% --------------------------------------------------------
fprintf('---------------------------------------------\n')
fprintf('  Force Plate Foot Sequences from Event Key  \n')
fprintf('---------------------------------------------\n')
fprintf('   FRAME NUMBERS           PLATE NUMBER      \n')         
fprintf('   start    end          rightF    leftF     \n')
fprintf('---------------------------------------------\n')

for i = 1:numSeq
    fprintf('%7d\t%7d\t\t||\t%5d\t%5d\n', ...
        C3Dkey.sequence.frames(i,1), C3Dkey.sequence.frames(i,2), ...
        C3Dkey.sequence.plates(i,1), C3Dkey.sequence.plates(i,2));
end



% ------------------------------------------------------------------------
% Prepare and extract all required analog parameters from the c3d file
% at the analog frame rate. This will come in handy for CoP and GRM(CoP)
% processing.
% 
% Forces are extracted in N
% Moments are extracted in Nmm and converted to Nm
% i.e. F(direction, frame, plate)    where X=1, Y=2, Z=3
% ------------------------------------------------------------------------
fprintf('\nExtracting GRF, CoP, GRMo, GRMx from %s...\n\n', C3Dkey.c3dFile);

if isfield(glab.FP, 'filterOrder')
    filterOrder = glab.FP.filterOrder;
else
    filterOrder = 4;        % default value
end

for i = 1:C3Dkey.numPlatesUsed      % plate loop
    for j = 1:3                     % direction loop
        
        % F_all(direction, frame, plate in stepping order), same for Mo_all
        F_all_raw(j,:,i) = plateSign(6*(i-1) + j) * ...
            double(getanalogchannel(itf, plateLabel{6*(i-1) + j}, ...
            C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end)));
        
        Mo_all_raw(j,:,i) = plateSign(6*(i-1) + j+3) * ...
            double(getanalogchannel(itf, plateLabel{6*(i-1) + j+3}, ...
            C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end)) / C3Dkey.divide_to_meters);
        
        
        if smoothAll.on
            F_all(j,:,i) = plateSign(6*(i-1) + j) * ...
                smooth(double(getanalogchannel(itf, plateLabel{6*(i-1) + j}, ...
                C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end))), ...
                filterFreq, C3Dkey.aFreq, filterOrder);

            Mo_all(j,:,i) = plateSign(6*(i-1) + j+3) * ...
                smooth(double(getanalogchannel(itf, plateLabel{6*(i-1) + j+3}, ...
                C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end)) / C3Dkey.divide_to_meters), ...
                filterFreq, C3Dkey.aFreq, filterOrder);
            
        elseif ~smoothAll.on           % assumes only one foot per plate
            smoothAll.startPriorViconFrames = 0;    % start filtering this many 
                                                    % Vicon frames before the event
            
            frm = C3Dkey.event.Aframes0(C3Dkey.sequence.eventIndex(i,1)) - smoothAll.startPriorViconFrames*C3Dkey.r : ...
                   C3Dkey.event.Aframes0(end);

            F_all(j,:,i) = zeros(size(F_all_raw(j,:,i)));
            F_all(j,frm,i) = plateSign(6*(i-1) + j) * ...
                smooth(double(getanalogchannel(itf, plateLabel{6*(i-1) + j}, ...
                C3Dkey.event.Vframes(C3Dkey.sequence.eventIndex(i,1)) - smoothAll.startPriorViconFrames, ...
                C3Dkey.event.Vframes(end))), ...
                filterFreq, C3Dkey.aFreq, filterOrder);

            Mo_all(j,:,i) = zeros(size(Mo_all_raw(j,:,i)));
            Mo_all(j,frm,i) = plateSign(6*(i-1) + j+3) * ...
                smooth(double(getanalogchannel(itf, plateLabel{6*(i-1) + j+3}, ...
                C3Dkey.event.Vframes(C3Dkey.sequence.eventIndex(i,1)) - smoothAll.startPriorViconFrames, ...
                C3Dkey.event.Vframes(end)) / C3Dkey.divide_to_meters), ...
                filterFreq, C3Dkey.aFreq, filterOrder);        
        end
        
        

    end
end

if filterFreq > 0
    fprintf('Smoothall.on = %d\nLow pass filtered (order=%d Butterworth) GRF, GRMo @ %d Hz\n', ...
        smoothAll.on, filterOrder, filterFreq);
else
    fprintf('No kinetic filtering performed\n')
end
closec3d(itf);




% --------------------------------------------------
% Calculate CoP & GRMx data from extracted variables
% --------------------------------------------------
% Finds the global CoP (VICON COORDS) for each plate in each direction
% Finds the GRM about CoP for each plate in each direction
% 
% referenced from: http://www.kwon3d.com/theory/grf/cop.html
% --------------------------------------------------------------------
warning off MATLAB:divideByZero
for i = 1:C3Dkey.numPlatesUsed       % Force Plate Loop (i)
      
    % Determine local force plate origin values (FP coordinates)
    % This comes from the individual callibration of each force plate
    % According to C3Dmanual (p81), The ORIGIN parameter contains
    % the vector from the origin of the FP coordinate system to the point 
    % at the geometric center of the working surface of the force platform. 
    % However the set of equations here assume the opposite: 
    % FP.localOriginsFP points from the geometric center of the FP to the
    % true origin and is expressed in the force plate coordinate system.
    % = FP_GEOMETRIC_CENTER -> FP_ORIGIN 
    % ---------------------------------------------------------------
    a = FP.localOriginsFP(i,x);
    b = FP.localOriginsFP(i,y);
    c = FP.localOriginsFP(i,z);
    
    
    % Determine vector from FP true origin to CoP in FORCE PLATE coordinates
    % vecA::  FP_ORIGIN -> CoP (FP COORDS)
    % -----------------------------------------------------------------
    vecA(z,:) = -ones(1, C3Dkey.numFrames.croppedA) * c;
    vecA(x,:) = ((-Mo_all(y,:,i) - (c * F_all(x,:,i))) ./ F_all(z,:,i)) + a;
    vecA(y,:) = (( Mo_all(x,:,i) - (c * F_all(y,:,i))) ./ F_all(z,:,i)) + b;

    
    % Convert now to VICON coordinate system
    % FPorig2CoP:: FP_ORIGIN -> CoP (VICON COORDS)
    % --------------------------------------------
    FPorig2CoP = coordChange(vecA, C3Dkey.transform.FPVICON{i});
    % A small check to see if the force plate origin is correct.
    % FPorig2CoP = zeros(3, C3Dkey.numFrames.croppedA);
    
    
    % Calculate CoP from VICON ORIGIN (all in VICON coordinates)
    % X_all:: VICON_ORIGIN -> FP_GEOMETRIC_CENTER -> FP_ORIGIN (VICON COORDS) 
    %          + FP_ORIGIN -> CoP                              (VICON COORDS)
    % (this will be converted to MODEL coordinates later on)
    % ----------------------------------------------------------
    for j = 1:3               % Coordinate Loop  (j)
        X_all(j,:,i) = ones(1, C3Dkey.numFrames.croppedA) ...
            * FP.viconOrig2FPorigGLOB(i,j) + FPorig2CoP(j,:);
    end
    
    
    % Calculate Moments about CoP (FP COORDS)
    % ---------------------------------------
    
    % Tx: Should be 0 by definition
    % Tx = Mxo -c.Fy - (Ycop - b)Fz
    Mx_all(x,:,i) = Mo_all(x,:,i) ...
        - ( (c * F_all(y,:,i)) ) ...
        - ( (vecA(y,:) - b) .* F_all(z,:,i) );
    
    % Ty: Should be 0 by definition
    % Ty = Myo +c.Fx + (Xcop - a)Fz
    Mx_all(y,:,i) = Mo_all(y,:,i) ...
        + ( (c * F_all(x,:,i)) ) ...
        + ( (vecA(x,:) - a) .* F_all(z,:,i) );
    
    % Tz: Free moment is nonzero
    % Tz = Mzo + (Ycop - b)Fx - (Xcop - a)Fy
    Mx_all(z,:,i) = Mo_all(z,:,i) ...
        + ( (vecA(y,:) - b) .* F_all(x,:,i) ) ...
        - ( (vecA(x,:) - a) .* F_all(y,:,i) );  
    

    % ========================================================================
    % Rigid Body transformations for GRF, GRMo and GRMx from FP to MODEL
    % ========================================================================
    % Currently at this stage, the directly extracted variables (GRF & GRMo)
    % are expressed in the FORCE PLATE coordinate system because that is how
    % they were originally recorded during the trial. GRMx was also calculated
    % based on data in the FORCE PLATE coordinate system.
    % C3Dkey.transform.FPMODEL{i} converts from the FP to MODEL coordinate system
    % for force plate i
    % ========================================================================

    % Convert from FP to MODEL coordinate system
    % ------------------------------------------
    F_all(:,:,i) = coordChange(F_all(:,:,i), C3Dkey.transform.FPMODEL{i});
    F_all_raw(:,:,i) = coordChange(F_all_raw(:,:,i), C3Dkey.transform.FPMODEL{i});
    Mo_all(:,:,i) = coordChange(Mo_all(:,:,i), C3Dkey.transform.FPMODEL{i});
    Mo_all_raw(:,:,i) = coordChange(Mo_all_raw(:,:,i), C3Dkey.transform.FPMODEL{i});
    Mx_all(:,:,i) = coordChange(Mx_all(:,:,i), C3Dkey.transform.FPMODEL{i});

end
warning on MATLAB:divideByZero



% -----------------------------------------------------------
% Extract GRF CoP, GRMo and GRMx data for each foot from the 
% c3d file according to the generated force plate sequences
% defined in the event key (from getEvents.m).
% -----------------------------------------------------------

for i = 1:numSeq
    startF = C3Dkey.sequence.frames(i,1);    % analog frame numbers
    endF = C3Dkey.sequence.frames(i,2);
    
    for foot = 1:2
        lines = (foot-1)*3 + 1 : (foot-1)*3 + 3;
        plate = C3Dkey.sequence.plates(i,foot);   % actual vicon labelled FP
        
        if plate == 0,
            GRF(lines, startF:endF)  = zeros(3, endF-startF+1);
            CoP(lines, startF:endF)  = zeros(3, endF-startF+1);
            GRMo(lines, startF:endF) = zeros(3, endF-startF+1);
            GRMx(lines, startF:endF) = zeros(3, endF-startF+1);
            
        else
            GRF(lines, startF:endF) = squeeze(F_all(1:3, startF:endF, C3Dkey.FP_order_inv(plate)));
            CoP(lines, startF:endF) = squeeze(X_all(1:3, startF:endF, C3Dkey.FP_order_inv(plate)));
            GRMo(lines, startF:endF) = squeeze(Mo_all(1:3, startF:endF, C3Dkey.FP_order_inv(plate)));
            GRMx(lines, startF:endF) = squeeze(Mx_all(1:3, startF:endF, C3Dkey.FP_order_inv(plate)));
        end
    end  
end



% Detect NaNs in CoP, GRMx (which could come about by 
% dividing by a zero vertical GRF)
% ------------------------------------------------------
for i = 1:2
    switch i
        case 1
            tmpData = CoP;
            tmpTxt = 'CoP';
            tmpLab = glab.X;
        case 2
            tmpData = GRMx;
            tmpTxt = 'GRMx';
            tmpLab = glab.Mx;
    end
    
    [xnan, ynan] = findNaN(tmpData);
    for j = 1:length(xnan),
        fprintf('NaN found on: %s at Matlab frame %d (Vicon: %d)--> %s(%d, %d)\n', ...
            tmpLab{ynan(j)}{1}, xnan(j), ...
            floor(C3Dkey.event.Vframes(1)+(xnan(j)/C3Dkey.r)), ...
            tmpTxt, ynan(j), xnan(j));
    end

    if isempty(xnan),
        fprintf('***** No NaNs found at any time in %s data :) *****\n', tmpTxt);
    else
        fprintf('%d NaNs found in %s data\nWARNING: OpenSim may not like this very much...\n\n', length(xnan), tmpTxt);
    end
end
%pause(3);  % time to visually inspect matlab output to check that there are no NaNs


% First attempt to fix rapid gradients in CoP
% (I have noticed that the first and last nonzero frame of the CoP spike
% perhaps due to the low vertical GRF value, so here I take the first and
% last non zero frame (frame n) and if gradients on either side of that
% point are opposite in sign, then set it equal to the average of the
% previous (n-1) and next (n+1) frame). Repeat in an inner/outer loop.
% ----------------------------------------------------------------------------
if fixCoP
    m = size(CoP,1);
    CoPorig = CoP;
    CopRenderOuter = 100;
    CopRenderInner = 3*C3Dkey.r;

    % For each START interval of CoP(a,b), set CoP(a-1) = CoP(a) --> (vicon frames).
    % zero order hold over the analog frame interval.
    f1 = 2;     % analog frames
    for i = 1:m
        p = findPotentialSpikePoints(CoP(i,:), 0);
        for j = p
           if CoP(i, j-1) == 0      % ensure we are at the START of the interval
               CoP(i, j-C3Dkey.r:j-1+f1) = CoP(i, j+f1+1);
           end
        end
    end
    
    
    fprintf('Render CoP ON: OuterLoop = %d, InnerLoop = %d\n', CopRenderOuter, CopRenderInner);
    for i = 1:CopRenderOuter
        for j = 1:CopRenderInner
    %         fprintf('\nCoP Render (Inner: %d\\%d), (Outer: %d\\%d)\n', ...
    %             j, CopRenderInner, i, CopRenderOuter);
            CoP = renderCoP(CoP, j-1, i);
        end
    end

    
    % For each END interval of CoP(a,b), set CoP(b+1) = CoP(b) --> (vicon frames).
    % zero order hold over the analog frame interval.
    f2 = 2;     % analog frames
    for i = 1:m
        p = findPotentialSpikePoints(CoP(i,:), 0);
        for j = p
           if CoP(i, j+1) == 0      % ensure we are at the END of the interval
               CoP(i, j-1-f2:j+C3Dkey.r) = CoP(i, j-f2-1);
           end
        end
    end

    
    % Plot old/new CoP after rendering (only if plottog ~= 0)
    if plottog ~= 0
        figure('Name', 'Original/Modified CoP', 'NumberTitle', 'off')
        maximize;
        coptitl = {'Right X CoP', 'Right Y CoP', 'Right Z CoP', 'Left X CoP', 'Left Y CoP', 'Left Z CoP'};

        for i = 1:m
           subplot(2,3,i)
           hold on
           title(sprintf('%s\nFP Coordinates', coptitl{i}))
           plot(CoPorig(i,:))
           plot(CoP(i,:), 'r')
        end

        if glab.storeInfo == 1
            saveas(gcf, sprintf('%sModified_CoP_%s.emf', glab.infoDirectory, ...
                C3Dkey.c3dFile), 'emf');
            saveas(gcf, sprintf('%sModified_COP_%s.fig', glab.infoDirectory, ...
                C3Dkey.c3dFile), 'fig');
        end
    end

end




% Verify kinetic data (Vicon coordinates)
% Removed this implementation (TD Sept 2010)
% ----------------------------------------------
% if plottog == 3 || plottog == 5,
%     % need to first convert the marker positions into the vicon
%     % coordinate system if they are given
%    
%     if isstruct(markersDyn),
%        markerpos = markersDyn.data(:, 2:end);
%        n = (length(markerpos(1,:))) / 3;
% 
%        % get transform inverse because markers here are in model coords
%        % and need to be plotted in Vicon coords
%        for k = 1:3,
%            tmp = markersDyn.dirVecVICMODEL(k);
%            dirVecInv(abs(tmp)) = sign(tmp)*k;
%        end
% 
%        for i = 1:n,
%            for j = 1:3,
%                markerposVicon(:, (i-1)*3 + j) = ...
%                    markerpos(:, (i-1)*3 + abs(dirVecInv(j))) * sign(dirVecInv(j));
%            end
%        end
%        markersDyn.data(:, 2:end) = markerposVicon;
%     end
%     
%     % convert temporary GRF from MODEL to VICON (CoP is already in VICON)
%     GRF = coordChange(GRF, C3Dkey.transform.MODELVIC);
% 
%     % now verify the data extraction all in vicon coordinate system
%     C3Dkey.coords = 'vicon';
%     verifyKinetics(numVer, C3Dkey, [C3Dkey.timeVec.analogFrame ; GRF], ...
%         [C3Dkey.timeVec.analogFrame ; CoP], FP, 1, markersDyn, [0 0 1]);
%     close all
% end





% During the calculation of the GLOBAL CoP, things were converted to the 
% VICON coordinate system. C3Dkey.transform.VICMODEL converts from the 
% VICON to MODEL coordinate system.
% ---------------------------------------------
CoP = coordChange(CoP, C3Dkey.transform.VICMODEL);





% ======================================================================
% Save, plot, and verify kinetic data (MODEL coordinates)
% ======================================================================

if plottog == 2 || plottog == 4,
    FP.corners = FP.corners_model;

    % verify the data extraction all in model coordinate system
    C3Dkey.coords = 'model';
    verifyKinetics(numVer, C3Dkey, [C3Dkey.timeVec.analogFrame ; GRF], ...
        [C3Dkey.timeVec.analogFrame ; CoP], FP, 1, markersDyn, [0 1 0]);
end
   
if plottog > 0 && plottog < 4,
    % plot the kinetic graphs
    fprintf('\nCreating Plots...\n')
    plotPlateKinetics(C3Dkey, F_all, F_all_raw, Mo_all, Mo_all_raw, ...
        filterFreq, glab.storeInfo, glab.infoDirectory, smoothAll);
    plotFootKinetics(C3Dkey, GRF, CoP, GRMo, GRMx, glab.storeInfo, glab.infoDirectory);  
end

timeVec = C3Dkey.timeVec;

if glab.storeInfo == 1,
    f1 = sprintf('%sKinetics_%s.mat', glab.infoDirectory, C3Dkey.c3dFile);
    f2 = sprintf('%sForcePlateData_%s.mat', glab.infoDirectory, C3Dkey.c3dFile);

    save(f1, 'timeVec', 'GRF', 'CoP', 'GRMo', 'GRMx', 'F_all', 'X_all', 'Mo_all', 'Mx_all')
    save(f2, 'FP')
end



% Save to *.mot file for use with OpenSim
% ---------------------------------------

% first offset the cop according to C3Dkey.offset
offset = [C3Dkey.offset/C3Dkey.divide_to_meters, ...
    C3Dkey.offset/C3Dkey.divide_to_meters];    % convert to m

for i = 1:6
    CoP_offset(i,:) = CoP(i,:) + offset(i);
end
fprintf('X CoP offset: %d mm\nY CoP offset: %d mm\nZ CoP offset: %d mm\n', ...
    C3Dkey.offset(1), C3Dkey.offset(2), C3Dkey.offset(3));

% now save to kinetics file
file = sprintf('%s_kinetics', C3Dkey.c3dFile);
bigM = [timeVec.Asec', GRF(1:3, :)', CoP_offset(1:3, :)', GRF(4:6, :)', CoP_offset(4:6, :)', GRMx'];
colnamesKinetics = {'time', 'ground_force_vx', 'ground_force_vy', 'ground_force_vz', ...
    'ground_force_px', 'ground_force_py', 'ground_force_pz', '1_ground_force_vx', ...
    '1_ground_force_vy', '1_ground_force_vz', '1_ground_force_px', '1_ground_force_py', ...
    '1_ground_force_pz', 'ground_torque_x', 'ground_torque_y', 'ground_torque_z', ...
    '1_ground_torque_x', '1_ground_torque_y', '1_ground_torque_z'};
generateMotFile(bigM, colnamesKinetics, sprintf('%s.mot', file));

% Generate kinetics XML file to be associated with this mot file
%generateKineticsXML(file);    % removed by Prasanna Sritharan, January 2016



% Save coordinates file *.mot for use with OpenSim
% (coordinates have zeros by default). If you have actual coordinate
% values, change them manually in the output mot file.
% ------------------------------------------------------------------

% the coordinates file is reported in vicon freq (in OpenSim) so the
% analog kinetics must be resampled at vicon frequency before saving.
vGRF = []; vCoP_offset = []; vGRMx = [];
if C3Dkey.r == 1,      % vicon freq = analog freq
    vGRF = GRF;
    vCoP_offset = CoP_offset;
    vGRMx = GRMx;
else
    for i = 1:length(GRF(1,:))    % scale time to video frames
        if mod(i, C3Dkey.r) == 1,
           vGRF = [vGRF,  GRF(:,i)];
           vCoP_offset = [vCoP_offset,  CoP_offset(:,i)];
           vGRMx = [vGRMx,  GRMx(:,i)];
        end
    end
end

colnamesJoints = {'time'};
datarows = length(timeVec.Vsec');
numJoints = length(glab.(glab.jointModel));

file = sprintf('%s_coordinates.mot', C3Dkey.c3dFile);
bigM = [timeVec.Vsec', zeros(datarows, numJoints), vGRF(1:3, :)', ...
    vCoP_offset(1:3, :)', vGRF(4:6, :)', vCoP_offset(4:6, :)', vGRMx'];
for i = 1:numJoints
    colnamesJoints(i+1) = glab.(glab.jointModel){i}(1);
end
colnamesTotal = [colnamesJoints, colnamesKinetics(2:end)];
generateMotFile(bigM, colnamesTotal, file);
fclose all;

% diary off













% ========================================================================
% SUBFUNCTION: plotPlateKinetics
% ========================================================================

function plotPlateKinetics(C3Dkey, F_all, F_all_raw, M_all, M_all_raw, ...
    filterFreq, savePlotAsFile, savePlotPath, smoothAll)

% do not change this value...
t = C3Dkey.timeVec.analogFrame;    % Set time axis 
                            % (timeVec.c3dAnalogFrame, timeVec.c3dVideoFrame,
                            %  timeVec.analogFrame, timeVec.videoFrame,
                            %  timeVec.Asec, timeVec.Vsec
                            %  timeVec.Apercent, timeVec.Vpercent)

lw = 2;                 % Set plot line width
titlsize = 14;          % Title font size
axissize = 12;          % Axis font size


% Plot extracted data
% -------------------
frames = C3Dkey.event.Aframes0;

% time vector is interpolated to analog frequency length
t = resamp2(t, length(C3Dkey.timeVec.analogFrame));
dirLabel = ['X', 'Y', 'Z'];
shadeCol = ['y' 'g'];


% Plot GRF/GRM PLATE Data
% -----------------------
for h = 1:2
    switch h
        case 1
            plotName = 'GRF';
            dataAllPlates = F_all;
            dataAllPlates_raw = F_all_raw;
        case 2
            plotName = 'GRMo';
            dataAllPlates = M_all;
            dataAllPlates_raw = M_all_raw;
    end
    
    figure('Name', sprintf('%s All Plates: %s', plotName, C3Dkey.c3dFile), ...
        'NumberTitle', 'off')
    maximize;
    
    
    for i = 1:C3Dkey.numPlatesUsed      % plate loop
        for j = 1:3                     % direction loop

            % note: subplot fills rows, then columns
            subplot(C3Dkey.numPlatesUsed, 3, (i-1)*3+j)
            hold on
            if filterFreq > 0
                plot(t, dataAllPlates_raw(j, :, C3Dkey.FP_order_inv(abs(C3Dkey.FP_order(i)))), ...
                    'LineWidth', lw, 'Color', [0.678 0.922 1.0])
            end
            
            plot(t, dataAllPlates(j, :, C3Dkey.FP_order_inv(abs(C3Dkey.FP_order(i)))), ...
                'b', 'LineWidth', lw)

            xlim([C3Dkey.timeVec.analogFrame(1), C3Dkey.timeVec.analogFrame(end)])
            ylim('manual')

            % shade areas used for each foot
            for k = 1:length(C3Dkey.sequence.plates(:,1))
                for l = 1:2
                    if C3Dkey.sequence.plates(k,l) == abs(C3Dkey.FP_order(i)),
                        ShadePlotForEmpahsisVert(C3Dkey.sequence.frames(k,:), shadeCol(l), 0.3);
                    end
                end
            end

            if ~smoothAll.on
               plot(C3Dkey.event.Aframes0(C3Dkey.sequence.eventIndex(i,1)) - smoothAll.startPriorViconFrames*C3Dkey.r, 0, 'rx') 
            end
            
            ylabel('N', 'FontSize', axissize)
            titl = sprintf('%s Plate %d   (FP%d) -- %s', plotName, i, ...
                abs(C3Dkey.FP_order(i)), dirLabel(j)); 
            
%             ['GRF Plate ', num2str(i), '  (FP', ...
%                 num2str(abs(C3Dkey.FP_order(i))), ') -- ', dirLabel(j)];
            title(sprintf('%s\nModel\n', titl), 'FontSize', titlsize) 

            % draw some dotted lines to mark labeled events
            plotEventLines(C3Dkey, 'Aframes0', 0);
            hline(0, 'k:', ''); 
        end
    end
    
    if savePlotAsFile == 1,
        saveas(gcf, sprintf('%sallPlates%s_%s.emf', savePlotPath, plotName, ...
            C3Dkey.c3dFile), 'emf');
        saveas(gcf, sprintf('%sallPlates%s_%s.fig', savePlotPath, plotName, ...
            C3Dkey.c3dFile), 'fig');
    end
    

    % Create textbox legend #1  (yellow: right foot)
    annotation(gcf,'textbox','String',{'Right Foot'},'FontSize',12,...
        'FitHeightToText','on',...
        'BackgroundColor',[1 1 0],...
        'Position',[0.016 0.9 0.06 0.068]);

    % Create textbox legend #2  (green: left foot)
    annotation(gcf,'textbox','String',{'Left Foot'},'FontSize',12,...
        'FitHeightToText','on',...
        'BackgroundColor',[0 1 0],...
        'Position',[0.016 0.8 0.06 0.068]);

    % Create textbox legend #3  (blue: filter frequency)
    annotation(gcf,'textbox','String',{sprintf('Filt Freq = %d Hz', filterFreq)}, ...
        'FontSize',12,...
        'FitHeightToText','on',...
        'BackgroundColor',[0.678 0.922 1.0],...
        'Position',[0.016 0.7 0.06 0.068]);

end





% ========================================================================
% SUBFUNCTION: generateKineticsXML
% ========================================================================

function generateKineticsXML(file)

forcename = {'ExternalForce_1', 'ExternalForce_2'};
bodyname = {'calcn_r', 'calcn_l'};
ind = [0 0 ; 6 3];
root = 'ForceSet';
data.ATTRIBUTE.name = file;

for i = 1:2
    data.objects.PrescribedForce(i).ATTRIBUTE.name = forcename{i};
    data.objects.PrescribedForce(i).body = bodyname{i};
    data.objects.PrescribedForce(i).pointIsGlobal = 'false';
    data.objects.PrescribedForce(i).forceIsGlobal = 'true';
    
    k=1;
    data.objects.PrescribedForce(i).FunctionSet(1).ATTRIBUTE.name = 'forceFunctions';
    for j = 1:3
        data.objects.PrescribedForce(i).FunctionSet(1).objects.GCVSpline(k).ATTRIBUTE.name = sprintf('#%d', j+ind(i,1));
        data.objects.PrescribedForce(i).FunctionSet(1).objects.GCVSpline(k).half_order = 0;
        data.objects.PrescribedForce(i).FunctionSet(1).objects.GCVSpline(k).error_variance = 0;
        k=k+1;
    end
    
    k=1;
    data.objects.PrescribedForce(i).FunctionSet(2).ATTRIBUTE.name = 'pointFunctions';
    for j = 4:6
        data.objects.PrescribedForce(i).FunctionSet(2).objects.GCVSpline(k).ATTRIBUTE.name = sprintf('#%d', j+ind(i,1));
        data.objects.PrescribedForce(i).FunctionSet(2).objects.GCVSpline(k).half_order = 0;
        data.objects.PrescribedForce(i).FunctionSet(2).objects.GCVSpline(k).error_variance = 0;
        k=k+1;
    end
    
    k=1;
    data.objects.PrescribedForce(i).FunctionSet(3).ATTRIBUTE.name = 'torqueFunctions';
    for j = 13:15
        data.objects.PrescribedForce(i).FunctionSet(3).objects.GCVSpline(k).ATTRIBUTE.name = sprintf('#%d', j+ind(i,2));
        data.objects.PrescribedForce(i).FunctionSet(3).objects.GCVSpline(k).half_order = 0;
        data.objects.PrescribedForce(i).FunctionSet(3).objects.GCVSpline(k).error_variance = 0;
        k=k+1;
    end
    
end
data.datafile = sprintf('%s.mot', file);

fileName = sprintf('%s.xml', file);
xmlString = xml_formatany(data, root);
fid = fopen(fileName, 'w');
if fid < 0
    fprintf('\nERROR: %s could not be opened for writing...\n\n', fileName);
    return
end
fprintf(fid, '%s\n', xmlString);
fclose(fid);
% xml_save(fileName, data, 'off');
fprintf('%s saved sucessfully...\n', fileName);






% ========================================================================
% SUBFUNCTION: plotFootKinetics
% ========================================================================

function plotFootKinetics(C3Dkey, GRF, CoP, GRMo, GRMx, savePlotAsFile, savePlotPath)

% if t is changed, make sure you change the plotEventLines input below.
t = C3Dkey.timeVec.c3dVideoFrame;    % Set time axis 
                            % (timeVec.c3dAnalogFrame, timeVec.c3dVideoFrame,
                            %  timeVec.analogFrame, timeVec.videoFrame,
                            %  timeVec.Asec, timeVec.Vsec
                            %  timeVec.Apercent, timeVec.Vpercent)

lw = 3;                 % Set plot line width
opt = 'b';              % Plotting options
titlsize = 14;          % Title font size
axissize = 12;          % Axis font size


% Plot extracted data
% -------------------
frames = C3Dkey.event.Aframes0;

% time vector is interpolated to analog frequency length
t = resamp2(t, length(C3Dkey.timeVec.analogFrame));


% Plot FOOT Data
% --------------

loadLabels;
for i = 1:6
    lab{i} = glab.S{i};
    lab{i+6} = glab.X{i};
    lab{i+12} = glab.Mo{i};
    lab{i+18} = glab.Mx{i};
end

for i = 1:4
    figure('Name', sprintf('%s: %s', glab.name{i}, C3Dkey.c3dFile), 'NumberTitle', 'off')
    maximize;
    
    switch i
        case 1
            data = GRF;
        case 2
            data = CoP;
        case 3
            data = GRMo;
        case 4
            data = GRMx;
    end
    
    for j = 1:6
        subplot(2,3,j)
        
        plot(t, data(j,:), opt, 'LineWidth', lw)    

        hold on     
        titl = lab{6*(i-1) + j}{1};
        title(sprintf('%s\nModel\n', titl), 'FontSize', titlsize) 
        xlabel(' ', 'FontSize', axissize)
        ylabel(lab{6*(i-1) + j}{2}, 'FontSize', axissize)
        xlim([t(1), t(end)])
                
        % scale vertical axis for X & Z GRMx to [-0.1 0.1]
        if i == 4 && (mod(j,3) == 0 || mod(j,3) == 1),
            ylim([-0.1 0.1])
        else
            ylim('manual')
        end
        
        % Now draw some dotted lines to mark labeled events... 
        plotEventLines(C3Dkey, 'Vframes', 0);
        hline(0, 'k', '');              
    end
    
    if savePlotAsFile == 1,  
        saveas(gcf, sprintf('%s%s_%s.emf', savePlotPath, glab.name{i}, ...
            C3Dkey.c3dFile), 'emf');
        saveas(gcf, sprintf('%s%s_%s.fig', savePlotPath, glab.name{i}, ...
            C3Dkey.c3dFile), 'fig');
    end
end



% ========================================================================
% SUBFUNCTION: plateDirVec
% ========================================================================
% To manipulate further the sign of the plate. Not used (always equal to 1)

function plateSign = plateDirVec(plateNum)

if sign(plateNum) > 0
    plateSign = [1 1 1 1 1 1];
else
    plateSign = [-1 -1 1 -1 -1 1];
end



% ========================================================================
% SUBFUNCTION: renderCoP
% ========================================================================
% Render the CoP

function CoP = renderCoP(CoP, shift, outerIter)

m = size(CoP,1);
for i = 1:m         % rows
    p = findPotentialSpikePoints(CoP(i,:), shift);
    
%     remove foot off spike in AP direction from potential spike points
%     if i == 2 || i == 5
%         p = p(1:2:length(p));
%     end
    
%     debug step
%     if i == 5
%         fprintf('shift = %d      ', shift); disp(p)
%         figure(1)
%         hold on
%         plot(CoP(i,:))
%         hold on
%     end
    
    if ~isempty(p)
        for j = 1:length(p)
            try
            if sign(CoP(i,p(j)) - CoP(i,p(j)-1)) * sign(CoP(i,p(j)+1) - CoP(i,p(j))) == -1
                if (i == 2 || i == 5) && outerIter > 1 && mod(j,2) == 0
                
%                 elseif (i == 2 || i == 5) && mod(j,2) == 1 && CoP(i,p(j)) - CoP(i,p(j)-1) < 0
                
                else
%                     oldPoint = CoP(i,p(j));
                    CoP(i,p(j)) = (CoP(i,p(j)+1) + CoP(i,p(j)-1)) / 2;
%                     fprintf('Modified CoP spike in FP coordinates CoP(%d,%d) from %f to %f\n', ...
%                         i, p(j), oldPoint, CoP(i,p(j)));
                end
            end
            catch
            end
        end
    end
end



% ========================================================================
% SUBFUNCTION: findPotentialSpikePoints
% ========================================================================
% Potential spike points occur at the first and last nonzero CoP of a region

function p = findPotentialSpikePoints(row, shift)

p = [];
nonZero = find(abs(row)>0);

if isempty(nonZero)
    p = [];
    return
end

if nonZero(1) ~= 1
    p = [p nonZero(1)+shift];
end

for i = 1:length(nonZero)-1
    if nonZero(i)+1 ~= nonZero(i+1)
        p = [p nonZero(i)-shift nonZero(i+1)+shift];
    end
end
if nonZero(end) ~= length(row)
    p = [p nonZero(end)-shift];
end

