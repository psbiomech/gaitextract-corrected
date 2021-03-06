% Get kinematic marker positions
% Tim Dorn
% October 2010
% 
% --------------------------------------------------------------------
% Usage: markers = getMarkers(C3Dkey, markerSetName, saveTRC*, filter*)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
% 
%           markerSetName: the label of markers contained in the marker set
%               (this must be defined in loadlabels.m as glab.[markerSetName])
% 
%           saveTRC*: toggle the trc saving on (1) or off(0). On is default
% 
%           filter*: optional low pass filter the marker positions
%                 if filter = 0, no filtering is done (default)
%                 if filter > 0, filtering is done at (filter) Hz
% 
% 
% Outputs:  markers.data:  the marker position data
%           markers.label: the marker position labels
%           markers.units: units of marker positions
%           markers.divide_to_meters: scale to divide to convert to meters
%           markers.SUCCESS: marker names that have been extracted successfully
%           markers.FAILED: marker names that have failed (do not exist)
%           markers.MISSINGMARKERS: marker names that have markers missing
% 
% 
% Notes: 
% 
% It is vital that the label names stored in the 'glab.markers' variable
% in loadLabels.m match the label names used in Vicon for extraction to be
% successful. Sometimes Vicon stores the marker labels as 
% [markerName:SubjectName]. This will result in an error since the strings
% do not match. To resolve this, ensure that �Include subject names in labels�
% check box is switched OFF in Vicon ? Trials ? Options.
% The data is extracted in the order defined in glab.markers.
% 
% Output is saved as a *.trc file used by OpenSim
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


function markers = getMarkers(C3Dkey, markerSetName, saveTRC, filter)

usage = 'Usage: markers = getMarkers(C3Dkey, markerSetName, saveTRC*, filter*)';

if nargin == 2,
    saveTRC = 1;
    filter = 0;
elseif nargin == 3,
    filter = 0;
elseif nargin ~= 4,
    disp(usage)
    return
end

if C3Dkey.allowed.markers == 0
   error('C3Dkey was generated with no direction input. Can not extract markers. Regenerate the C3Dkey with the direction variable and try again...\n'); 
end



% Load c3d files
% --------------
itf = c3dserver();
openc3d(itf, 0, C3Dkey.c3dFileFULL);

framesTotal = nframes(itf);
C_index = itf.GetParameterIndex('TRIAL', 'CAMERA_RATE');
if C_index<0, C_index = itf.GetParameterIndex('POINT', 'RATE'); end      % non-Vicon-generated files may not have the TRIAL group
MOTIONfreq = double(itf.GetParameterValue(C_index, 0));

if isstruct(C3Dkey),
    startF = C3Dkey.event.Vframes(1);
    endF = C3Dkey.event.Vframes(end);
    markers.a2vRatio = C3Dkey.r;
    markers.c3dFile = C3Dkey.c3dFile;
    transformVICMODEL = C3Dkey.transform.VICMODEL;
    
else
    startF = itf.GetVideoFrame(0);
    endF = itf.GetVideoFrame(1);
    transformVICMODEL = [-2, 3, -1];             % manual override here
end

markers.transformVICMODEL = transformVICMODEL;

fprintf('Processing C3D File: %s  (%d frames --> #%d - #%d)\n', ...
    C3Dkey.c3dFile, framesTotal, startF, endF);
fprintf('Motion Capture Frequency: %d Hz\nUsing marker set: %s\n', ...
    MOTIONfreq, markerSetName);



% Extract marker position data from C3D files (VICON coordinates) and 
% transform into the MODEL coordinate system using C3Dkey.transform.VICON
% -----------------------------------------------------------------------
loadLabels;
markerset = glab.(markerSetName);
markersetActualExtracted = glab.(markerSetName);
 
numMarkers = length(markerset);
var = [];
s = 0;
success = [];
frames = (startF:endF)';
markers.SUCCESS = [];
markers.FAILED = [];
markers.MISSINGTIMEFRAMES = [];

for j = 1:numMarkers
    try
        fprintf('[%s] --> ', C3Dkey.c3dFile);
        tmp = get3dtarget(itf, markerset{j}, 0, startF, endF) / C3Dkey.divide_to_mm;
        tmp = coordChange(tmp', transformVICMODEL)';      
        var = [var, tmp];
        fprintf('... DONE\n');
        s = s+1;
        success = [success j];
        markers.SUCCESS = [markers.SUCCESS, markerset(j)];
    catch
        fprintf('... ERROR\nMarker [%s] not found in markerset [%s]\n', ...
            markerset{j}, markerSetName);
        markers.FAILED = [markers.FAILED, markerset(j)];
        
        % remove the markerset entry from markersetActualExtracted
        markersetActualExtracted = removeEntryFromCell(markersetActualExtracted, markerset{j});
    end
end
closec3d(itf);



% Organise data into an output form
% ---------------------------------
d = {'-X', '-Y', '-Z'};
markers.data = double([frames, var]);
markers.label{1} = 'FRAME';
for j = 1:s
    for k = 1:3
        markers.label{3*(j-1) + k+1} = [markerset{success(j)}, d{k}];
    end
end

[xnan, ynan] = findNaN(markers.data);
for i = 1:3:length(xnan),
    mkr = strtok(markers.label{xnan(i)}, '-');
    fprintf('NaN found on marker: %s at frame %d --> output(%d, %d)\n', ...
        mkr, frames(ynan(i)), xnan(i), ynan(i));
    tmp1.marker = mkr;
    tmp1.time = [xnan(i), ynan(i)];
    tmp1.frame = frames(ynan(i));
    markers.MISSINGTIMEFRAMES = [markers.MISSINGTIMEFRAMES tmp1];
end

if isempty(xnan),
    fprintf('\n***** No missing markers found at any time :) *****\n');
else
    fprintf('\n%d marker-frames missing from trial\nWARNING: OpenSim may not like this very much...\n', length(xnan)/3);
end

f = numMarkers - s;
fprintf('%d MARKERS IN MARKERSET\n%d SUCCESSFUL MARKER EXTRACTIONS\n%d FAILED MARKER EXTRACTIONS\n', ...
    numMarkers, s, f);

markers.units = C3Dkey.units;
markers.divide_to_meters = C3Dkey.divide_to_meters;



% Optional filtering of the marker trajectories
% ---------------------------------------------
if filter > 0,
    [m,n] = size(var);
    var = var';  
    fprintf('Filtering kinematic marker positions at %d Hz...\n', filter)
    var = smooth1(var', filter*ones(1,m), MOTIONfreq)';
    var = var';
    markers.data = [frames, var];
end



% Apply offset to adjust the position on the ground
% -------------------------------------------------
for i = 1:s
    for j = 1:3
        % include offset here from getEvents.m
        markers.data_offset(:, (i-1)*3 + j + 1) = ...
            markers.data(:, (i-1)*3 + j + 1) + C3Dkey.offset(j);
    end
end
markers.data_offset(:,1) = frames;
fprintf('X marker offset: %d mm\nY marker offset: %d mm\nZ marker offset: %d mm\n', ...
    C3Dkey.offset(1), C3Dkey.offset(2), C3Dkey.offset(3));



% Save to *.trc file for use with OpenSim
% ---------------------------------------
if isstruct(C3Dkey) && saveTRC > 0,
    generateTrcFile(C3Dkey, markers.data_offset, markersetActualExtracted, saveTRC);
end
fclose all;


%pause(3)


% Plot markers against OpenSim tutorial (optional)
% ------------------------------------------------
% figure('Name', sprintf('%s platform placement', C3Dkey.c3dFile, 'NumberTitle', 'off'))
% plotmarkers(markers.data_offset(1, 2:end));   % blue with red outline
% load sampleMarkersOnPlatform
% hold on
% plotmarkers(sampleMarkersOnPlatform, [1,0,0], [0,0,1]);    % red with blue outline


