% Generate a marker *.trc file readable by OpenSim
% 
% Tim Dorn
% June 2009
% 
% --------------------------------------------------------------------
% Usage: generateTrcFile(C3Dkey, markerpos, markerset, normTime)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
%           markerpos = array of marker positions
%                   for M markers: should contain 1+3M columns
%                   (time + XYZ of each marker)
%           markerset = cell array of strings containing the names of markers
%                   e.g. markerset = {'M1', 'M2', 'M3'};
%           normTime = flag to normalize time between 0 and 100% 
%                      1 = original time, 2 = normalized time
% 
% Outputs:  output trc file
% 
% 
% Important note: the opensim convention is to output the GRF and CoP from
% start time to end time (columns) and RIGHT foot forces/positions then
% LEFT foot forces/positions (rows).
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

function generateTrcFile(C3Dkey, markerpos, markerset, normTime)

if nargin == 3
    normTime = 1;
end

PathFileType = 4;
name = sprintf('%s.trc', C3Dkey.c3dFile);
datatype = '(X/Y/Z)';
DataRate = C3Dkey.vFreq;
CameraRate = C3Dkey.vFreq;
NumFrames = length(markerpos(:,1));
NumMarkers = length(markerset);
Units = 'mm';
OrigDataRate = C3Dkey.vFreq;
OrigDataStartFrame = C3Dkey.event.Vframes(1);
OrigNumFrames = NumFrames;

% Adjust frames & times so that the effective trial period
% starts at frame 1 at time 0
markerpos(:,1) = markerpos(:,1) - C3Dkey.event.Vframes(1) + 1;
if normTime == 2        % normalize time
    t = C3Dkey.timeVec.Vpercent';
else                    % original time
    t = C3Dkey.timeVec.Vsec'; 
end
frame = markerpos(:,1);
markerpos = [t, markerpos(:,2:end)];



% TRC File Header
% ---------------

fid = fopen(name, 'w');
if fid < 0
    fprintf('\nERROR: %s could not be opened for writing...\n\n', name);
    return
end
fprintf(fid, 'PathFileType\t%d\t%s\t%s\t\n', PathFileType, datatype, name);
fprintf(fid, 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid, '%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', ...
    DataRate, CameraRate, NumFrames, NumMarkers, Units, OrigDataRate, OrigDataStartFrame, OrigNumFrames);
fprintf(fid, 'Frame#\tTime\t');


% TRC File Body
% -------------

for i = 1:NumMarkers
    fprintf(fid, '%s\t\t\t', markerset{i});
end
fprintf(fid, '\n\t\t');

for i = 1:NumMarkers
    fprintf(fid, 'X%d\tY%d\tZ%d\t', i, i, i);
end
fprintf(fid, '\n\n');

% marker position values
for i = 1:NumFrames
    fprintf(fid, '%d\t', frame(i));
    fprintf(fid, '%.5f\t', markerpos(i,:));
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Saved (tab delimited) marker positions to: %s\n', name);

