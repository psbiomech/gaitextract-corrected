% Batch Process multiple stride EMG signals from a single C3D file
% Tim Dorn
% July 2009
% 
% --------------------------------------------------------------------
% Usage: [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)
% --------------------------------------------------------------------
% 
% Inputs:   c3dFile = the name of the c3d file
% 
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as a cell: glab.[emgSetName])
% 
%           emgProcessTasks: the EMG processing options in the order of execution
%               (this must be defined in loadlabels.m as a cell: glab.[emgProcessTasks])
% 
%           interval4Time = the interval number to set for the time vector
% 
% 
% Outputs:  eVecGlob = structure of processed EMG
%           EMGVecGlob = structure of raw EMG
% 
% 
% Notes:    Ensure that events are places in the c3d file at the start & end
%           of each interval. e.g. 4 events == 3 intervals. If we want to set
%           the time vector to be the middle interval (between events 2&3), then
%           interval4Time = 2.
% 
%           e.g. multipleEMGprocess('myfile.c3d', 'emgset', 2)
% 
%           The output file will contain the time column (from interval4Time)
%           and the time normalized AVERAGE emg data over all intervals
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

function [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)

usage = 'Usage: [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)';

if nargin ~= 4
    disp(usage)
    return
end

C3DkeyEMG = getEvents(c3dFile);
C3DkeyEMG.interval4Time = interval4Time;
loadLabels;
if ischar(emgProcessTasks)
    emgProcessTasks = glab.(emgProcessTasks);
end


% look for 'plot' and change value to C3Dkey to change the plotting option
l = length(emgProcessTasks);
ind = getIndex(emgProcessTasks, 'plot');
if ind == 0
    emgProcessTasks{l+1} = 'plot';
    emgProcessTasks{l+2} = 'C3Dkey';
else
    emgProcessTasks{ind+1} = 'C3Dkey';
end
    

% Save EMG
[eVecGlob, EMGVecGlob] = batchEMGprocess(C3DkeyEMG, emgSetName, emgProcessTasks, 'AverageEMG');


