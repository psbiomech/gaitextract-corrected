% Batch Process EMG signals
% Tim Dorn
% November 2008
% 
% --------------------------------------------------------------------
% Usage: [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix*)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey = key of dynamic C3D file (from getEvents.m)
% 
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as a cell: glab.[emgSetName])
% 
%           emgProcessTasks: the EMG processing options in the order of execution
%               (this must be defined in loadlabels.m as a cell: glab.[emgProcessTasks])
% 
%           fileSuffix* = suffix of mot file that is saved    
%                   if this is not included, or empty, then file is not saved
% 
% 
% Outputs:  eVecGlob = structure of processed EMG
%           EMGVecGlob = structure of raw EMG
% 
% 
% Notes:    N/A
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

function [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix)

usage = 'Usage: [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix*)';
    
if nargin == 3,
    fileSuffix = [];

elseif nargin ~= 4,
    disp(usage)
    return
end

loadLabels;
EMGset = glab.(emgSetName);
for i = 1:length(EMGset)
    musc_name_orig = EMGset{i}{1};
    musc_name = removeSpaces(strtok(musc_name_orig, '.()'));
end

if ischar(emgProcessTasks)
    emgProcessTasks = glab.(emgProcessTasks);
end

ind = getIndex(emgProcessTasks, 'plot');
if ind > 0
    if strcmpi(emgProcessTasks{ind+1}, 'C3Dkey')
        emgProcessTasks{ind+1} = C3Dkey;
    end
    ind = getIndex(emgProcessTasks, 'vertlines');
    if strcmpi(emgProcessTasks{ind+1}, 'C3Dkey')
        emgProcessTasks{ind+1} = C3Dkey;
    end
end


% ---------------
% Fixed Variables
% ---------------

% Get Dynamic EMG Values
bigM = [];
p = 0;
colnames = {'time'};
EMGVecGlob = extractRawEMG(C3Dkey, emgSetName);
names = fieldnames(EMGVecGlob);
numMuscles = length(names);
file = sprintf('%s_%s.mot', C3Dkey.c3dFile, fileSuffix);



% Process all EMG channels
% ------------------------
for i = 1:numMuscles
    l = length(glab.(emgSetName){i});
    musc_name = names{i};
    
    if l < 3
        nM = 1;
        musc_name_in_model = {musc_name};
    else
        nM = l-2;
        musc_name_in_model = glab.(emgSetName){i}(3:3+nM-1);
    end
    
    EMGVecTmp = EMGVecGlob.(musc_name);
    
    for j = 1:nM
        
        % Process EMGVec(t) -> Processed EMG e(t)
        % ---------------------------------------
        if ~isempty(EMGVecTmp.data)    % only process data that we were able to extract
            eVecTmp = processEMG(EMGVecTmp, emgProcessTasks);
        
            % Set variables
            eVecGlob.(musc_name) = eVecTmp;
            if isfield(eVecTmp, 'SplitStrides')     % multiple cycles in a single C3Dfile
                                                    % uses AVERAGE values here
                colnames = [colnames, {sprintf('%s__Raw', musc_name_in_model{j})}, ...
                    {sprintf('%s__Processed', musc_name_in_model{j})}];
                bigM = [bigM, eVecGlob.(musc_name).SplitStrides.rawAvr, eVecGlob.(musc_name).SplitStrides.processedAvr];
                timeAvr = eVecGlob.(musc_name).SplitStrides.timeAvr1;
                p = 1;

            else        % normal use: single cycle in a single C3Dfile
                colnames = [colnames, {sprintf('%s__Raw', musc_name_in_model{j})}, ...
                    {sprintf('%s__Processed', musc_name_in_model{j})}];
                bigM = [bigM, EMGVecGlob.(musc_name).data, eVecGlob.(musc_name).data];
            end
        end
    end
end


if findstr(emgSetName, 'right')
    testLeg = 'r';
    fprintf('Test Leg: RIGHT\n')
elseif findstr(emgSetName, 'left')
    testLeg = 'l';
    fprintf('Test Leg: LEFT\n')
end
if exist('testLeg', 'var')
    fid = fopen('testLeg.data', 'wb');
    fwrite(fid, testLeg);
    fclose(fid);
end


if ~isempty(fileSuffix)
    if p == 0
        bigM = [C3Dkey.timeVec.Asec', bigM];       % add time column
        fprintf('Saving EMG values\n')
    else
        bigM = [timeAvr, bigM];       % add time column
        fprintf('Saving AVERAGE EMG values\n')
    end
    generateMotFile(bigM, colnames, file);
end

