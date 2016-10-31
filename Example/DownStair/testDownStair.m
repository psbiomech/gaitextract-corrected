% Example: Normal Level Walking Gait
% Tim Dorn
% -----------------------------------------------------------------------
% Note: This file can act as a template for extracting any experimental
% data from a motion capture system. loadlabels.m needs to be configured
% for the specific laboratory (see the user manual).
% -----------------------------------------------------------------------


%% Extract Data From Raw C3D Files
% ================================
clear all
close all

fprintf('\nAnswer YES when asked about the static trial.\n');


% File Name Descriptors
% ---------------------
c3dFileStatic = 'example2_static.c3d';
c3dFileDynamic = 'example2_downstair.c3d';


% Extract Event Keys
% ------------------
keySta = getEvents(c3dFileStatic, 2);
keyDyn = getEvents(c3dFileDynamic, -2, [5 4 1]);


% Extract Markers
% ---------------
markersSta = getMarkers(keySta, 'markersStaticArms');
markersDyn = getMarkers(keyDyn, 'markersDynamicArms');


% Extract Kinetics
% ----------------
opt = 2;        % Option to display all kinetic plots
filtFreq = 40;  % Kinetic low pass filter frequency (4th order butterworth)
[GRF, CoP, GRMo, GRMx] = getKinetics(keyDyn, opt, filtFreq, markersDyn);


% Extract EMG
% -----------
[eVecGlob, EMGVecGlob] = batchEMGprocess(keyDyn, 'quadsEMG', 'EMGprocessFILTERS', 'myEMGsetfile');
fprintf('\nLets plot some EMG signals...\n');
extractMotFile('file', 'example2_downstair_myEMGsetfile.mot');


% Write XML Files for Opensim
% ----------------------------------------------------------------
% Note that these XML files are only templates. They may need to be fine
% tuned in terms of the paths of the models you are using, and any
% additional settings that the analyses offer. The XML files can be
% modified in Notepad++ or for the more advanced Matlab users, you can go
% into the writeXML.m file and modify some of the default XML output
% settings.
% ----------------------------------------------------------------
writeXML('scale', keySta);
writeXML('ik', keyDyn);
writeXML('id_ik', keyDyn);

