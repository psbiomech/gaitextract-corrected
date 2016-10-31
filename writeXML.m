% Write XML files for use in OpenSim.
% Tim Dorn
% Last Updated: Sept 2010
% 
% --------------------------------------------------------------------
% Usage: [xmlString, fileName] = writeXML(type, C3Dkey*)
% --------------------------------------------------------------------
% 
% Inputs:  type: the type of XML setup file to be created (case sensitive)
%               type = 'scale' --> create scale XML file
%               type = 'ik' --> create inverse kinematics XML setup file
%               type = 'rra' --> create rra XML setup file
%               type = 'id_ik' --> create inverse dynamics XML setup file
%               type = 'id_rra' --> create inverse dynamics XML setup file (using rra kinematics)
%               type = 'so' --> create static optimization XML setup file
%               type = 'jr' --> create joint reaction XML setup file
%               type = 'ma' --> create muscle analysis XML setup file
%               type = 'pi' --> create pseudo inverse GRF decomposition XML setup file
%               type = 'vis' --> create visualization tool XML setup file
%
%           C3Dkey*: the C3D key structure from getEvents
%                   (if not given, default parameters are used. These
%                    can be then modified manually using an XML editor)
% 
%           saveFile: 0 = don't save file, 1 = save file (default)
% 
%           
% 
% Outputs:  xmlString: the generated XML string
%           fileName: the generated file name
% 
%           the XML file is saved in the current working directory as
%           [C3Dkey.c3dFile]_Setup_[type].xml
% 
% 
% Notes:    The output files are templates only. They may need to be
%           altered to conform to the precise simulation (i.e. paths, 
%           filtering frequencies, input files, etc).
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


function [xmlString, fileName] = writeXML(type, C3Dkey, saveFile)

usage = 'Usage: [xmlString, fileName] = writeXML(type, C3Dkey*)';
types = 'Valid Types: ''scale'', ''ik'', ''id'', ''so'', ''jr'', ''ma'', ''pi'', ''vis''';

if nargin == 1,
    C3Dkey = [];
    saveFile = 1;
    
elseif nargin == 2,
    saveFile = 1;
   
elseif nargin ~= 3,
    fprintf('%s\n%s\n', usage, types);
    return
end


if isempty(C3Dkey)
    % use default C3Dkey for manual editing
    C3Dkey.c3dFile = 'C3DFILE';
    C3Dkey.name = 'SUBJECTNAME';
    C3Dkey.timeVec.Asec = [0 9999];
end



% =============================
% Simulation XML Setting Files
% =============================
F.C3DFILE = C3Dkey.c3dFile;      % c3d filename without extension
offsetTime = 0.001;

% Model Files
% -----------
modelPrefix = '..\MU2392_arms';
F.unscaledModelFile = 'c:\!PROJECT\Sim_Models\PFJOA\!Template\GenericModel\PFJOA_Generic.osim';
% F.unscaledModelFile = 'c:\!PROJECT\Sim_Models\MU2392_Arms\MU2392_genericArms_noWrap_Mod2.osim';
F.scaledModelFile = sprintf('..\\%s_SCALED.osim', C3Dkey.name);
F.rraModelFile = sprintf('%s_SCALED_RRA1.osim', F.C3DFILE);
F.modModelFile = sprintf('%s_SCALED_MOD1.osim', F.C3DFILE);


% Result Directories
% ------------------
F.ID_Results = '.\ID_Results';
F.RRA_Results = '.\RRA_Results';
F.SO_Results = '.\SO_Results';
F.JR_Results = '.\JR_Results';
F.MA_Results = '.\MA_Results';
F.PI_Results = '.\PI_Results';



% Experimental / Model Output Data Files
% --------------------------------------
F.manualScaleSetFile = sprintf('%s_Scale_ScaleSet.xml', F.C3DFILE);
F.markerFile = sprintf('%s.trc', F.C3DFILE);
F.coordinateFile = sprintf('%s_coordinates.mot', F.C3DFILE);
F.kineticsFile = sprintf('%s_kinetics.mot', F.C3DFILE);
F.IKFile = sprintf('%s_ik_arms.mot', F.C3DFILE);
F.IDFile = sprintf('%s\\%s_INVDYN_RRA_arms_InverseDynamics.sto', F.ID_Results, F.C3DFILE);
F.GRFFile = sprintf('%s_kinetics.xml', F.C3DFILE);
F.RRAcoordinateFile = sprintf('%s\\%s_RRA2_arms_Kinematics_q.mot', F.RRA_Results, F.C3DFILE);
F.SOforces = sprintf('%s\\%s_SO_arms_StaticOptimization_force.sto', F.SO_Results, F.C3DFILE);
F.markerSetFile = sprintf('%s_Scale_MarkerSet.xml', modelPrefix);
F.measurementSetFile = sprintf('%s_Scale_MeasurementSet.xml', modelPrefix);
F.taskSetFile.SCALE = sprintf('%s_Scale_Tasks.xml', modelPrefix);
F.taskSetFile.IK = sprintf('%s_IK_Tasks.xml', modelPrefix);


% Other Settings
% --------------
S.kinematics_filter_freq = 7;        % Hz
S.stepInt = 5;
S.output_precision = 20;
S.optIK.optimizer_algorithm = 'cfsqp';



% ====================================================================
% Based on the type of XML file to create, the data structure is built
% and saved to an xml file
% ====================================================================
timeRange = [C3Dkey.timeVec.Vsec(1), C3Dkey.timeVec.Vsec(end)];

switch type
    
    % SCALE XML FILE (SCALE)
    % ----------------------
    case 'scale',
        root = 'ScaleTool';
        
        % ScaleTool
        data.ATTRIBUTE.name = sprintf('%s_SCALED', C3Dkey.name);
        if isfield(C3Dkey, 'mass')
            data.mass = C3Dkey.mass;
        else
            data.mass = '9999999';
        end
        data.height = '9999999';
        data.age = '9999999';
        data.notes = 'enter notes';
        
        % ScaleTool -> GenericModelMaker
        data.GenericModelMaker.ATTRIBUTE.name = '';
        data.GenericModelMaker.model_file = F.unscaledModelFile;
        data.GenericModelMaker.marker_set_file = F.markerSetFile;
        
        % ScaleTool -> ModelScaler
        data.ModelScaler.ATTRIBUTE.name = '';
        data.ModelScaler.apply = 'true';
        data.ModelScaler.scaling_order = 'measurements manualScale';
        
%         the following line is used for manualScaling (if necessary)
%         data.ModelScaler.ScaleSet.ATTRIBUTE.file = F.manualScaleSetFile;

        data.ModelScaler.MeasurementSet.ATTRIBUTE.file = F.measurementSetFile;
        data.ModelScaler.marker_file = F.markerFile;
        data.ModelScaler.time_range = num2str(timeRange);
        data.ModelScaler.preserve_mass_distribution = 'true';
%         data.ModelScaler.output_model_file = F.scaledModelFile;         % output model without markers
        data.ModelScaler.output_scale_file = sprintf('%s_scaleSet_applied.xml', C3Dkey.name);
        
        % ScaleTool -> MarkerPlacer
        data.MarkerPlacer.ATTRIBUTE.name = '';
        data.MarkerPlacer.apply = 'true';
        data.MarkerPlacer.optimizer_algorithm = S.optIK.optimizer_algorithm;
        data.MarkerPlacer.IKTaskSet.ATTRIBUTE.file = F.taskSetFile.SCALE;
        data.MarkerPlacer.marker_file = F.markerFile;
        data.MarkerPlacer.coordinate_file = F.coordinateFile;
        data.MarkerPlacer.time_range = num2str(timeRange);
        data.MarkerPlacer.output_model_file = F.scaledModelFile;          % output model with markers
        data.MarkerPlacer.output_motion_file = sprintf('%s_static_output.mot', C3Dkey.name);
        
        
        
        
    % INVERSE KINEMATICS XML FILE (IK)
    % --------------------------------
    case 'ik',
        
        % IKTool
        root = 'IKTool';
        data.ATTRIBUTE.name = F.C3DFILE;
        data.model_file = F.scaledModelFile;
        data.optimizer_algorithm = S.optIK.optimizer_algorithm;
        data.IKTaskSet.ATTRIBUTE.file = F.taskSetFile.IK;
        
        % IKTool -> IKTrialSet
        data.IKTrialSet.ATTRIBUTE.name = '';
        data.IKTrialSet.objects.IKTrial.ATTRIBUTE.name = F.C3DFILE;
        data.IKTrialSet.objects.IKTrial.marker_file = F.markerFile;
        data.IKTrialSet.objects.IKTrial.coordinate_file = F.coordinateFile;
        data.IKTrialSet.objects.IKTrial.time_range = num2str(timeRange);
        data.IKTrialSet.objects.IKTrial.output_motion_file = F.IKFile;
%         data.IKTrialSet.objects.IKTrial.include_markers = 'true';
    
     

    % RRA XML FILE (RRA)
    % -------------------
    case 'rra',
        
        % CMCTool
        root = 'CMCTool';
        data.ATTRIBUTE.name = sprintf('%s_RRA2_arms', F.C3DFILE);
        data.model_file = F.scaledModelFile;
        data.replace_force_set = 'true';
        data.force_set_files = sprintf('%s_RRA_Forces.xml', F.C3DFILE);
        data.task_set_file = sprintf('%s_RRA_Tasks.xml', F.C3DFILE);
        data.constraints_file = sprintf('%s_RRA_ControlConstraints.xml', F.C3DFILE);
        data.results_directory = F.RRA_Results;
        data.output_precision = S.output_precision;
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimesRRA(data, timeRange);
        data.maximum_number_of_integrator_steps = 20000;
        data.maximum_integrator_step_size = 1e-4;
        data.integrator_error_tolerance = 1e-6;
        data.adjust_com_to_reduce_residuals = 'true';
        data.adjusted_com_body = 'torso';
        data.output_model_file = F.rraModelFile;
        data.desired_kinematics_file = F.IKFile;
        data.lowpass_cutoff_frequency = S.kinematics_filter_freq;
        data = addExternalLoadSettingsRRA(data, S, F);
        data.use_fast_optimization_target = 'false';
        data.optimizer_derivative_dx = 1e-4;
        data.optimizer_convergence_criterion = 1e-8;
        data.optimizer_max_iterations = 2000;
        data.optimizer_print_level = 0;
        data.optimizer_algorithm = S.optIK.optimizer_algorithm;
        data.cmc_time_window = 0.001;
        data.use_curvature_filter = 'false';
        data.compute_average_residuals = 'true';
        
        % AnalyzeTool -> AnalysisSet
        data.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addKinematicsAnalysis(data);
        data = addActuationAnalysis(data);
        data = addBodyKinematicsAnalysis(data);
        
       

    % INVERSE DYNAMICS XML FILE (ID)
    % using inverse kinematics solution
    % ------------------------------
    case 'id_ik',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_INVDYN_IK_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.scaledModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = ' ';
        data.AnalyzeTool.results_directory = F.ID_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.IKFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addInverseDynamicsAnalysis(data);
        
        
        
    % INVERSE DYNAMICS XML FILE (ID_RRA)
    % using rra solution
    % ------------------------------
    case 'id_rra',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_INVDYN_RRA_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.modModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = ' ';
        data.AnalyzeTool.results_directory = F.ID_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.RRAcoordinateFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addInverseDynamicsAnalysis(data);
        data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.step_interval = 5;
  
        
        
    % STATIC OPTIMIZATION XML FILE (SO)
    % ---------------------------------
    case 'so',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_SO_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.scaledModelFile;    % F.modModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = 'Reserve_Actuators.xml';
%         data.AnalyzeTool.force_set_files = '..\..\..\Reserve_Actuators.xml';
%         data.AnalyzeTool.force_set_files = sprintf('%s_Reserve_Actuators.xml', F.C3DFILE);
        data.AnalyzeTool.results_directory = F.SO_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        data.AnalyzeTool.solve_for_equilibrium_for_auxiliary_states = 'false';
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.IKFile;      % F.RRAcoordinateFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addStaticOptimizationAnalysis(data, 1);
        % data = addStaticOptimizationAnalysis(data, S.stepInt);
        
        
    % JOINT REACTION XML FILE (JR)
    % -------------------------------
    case 'jr',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_JR_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.modModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = sprintf('%s_Reserve_Actuators.xml', F.C3DFILE);
        data.AnalyzeTool.results_directory = F.JR_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        data.AnalyzeTool.solve_for_equilibrium_for_auxiliary_states = 'false';
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.RRAcoordinateFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addJointReactionAnalysis(data, S.stepInt);
        data.AnalyzeTool.AnalysisSet.objects.JointReaction.forces_file = ...
            sprintf('%s\\%s_SO_arms_StaticOptimization_force.sto', F.SO_Results, F.C3DFILE);
                
        
        
    % MUSCLE ANALYSIS XML FILE (MA)
    % -----------------------------
    case 'ma',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.modModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = ' ';
        data.AnalyzeTool.results_directory = F.MA_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        data.AnalyzeTool.solve_for_equilibrium_for_auxiliary_states = 'true';
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.RRAcoordinateFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        % data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses'; 
        data = addMuscleAnalysis(data, S.stepInt);
       
     
        
    % PSEUDO INVERSE GRF DECOMPOSITION XML FILE (PI)
    % ----------------------------------------------
    case 'pi',
        
        % AnalyzeTool
        root = 'OpenSimDocument';
        data.ATTRIBUTE.Version = '20001';
        data.AnalyzeTool.ATTRIBUTE.name = sprintf('%s_PI_arms', F.C3DFILE);
        data.AnalyzeTool.model_file = F.modModelFile;
        data.AnalyzeTool.replace_force_set = 'false';
        data.AnalyzeTool.force_set_files = sprintf('%s_Reserve_Actuators.xml', F.C3DFILE);
        data.AnalyzeTool.results_directory = F.PI_Results;
        data.AnalyzeTool.output_precision = S.output_precision;
        data.AnalyzeTool.solve_for_equilibrium_for_auxiliary_states = 'false';
        timeRange(1) = timeRange(1) + offsetTime;
        timeRange(2) = timeRange(2) - offsetTime;
        data = addInitialFinalTimes(data, timeRange);
        data.AnalyzeTool.coordinates_file = F.RRAcoordinateFile;
        data.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addPseudoInverseAnalysis(data, 1, F);
        

         
    % VISUALIZATION TOOL XML FILE (VIS)
    % ---------------------------------
    case 'vis',
        
        % AnalyzeTool
        root = 'VisualSetup';
        
        data(1).VisualTask.ATTRIBUTE.name = 'IK_MARKERS';
        data(1).VisualTask.AssignmentSet = ' ';
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.alpha.ChannelMap.RegisteredMaps.inputValues = [0 1];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.alpha.ChannelMap.RegisteredMaps.channelValues = [0.5 1];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.red.ChannelMap.RegisteredMaps.inputValues = [0.005 0.032];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.red.ChannelMap.RegisteredMaps.channelValues = [0 1];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.green.ChannelMap.RegisteredMaps.inputValues = [0.005 0.032];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.green.ChannelMap.RegisteredMaps.channelValues = [0 1];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.blue.ChannelMap.RegisteredMaps.inputValues = [0.005 0.032];
        data(1).VisualTask.PerformersDB.Streamlines.Colormap.blue.ChannelMap.RegisteredMaps.channelValues = [0 1];
        
        data(2).VisualTask.ATTRIBUTE.name = 'ID_MOMENT';
        data(2).VisualTask.AssignmentSet = ' ';
        data(2).VisualTask.PerformersDB.ForceMomentViewer.data_file = F.IDFile;
        data(2).VisualTask.PerformersDB.ForceMomentViewer.scale_factor = 0.005;
       
   
        
    % UNKNOWN XML TYPE
    % ----------------
    otherwise
        fprintf('Invalid XML type. The type is case sensitve\n%s', types);
end



% Save the data structure to an xml file
% --------------------------------------
if saveFile > 0
    fileName = sprintf('%s_Setup_%s.xml', F.C3DFILE, upper(type));
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
end








% ========================================================================
% SUBFUNCTION: data = addExternalLoadSettings(data, S, F)
% ========================================================================
% Loads the data structure with "experimental" external load information
% Used for ID, RRA1, RRA2, CMC and FORWARD
% Note: If kinetics are modified in RRA for example, run this function
%       and then overwrite data.external_loads_file with the new source

function data = addExternalLoadSettings(data, S, F)

data.AnalyzeTool.external_loads_file = F.GRFFile;   % xml file, so we don't need external_loads_body1, external_loads_body2 fields 
data.AnalyzeTool.external_loads_model_kinematics_file = F.IKFile;
% data.AnalyzeTool.external_loads_body1 = 'calcn_r';      % Right foot first
% data.AnalyzeTool.external_loads_body2 = 'calcn_l';      % Left foot second
data.AnalyzeTool.lowpass_cutoff_frequency_for_load_kinematics = S.kinematics_filter_freq;



% ========================================================================
% SUBFUNCTION: data = addExternalLoadSettingsRRA(data, S, F)
% ========================================================================
% Loads the data structure with "experimental" external load information
% Used for ID, RRA1, RRA2, CMC and FORWARD
% Note: If kinetics are modified in RRA for example, run this function
%       and then overwrite data.external_loads_file with the new source

function data = addExternalLoadSettingsRRA(data, S, F)

data.external_loads_file = F.GRFFile;               % xml file, so we don't need external_loads_body1, external_loads_body2 fields 
data.external_loads_model_kinematics_file = F.IKFile;
% data.external_loads_body1 = 'calcn_r';      % Right foot first
% data.external_loads_body2 = 'calcn_l';      % Left foot second
data.lowpass_cutoff_frequency_for_load_kinematics = S.kinematics_filter_freq;



% ========================================================================
% SUBFUNCTION: data = addInitialFinalTimes(data, timeRange)
% ========================================================================
% Loads the data structure with initial and final time information
% Used for ID, RRA1, RRA2, CMC and FORWARD

function data = addInitialFinalTimes(data, timeRange)

data.AnalyzeTool.initial_time = timeRange(1);
data.AnalyzeTool.final_time = timeRange(end);



% ========================================================================
% SUBFUNCTION: data = addInitialFinalTimesRRA(data, timeRange)
% ========================================================================
% Loads the data structure with initial and final time information
% Used for ID, RRA1, RRA2, CMC and FORWARD

function data = addInitialFinalTimesRRA(data, timeRange)

data.initial_time = timeRange(1);
data.final_time = timeRange(end);



% ========================================================================
% SUBFUNCTION: data = addInverseDynamicsAnalysis(data)
% ========================================================================
% Adds an INVERSE DYNAMICS analysis

function data = addInverseDynamicsAnalysis(data)

data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.ATTRIBUTE.name = 'InverseDynamics';
data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.on = 'true';
data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.in_degrees = 'true';
data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.step_interval = 1;
data.AnalyzeTool.AnalysisSet.objects.InverseDynamics.use_model_force_set = 'false';         % > OS2.0



% ========================================================================
% SUBFUNCTION: data = addStaticOptimizationAnalysis(data)
% ========================================================================
% Adds a STATIC OPTIMIZATION analysis

function data = addStaticOptimizationAnalysis(data, step_interval)

data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.ATTRIBUTE.name = 'StaticOptimization';
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.on = 'true';
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.in_degrees = 'true';
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.step_interval = step_interval;
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.use_model_force_set = 'true';         % > OS2.0
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.activation_exponent = 2.0;
data.AnalyzeTool.AnalysisSet.objects.StaticOptimization.use_muscle_physiology = 'true';



% ========================================================================
% SUBFUNCTION: data = addJointReactionAnalysis(data)
% ========================================================================
% Adds a JOINT REACTION analysis

function data = addJointReactionAnalysis(data, step_interval)

data.AnalyzeTool.AnalysisSet.objects.JointReaction.ATTRIBUTE.name = 'JointReaction';
data.AnalyzeTool.AnalysisSet.objects.JointReaction.on = 'true';
data.AnalyzeTool.AnalysisSet.objects.JointReaction.in_degrees = 'true';
data.AnalyzeTool.AnalysisSet.objects.JointReaction.step_interval = step_interval;
data.AnalyzeTool.AnalysisSet.objects.JointReaction.joint_names = 'all';
data.AnalyzeTool.AnalysisSet.objects.JointReaction.apply_on_bodies = 'child';
data.AnalyzeTool.AnalysisSet.objects.JointReaction.express_in_frame = 'child';



% ========================================================================
% SUBFUNCTION: data = addMomentArmsAnalysis(data)
% ========================================================================
% Adds a MOMENT ARM analysis

function data = addMuscleAnalysis(data, step_interval)

data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.ATTRIBUTE.name = 'MuscleAnalysis';
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.on = 'true';
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.in_degrees = 'true';
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.step_interval = step_interval;
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.muscle_list = 'all';
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.moment_arm_coordinate_list = ...
    'hip_flexion_r hip_adduction_r hip_rotation_r knee_angle_r ankle_angle_r subtalar_angle_r hip_flexion_l hip_adduction_l hip_rotation_l knee_angle_l ankle_angle_l subtalar_angle_l';
data.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.compute_moments = 'false';



% ========================================================================
% SUBFUNCTION: data = addPseudoInverseAnalysis(data)
% ========================================================================
% Adds a PSEUDOINVERSE FORCE DECOMPOSITION analysis

function data = addPseudoInverseAnalysis(data, step_interval, F)

data.AnalyzeTool.AnalysisSet.objects.IndAccPI.ATTRIBUTE.name = 'IndAccPI';
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.on = 'true';
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.step_interval = step_interval;
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.in_degrees = 'true';
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.kinetics_file = F.kineticsFile;
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.forces_file = F.SOforces;
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.force_threshold = 15;
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.weights = [1000 100 1];
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.footpoint_markers = ...
        'r_fp_heel1 r_fp_heel2 r_fp_mt1 r_fp_mt2 r_fp_toe l_fp_heel1 l_fp_heel2 l_fp_mt1 l_fp_mt2 l_fp_toe';
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.compute_potentials_only = 'false';
% data.AnalyzeTool.AnalysisSet.objects.IndAccPI.contact_point_type = '5points';
% data.AnalyzeTool.AnalysisSet.objects.IndAccPI.cop_weights = [1 1 1 1 1 1];
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.output_grf_contribution = 2;
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.acceleration_match = 'false';
data.AnalyzeTool.AnalysisSet.objects.IndAccPI.debug_output = 0;



% ========================================================================
% SUBFUNCTION: data = addKinematicsAnalysis(data)
% ========================================================================
% Adds a KINEMATICS analysis

function data = addKinematicsAnalysis(data)

data.AnalysisSet.objects.Kinematics.ATTRIBUTE.name = 'Kinematics';
data.AnalysisSet.objects.Kinematics.on = 'true';
data.AnalysisSet.objects.Kinematics.step_interval = 10;
data.AnalysisSet.objects.Kinematics.in_degrees = 'true';
data.AnalysisSet.objects.Kinematics.coordinates = 'all';



% ========================================================================
% SUBFUNCTION: data = addActuationAnalysis(data)
% ========================================================================
% Adds an ACTUATION analysis

function data = addActuationAnalysis(data)

data.AnalysisSet.objects.Actuation.ATTRIBUTE.name = 'Actuation';
data.AnalysisSet.objects.Actuation.on = 'true';
data.AnalysisSet.objects.Actuation.step_interval = 10;
data.AnalysisSet.objects.Actuation.in_degrees = 'true';



% ========================================================================
% SUBFUNCTION: data = addBodyKinematicsAnalysis(data)
% ========================================================================
% Adds a BODY KINEMATICS analysis

function data = addBodyKinematicsAnalysis(data)

data.AnalysisSet.objects.BodyKinematics.ATTRIBUTE.name = 'BodyKinematics';
data.AnalysisSet.objects.BodyKinematics.on = 'true';
data.AnalysisSet.objects.BodyKinematics.step_interval = 10;
data.AnalysisSet.objects.BodyKinematics.in_degrees = 'true';
data.AnalysisSet.objects.BodyKinematics.express_results_in_body_local_frame = 'false';

