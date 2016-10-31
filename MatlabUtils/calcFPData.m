% ========================================================================
% FUNCTION: calcFPData
% Tim Dorn
% Last Updated: September 2010
% ========================================================================
% Extract force plate data from C3D file
% (all FP data is represented in the Vicon coordinate system unless
% otherwise stated)
% -------------------------------------------------------------------------
% 
% from FP origin (transducer) to center of plate (surface) in FP coords:
%               FP.localOriginsFP(PlateNumber, Coord)
% 
% from FP origin (transducer) to center of plate (surface) in VICON coords:
%               FP.localOriginsGLOB(PlateNumber, Coord)
% 
% from VICON origin to forceplate corners in VICON coords: 
%               FP.corners(CornerNumber, Coord, PlateNumber)
% 
% from VICON origin to forceplate center (surface) in VICON coords:
%               FP.viconOrig2FPcenterSurfaceGLOB(PlateNumber, Coord)
% 
% from VICON origin to forceplate origin (transducer) in VICON coords:
%               FP.viconOrig2FPorigGLOB(PlateNumber, Coord)
% 
%     where X = 1, Y = 2, Z = 3 in the relevant coordinate system
% 
%     where the plate number is the ORDER of stepping, NOT the plate
%     identity number. i.e the first plate to be stepped on is plate 1, the
%     2nd plate to be stepped on is plate 2, etc. See C3Dkey.FP_order for plate
%     ordering details...
% 
%     The corners need to be labelled sequentially (either CW or ACW)
% 
%     opt == 0 (default): calculate everything (used in getKinetics.m)
%     opt == 1          : omit FP.localOriginsGLOB and FP.viconOrig2FPorigGLOB (used in getEvents.m)
% 
% -------------------------------------------------------------------------

function FP = calcFPData(itf, C3Dkey, opt)

if nargin == 2
    opt = 0;
elseif nargin ~= 3
    error ('Usage: FP = calcFPData(itf, C3Dkey, opt*)\n\n');
end

FP.order = abs(C3Dkey.FP_order);
   
x = 1; y = 2; z = 3;
Corner_index = itf.GetParameterIndex('FORCE_PLATFORM', 'CORNERS');
Origin_index = itf.GetParameterIndex('FORCE_PLATFORM', 'ORIGIN');
Type_index = itf.GetParameterIndex('FORCE_PLATFORM', 'TYPE');

for i = 1:length(FP.order)      % Force Plate Loop (i)
    
    plateNum = FP.order(i);
    FP.type = double(itf.GetParameterValue(Type_index, plateNum-1));
    
    % Currently there is only support for type 2 force plates (6 channel outputs)
    if FP.type ~= 2
        error(sprintf('Force Plate %d == TYPE %d\nTOOLBOX CURRENTLY ONLY SUPPORTS TYPE 2 FORCE PLATES\n (3 Force & 3 Moment analog channel outputs). i.e. AMTI, BERTEC\n', plateNum, FP.type))
    end
    
    for j = 1:4                 % Corner Loop for extracting corners
        for k = 1:3             % Direction Loop for extracting corners
            FP.corners(j,k,i) = double(itf.GetParameterValue(Corner_index, ...
                12*(plateNum-1)+3*(j-1)+k-1)) / C3Dkey.divide_to_meters;
        end
    end
    
    for k = 1:3                 % Direction loop to extract local origins
        FP.localOriginsFP(i,k) = double(itf.GetParameterValue(Origin_index, ...
            3*(plateNum-1)+k-1)) / C3Dkey.divide_to_meters;
    end
    
    % Z value for Force Plate ORIGIN should be POSITIVE here because it is
    % defined as the vector FROM the geometric center of the force plate
    % (defined by the centroid of the CORNERS) to the FP TRUE ORIGIN. This
    % value is expressed in the FP coordinate system and used to calculate 
    % the correct CoP based on http://www.kwon3d.com/theory/grf/cop.html. 
    % Since the vector is downward and the FP Z-axis is downward, the Z
    % value of the ORIGIN must be negative. Note that in the c3dFormat, 
    % typically this vector is specified in reverse -> therefore, we do a 
    % check here and flip the vector if necessary.
    if FP.localOriginsFP(i,3) < 0
        FP.localOriginsFP(i,1) = -1*FP.localOriginsFP(i,1);   % reverse sign for X Origin
        FP.localOriginsFP(i,2) = -1*FP.localOriginsFP(i,2);   % reverse sign for Y Origin
        FP.localOriginsFP(i,3) = -1*FP.localOriginsFP(i,3);   % reverse sign for Z Origin
        fprintf('Note: FP%d ORIGIN vector has been negated.\n', plateNum);
    end
    
    % Get force plate surface centers (global vicon coordinates)
    % corners must be labelled sequentially
    % (bug fixed TD June 2009)
    FP.viconOrig2FPcenterSurfaceGLOB(i,x) = (FP.corners(3,x,i) - ...
        FP.corners(1,x,i))/2 + FP.corners(1,x,i);
    FP.viconOrig2FPcenterSurfaceGLOB(i,y) = (FP.corners(3,y,i) - ...
        FP.corners(1,y,i))/2 + FP.corners(1,y,i);
    FP.viconOrig2FPcenterSurfaceGLOB(i,z) = mean(FP.corners(:,z,i));
    
    % Convert the vector from FP true origin --> FP top surface center
    % from FP coordinate system to VICON global coordinate system
    if opt == 0,
        FP.localOriginsGLOB(i,:) = coordChange(FP.localOriginsFP(i,:)', C3Dkey.transform.FPVICON{i})';
    end
    
end


% Vector addition from VICON ORIGIN --> CENTER OF PLATE SURFACE
% and CENTER OF PLATE SURFACE --> PLATE TRUE ORIGIN (given by FP.localOriginsGLOB)
% (all represented in global VICON coordinate system)
if opt == 0,
    FP.viconOrig2FPorigGLOB = FP.viconOrig2FPcenterSurfaceGLOB + FP.localOriginsGLOB;
end

% Convert the force plate corner data into the model coordinate system
% (it may be usful later on)
cornersNew = [];
for i = 1:C3Dkey.numPlatesUsed
    cornersNew(:,:,i) = coordChange(FP.corners(:,:,i)', C3Dkey.transform.VICMODEL)';
end
FP.corners_model = cornersNew;


