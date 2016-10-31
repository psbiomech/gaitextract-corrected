% Change coordinate system between Vicon / Force Plate / 
% OR any other arbitrary coordinate system
% Tim Dorn
% Updatae: September 2010
% 
% --------------------------------------------------------------------
% Usage: newSysData = coordChange(oldSysData, transform)
% --------------------------------------------------------------------
% 
% Inputs:   oldSysData = matrix of a multiple of 3 ROWS in old coord sys
%                       e.g. [x, y, z] OR [x1, y1, z1, x2, y2, z2]
% 
%           transform = transformation vector [x, y, z] form
%                   OR transformation matrix [xx xy xz ; yx yy yz ; zx xy zz]
% 
%
% Outputs:  newSysData = matrix of a multiple of 3 ROWS in new coord sys
%                       e.g. [x, y, z] OR [x1, y1, z1, x2, y2, z2] 
% 
% 
% Notes
% -----
% 
% newSystem = sign(transform) * oldSysData(|transform|)
% 
% No time vectors should be in the data set as inputs. 
% No time vectors are returned as output.
% 
% --------------------------------------------------------------------

function newSysData = coordChange(oldSysData, transform)

usage = 'Usage: newSysData = coordChange(oldSysData, transform)';

if nargin ~= 2,
    disp(usage)
    return
end



% Check input data for correctness
% --------------------------------
if isvector(transform) && length(transform) == 3  % a direction vector
    transform = make3DTransform(transform);
elseif sum(size(transform)==[3 3]) ~= 2   % not a 3x3 matrix
    error('Incorrect transform input parameters...\n%s', usage)
end

% Perform the transformation
% --------------------------------
newSysData = [];
[m,n] = size(oldSysData);
if mod(m,3) == 0,
    r = floor(m/3);
    for i = 1:r
        tmp = oldSysData(3*(i-1)+1:3*(i-1)+3,:);
        
        % Transformation done HERE
        newSysData = [newSysData ; (transform * tmp)];
    end
else
    error('Input data must have a multiple of 3 rows i.e. (x,y,z) / (x,y,z,x,y,z)...\n%s', usage)
end


