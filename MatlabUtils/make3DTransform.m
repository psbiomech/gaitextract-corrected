% ========================================================================
% FUNCTION: make3DTransform
% mat = make3DTransform(vec)
% Tim Dorn
% Sept 2010
% ========================================================================
function mat = make3DTransform(vec)

l = length(vec);
mat = zeros(l,l);
for i = 1:length(vec)
    mat(i,abs(vec(i))) = 1*sign(vec(i));
end

