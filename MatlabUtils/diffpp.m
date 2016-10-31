function ppdX = diffpp(ppX)
% Differentiate a function defined as a piecewise polynomial (pp).
% Returns a piecewise polynomial of the differentiated input function.
% USAGE = ppdX = diffpp(ppX)
%
% Ajay Seth
% June 25, 2004

[breaks, coefs, np, op, nd] = unmkpp(ppX);

% breaks are points defining the np intervals
% np is lengths(breaks)-1
% op is the pp order+1 which defines the number of coefs per pp
% nd is the number of curves being fit on the same abscissa

% coefs are arranged so that the for interval the coefficients of 
% of the nd curves appear as rows => np*nd rows and op columns
% to differntiate the pp you get reduce the order, (eliminate the last
% column) and multiply each column by its associated order in descending
% order (op-1:-1:1);
dcoefs = coefs(:,1:op-1).*(ones(np*nd,1)*(op-1:-1:1));

% reassemble the spline by calling MATLAB's pp maker function
ppdX = mkpp(breaks, dcoefs, nd);
        