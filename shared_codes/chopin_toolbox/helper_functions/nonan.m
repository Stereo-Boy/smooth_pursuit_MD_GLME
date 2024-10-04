function y = nonan(x,dim)
% remove nan in x
% if dim = 1, group data along lines and remove lines where one of the columns is a nan

if ~exist('dim','var'); dim = 0; end

if dim>0
    y = x(any(isnan(x)')'==0,:);
else
    y = x(isnan(x)==0);
end