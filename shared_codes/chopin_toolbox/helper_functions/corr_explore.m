function [coef, pval] = corr_explore(x,y,type,header,figH,subfigH,xname,yname)
% This function explores the correlation between x and y using type
% ('Spearman', 'Pearson' or 'Kendall') correlation, plotting the results
% and formating the stats on command window using corr_format
% Header is printed in front and can be let empty like type.
% It plots in figure handle fig and subplot of handle subfigH and names the axes
% with xname and yname 
%
% Ex: fig=figure('Color', 'w'); h1=subplot(2,3,1); 
% corr_explore(x,y, 'Spearman','action group',fig,h1,'MOT post', 'Dual nback slope')

if ~exist('type','var')||isempty(type)||(strcmpi(type,'Spearman')==0&&strcmpi(type,'Pearson')==0&&strcmpi(type,'Kendall')==0)
    type = 'Spearman'; %default
end
if ~exist('header','var')||isempty(header)
    header = '';
end
if ~exist('xname','var')||isempty(xname)
    xname = 'x';
end
if ~exist('yname','var')||isempty(yname)
    yname = 'y';
end

if numel(size(x))>2||numel(size(y))>2; warning('corr_explore: Data should be vectors - not sure what will happen here...'); end

% we put x and y as columns
if size(x,2)>1 && size(x,1)==1; x=x'; end
if size(y,2)>1 && size(y,1)==1; y=y'; end

% we remove lines with nans
data = nonan([x,y],1);
x=data(:,1); y=data(:,2);

figure(figH)
subplot(subfigH)
    hold on
    scatter(x,y);
    ylabel(yname);
    xlabel(xname);
   [b] = robustfit(x, y, 'bisquare');
    xx = sort(x);
    Y_hat = b(1)+xx.*b(2);
    line(xx,Y_hat, 'LineWidth', 2,'Color', 'black');
    [coef, pval] = corr_format(x, y, type, header);
    txt1 = ['r = ', num2str(coef,'%.3f'), ' p = ', num2str(pval,'%.3f') ];
    fontsz = 10;
    text(0.2,0.2, txt1, 'Units', 'Normalized', 'Fontsize', fontsz);
    title(subfigH,header);
    set(gca,'TickDir','out');
end