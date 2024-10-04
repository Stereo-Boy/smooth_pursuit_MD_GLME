function [coef, pval] = corr_format(x, y, type, header, nanIgnore)
% Find correlation between x and y using type ('Spearman', 'Pearson' or
% 'Kendall') correlations AND format it nicely on the command window. 
% Header is printed in front and can be omitted
% If nanIgnore == 1 (default), then remove nan lines
% 
% It also reads the effet size Cohen's f2 and its interpretation following Cohen’s (1988)
% guidelines, f2? 0.02, f2? 0.15, and f2 ? 0.35 represent small, medium, and large effect sizes
% typical use: corr_format(x, y, 'Spearman', 'x vs. y');

% default values
if ~exist('type','var')||isempty(type)||(strcmpi(type,'Spearman')==0&&strcmpi(type,'Pearson')==0&&strcmpi(type,'Kendall')==0);    type = 'Spearman'; end
if ~exist('header','var');    header = '';  end
if ~exist('nanIgnore','var');    nanIgnore = 1;  end

% orient data in columns
if size(x,2)>size(x,1); x = x'; end
if size(y,2)>size(y,1); y = y'; end

% remove lines with nan if necessary
data = [x, y];
data(any(isnan(data),2),:) = []; 
x = data(:,1);  y = data(:,2);

%calculating r
[coef, pval] = corr(x, y,'type',type);

%calculating coefficient of determination r2 (from Pearson!)
r = corr(x, y,'type','Pearson');
if abs(r)<0.1
    r_effect = 'none';
elseif abs(r)<0.3
    r_effect = 'small';
elseif abs(r)<0.5
    r_effect = 'medium';
else
    r_effect = 'large';
end
r2 = r.^2;

%calculating effect size Cohen's f2
f2 = r2/(1-r2);
if f2<0.02
    f2_effect = 'none';
elseif f2<0.15
    f2_effect = 'small';
elseif f2<0.35
    f2_effect = 'medium';
else
    f2_effect = 'large';
end
if pval<0.0001
    dispi(header,' - ',type,' correlation: r = ', num2str(coef,'%.2f'),', p = ',...
        num2str(pval,'%.6f'),', R2 = ',num2str(r2,'%.2f'),' (',r_effect,'), Effect size Cohen''s f2 = ',num2str(f2,'%.2f'),' (',f2_effect,')')
else
    dispi(header,' - ',type,' correlation: r = ', num2str(coef,'%.2f'),', p = ',...
        num2str(pval,'%.4f'),', R2 = ',num2str(r2,'%.2f'),' (',r_effect,'), Effect size Cohen''s f2 = ',num2str(f2,'%.2f'),' (',f2_effect,')')
end

if strcmp(type,'Pearson') %check for normality
   SH_stats1 = shapiro(x);        P1 = SH_stats1{2,7};  w1 = SH_stats1{2,5}; SH_stats2 = shapiro(y);        P2 = SH_stats2{2,7};  w2 = SH_stats2{2,5};
    dispi('Shapiro-Wilk test for normality (alpha 5%):  W = ',sprintf('%.2f',w1),', p = ',sprintf('%.4f',P1), ' and W = ',sprintf('%.2f',w2),', p = ',sprintf('%.4f',P2)); 
end
end

