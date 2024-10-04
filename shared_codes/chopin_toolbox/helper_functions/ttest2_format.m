function [h,p,ci,stats] = ttest2_format(x,y,header, varargin)
% Does a t-test for 2 independant samples and formats the results, printing
% the header at the beginning on the command window
% Issue an effect size (Hedge's g) and interpret it, and an interpreted bayes Factor
% All varargin are transfered to ttest2, following the same rules.
% To use the Bayes Factor, this function requires the BayesFactor library https://klabhub.github.io/bayesFactor/
% Install it to the path so that we can use it here. It follows Rouder, J. N., Morey, R. D., Speckman, P. L. & Province, J. M. Default Bayes factors for ANOVA designs. J. Math. Psychol. 56, 356?374 (2012).
% To test for the normality of the data, we provide a test: you can choose which one by adding 'normality' and either
% 'kolmogorov' or 'shapiro' respectively for the Kolmogorov-Smirnov tests and the Shapiro-Wilk test (default).
% If the differences are not normally distributed, we also add a non-parametric test. 
% By default, we are using a Welch test assuming unequal variance with Satterthwaite''s approximation for correcting degrees of freedom'.
% Typical use: ttest2_format(data(cond==1),data(cond==1),'cond 1 - cond 2')
if any(strcmpi(varargin,'alpha')); i=find(strcmpi(varargin,'alpha')); alpha = varargin{i+1}; else; alpha = 0.05; end
if any(strcmpi(varargin,'dim')); i=find(strcmpi(varargin,'dim')); dim = varargin{i+1}; else; dim = []; end
if any(strcmpi(varargin,'tail')); i=find(strcmpi(varargin,'tail')); tail = varargin{i+1}; else; tail = 'both'; end
if any(strcmpi(varargin,'vartype')); i=find(strcmpi(varargin,'vartype')); vartype = varargin{i+1}; else; vartype = 'unequal'; end
if any(strcmpi(varargin,'normality')); i=find(strcmpi(varargin,'normality')); norm_test = varargin{i+1}; else; norm_test = 'shapiro'; end
if ~exist('header','var'); header = ''; end

%reshaping x and y so that row data are subjects
if size(x,2)>size(x,1); x=x'; end
if size(y,2)>size(y,1); y=y'; end

%computing t-test
[h,p,ci,stats] = ttest2(x,y,'alpha',alpha,'dim',dim,'tail',tail,'vartype',vartype);

%compting Hedge's g and interpreting its strength
% Note: according to Cohen and Sawilowsky:
%
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
d = computeCohen_d(x,y, 'independent'); % this is actually an Hedge's g because it uses sample estimate of standard deviation, not population estimate
if abs(d)<0.2
    d_effect = 'none';
elseif abs(d)<0.5
    d_effect = 'small';
elseif abs(d)<0.8
    d_effect = 'medium';
else
    d_effect = 'large';
end

[K, k_effect] = bayes_factor_ttest2(x,y,0);

if strcmp(tail,'both'); tail_add=''; else;    tail_add = 'one-tailed ';  end
if p<0.0001
    dispi(header, ' - ', 'Independent-sample Welch ',tail_add,'t-test t(',stats.df,') = ',sprintf('%.2f',stats.tstat),', p = ',sprintf('%.6f',p),...
    ', Hedge''s g = ',sprintf('%.2f',d),' (',d_effect,'), BF10 = ',sprintf('%.2f',K),' (',k_effect,')');
else
    dispi(header, ' - ', 'Independent-sample Welch ',tail_add,'t-test t(',stats.df,') = ',sprintf('%.2f',stats.tstat),', p = ',sprintf('%.4f',p),...
    ', Hedge''s g = ',sprintf('%.2f',d),' (',d_effect,'), BF10 = ',sprintf('%.2f',K),' (',k_effect,')');
end
%if strcmpi(vartype,'unequal'); disp('Unequal variance: using Welch test with Satterthwaite''s approximation for correcting degrees of freedom') ;end

%providing a normality test
switch norm_test
    case 'kolmogorov'
        %detecting and removing nan values to be able to run Kolmogorov-Smirnov tests for normality
        x2=x; y2=y; 
        if any(isnan(x))||any(isnan(y))
            tot=sum(isnan(x))+sum(isnan(y));
            dispi('Found a total of ',tot,' nan values - removing them for KS test');
            x2(isnan(x))=[]; y2(isnan(y))=[];
        end
        [H1, P1, KSstat1] = kstest(zscore(x2)); [H2, P2, KSstat2] = kstest(zscore(y2));
        dispi('Kolmogorov-Smirnov tests for normality:  KS = ',sprintf('%.2f',KSstat1),', p = ',sprintf('%.4f',P1), ...
            ' and KS = ',sprintf('%.2f',KSstat2),', p = ',sprintf('%.4f',P2));
    case 'shapiro'
        SH_stats_x = shapiro(x); SH_stats_y = shapiro(y);
        P1 = SH_stats_x{2,7}; H1 = (P1<0.05); w1 = SH_stats_x{2,5};
        P2 = SH_stats_y{2,7}; H2 = (P2<0.05); w2 = SH_stats_y{2,5};
        dispi('Shapiro-Wilk test for normality (alpha 5%):  W = ',sprintf('%.2f',w1),', p = ',sprintf('%.4f',P1),...
            ' and W = ',sprintf('%.2f',w2),', p = ',sprintf('%.4f',P2));
    otherwise
        error(['Normality test requested not found: ',norm_test]);
end

if H1==1 || H2==1
    [pp, hh, statmww] = ranksum(x,y,'alpha',alpha,'tail',tail);
    if pp<0.0001
        dispi('Non-normal distributions: practicing a non-parametric ',tail_add,'Wilcoxon-Mann–Whitney (signed-ranks) test too: U = ',sprintf('%.0f',statmww.ranksum),', p = ',sprintf('%.6f',pp))
    else
        dispi('Non-normal distributions: practicing a non-parametric ',tail_add,'Wilcoxon-Mann–Whitney (signed-ranks) test too: U = ',sprintf('%.0f',statmww.ranksum),', p = ',sprintf('%.4f',pp))   
    end
    
    %replace output stats with non-parametric
    h = hh; p = pp; stats = statmww; ci = [];
else
    disp('Distributions do not differ from normality.')
end

