function [K, k_effect] = bayes_factor_ttest2(x, y, format)
% Bayes Factor for the equivalent of an INDEPENDANT group t-test
% k_effect is a string interpreting the strength of the effect
% if format = 1, disp the stats in a formatted way
% K = BF10

%computing Bayes Factor for independent t-test (anova with one factor) and interpreting its strength
if exist(which('bf.anova'),'file')
    %removes nan first
    if any(isnan(x))==1; x(isnan(x)==1)=[]; y(isnan(x)==1)=[]; end
    if any(isnan(y))==1; y(isnan(y)==1)=[]; x(isnan(y)==1)=[]; end
    data=cell2table([num2cell(x), cellstr(num2str(ones(size(x)))); num2cell(y), cellstr(num2str(2.*ones(size(y))))],'VariableNames',{'VD','factor'});
    K = bf.anova(data,'VD~factor');
    %interpretation follows Jeffreys, Harold (1998) [1961]. The Theory of Probability (3rd ed.). Oxford, England. p. 432. ISBN 9780191589676.
    if K<1
        k_effect = 'none';
    elseif K<5
        k_effect = 'barely worth mentioning';
    elseif K<10
        k_effect = 'substantial';
    elseif K<10^(3/2)
        k_effect = 'strong';
    elseif K<100
        k_effect='very strong';
    else
        k_effect='decisive';
    end
else
    K = nan;
    k_effect = 'N/A';
end

if format==1
    dispi('Bayes Factor K = ',sprintf('%.2f',K),' (',k_effect,')');
end
