function [bf10, k_effect] = bayes_factor_ttest(x, y, format)
% Bayes Factor for the equivalent of a PAIRED group t-test 2
% k_effect is a string interpreting the strength of the effect
% if format = 1, disp the stats in a formatted way
% K = BF10

%computing Bayes Factor for independent t-test (anova with one factor) and interpreting its strength
if exist(which('bf.ttest2'),'file')
    %removes nan first
  %  if any(isnan(x))==1; x(isnan(x)==1)=[]; y(isnan(x)==1)=[]; end
  %  if any(isnan(y))==1; y(isnan(y)==1)=[]; x(isnan(y)==1)=[]; end
    [bf10,p,CI,stats] = bf.ttest2(x,y);
    %interpretation follows Jeffreys, Harold (1998) [1961]. The Theory of Probability (3rd ed.). Oxford, England. p. 432. ISBN 9780191589676.
    if bf10<1
        k_effect = 'none';
    elseif bf10<5
        k_effect = 'barely worth mentioning';
    elseif bf10<10
        k_effect = 'substantial';
    elseif bf10<10^(3/2)
        k_effect = 'strong';
    elseif bf10<100
        k_effect='very strong';
    else
        k_effect='decisive';
    end
else
    bf10 = nan;
    k_effect = 'N/A';
end

if format==1
    dispi('Bayes Factor BF10 = ',sprintf('%.2f',bf10),' (',k_effect,')');
end
