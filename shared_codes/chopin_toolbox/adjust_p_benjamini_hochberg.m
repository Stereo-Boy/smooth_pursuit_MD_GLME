function table = adjust_p_benjamini_hochberg(table, additional_comparisons, model)
% Adjust p values and hypothesis test in the stat table provided by the function display_model, and display it
% table is the stat table provided by the function display_model and should at least contain the column pValue
% optional additional_comparisons: a number reflecting how many other comparisons have been already tested, 
% as part of a primary hypothesis test structure for example
% model is an optional structure that would contain the field model.nb_tests, the nb of tests to correct the table for (if not provided, 
% default is the number of lines in table)
%
% typical use:
% model.p_adjust_method = 'none';
% table_stats_1 = display_model(mdls1{1});
% table_stats_2 = display_model(mdls2{1});
% adjust_p_benjamini_hochberg([table_stats_1; table_stats_2])
if ~exist('model','var') || isfield(model,'nb_tests')==0; model.nb_tests = size(table,1); end
if ~exist('table','var') || ~ismember('pValue',table.Properties.VariableNames); disp('[adjust_p_benjamini_hochberg: incorrect table provided - exiting.]'); return; end
if ~exist('additional_comparisons','var'); additional_comparisons = 0; end

% adjusting the nb of comparisons to reflect potential past comparisons
model.nb_tests =  model.nb_tests + additional_comparisons;
 
disp('Adjustment for multiple comparisons: method of Benjamini-Hochberg')

% reorder the factors by increasing p values
[~,idx]=sort(table.pValue);
table = table(idx,:);

% define thresholds
alpha = 0.05;
H0_reject = zeros(size(table,1),1); % hypotheses rejection (1)
adj_alpha = ones(size(table,1),1).*alpha; % values of adjusted alphas

first_no_reject = 0; % as soon as one hypothesis is not rejected, the larger p-values are also rejected
for i=1:size(table,1)
    % apply p adjustment
    adj_alpha(i) = i.*alpha./model.nb_tests;
    if table.pValue(i)<=adj_alpha(i) && first_no_reject == 0
        H0_reject(i) = 1;
    else
        first_no_reject = 1;
    end
end
table.adj_alpha = adj_alpha;
table.adj_pValue = min(1,table.pValue.*alpha./adj_alpha);
table.H0_reject = H0_reject;

% show the results
disp(table)
for j=1:size(table,1)
    if table.H0_reject(j) == 1 %reject
        result = 'Significant';
    else
        result = 'No significant';
    end
    if ismember('Name',table.Properties.VariableNames) && ismember('DF',table.Properties.VariableNames) && ismember('tStat',table.Properties.VariableNames)
        dispi(result,' effect of ',table.Name{j},' (t(',table.DF(j),') = ',round(table.tStat(j),2),', adjusted p = ',round(table.adj_pValue(j),4),')')
    end
end