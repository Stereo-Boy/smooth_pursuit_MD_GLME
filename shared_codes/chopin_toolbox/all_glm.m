function mdls = all_glm(model, verbose)
% This function tries to apply all combinations of factors to determine the best GLM (Generalized Linear Model) using model comparison.
% It needs a structure (model) containing the following fields:
%   model.solid_factors = {'name_of_factor_or_list'}; % a factor or a list of factors that are always included in the model (use '' for none)
%   model.liquid_factors =   {'name_factor1','name_factor2','name_factor1:name_factor2'}; %  a list of possible factors to be included, that can be removed if needed, and the interactions terms to explore
%   model.max_nb_factors = 5; % the maximal nb of factors to explore in the model - as a rule of thumb, you need ~10 datapoints for each
%   model.dv = 'dependent_variable'; % the name of the dependent variable
%   model.distribution = 'normal'; % its distribution among poisson, normal, gamma, inverse gaussian, binomial
%   model.data = data; % a table with the data
%   model.links = {'log', 'reciprocal','identity','-2','-3','probit','logit', 'loglog','comploglog'}; %  a list of possible link functions
%   model.glme = 0; % whether to use a GLM (0) or a GLME (1).
%   model.p_adjust_method = 'benjamini-hochberg'; % optional - method for adjusting p for multiple comparison - default 'none', other options are 'benjamini-hochberg', 'bonferroni' - optionnaly, you can enter the nb of hypothesis tests in model.nb_tests - actually used in display_model only
%   model.nb_tests = 3; % optional - nb of hypothesis tests used for p adjustment - if not specified, the nb of tests is the nb of factors in the best model - actually used in display_model only
% Example of use:
%           model.solid_factors = {'meditation'}; %keep these between {} - intercept is included by default, use '-1' as a solid factor to remove intercept
%           model.liquid_factors = {'music','sport','expect','music:meditation','sport:meditation','expect:meditation'}; %keep these between {}
%           model.data = data;
%           model.max_nb_factors = 5;
%           model.warning_off = 1;
%           model.dv = 'initial_work_mem';
%           model.distribution = 'normal';
%           model.links = { 'identity'}; %   model.links = {'log', 'reciprocal','identity'};
%           model.exclude = [8,12];
%           model.glme = 1; %default 0
%     mdls = all_glm(model);
%     display_model(mdls{1}, model)
%     h=subplot(1,4,4); plot_group_effect(data.initial_work_mem, data.meditation, h, 'Meditation group', 'initial working memory performance', {'Meditators','Non-meditators'})
%     saveas(gcf,fullfile(figure_path,'working_memory_initial_glm.png')); 

if isfield(model,'warning_off') || model.warning_off==0; end
if isfield(model,'exclude') ; exclude = 1; else; exclude = 0; end
if ~isfield(model,'glme') ; model.glme = 0; end
if ~isfield(model,'p_adjust_method'); model.p_adjust_method = 'none'; end
if ~exist('verbose','var')||isempty(verbose); verbose = 1; end

if exclude % here I prefer to exclude the observations, rather than using the Exclude option in fitglm, otherwise, the excluded data are then wrongly reincorporated in the diagnostic plots.
   model.data(model.exclude,:) = []; 
end
rot_fact_nb = model.max_nb_factors - numel(model.solid_factors);
skip=0;
if any(strcmp(model.solid_factors,'-1')) % detected -1 as a solid factor to remove intercept
    formula_start = [model.dv,' ~ -1'];    %removing intercept
    model.solid_factors(strcmp(model.solid_factors,'-1')) = []; % removing -1 as a factor
else
    formula_start = [model.dv,' ~ 1'];
end
for ii=1:numel(model.solid_factors)
    if ~isempty(model.solid_factors{ii})
        formula_start = [formula_start,' + ',model.solid_factors{ii}];
    end
end

if rot_fact_nb>0 % generates a list of models with various liquid factors to test
    list_models = cell(1,1);
    n=1;
    for i=1:min(rot_fact_nb,numel(model.liquid_factors)) %correct max nb to nb of liquid factors
        list_models_nb = nchoosek(1:numel(model.liquid_factors),i);
        for j=1:size(list_models_nb,1) % each line is a permutation of factors
            list_models(n+(j-1),1) = {model.liquid_factors(list_models_nb(j,:))};
        end
        n=n+size(list_models_nb,1);
    end
else
    skip = 1;
    if model.glme == 1 % This intend at finding all (1|factor) to move them at the end of the formula to avoid a weird bug with fitglme
       formula_start = moveSubstringToEnd(formula_start);
    end
    formulas = formula_start;
end

% generates the model formulas to test
if skip==0
    formulas = cell(size(list_models,1),1);
    for i=1:size(list_models,1)
        formula = formula_start;
        this_mdl = list_models{i};
        for j=1:numel(this_mdl)
            if ~isempty(this_mdl{j})
                formula = [formula,' + ', this_mdl{j}];
            end
        end
        if model.glme == 1 % This intend at finding all (1|factor) to move them at the end of the formula to avoid a weird bug with fitglme
            formula = moveSubstringToEnd(formula);
        end
        formulas{i}= formula;
    end
end

% generates the possible link functions to test
model.links = cellfun(@get_distr,model.links,'UniformOutput',false);

% run all the models formulas into glm
mdls = cell(size(formulas,1)*numel(model.links),1);
mdl_formulas = mdls; mdl_links = mdls; mdl_aiccs = zeros(numel(mdls),1); mdl_r2_adj = zeros(numel(mdls),1); mdl_r2 = mdl_r2_adj; norm_res = mdls;
try
if verbose==1
    if model.glme==0
        dispi('Running ',numel(formulas).*numel(model.links),' GLMs...');
    else
        dispi('Running ',numel(formulas).*numel(model.links),' GLMEs...');
    end
end
for i=1:size(formulas,1)
    for j=1:numel(model.links)
        idx = numel(model.links)*(i-1)+j;
        if size(formulas,1)==1
            formula = formulas;
        else
            formula = formulas{i};
        end
        if model.glme==0
            mdl = fitglm(model.data,formula,'Distribution',model.distribution,'Link',model.links{j});
        else
            mdl = fitglme(model.data,formula,'Distribution',model.distribution,'Link',model.links{j});
        end
        mdls{idx} = mdl;
        mdl_formulas{idx} = formula;
        mdl_links{idx} = model.links{j};
        if model.glme==0
            mdl_aiccs(idx) = mdl.ModelCriterion.AICc;
        else
            mdl_aiccs(idx) = mdl.ModelCriterion.AIC; % AICc is not available for GLMEs, let's use AIC
        end
        mdl_r2_adj(idx) = round(100.*mdl.Rsquared.Adjusted,1);
        mdl_r2(idx) = round(100.*mdl.Rsquared.Ordinary,1);
        if sum(isnan(mdl.Residuals.raw))==numel(mdl.Residuals.raw)
            %all values are nan, something went wrong
            norm_res{idx} = 'no';
        else
            H = kstest((mdl.Residuals.raw-nanmean(mdl.Residuals.raw))./nanstd(mdl.Residuals.raw));
            if H==1; norm_res{idx} = 'no'; else; norm_res{idx} = 'yes'; end
        end
    end
end

if verbose==1
    dispi('We tested ',numel(mdl_aiccs),' models.')
end
if model.glme==0
    models = sortrows(table((1:numel(mdl_aiccs))',mdl_formulas,mdl_links,mdl_aiccs,mdl_r2_adj,mdl_r2,norm_res,'VariableNames',{'Rank','formula','link','AICc','adj.R2(%)','R2(%)','norm.res.'}),'AICc');
else
    models = sortrows(table((1:numel(mdl_aiccs))',mdl_formulas,mdl_links,mdl_aiccs,mdl_r2_adj,mdl_r2,norm_res,'VariableNames',{'Rank','formula','link','AIC','adj.R2(%)','R2(%)','norm.res.'}),'AIC');
end
mdls = mdls(models.Rank); %reorder mdls so it is in the same order as models
models.Rank = (1:numel(mdl_aiccs))'; %make their rank increase too
if verbose==1
    disp(models)
    disp(' ------------------------------------------------------------------------------- ')
end
catch err
    disp('Error caught: for debugging, write rethrow(err)')
    keyboard
end
end

function link = get_distr(linkn)
    switch linkn 
        case{'reciprocal','Reciprocal','inverse'}
            link = 'reciprocal';
        case{'identity','Identity'}
            link = 'identity';
        case{'log', 'Log','log10','ln','Log10','Ln'}
            link = 'log';
        case{'logit' , 'Logit'} 
            link = 'logit';
        case{'Probit' , 'probit'} 
            link = 'probit';
        case{'-2','inverse_square','reciprocal_square'}
            link = -2;
        case{'Loglog','loglog'}
            link = 'loglog';
    end
end

function modifiedString = moveSubstringToEnd(inputString)
    % Regular expression pattern to match '+ (' and everything until ')'
    % This intend at finding all (1|factor) to move them at the end of the formula to avoid a weird bug with fitglme
    pattern = '\+\s*\([^)]*\)';
    
    % Find all matches in the input string
    matches = regexp(inputString, pattern, 'match');
    % Replace all matches with an empty string in the original string
    modifiedString = regexprep(inputString, pattern, '');
    
    % Append the matched substrings at the end
    for i = 1:length(matches)
        modifiedString = [modifiedString ' ' matches{i}];
    end
end

