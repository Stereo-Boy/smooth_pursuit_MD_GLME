function effect_sizes_table = effect_sizes(mdl, model)
% Calculate local Cohen's f2 (f squared) for mixed effect model.
% For each continuous predictor, we follow the rational in Selya et al. (2012) - A Practical Guide to Calculating Cohents f2, a Measure of Local Effect Size, from PROC MIXED.
% It will output the local effect sizes for each continous factor in the model.
%
% Inputs:
% mdl is one particular cell of the output from all_glm
% model is a model structure
%
% Example of usage:
% model.data = data;
% model.dv = 'Time';
% model.distribution = 'normal';
% model.max_nb_factors = 10;
% model.solid_factors = {''}; %keep these between {}
% model.liquid_factors = {'(1|Participant)','(1|task)','trial','load','stereo','ageGroup','ageGroup:stereo','load:stereo','ageGroup:load','load:ageGroup:stereo'}; %keep these between {}
% model.links = {'log'}; %{'log','identity'}; % identity was also tested but consistently underperformed
% mdls = all_glm(model);
% effect_sizes(mdls{1}, model)

try

effect_sizes_table = [];

% first check that dependent variable is continuous, and if yes, continue, otherwise just stop
if iscategorical(mdl.ResponseName)
    disp('Dependent variable is not continuous so we cannot easily extract a local effect size.')
    return
end

f2s = []; f2sizes = {}; types = {};
% save a list of all factors in the model and redefine a model structure accordingly
all_factors_in_model = split_factors(mdl.Formula);
model.liquid_factors = {};
model.p_adjust_method = 'none';

% find the list of factors in the mdl structure and continue redefininge a model structure by choosing a unique distribution and link function
if model.glme
    mdl_factor_list = mdl.Coefficients.Name; %GLME
    model.distribution = mdl.Distribution; % GLME case
else
    mdl_factor_list = mdl.Coefficients.Row; %GLM
    model.distribution = mdl.Distribution.Name; % GLM case
end
model.links = {mdl.Link.Name}; 

% loop through that list of factors to extract effect sizes depending on whether they are interaction terms or not, and if not, continous or categorical
for f = 1:numel(mdl_factor_list)
    factor = mdl_factor_list{f};
    if strcmp(factor,'(Intercept)') % skip intercept factor
        continue
    end
    underscores = strsplit(factor,'_');    % check for categorical terms (be sure to have initially clean factor names of any _ )
    if numel(underscores)==1       % continuous factor only - calculate local Cohen's f2
        model.solid_factors = all_factors_in_model(~cellfun(@(x) strcmp(x,factor), all_factors_in_model)); % select all factors but the one of interest
        model.max_nb_factors = numel(model.solid_factors);
        mdls = all_glm(model,0); % run the model with verbose off
        f2s(f) = round((mdl.Rsquared.Ordinary - mdls{1}.Rsquared.Ordinary)/(1-mdl.Rsquared.Ordinary),2);
        types(f) = {'Cohen s f2'};
    else % categorical case (or at least one categorical in an interaction): calculate cohen's d
        % colons = strsplit(factor,':');
        % % interaction case
        % if numel(colons)>1 
        % 
        % end
        idx = find(strcmp(mdl_factor_list,factor));
        if numel(idx)==0 % problem detected
                warning('We are not able to find the factor of interest in the list of factors - check code.');
                f2s(f) = nan;
                types(f) = {'None'};
        else % all good, calculate Cohen's d
            t = mdl.Coefficients.tStat(idx); % read the t-stat
            if model.glme
                df = mdl.Coefficients.DF(idx); % degrees of freedom - GLME
            else
                df = mdl.DFE; % GLM
            end
            f2s(f) = round(2*abs(t)/sqrt(df),2); % this is actually a cohen's d using formula from Rosenthal and Rosnow, 1991
            types(f) = {'Cohen s d'};
        end
    end
    f2sizes(f) = {interpret_f2(f2s(f),types{f})};
    dispi('Local effect size for ',factor,' : ',types{f},' = ',f2s(f),' (',f2sizes{f},')')
end
effect_sizes_table = table(mdl_factor_list,types',f2s',f2sizes','VariableNames',{'Factor','Type','ES','Interpretation'});

catch err
    disp('Error caught: for debugging, write rethrow(err)')
    keyboard
end
end

function new_terms = split_factors(formula) 
    try
    % Split the formula into terms and dv and deal with special cases like * interaction reintroduced elsewhere 
        % first convert the model class in string
        formula = char(formula);
    
         % remove '~' and leading/trailing spaces
        terms = strsplit(formula, '~');
        terms = strtrim(terms{2});
    
        % split between + terms
        terms = strsplit(terms, '+'); 
        
        % remove the term '1'
        terms = terms(~strcmp(strtrim(terms), '1'));
        
        % remove empty spaces
        terms = cellfun(@(x) strrep(x, ' ', ''), terms, 'UniformOutput', false);
        
        % remove empty cells if any
        terms(cellfun(@isempty, terms))=[];
    
        % split x*y interaction factors (fit_glm reintroduces this notation in the formula)
        new_terms = {};
        for i = 1:numel(terms)
            celli = terms{i}; 
            if contains(celli, '*')  
              % Split the cell at the asterisk
              [parts1, ~] = strsplit(celli, '*');
              % Create three new cells from the split parts
              new_terms(end+1) = parts1(1); 
              new_terms(end+1) = parts1(2); 
              new_terms(end+1) = {[parts1{1}, ':', parts1{2}]};
            else
              new_terms{i} = celli;
            end
        end
    
        % remove potential doublon factors
        new_terms = unique(new_terms);
    
    catch err
        disp('Error caught: for debugging, write rethrow(err)')
        keyboard
    end
end

function f2size = interpret_f2(f2,type)
% finding effect size using guidelines for interpretation of f2 indicating that 0.02 is a small effect, 0.15 is a medium effect, and 0.35 is a large effect (Cohen 1992)
    % % removes parentheses in potential (X^2) factor
    % factor(strfind(factor,'(')) = []; 
    % factor(strfind(factor,')')) = []; 
f2size = 'dubious';
    switch type
        case {'Cohen s f2'}
            if f2>=0.35
                f2size = 'large';
            elseif f2>=0.15
                f2size = 'medium';
            elseif f2>=0.02
                f2size = 'small';
            end
        case {'Cohen s d'} % using Cohen's (1992) guidelines expanded by Sawilowsky (2009)
            if f2>=2
                f2size = 'huge';
            elseif f2>=1.2
                f2size = 'very large';
            elseif f2>=0.8
                f2size = 'large';
            elseif f2>=0.5
                f2size = 'medium';
            elseif f2>=0.2
                f2size = 'small';
            elseif f2>=0.01
                f2size = 'very small';
            end
        case {'None'}   % not calculated
            f2size = 'N/A'; 
        otherwise
            warning('unrecognized effect size type')
            f2size = 'N/A';
    end
end