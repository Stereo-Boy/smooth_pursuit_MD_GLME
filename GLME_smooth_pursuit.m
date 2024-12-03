%% Smooth pursuit in macular degeneration
% * This functions carries on the model statistical analysis for that project.
% * Eye movement data, visibility and d' have been preprocessed and gathered in the master datafile.
% * runs the GLM/GLME with the chopin toolbox.
% * plots various graphs
% * displays statistic conclusion messages
% * saves figures from statistical analysis in figures folder
% * Inputs: None. The function must run from the analysis folder with the master datafile being located in data. 
% * Functions used in this code should all be located in shared_codes.
% * Outputs: None. But generates figures (folder figures), statistics (command window) and generates a joint markdown output in html by pressing Publish.

function GLME_smooth_pursuit
try
 close all;
 clc;
%% Preprocessing
% Creates the correct path structure
analysis_path = fileparts(mfilename('fullpath')); % root path
addpath(genpath(fullfile(analysis_path,'shared_codes'))) % path with necessary codes
data_path = fullfile(analysis_path,'data'); % path with the masterfile with data for analysis
figure_path = fullfile(analysis_path,'figures');  % path for saving the figures

% Reads data files
file = 'data_smooth_pursuit_MD.xlsx';
data = readtable(fullfile(data_path,file),'Sheet','data'); % read file
dispi('Loaded file ',file);

% exclude the 3 instances for which participants did not follow pursuit instructions - they confirm not following instuctions and ended up with eccentricities >20
data(data.eccentricity>20,:) = [];
data(isnan(data.fraction_occluded),:) = []; % also remove a line with NaN data

% Normalize variable formats
data.subject_ID = categorical(data.subject_ID);
data.status = categorical(data.status);
data.condition = categorical(data.condition);
data = clean_names(data);

% show a sample of the data
disp('Data sample and size:')
disp(data(1:10,:))
dispi(size(data,1),' lines of data found.')

%% Primary hypotheses on d' - fixation data only
% * Explore how the data can be explained by oculomotor noise, fraction occluded and eccentricity
% * Hypotheses: 
% * more eccentric stimuli result in worse d'
% * more fraction occluded for the stimuli result in worse d'
% * more oculomotor noise result in worse d'

% restrain data to fixation
dataFixation = data(data.condition=='fix',:);
disp(dataFixation(1:10,:))
dispi(size(dataFixation,1),' lines of data found.')

% verify that distribution did not change
check_distrib_indep(dataFixation.dprime(dataFixation.status=='control'),dataFixation.dprime(dataFixation.status=='patient'),'d''');
disp('Data are not distributed differently than normal')
saveas(gcf,fullfile(figure_path,'fixation_dprime_distribution.fig')); snapnow;

% Only apply the collinearity analysis on continuous predictors
% Use Kendall correlations that are more powerful when using small samples
disp('Collinearity')
disp('Lines and columns are ordered as eccentricity,fractionOccluded and OMN.')
[Tau,PValue] = corrplot2(dataFixation(:,{'eccentricity','fractionOccluded','omn'}),'type','Kendall')
disp('No particular worries for collinearity given all correlations are well <0.8.')
saveas(gcf,fullfile(figure_path,'fixation_factor_collinearity.fig')); snapnow; close;

% define a model structure
% a table with the data, here called data
model.data = dataFixation;
% the name of the dependent variable in the data structure, here it is data.initial_work_mem
model.dv = 'dprime';
% its distribution among poisson, normal, gamma, inverse gaussian, binomial as previously determined
model.distribution = 'normal';
% the maximal nb of factors to explore in the model
model.max_nb_factors = 4; % n=10 so 100 data are used to estimate subjects' random effects (1 factor). 16 data can be used to estimate 2 more factors - we add 1 so that we can also look at the effect of non-significant factors.
model.solid_factors = {'(1|subjectID)'}; %keep these between {}
% a list of possible factors to be included, that can be removed if needed, and the interactions terms to explore
model.liquid_factors = {'eccentricity','omn','fractionOccluded','fractionOccluded:omn','eccentricity:omn','fractionOccluded:eccentricity'}; %keep these between {}
% a list of potential model links
model.links = {'identity'};
% outliers/subjects to be removed - can be left empty
model.exclude = []; 
% no warnings if 1 - careful with that option
model.warning_off = 0;
% whether to use a GLM (0) or a GLME (1).
model.glme = 1;

% run the model
mdls = all_glm(model);

% display diagnostics and results
display_model(mdls{1}, model) %plot model ranked 1 - you can select any other models by rank according to the results on the various indicators provided
saveas(gcf,fullfile(figure_path,'fixation_dprime_GLME_diagnostics.fig')); snapnow; close;
figure('Color', 'w', 'units','normalized','outerposition',[0.1 0.1 0.8 0.6]);
h=subplot(1,3,1);       plot_covariate_effect(dataFixation.dprime, dataFixation.fractionOccluded.*100, h, 'fraction occluded (%)', 'd prime', 0, 0, mdls{1},1,model); 
h=subplot(1,3,2);       plot_covariate_effect(dataFixation.dprime, dataFixation.eccentricity, h, 'eccentricity (deg)', 'd prime', 0, 0, mdls{1},1,model);       title('Fixation condition only');
h=subplot(1,3,3);       plot_interaction(dataFixation.dprime, (dataFixation.eccentricity<=median(dataFixation.eccentricity)),dataFixation.fractionOccluded.*100, h, 'fraction occluded (%)','d prime', {'eccentricity > median','eccentricity <= median'},mdls{1}, 1, model); % median is 2.36
saveas(gcf,fullfile(figure_path,'fixation_dprime_GLME_effects.fig')); snapnow; close;
disp('Interpretation: there is no good model of the data, other than an intercept model. The first 3 best gave no significant factors.')
disp('Figures are only here for comparison.')
disp('To get (non-significant) p-values, we run the model #41.')
display_model(mdls{41}, model)

%% Same analysis restrained to pursuit data only
% * restrain data to fixation
dataPursuit = data(data.condition=='pursuit',:);
disp(dataPursuit(1:10,:))
dispi(size(dataPursuit,1),' lines of data found.')

% verify that distribution did not change
check_distrib_indep(dataPursuit.dprime(dataPursuit.status=='control'),dataPursuit.dprime(dataPursuit.status=='patient'),'d''');
disp('Data are not distributed differently than normal')
saveas(gcf,fullfile(figure_path,'pursuit_dprime_distribution.fig')); snapnow;

% Only apply the collinearity analysis on continuous predictors
% Use Kendall correlations that are more powerful when using small samples
disp('Collinearity')
disp('Lines and columns are ordered as eccentricity,fractionOccluded and OMN.')
[Tau,PValue] = corrplot2(dataPursuit(:,{'eccentricity','fractionOccluded','omn'}),'type','Kendall')
disp('No particular worries for collinearity given all correlations are well <0.8.')
saveas(gcf,fullfile(figure_path,'pursuit_factor_collinearity.fig')); snapnow; close;

% define a model structure
% a table with the data, here called data
model.data = dataPursuit;
% the name of the dependent variable in the data structure, here it is data.initial_work_mem
model.dv = 'dprime';
% its distribution among poisson, normal, gamma, inverse gaussian, binomial as previously determined
model.distribution = 'normal';
% the maximal nb of factors to explore in the model
model.max_nb_factors = 4; % n=10 so 100 data are used to estimate subjects' random effects (1 factor). 16 data can be used to estimate 2 more factors - we add 1 so that we can also look at the effect of non-significant factors.
model.solid_factors = {'(1|subjectID)'}; %keep these between {}
% a list of possible factors to be included, that can be removed if needed, and the interactions terms to explore
model.liquid_factors = {'eccentricity','omn','fractionOccluded','fractionOccluded:omn','eccentricity:omn','fractionOccluded:eccentricity'}; %keep these between {}
% a list of potential model links
model.links = {'identity'};
% outliers/subjects to be removed - can be left empty
model.exclude = []; 
% no warnings if 1 - careful with that option
model.warning_off = 0;
% whether to use a GLM (0) or a GLME (1).
model.glme = 1;

% run the model
mdls = all_glm(model);

% display diagnostics and results
display_model(mdls{1}, model) %plot model ranked 1 - you can select any other models by rank according to the results on the various indicators provided
saveas(gcf,fullfile(figure_path,'pursuit_dprime_GLME_diagnostics.fig')); snapnow; close;
figure('Color', 'w', 'units','normalized','outerposition',[0.1 0.1 0.8 0.6]);
h=subplot(1,3,1);       plot_covariate_effect(dataPursuit.dprime, dataPursuit.fractionOccluded.*100, h, 'fraction occluded (%)', 'd prime', 0, 0, mdls{1},1,model); 
h=subplot(1,3,2);       plot_covariate_effect(dataPursuit.dprime, dataPursuit.eccentricity, h, 'eccentricity (deg)', 'd prime', 0, 0, mdls{1},1,model);       title('Pursuit condition only');
h=subplot(1,3,3);       plot_interaction(dataPursuit.dprime, (dataPursuit.fractionOccluded<=median(dataPursuit.fractionOccluded)),dataPursuit.eccentricity, h, 'eccentricity','d prime', {'fraction occluded > median','fraction occluded <= median'},mdls{1}, 1, model); % median is 19.7%
saveas(gcf,fullfile(figure_path,'pursuit_dprime_GLME_effects.fig')); snapnow; close;
disp('The first model is undeniably the best fit for the data with best AIC but also best adj. R2 and best R2.')
disp('Interpretation: there is mostly a large effect of fraction occluded. There is also a small effect of eccentricity and a small interaction effect between them.')
disp('Larger eccentricities and fraction occluded are associated with worse d'' results, and interacting, with the effect of fraction occluded decreasing with larger eccentricities.')
disp('When running the first model that includes OMN to see the effect of oculo-motor noise, we find that its effect is dubious and non-significant.')

%% Same analysis for all conditions taken together
% * Additionnal analysis not included in the report.

% explore which distribution is correct
check_distrib_indep(data.dprime(data.status=='control'),data.dprime(data.status=='patient'),'d''');
disp('Data are not distributed differently than normal.')
saveas(gcf,fullfile(figure_path,'dprime_distribution.fig')); snapnow;

% Only apply the collinearity analysis on continuous predictors
% Use Kendall correlations that are more powerful when using small samples
disp('Collinearity')
disp('Lines and columns are ordered as eccentricity,fractionOccluded and OMN.')
[Tau,PValue] = corrplot2(data(:,{'eccentricity','fractionOccluded','omn'}),'type','Kendall')
disp('No particular worries for collinearity given all correlations are well <0.8.')
saveas(gcf,fullfile(figure_path,'factor_collinearity.fig')); snapnow; close;

% define a model structure
% a table with the data, here called data
model.data = data;
% the name of the dependent variable in the data structure, here it is data.initial_work_mem
model.dv = 'dprime';
% its distribution among poisson, normal, gamma, inverse gaussian, binomial as previously determined
model.distribution = 'normal';
% the maximal nb of factors to explore in the model
model.max_nb_factors = 4; % n=10 so 100 data are used to estimate subjects' random effects (1 factor). 16 data can be used to estimate 2 more factors - we add 1 so that we can also look at the effect of non-significant factors.
model.solid_factors = {'(1|subjectID)'}; %keep these between {}
% a list of possible factors to be included, that can be removed if needed, and the interactions terms to explore
model.liquid_factors = {'eccentricity','omn','fractionOccluded','fractionOccluded:omn','eccentricity:omn','fractionOccluded:eccentricity'}; %keep these between {}
% a list of potential model links
model.links = {'identity'};
% outliers/subjects to be removed - can be left empty
model.exclude = [];
% no warnings if 1 - careful with that option
model.warning_off = 0;
% whether to use a GLM (0) or a GLME (1).
model.glme = 1;

% run the model
mdls = all_glm(model);

% display diagnostics and results
display_model(mdls{1}, model) %plot model ranked 1 - you can select any other models by rank according to the results on the various indicators provided
saveas(gcf,fullfile(figure_path,'dprime_GLME_diagnostics.fig')); snapnow; close;
figure('Color', 'w', 'units','normalized','outerposition',[0.1 0.1 0.8 0.6]);
h=subplot(1,3,1);       plot_covariate_effect(data.dprime, data.fractionOccluded.*100, h, 'fraction occluded (%)', 'd prime', 0, 0, mdls{1},1,model);
h=subplot(1,3,2);       plot_covariate_effect(data.dprime, data.eccentricity, h, 'eccentricity (deg)', 'd prime', 0, 0, mdls{1},1,model);
h=subplot(1,3,3);       plot_interaction(data.dprime, (data.fractionOccluded<=median(data.fractionOccluded)),data.eccentricity, h, 'eccentricity','d prime', {'fraction occluded > median','fraction occluded <= median'},mdls{1}, 1, model); % median is 12.8%
disp('The first model is undeniably the best fit for the data with best AIC but also best adj. R2 and best R2.')
disp('Interpretation: there is mostly a medium effect of fraction occluded. There is also a small effect of eccentricity and a small interaction effect between them.')
disp('Larger eccentricities and fraction occluded are associated with worse d'' results, and interacting, with the effect of fraction occluded decreasing with larger eccentricities.')
disp('When running the first model that includes OMN to see the effect of oculo-motor noise, we find that its effect is dubious and non-significant.')
saveas(gcf,fullfile(figure_path,'dprime_GLME_effects.fig')); snapnow; close;

%% debugging
catch err
   keyboard 
end
end

