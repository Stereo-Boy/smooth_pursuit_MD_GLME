%% Smooth pursuit in macular degeneration
% * This functions carries on the model statistical analysis for that project.
% * Eye movement data, visibility and d' have been preprocessed and gathered in the master datafile.
% * runs the GLM/GLME with the chopin toolbox.
% * plots various graphs
% * displays statistic conclusion messages
% * saves figures from statistical analysis in figures folder

%% Inputs
% None. The function must run from the analysis folder with the master datafile being located in data. 
% Functions used in this code should all be located in shared_codes.

%% Outputs
% None. But generates figures (folder figures), statistics (command window) and generates a joint markdown output in html 
% by pressing Publish.

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

% Normalize variable formats
data.subject_ID = categorical(data.subject_ID);
data.status = categorical(data.status);
data = clean_names(data);

% show a sample of the data
disp('Data sample and size:')
disp(data(1:10,:))
dispi(size(data,1),' lines of data found.')

%% Primary hypotheses on d'
% Explore how the data can be explained by oculomotor noise, fraction occluded and eccentricity
% Hypotheses: 
% * more eccentric stimuli result in worse d'
% * more fraction occluded for the stimuli result in worse d'
% * more oculomotor noise result in worse d'

% explore which distribution is correct
check_distrib_indep(data.dprime(data.status=='control'),data.dprime(data.status=='patient'),'d''');
disp('Data are not distributed differently than normal')
saveas(gcf,fullfile(figure_path,'dprime_distribution.fig')); snapnow;

% Only apply the collinearity analysis on continuous predictors
% Use Kendall correlations that are more powerful when using small samples
disp('Collinearity')
disp('Lines and columns are ordered as eccentricity,fractionOccluded and ocm.')
[Tau,PValue] = corrplot2(data(:,{'eccentricity','fractionOccluded','ocm'}),'type','Kendall');
disp('No particular worries for collinearity given all correlations are <0.8 but keep in mind the significant correlation between eccentricity and ocm.')
saveas(gcf,fullfile(figure_path,'factor_collinearity.fig')); snapnow;

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
model.liquid_factors = {'eccentricity','ocm','fractionOccluded','fractionOccluded:ocm','eccentricity:ocm','fractionOccluded:eccentricity'}; %keep these between {}
% a list of potential model links
model.links = {'identity'};
% outliers/subjects to be removed - can be left empty
model.exclude = [98,107,110]; % exclude the 3 instances for which participants did not follow pursuit instructions
% no warnings if 1 - careful with that option
model.warning_off = 0;
% whether to use a GLM (0) or a GLME (1).
model.glme = 1;

% run the model
mdls = all_glm(model);

% display diagnostics and results
display_model(mdls{1}, model) %plot model ranked 1 - you can select any other models by rank according to the results on the various indicators provided
saveas(gcf,fullfile(figure_path,'dprime_GLME_diagnostics.fig')); snapnow;
h=subplot(1,3,1);       plot_covariate_effect(data.dprime, data.fractionOccluded.*100, h, 'fraction occluded (%)', 'd prime', 0, 0, mdls{1},1,model);
h=subplot(1,3,2);       plot_covariate_effect(data.dprime, data.eccentricity, h, 'eccentricity (deg)', 'd prime', 0, 0, mdls{1},1,model);
h=subplot(1,3,3); plot_interaction(data.dprime, (data.eccentricity<=median(data.eccentricity)),data.fractionOccluded.*100, h, 'fraction occluded (%)','d prime', {'eccentricity > median','eccentricity <= median'},mdls{1}, 1, model); % median is 2.36
disp('The first model is undeniably the best fit for the data with best AIC but also best adj. R2 and best R2.')
disp('Interpretation: there is mostly a medium effect of fraction occluded. There is also a small effect of eccentricity and a small interaction effect between them.')
disp('Larger eccentricities and fraction occluded are associated with worse d'' results, but for higher fraction occluded, eccentricity effects decrease.')
disp('In lay terms, when it is occluded, it does not matter that it is blurry.')
disp('When running model #3 to see the effect of oculo-motor noise, its effect is small and non-significant.')
saveas(gcf,fullfile(figure_path,'dprime_GLME_effects.fig')); snapnow;

%% debugging
catch err
   keyboard 
end
end

