
function glm_res = run_me(data,session,celln,varargin)
%%% Description of run_me

% This script is segmented into several parts. 
% 1. data is loaded. 
% 2. 7 LN models are fit to the cell's spike train. 
% Each model uses information (below) to predict a section of the
% spike train: 
%   - position
%   - head direction
%   - running speed 
%   - egocentric bearing
% 3. Model fitting and model performance is computed through
%    10-fold cross-validation, and the minimization procedure is carried out
%    through fminunc. 
% 4. forward-search procedure is implemented to find the simplest 'best' model 
%    describing this spike train. 
% 
% Input: 
%  - data: processed data file from ephys_tools. 
%  - session: session number dictates which event files are used to index
%  spikes/timestamps.
%  - celln: cell number used as index for data.Spikes. 
% Output:
%  - glm_res: structure that contains train/test data sets, parameters
%  for each LN model, and the selected model.  
%  - smooth_fr: Firing rate smoothed with guassian filter for plotting
%  purposes. 
%

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017
% Adapted for ephys_tools RH 2019 
% Added EBC predictors LB 02/2020
% Simplified for Place, HD, Ego & made into function LB 09/2020 

p = inputParser;
addOptional(p,'numFolds',10,@isnumeric)
addOptional(p,'n_pos_bins',20,@isnumeric)
addOptional(p,'n_circ_bins',18,@isnumeric)
addOptional(p,'n_speed_bins',10,@isnumeric)
parse(p,varargin{:})

% Initialize model hyperparameters
model_params = struct;
model_params.numFolds = p.Results.numFolds;
model_params.pos_bins = p.Results.n_pos_bins;
model_params.hd_bins = p.Results.n_circ_bins;
model_params.ego_bins = p.Results.n_circ_bins;
model_params.speed_bins = p.Results.n_speed_bins;
model_params.numModels = 7;

% boxSize = length (in cm) of one side of the square box
boxSize=data.maze_size_cm(session);

frames=data.frames(data.frames(:,1)>data.events(1,session) &...
    data.frames(:,1)<data.events(2,session),:); 

% convert xy to cm 
frames(:,2)=rescale(frames(:,2),0,boxSize);
frames(:,3)=rescale(frames(:,3),0,boxSize);

% post = vector of time (seconds) at every 33 ms time bin
post=frames(:,1);

% spiketrain = vector of the # of spikes in each 33 ms time bin
spiketrain=histcounts(data.Spikes{celln}(data.Spikes{celln}>data.events(1,session) &...
    data.Spikes{celln}<data.events(2,session)),...
    linspace(post(1),post(end),length(post)+1))';

% posx_c = x-position in middle of LEDs
posx_c = frames(:,2);

% posy_c = y-position in middle of LEDs
posy_c = frames(:,3);
pos = [posx_c posy_c];

% sampleRate = sampling rate of neural data
sampleRate = data.samplerate;

% head direction
hd = deg2rad(frames(:,4));

%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')

[glm_res,smooth_fr] = fit_all_ln_models(pos,hd,spiketrain,post,model_params,boxSize);


%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
selected_model = select_best_model(glm_res.test,model_params);

glm_res.best_model = selected_model; 
glm_res.smooth_fr = smooth_fr;

end


%% Local Functions 

function [glm_res,smooth_fr] = fit_all_ln_models(pos,hd,spiketrain,post,model_params,boxSize)
%% Description
% The model: r = exp(W*theta), where r is the predicted # of spikes, W is a
% matrix of one-hot vectors describing variable (P, H, E, or S) values, and
% theta is the learned vector of parameters.
glm_res = struct;
glm_res.model_type = {'Place-HD-Ego';'Place-HD';'Place-Ego'; 'HD-Ego';'Place';'HD';'Ego'};

%% Compute the position, head direction, speed, and egocentric bearing matrices

% initialize the number of bins that position, head direction, speed, and
% egocentric bearing will be divided into
n_pos_bins = model_params.pos_bins;
n_dir_bins = model_params.hd_bins;
n_speed_bins = model_params.speed_bins;
n_ego_bins = model_params.ego_bins;

% compute position matrix
[posgrid, ~] = pos_map(pos, n_pos_bins, boxSize);

% compute head direction matrix
[hdgrid,~,~] = hd_map(hd,n_dir_bins);

% compute speed matrix
[~,~,speed] = speed_map(pos(:,1),pos(:,2),n_speed_bins);

% compute EBC matrices (bearing and distance)
[egogrid, ~, ~, ~, ~, ~,~] = ebc_map(pos,hd, n_ego_bins);

% remove times when the animal ran > 100 cm/s (these data points may contain artifacts)
too_fast = find(speed >= 100);
posgrid(too_fast,:) = []; hdgrid(too_fast,:) = []; 
spiketrain(too_fast) = []; 
egogrid(too_fast,:) = [];


%% Fit all 7 LN models
% Initialize saved variables
testFit = cell(model_params.numModels,1);
trainFit = cell(model_params.numModels,1);
param = cell(model_params.numModels,1);
A = cell(model_params.numModels,1);
modelType = cell(model_params.numModels,1);

% ALL VARIABLES
A{1} = [ posgrid hdgrid egogrid]; modelType{1} = [1 1 1];

% TWO VARIBALE 
A{2} = [ posgrid hdgrid]; modelType{2} = [1 1 0];
A{3} = [ posgrid egogrid]; modelType{3} = [1 0 1];
A{4} = [ hdgrid egogrid]; modelType{4} = [0 1 1];

% ONE VARIABLE 
A{5} = posgrid; modelType{5} = [1 0 0];
A{6} = hdgrid; modelType{6} = [0 1 0];
A{7} = egogrid; modelType{7} = [0 0 1];

% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
dt = mean(diff(post)); fr = spiketrain/dt;
smooth_fr = conv(fr,filter,'same');
num_models = model_params.numModels;
num_folds = model_params.numFolds;

parfor n = 1:num_models
    fprintf('\t- Fitting model %d of %d\n', n, num_models);
    [testFit{n},trainFit{n},param{n}] = fit_model(A{n},dt,spiketrain,filter,modelType{n},num_folds);
end

% save model fits 
glm_res.train = cell2mat(trainFit);
glm_res.test = cell2mat(testFit);
glm_res.param = param;

end

function [testFit,trainFit,param_mean] = fit_model(A,dt,spiketrain,filter,modelType,numFolds)

%% Description
% This code will section the data into 10 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 9 sections, and test the model performance on the
% remaining section. This procedure will be repeated 10 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.


%% Initialize matrices and section the data for k-fold cross-validation

[~,numCol] = size(A);
sections = numFolds*5;

% divide the data up into 5*num_folds pieces
edges = round(linspace(1,numel(spiketrain)+1,sections+1));

% initialize matrices
testFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of test data
trainFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of train data
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
for k = 1:numFolds
    fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
    
    % get test data from edges - each test data chunk comes from entire session
    test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1]   ;
    
    test_spikes = spiketrain(test_ind); %test spiking
    smooth_spikes_test = conv(test_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_test = smooth_spikes_test./dt;
    test_A = A(test_ind,:);
    
    % training data
    train_ind = setdiff(1:numel(spiketrain),test_ind);
    train_spikes = spiketrain(train_ind);
    smooth_spikes_train = conv(train_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_train = smooth_spikes_train./dt;
    train_A = A(train_ind,:);
    
    opts = optimset('Gradobj','on','Display','off'); % removed ,'Hessian','on'
    
    data{1} = train_A; data{2} = train_spikes;
    if k == 1
        init_param = 1e-3*randn(numCol, 1);
    else
        init_param = param;
    end
    [param] = fminunc(@(param) ln_poisson_model(param,data,modelType),init_param,opts);
    
    %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_test = exp(test_A * param)/dt;
    smooth_fr_hat_test = conv(fr_hat_test,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
    sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % compute correlation
    correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson');
    
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
    r = exp(test_A * param); n = test_spikes; meanFR_test = nanmean(test_spikes); 
    
    log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
    log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
    log_llh_test = (-log_llh_test_model + log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test;
    
    % compute MSE
    mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
    
    % fill in all the relevant values for the test fit cases
    testFit(k,:) = [varExplain_test correlation_test log_llh_test mse_test sum(n) numel(test_ind) ];
    
    %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_train = exp(train_A * param)/dt;
    smooth_fr_hat_train = conv(fr_hat_train,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_train-smooth_fr_train).^2);
    sst = sum((smooth_fr_train-mean(smooth_fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % compute correlation
    correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson');
    
    % compute log-likelihood
    r_train = exp(train_A * param); n_train = train_spikes; meanFR_train = nanmean(train_spikes);   
    log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train;
    
    % compute MSE
    mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
    
    trainFit(k,:) = [varExplain_train correlation_train log_llh_train mse_train sum(n_train) numel(train_ind)];

    % save the parameters
    paramMat(k,:) = param;

end

param_mean = nanmean(paramMat);

end

function selected_model = select_best_model(testFit,model_params)
testFit_mat = testFit;
LLH_values = reshape(testFit_mat(:,3),model_params.numFolds,model_params.numModels);

% find the best single model
singleModels = 5:7;
[~,top1] = max(nanmean(LLH_values(:,singleModels))); top1 = top1 + singleModels(1)-1;

% find the best double model that includes the single model
if top1 == 5 % P-> PH, PE
    [~,top2] = max(nanmean(LLH_values(:,[2 3])));
    vec = [2 3]; top2 = vec(top2);
elseif top1 == 6 % H -> HP, HE
    [~,top2] = max(nanmean(LLH_values(:,[2 4])));
    vec = [2 4]; top2 = vec(top2);
elseif top1 == 7 % E -> PE, HE
    [~,top2] = max(nanmean(LLH_values(:,[3 4])));
    vec = [3 4]; top2 = vec(top2);
end

% full model PHE
top3 = 1; 

% Grab logliklihood values for top models
LLH1 = LLH_values(:,top1); LLH2 = LLH_values(:,top2);LLH3 = LLH_values(:,top3);

% compare increasingly complicated models 
[p_llh_12,~] = signrank(LLH2,LLH1,'tail','right');
[p_llh_23,~] = signrank(LLH3,LLH2,'tail','right');


if p_llh_12 < 0.05 % double model is sig. better
    if p_llh_23 < 0.05  % triple model is sig. better
        selected_model = top3; %full model
    else
        selected_model = top2; %double model
    end
else
    selected_model = top1; %single model
end

% re-set if selected model is not above baseline
pval_baseline = signrank(LLH_values(:,selected_model),[],'tail','right');

if pval_baseline > 0.05
    selected_model = NaN;
end
end
