function plot_hdVel_glm_res(data,session,glm_res,varargin)

p = inputParser;
addOptional(p,'n_pos_bins',20,@isnumeric)
addOptional(p,'n_circ_bins',18,@isnumeric)
addOptional(p,'n_speed_bins',10,@isnumeric)
addOptional(p,'n_angvel_bins',20,@isnumeric)
parse(p,varargin{:})


% Initialize model hyperparameters
n_pos_bins = p.Results.n_pos_bins;
n_dir_bins = p.Results.n_circ_bins;
n_speed_bins = p.Results.n_speed_bins;
n_angvel_bins = p.Results.n_angvel_bins;

% boxSize = length (in cm) of one side of the square box
boxSize=data.maze_size_cm(session);
frames=data.frames(data.frames(:,1)>data.events(1,session) &...
    data.frames(:,1)<data.events(2,session),:); 

% convert xy to cm 
frames(:,2)=rescale(frames(:,2),0,boxSize);
frames(:,3)=rescale(frames(:,3),0,boxSize);

% Grab position and head direction 
posx_c = frames(:,2);
posy_c = frames(:,3);
direction = deg2rad(frames(:,4));
speed = frames(:,5);
angvel = insta_angvel(direction,data.samplerate,3);

% take out times when the animal ran >= 50 cm/s
too_fast = find(speed >= 100);
posx_c(too_fast) = []; posy_c(too_fast) = []; 
direction(too_fast) = [];
speed(too_fast) = []; angvel(too_fast) = [];

% compute tuning curves for position, head direction, speed, and theta phase
[pos_curve] = compute_2d_tuning_curve(posx_c,posy_c,glm_res.smooth_fr,n_pos_bins,0,boxSize);
[hd_curve] = compute_1d_tuning_curve(direction,glm_res.smooth_fr,n_dir_bins,0,2*pi);
[speed_curve] = compute_1d_tuning_curve(speed,glm_res.smooth_fr,n_speed_bins,0,50);
[angvel_curve] = compute_1d_tuning_curve(angvel,glm_res.smooth_fr,n_angvel_bins,0,1000);

% create x-axis vectors
hd_vector = 2*pi/n_dir_bins/2:2*pi/n_dir_bins:2*pi - 2*pi/n_dir_bins/2;
speed_vector = 2.5:50/n_speed_bins:47.5;
angvel_vector = 0:1000/n_angvel_bins:1000-20;

% Get predicted tuning curves from models
[pos_response,hd_response,speed_response,angvel_response] = response_profile(glm_res.param);


%% plot the tuning curves
fig=figure;fig.Color=[1 1 1];
% Position
subplot(3,4,1)
imagesc(pos_curve); colorbar
ylabel('Tuning Curves')
title('Position')
% Head Direction
subplot(3,4,2)
plot(hd_vector,hd_curve,'k','linewidth',3)
box off
axis([0 2*pi -inf inf])
xlabel('direction angle')
title('Head direction')
% Linear Velocity
subplot(3,4,3)
plot(speed_vector,speed_curve,'k','linewidth',3)
box off
xlabel('Running speed')
axis([0 50 -inf inf])
title('Linear Velocity')
% Angular Velocity
subplot(3,4,4)
plot(angvel_vector,angvel_curve,'k','linewidth',3)
xlabel('Angular Velocity')
axis([0 1000 -inf inf])
box off
title('Angular Velcoity')

% plot the model-derived response profiles
% Predicted Position
subplot(3,4,5)
imagesc(reshape(pos_response,20,20)); 
colorbar
ylabel('LN Model Profiles')
% Predicted Heading
subplot(3,4,6)
plot(hd_vector,hd_response,'k','linewidth',3)
xlabel('direction angle')
axis([0 2*pi -inf inf])
box off
% Predicted Linear Velocity
subplot(3,4,7)
plot(speed_vector,speed_response,'k','linewidth',3)
xlabel('Linear Velocity')
axis([0 50 -inf inf])
box off
% Predicted Angular Velocity
subplot(3,4,8)
plot(angvel_vector,angvel_response,'k','linewidth',3)
xlabel('Angular Velocity')
axis([0 1000 -inf inf])
box off

%% make the axes match
subplot(3,4,1)
caxis([min(min(pos_response),min(pos_curve(:))) max(max(pos_response),max(pos_curve(:)))])
subplot(3,4,5)
caxis([min(min(pos_response),min(pos_curve(:))) max(max(pos_response),max(pos_curve(:)))])

subplot(3,4,2)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])
subplot(3,4,6)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])

subplot(3,4,3)
axis([0 50 min(min(speed_response),min(speed_curve)) max(max(speed_response),max(speed_curve))])
subplot(3,4,7)
axis([0 50 min(min(speed_response),min(speed_curve)) max(max(speed_response),max(speed_curve))])

subplot(3,4,4)
axis([0 1000 min(min(angvel_response),min(angvel_curve)) max(max(angvel_response),max(angvel_curve))])
subplot(3,4,8)
axis([0 1000 min(min(angvel_response),min(angvel_curve)) max(max(angvel_response),max(angvel_curve))])

% Loglikelihood values
LLH_values = reshape(glm_res.test(:,3),10,15);
LLH_increase_mean = mean(LLH_values);
LLH_increase_sem = std(LLH_values)/sqrt(10);

subplot(3,4,9:12)
errorbar(LLH_increase_mean,LLH_increase_sem,'ok','linewidth',3)
hold on
plot(glm_res.best_model,LLH_increase_mean(glm_res.best_model),'.r','markersize',25)
plot(0.5:15.5,zeros(16,1),'--b','linewidth',2)
hold off
box off
grid on
set(gca,'fontsize',12)
set(gca,'XLim',[0 16]); set(gca,'XTick',1:15)
set(gca,'XTickLabel',glm_res.model_type);
legend('Model performance','Selected model','Baseline')
ylabel('Log-Likelihood (bits/spike)')

end


function [tuning_curve] = compute_1d_tuning_curve(variable,fr,numBin,minVal,maxVal)

%bin it
var_vec = linspace(minVal,maxVal,numBin+1);
tuning_curve = nan(numBin,1);

% compute mean fr for each bin
for n = 1:numBin
    tuning_curve(n) = nanmean(fr(variable >= var_vec(n) & variable < var_vec(n+1)));
    
    if n == numBin
        tuning_curve(n) = nanmean(fr(variable >= var_vec(n) & variable <= var_vec(n+1)));
    end
    
end


end

function [tuning_curve] = compute_2d_tuning_curve(variable_x,variable_y,fr,numBin,minVal,maxVal)

% this assumes that the 2d environment is a square box, and that the
% variable is recorded along the x- and y-axes

%% define the axes and initialize variables

xAxis = linspace(minVal,maxVal,numBin+1);
yAxis = linspace(minVal,maxVal,numBin+1);

% initialize 
tuning_curve = zeros(numBin,numBin);

%% fill out the tuning curve

% find the mean firing rate in each position bin
for i  = 1:numBin
    start_x = xAxis(i); stop_x = xAxis(i+1);
    % find the times the animal was in the bin
    if i == numBin
        x_ind = find(variable_x >= start_x & variable_x <= stop_x);
    else
        x_ind = find(variable_x >= start_x & variable_x < stop_x);
    end
    for j = 1:numBin
        
        start_y = yAxis(j); stop_y = yAxis(j+1);
        
        if j == numBin
            y_ind = find(variable_y >= start_y & variable_y <= stop_y);
        else
            y_ind = find(variable_y >= start_y & variable_y < stop_y);
        end
        
        ind = intersect(x_ind,y_ind);
        
        % fill in rate map
        tuning_curve(numBin+1 - j,i) = mean(fr(ind));
    end
end


% smooth with Gaussian
H = fspecial('gaussian'); % using default values - size=[3 3] and sigma=0.5
% tuning_curve = imfilter(tuning_curve,H);
tuning_curve = nanconv(tuning_curve,H, 'nanout');

end

function[pos_response,hd_response,speed_response,angvel_response] = response_profile(param)
% pull out the parameter values
% show parameters from the full model
param_full_model = param{1};

n_pos_bins = sqrt(length(param{12}));
n_dir_bins = length(param{13});
n_speed_bins = length(param{14});
n_angvel_bins = length(param{15});

% pull out the parameter values
pos_param = param_full_model(1:n_pos_bins^2);
hd_param = param_full_model(n_pos_bins^2+1:n_pos_bins^2+n_dir_bins);
speed_param = param_full_model(n_pos_bins^2+n_dir_bins+1:n_pos_bins^2+n_dir_bins+n_speed_bins);
angvel_param = param_full_model(numel(param_full_model)-n_angvel_bins+1:numel(param_full_model));


% compute the scale factors
% NOTE: technically, to compute the precise scale factor, the expectation
% of each parameter should be calculated, not the mean.
scale_factor_pos = mean(exp(speed_param))*mean(exp(hd_param))*mean(exp(angvel_param))*30;
scale_factor_hd = mean(exp(speed_param))*mean(exp(pos_param))*mean(exp(angvel_param))*30;
scale_factor_spd = mean(exp(pos_param))*mean(exp(hd_param))*mean(exp(angvel_param))*30;
scale_factor_angvel = mean(exp(speed_param))*mean(exp(hd_param))*mean(exp(pos_param))*30;

% compute the model-derived response profiles
pos_response = scale_factor_pos*exp(pos_param);
hd_response = scale_factor_hd*exp(hd_param);
speed_response = scale_factor_spd*exp(speed_param);
angvel_response = scale_factor_angvel*exp(angvel_param);

end