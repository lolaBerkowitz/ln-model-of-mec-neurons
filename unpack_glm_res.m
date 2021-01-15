function results = unpack_glm_res(datadir)

file_info = struct2table(dir(fullfile([datadir,filesep,'*.mat'])));
files = file_info.name;

results = table(zeros(size(files,1),1),cell(size(files,1),1),'VariableNames',{'best_model','file_name'});
for i = 1:length(files)
    % load glm_res
    load([file_info.folder{i},filesep,files{i}])
    
    % for cells that didn't meet classification
    if isnan(glm_res.best_model)
        results.best_model(i) = glm_res.best_model;
        results.file_name{i} = [file_info.folder{i},filesep,files{i}];
        continue
    end
    
    % grab selected model
    results.best_model(i) = glm_res.best_model;
    results.file_name{i} = [file_info.folder{i},filesep,files{i}];
    
    % Create plot for cell and save
    plot_hdVel_glm_res(data,session,glm_res)
%     
end
end

function response = response_profile(param)
% pull out the parameter values
pos_param = param{5};

% compute the scale factors
% NOTE: technically, to compute the precise scale factor, the expectation
% of each parameter should be calculated, not the mean.
scale_factor = mean(exp(ego_param))*mean(exp(hd_param))*33;

% compute the model-derived response profiles
response = scale_factor_pos*exp(pos_param);


end