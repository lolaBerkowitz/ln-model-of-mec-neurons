function results = unpack_glm_res(resdir,datadir,save_plots)
% pull file info from results directory
file_info = struct2table(dir(fullfile([resdir,filesep,'*.mat'])));
files = file_info.name;

% pull data, tetrode, cell from res name.
proc_data_name = extractBetween(files,'results_','_T');
tetrode = extractBetween(files,'TT','_');
cell = extractBetween(extractAfter(files,'TT'),'_','.');

% initialize the results
results = table;
results.file_name{1} = zeros(0);
results.tet{1} = zeros(0);
results.cell{1} = zeros(0);
results.best_model{1} = zeros(0);
results.best_model_name{1} = zeros(0);

% loop through each results file and pull data
for i = 1:length(files)
    % load glm_res
    load([file_info.folder{i},filesep,files{i}])
    data = load([datadir,filesep,proc_data_name{i},'.mat'],'frames','spikesID','maze_size_cm','samplerate','events','Spikes');
    % for cells that didn't meet classification
    if ~isstruct(glm_res)
        results.best_model(i) = {NaN};
        results.best_model_name{i} = NaN;
        results.file_name{i} = [proc_data_name{i},'.mat'];
        results.tet{i} = ['TT',tetrode{i},'.mat'];
        results.cell{i} = cell{i};
        continue
    end
    
    if isnan(glm_res.best_model)
        results.best_model(i) = {glm_res.best_model};
        results.best_model_name{i} = NaN;
        results.file_name{i} = [proc_data_name{i},'.mat'];
        results.tet{i} = ['TT',tetrode{i},'.mat'];
        results.cell{i} = cell{i};
        continue
    end
    
    % grab selected model
    results.best_model(i) = {glm_res.best_model};
    results.best_model_name{i} = glm_res.model_type{glm_res.best_model};
    results.file_name{i} = [proc_data_name{i},'.mat'];
    results.tet{i} = ['TT',tetrode{i},'.mat'];
    results.cell{i} = cell{i};
    
    % Create plot for cell and save
    %     plot_hdVel_glm_res(data,1,glm_res,str2num(cell{i}))
    %     saveas(gca,[save_plots,filesep,proc_data_name{i},'_','tetrode_',num2str(tetrode{i}),'_cell_',num2str(cell{i}),'.png'],'png')
    %     close all
    %
end

writetable(results,[save_path,'results.csv'])
end

