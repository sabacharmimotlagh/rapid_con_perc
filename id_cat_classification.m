function [results] = id_cat_classification(num_of_participant, varargin)

%% parse inputs
p = inputParser;
addParameter(p, 'results_filename', 'results_cat.mat', @ischar);
addParameter(p, 'trial_directory', pwd, @ischar);
addParameter(p, 'number_of_repetitions', 100, @isnumeric);
addParameter(p, 'number_of_permutations', 100, @isnumeric);
parse(p, varargin{:});

results_filename = p.Results.results_filename;
trial_directory = p.Results.trial_directory;
num_of_rep = p.Results.number_of_repetitions;
num_of_permut = p.Results.number_of_permutations;

% Check if results file exists and load it
if isfile(results_filename)
    load(results_filename, 'results');
else
    results = struct('id_cat_timeseries', {{}});
end

%% loading epoched data for specific participant excluding bad trials
% navigate to trials directory
d = dir(trial_directory);
cd(trial_directory)

% this has bad trials' names for each subject
load(d(3).name);


counter_both_c = 0;
counter_category_c = 0;
counter_both_i = 0;

% loading trials
% data is 3D (variables(channels), timepoints, observations(trials))

for f=5:length(d)
    name = regexp(d(f).name,'_','split');
    if ~any(strcmp(BadTrials,d(f).name))
        if(strcmp([name{1,2},'_',name{1,3}],'Both_correct'))
            counter_both_c = counter_both_c + 1;
            load(d(f).name)
            data_both_c(:,:,counter_both_c) = F(ChannelFlag(1:64) == 1,:);
            data_both_c_label{counter_both_c,1} = 'Both_correct';
        elseif(strcmp([name{1,2},'_',name{1,3}],'Category_correct'))
            counter_category_c = counter_category_c + 1;
            load(d(f).name)
            data_category_c(:,:,counter_category_c) = F(ChannelFlag(1:64) == 1,:);
            data_category_c_label{counter_category_c,1} = 'Category_correct';
        elseif(strcmp([name{1,2},'_',name{1,3}],'Both_incorrect'))
            counter_both_i = counter_both_i + 1;
            load(d(f).name)
            data_both_i(:,:,counter_both_i) = F(ChannelFlag(1:64) == 1,:);
            data_both_i_label{counter_both_i,1} = 'Both_incorrect';
        end
    end
end

%% bootstraping on data to overcome imbalance in data
num_of_samples = min([size(data_both_c,3),size(data_category_c,3),size(data_both_i,3)]);

for rep=1:num_of_rep
    fprintf('Repetition number: %d \r', rep);
    % finding the minority class to subsample from other two classes
    if size(data_both_c,3) == num_of_samples
        n = randperm(size(data_category_c,3),num_of_samples);   
        data_category_c_s = data_category_c(:,:,n);   
        data_category_c_label_s = data_category_c_label(n,1); 
        n2 = randperm(size(data_both_i,3),num_of_samples);   
        data_both_i_s = data_both_i(:,:,n2);   
        data_both_i_label_s = data_both_i_label(n2,1); 
        Data = cat(data_both_c,3,data_category_c_s,data_both_i_s);
        Labels = [data_both_c_label; data_category_c_label_s; data_both_i_label_s];
    elseif size(data_category_c,3) == num_of_samples
        n = randperm(size(data_both_c,3),num_of_samples);   
        data_both_c_s = data_both_c(:,:,n);   
        data_both_c_label_s = data_both_c_label(n,1); 
        n2 = randperm(size(data_both_i,3),num_of_samples);   
        data_both_i_s = data_both_i(:,:,n2);   
        data_both_i_label_s = data_both_i_label(n2,1); 
        Data = cat(3,data_both_c_s,data_category_c,data_both_i_s);
        Labels = [data_both_c_label_s; data_category_c_label; data_both_i_label_s];
    else
        n = randperm(size(data_both_c,3),num_of_samples);   
        data_both_c_s = data_both_c(:,:,n);   
        data_both_c_label_s = data_both_c_label(n,1); 
        n2 = randperm(size(data_category_c,3),num_of_samples);   
        data_category_c_s = data_category_c(:,:,n2);   
        data_category_c_label_s = data_category_c_label(n2,1); 
        Data = cat(3,data_both_c_s,data_category_c_s,data_both_i);
        Labels = [data_both_c_label_s; data_category_c_label_s; data_both_i_label];
    end
   
    d_singleimage{rep} = fl_decodesvm(Data,Labels','numpermutation',num_of_permut,'verbose',10,'kfold',10);
    decoding(rep,:,:) = d_singleimage{1,rep}.d;
end

% average the decodings
avg_decoding = squeeze(mean(decoding,1));

%% save results
results.id_cat_timeseries{num_of_participant} = avg_decoding;
save(results_filename,'results');

end