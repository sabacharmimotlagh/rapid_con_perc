function [results] = pairwise_classification(num_of_participant, varargin)

%% parse inputs
p = inputParser;
addParameter(p, 'results_filename', 'results.mat', @ischar);
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
    results = struct('recog_unrecog_timeseries', {{}}, 'excluded_stimuli', {{}});
end

%% loading epoched data for specific participant excluding bad trials
% navigate to trials directory
d = dir(trial_directory);
cd(trial_directory)

% this has bad trials' names for each subject
load(d(3).name);
    
% loading recognized and unrecognized trials for each stimulus
% data is 3D (variables(channels) x timepoints x observations(trials))
data = cell(1,80);
label = cell(1,80);

for f=5:length(d)
    name = regexp(d(f).name,'_','split');
    if ~any(strcmp(BadTrials,d(f).name))    % if the trial was not a bad trial, add the trial to the data
        load(d(f).name);
        data{1,str2double(name{1,2})} = cat(3,data{1,str2double(name{1,2})},reshape(F(ChannelFlag(1:64) == 1,:),[size(F(ChannelFlag(1:64) == 1,:)), 1]));
        label{1,str2double(name{1,2})} = cat(2, label{1,str2double(name{1,2})}, [str2double(name{1,2})]);
    end
end

%% finding conditions with less than 5 trials so we can remove those conditions
[~, ~, nTrials] = cellfun(@size, data); % finding number of trials for each condition
excluded = [];

% setting the threshold for deleting the targets
for rep=1:(length(nTrials)/2)
    if min(nTrials(rep), nTrials(rep+40))<5
        excluded = [excluded, rep, rep+40];
    end
end

% removing excluded stimuli from data cell
data(excluded) = [];
label(excluded) = [];

results.excluded_stimuli{num_of_participant} = excluded;

%% perform bootstraping on data to overcome imbalance in data
% the minimum number of trials among all conditions is the number of samples to ensure consistent signal to noise ratio
num_of_samples =  min(cellfun(@(c) size(c,3), data));

for rep=1:num_of_rep

    Data = [];
    Labels = [];
    
    for i=1:length(data)
        % randomly choose the number of samples from trials of all conditions and store the results in Data and Labels
        n = randperm(size(data{1,i}, 3), num_of_samples); 
        Data = cat(3,Data,data{1,i}(:, :, n));
        Labels = cat(2,Labels,label{1,i}(:, n));
    end
    
    Labels = sprintfc('%02d',Labels);

    fprintf('Repetition number %d \r', rep);
    d_singleimage{rep} = fl_decodesvm(Data,Labels,'numpermutation',num_of_permut,'verbose',10,'kfold',5);
    decoding(rep,:,:) = d_singleimage{rep}.d;
end

avg_decoding = squeeze(mean(decoding,1));

%% save the single subjects result in the average matrix
results.recog_unrecog_timeseries{num_of_participant} = avg_decoding;
save('results.mat','results');
