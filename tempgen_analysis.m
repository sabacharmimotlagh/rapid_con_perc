function [] = tempgen_analysis(num_of_participant, varargin)

%% parse inputs
p = inputParser;
addParameter(p, 'trial_directory', pwd, @ischar);
addParameter(p, 'save_directory', pwd, @ischar);
addParameter(p, 'number_of_repetitions', 100, @isnumeric);
addParameter(p, 'number_of_permutations', 100, @isnumeric);
parse(p, varargin{:});

trial_directory = p.Results.trial_directory;
save_directory = p.Results.save_directory;
num_of_rep = p.Results.number_of_repetitions;
num_of_permut = p.Results.number_of_permutations;

recognized_tempgen_dir = fullfile(save_directory, tempgen_recognized);
unrecognized_tempgen_dir = fullfile(save_directory, tempgen_recognized);

% Check if the directory already exists
if ~exist(recognized_tempgen_dir, 'dir')
    % Create the directory
    mkdir(recognized_tempgen_dir);
end

if ~exist(unrecognized_tempgen_dir, 'dir')
    mkdir(unrecognized_tempgen_dir);
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
        % sownsampling to reduce the computational complexity and cost
        F_ds = downsample(F',2);
        F_downsample = F_ds';
        data{1,str2num(name{1,2})} = cat(3,data{1,str2num(name{1,2})},reshape(F_downsample(ChannelFlag(1:64) == 1,:),[size(F_downsample(ChannelFlag(1:64) == 1,:)), 1]));
        label{1,str2num(name{1,2})} = cat(2, label{1,str2num(name{1,2})}, [str2num(name{1,2})]);
    end
end


% separating recognized and unrecognized cases only for temporal generalization method to reduce computational complexity
data_recog = data(1:40);
label_recog = label(1:40);
    
data_unrecog = data(41:end);
label_unrecog = label(41:end);

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
data_recog(excluded(excluded<=40)) = [];
label_recog(excluded(excluded<=40)) = [];

data_unrecog(excluded(excluded>40)-40) = [];
label_unrecog(excluded(excluded>40)-40) = [];

%% perform bootstraping on data to overcome imbalance in data
num_of_samples_r =  min(cellfun(@(c) size(c,3), data_recog));
num_of_samples_u =  min(cellfun(@(c) size(c,3), data_unrecog));

% for recognized conditions
for rep=1:num_of_rep
   
    filename = sprintf('Temp_gen_r_participant%d_repetition%d', num_of_participant, rep);
    Data = [];
    Labels = [];
    
    for i=1:length(data_recog)
        n = randperm(size(data_recog{1,i},3),num_of_samples_r); 
        Data = cat(3,Data,data_recog{1,i}(:,:,n));
        Labels = cat(2,Labels,label_recog{1,i}(:,n));
    end
    
    Labels = sprintfc('%02d',Labels);

    fprintf('Repetition number %d \r', rep);
    d_singleimage = fl_decodesvm(Data,Labels,'numpermutation',50,'verbose',10,'kfold',5, 'method', 'temporalgen');
    decoding = d_singleimage.d;

    cd(recognized_tempgen_dir);
    save([filename '.mat'],'decoding');
end

% for unrecognized conditions
for rep=1:num_of_rep
    
    filename = sprintf('Temp_gen_u_participant%d_repetition%d', num_of_participant, i);

    Data = [];
    Labels = [];
    
    for i=1:length(data_recog)
        n = randperm(size(data_unrecog{1,i},3),num_of_samples_u); 
        Data = cat(3,Data,data_unrecog{1,i}(:,:,n));
        Labels = cat(2,Labels,label_unrecog{1,i}(:,n));
    end
    
    Labels = sprintfc('%02d',Labels);

    fprintf('repetition number %d \r', i);
    d_singleimage = fl_decodesvm(Data,Labels,'numpermutation',50,'verbose',10,'kfold',5, 'method', 'temporalgen');
    decoding = d_singleimage.d;

    cd(unrecognized_tempgen_dir);
    save([filename '.mat'],'decoding');
end