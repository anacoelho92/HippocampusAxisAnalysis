% some options
orden = 2; % order for Bezier function
decimate = 5; % decimation in Bezier function
optim = 1; % 1: use matlab's internal fmninunc, 0: use lbfsg in cluster 
hemi='rh'; % lh: left hemisphere, rh: right hemisphere
structure='Hippo'; % name of structure to compute centerline

dir='/Users/anacoelho/Documents/Axis_analysis/hippocampus_M1/';

clear Q master
clear subjID master
clear subj_list

% Open file with subjects ids
flogs = fopen([dir 'subjects.txt']);
logno = textscan(flogs,'%s');
fclose(flogs);

% Create array with paths to subjects folders

for index_log=1:size(logno{1},1)
    %subj = cell2mat(logno{1}(index_log));
    %subjID(index_log,:) = subj;
    %subj = [dir logno{1}(index_log)];
    subjID(index_log,:) = logno{1}(index_log);
    subj = strcat(dir,logno{1}(index_log));
    Q(index_log,:) = subj;
    
end

% Compute centerline for each subject
for master = 1:numel(Q(:,1))

    fprintf(['*** Working on subject ', subjID{master,:}, ' ... ***\n']);
    
    % load volume
    M = MRIread([deblank(Q{master,:}),'/mri/',hemi,'.aseg',structure,'.mgz'], 0);
    %M = MRIread([deblank(Q{master,:}),'/hippocampus_',hemi,'.nii.gz'], 0);
    
    % compute centerline
    compute_centerline(dir,M,structure,orden,decimate,optim,subjID{master,:},hemi);
    
end