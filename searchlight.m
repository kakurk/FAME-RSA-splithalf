%% Searchlight
%
% Load beta image data from one subject, apply 'vt' mask, compute difference
% of (fisher-transformed) between on- and off diagonal split-half
% correlation values.

if isunix % if we are on Hammer, a unix system
    addpath(genpath('/gpfs/group/n/nad12/RSA/Scripts/CoSMoMVPA-master'))
else % if not on unix, assume we are on Anvil
    addpath(genpath('S:\nad12\CoSMoMVPA-master'))
end

% configure the spm_jobmanager
spm_jobman('initcfg')

%% Set analysis parameters
subjects   = {'18y404'  '18y566'  '20y297' '20y415'  '20y439'}; % '20y396' <-- this subjects doesn't have a model run
study_path = '/gpfs/group/n/nad12/RSA/Analysis_ret/FAMEret8RSA_hrf'; % path on hammer

for ss = 1:length(subjects)
 
    % Path to this subjects analysis folder
    data_path  = fullfile(study_path, subjects{ss});
    
    % Edit the SPM.mat file to use paths here on Hammer
    if isunix % only execute if we are on a Unix system (i.e., Hammer)
        spm_changepath(fullfile(data_path, 'SPM.mat'), 'S:\nad12\FAME8', '/gpfs/group/n/nad12/RSA')
        spm_changepath(fullfile(data_path, 'SPM.mat'), '\', '/')
    end

    % Remove Existing Contrasts
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(data_path, 'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess = {};
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch)
    
    % Odd Even Contrast
    oddeven_contrasts(subjects{ss}, study_path, 'HREC', '.*HREC.*');
    oddeven_contrasts(subjects{ss}, study_path, 'HFAM', '.*HFAM.*');
    oddeven_contrasts(subjects{ss}, study_path, 'FAREC', '.*FAREC.*');
    oddeven_contrasts(subjects{ss}, study_path, 'FAFAM', '.*FAFAM.*');
    
    % Gather Odd and Even Runs spmT files into a single nii file
    odd_Ts  = spm_select('FPList', data_path, 'spmT_000[1357]');
    even_Ts = spm_select('FPList', data_path, 'spmT_000[2468]');
    spm_file_merge(odd_Ts, fullfile(data_path, 'glm_T_stats_odd.nii'));
    spm_file_merge(even_Ts, fullfile(data_path, 'glm_T_stats_even.nii'));
   
    % file locations for both halves
    Odd  = fullfile(data_path, 'glm_T_stats_odd.nii');
    Even = fullfile(data_path, 'glm_T_stats_even.nii');

    % Path to this subject's whole brain mask
    mask_fn    = fullfile(data_path, 'mask.img');
    
    % load two halves as CoSMoMVPA dataset structs.
    % Chunks = Runs  Targets = trial type conditions
    Odd_ds  = cosmo_fmri_dataset(Odd,'mask',mask_fn,...
                            'targets',1:4, ...
                            'chunks', 1);

    Even_ds  = cosmo_fmri_dataset(Even,'mask',mask_fn,...
                            'targets',1:4, ...
                            'chunks', 2);
                        
    labels_odd  = {'HREC';'HFAM';'FAREC';'FAFAM'};
    labels_even = {'HREC';'HFAM';'FAREC';'FAFAM'};
    Odd_ds.sa.labels  = labels_odd;
    Even_ds.sa.labels = labels_even;

    % Combine files at encoding and retrieval to create two files (i.e.,
    % stacking)
    % make sure all ds_* changed from here on
    ds_all = cosmo_stack({Odd_ds, Even_ds});
    
    % remove constant features
    ds_all=cosmo_remove_useless_data(ds_all);

    % cosmo checks to make sure data in right format
    cosmo_check_dataset(ds_all);

    % print dataset
    fprintf('Dataset input:\n');
    cosmo_disp(ds_all);

    % Use cosmo_correlation_measure.
    % This measure returns, by default, a split-half correlation measure
    % based on the difference of mean correlations for matching and
    % non-matching conditions (a la Haxby 2001).
    measure = @cosmo_correlation_measure;

    % define spherical neighborhood with radius of 3 voxels
    radius  = 3; % voxels
    nbrhood = cosmo_spherical_neighborhood(ds_all,'radius',radius);

    % Run the searchlight with a 3 voxel radius
    corr_results = cosmo_searchlight(ds_all, nbrhood, measure);

    % print output
    fprintf('Dataset output:\n');
    cosmo_disp(corr_results);

%     % Plot the output
%     cosmo_plot_slices(corr_results);

    % Define output location
    output_fn = fullfile(data_path, [subjects{ss} 'corr_searchlight.nii']);

    % Store results to disc
    cosmo_map2fmri(corr_results, output_fn);
    
    % display results
    spm_image('Display', output_fn)
    
end