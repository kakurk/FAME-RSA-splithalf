function oddeven_contrasts(subject, studypath, cond, express)

% Initalize SPM; load this subjects SPM.mat
SPM    = [];
SPMmat = fullfile(studypath, subject, 'SPM.mat');
load(SPMmat)

% Find all of the betas for this condition using the regular expression
matches = regexp(SPM.xX.name, express);
matches = find(~cellfun('isempty', matches));

%-- Calculate Contrast Vectors

cont_vector = zeros(2,length(SPM.xX.name));

% Odd Runs
cont_name{1} = [cond '_Odd'];
cont_vector(1, matches(1:2:length(matches))) = 1/length(matches(1:2:length(matches)));

% Even Runs
cont_name{2} = [cond '_Even'];
cont_vector(2, matches(2:2:length(matches))) = 1/length(matches(2:2:length(matches)));

% Set the conmanager parameters calculate odd and even run betas
matlabbatch = set_conmanger(SPMmat, cont_name, cont_vector);

% Run through spm_jobman
spm_jobman('run', matlabbatch)

% set_conmanager subfunction

function matlabbatch = set_conmanger(fullpath2SPM, cont_name, cont_vec)
    matlabbatch{1}.spm.stats.con.spmmat = {fullpath2SPM};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = cont_name{1};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = cont_vec(1,:);
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = cont_name{2};
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = cont_vec(2,:);
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
end
       

end