function model_per_job

hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(hpc_job_number), error('Problem with array assigment'); end
fprintf('job number: %i \n', hpc_job_number);

i_model = numJob;
fit_atheo_model(i_model, 1)

end