function data = load_subject_data(result_folder, sub_slc, file_str)
    % List files in the directory
    files = dir(fullfile(result_folder, file_str));
    
    % Initialize cell arrays for data and pred
    data = cell(1, numel(sub_slc));
    
    for ss = 1:numel(sub_slc)
        % Convert the current element to a two-digit string
        sub_str = sprintf('%02d', sub_slc(ss));
        
        % Find the file that matches 'sub-XX'
        file_name = '';
        for file = files'
            if contains(file.name, ['sub-', sub_str])
                file_name = file.name;
                break;
            end
        end
        
        % Load the data if the file was found
        if ~isempty(file_name)
            data{ss} = load(fullfile(result_folder, file_name));
        else
            error('File for sub-%s not found.', sub_str);
        end
    end
end