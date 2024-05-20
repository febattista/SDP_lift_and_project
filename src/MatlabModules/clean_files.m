function termcode = clean_files(model_dir, filename, splitted)
% CLEAN_FILES Clean up files previously extracted from the zipFile.
    termcode = true;
    if (splitted)
        num_pieces = 5;
        for k=1:num_pieces
            delete(fullfile(model_dir, strcat(filename, sprintf('_%d.mat', k))));
        end
    else 
        delete(fullfile(model_dir, strcat(filename, '.mat')));
    end
end
