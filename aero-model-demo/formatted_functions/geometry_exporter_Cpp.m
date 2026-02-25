close all
clear
clc

%% Load the geometry file

% Select the printing layout
layout = "horizontal";   % "vertical" or "horizontal"

% Select and load the configuration
config = 1;
file_names = {'3cat8_geometry_bare.mat', '3cat8_geometry_panel.mat', '3cat8_geometry_antenna.mat'};
file_name = char(file_names{config});
S = load(file_name);

% Extract the configuration name
conf_name = erase(file_name, '3cat8_geometry');
conf_name = erase(conf_name, '.mat');

%% Open output file

out_name = ['geometry' conf_name '.txt'];
fid = fopen(out_name,'w');

if fid == -1
    error('Cannot open output file');
end

%% Write arrays

vars = fieldnames(S);

for v = 1:numel(vars)

    name = vars{v};
    data = S.(name);

    if ~isnumeric(data)
        continue
    end

    [r,c] = size(data);
    cpp_name = [name conf_name];

    % --- HEADER ---
    if layout == "horizontal"
        fprintf(fid,'const double %s[%d][%d] = {{', cpp_name, r, c);
    else
        fprintf(fid,'const double %s[%d][%d] = {\n', cpp_name, r, c);
    end

    % --- DATA ---
    for i = 1:r

        % Start row
        if layout == "vertical"
            fprintf(fid,'    {');
        elseif i > 1
            fprintf(fid,' {');
        end

        % Print row values
        for j = 1:c
            if j < c
                fprintf(fid,'%.16g, ', data(i,j));
            else
                fprintf(fid,'%.16g', data(i,j));
            end
        end

        % End row
        if i < r
            if layout == "vertical"
                fprintf(fid,'},\n');
            else
                fprintf(fid,'},');
            end
        else
            fprintf(fid,'}');
        end
    end

    % --- FOOTER ---
    if layout == "horizontal"
        fprintf(fid,'};\n\n');
    else
        fprintf(fid,'\n};\n\n');
    end

end

%% Close file
fclose(fid);

fprintf('C++ arrays written to: %s\n', out_name);