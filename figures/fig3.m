% fig 3: illustration of recalibration models

% Define the proportions for rows and columns
rowProportions = [0.2, 0.2, 0.35];
colProportions = [0.2, 0.35, 0.2];

% Convert proportions to normalized units (sum should be 1)
rowProportions = rowProportions / sum(rowProportions);
colProportions = colProportions / sum(colProportions);

% Define the figure size
figureWidth = 420;
figureHeight = 500;

% Create the figure
figure('Position', [0, 0, figureWidth, figureHeight]);

% Loop over rows and columns to create subplots
for i = 1:3
    for j = 1:3
        % Calculate position for each subplot
        left = sum(colProportions(1:j-1));
        bottom = sum(rowProportions(1:3-i));
        width = colProportions(j);
        height = rowProportions(i);
        
        % Create the subplot with specified position
        subplot('Position', [left, bottom, width, height]);
        
        % Example content: display the subplot indices
        text(0.5, 0.5, sprintf('Subplot (%d, %d)', i, j), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end