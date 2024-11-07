
function MRI_GUI
    % Create a figure for the GUI
    hFig = figure('Name', 'MRI Simulation GUI', 'Position', [100 100 800 600]);
    
    % Persistent variables for phantom and k-space data
    persistent k_space phantom_image

    % Dropdown menu for selecting acquisition scheme
    uicontrol('Style', 'text', 'Position', [20 550 100 20], 'String', 'Acquisition Scheme:');
    scheme_menu = uicontrol('Style', 'popupmenu', 'Position', [130 550 100 20], 'String', {'Cartesian', 'Radial', 'Spiral'});
    
    % Button to generate phantom and k-space
    uicontrol('Style', 'pushbutton', 'Position', [20 500 100 30], 'String', 'Generate Phantom', 'Callback', @(~,~) generateAndDisplayPhantom);
    
    % Button to sample and reconstruct the image
    uicontrol('Style', 'pushbutton', 'Position', [20 450 100 30], 'String', 'Reconstruct Image', 'Callback', @(~,~) reconstructImageCallback);
    
    % Axes for displaying images
    hAx1 = axes('Position', [0.35 0.5 0.6 0.4]);
    hAx2 = axes('Position', [0.35 0.05 0.6 0.4]);

    % Matrix size and structure size input fields
    uicontrol('Style', 'text', 'Position', [20 400 100 20], 'String', 'Matrix Size:');
    matrix_size_input = uicontrol('Style', 'edit', 'Position', [130 400 100 20], 'String', '128');
    
    uicontrol('Style', 'text', 'Position', [20 370 100 20], 'String', 'Structure Size:');
    structure_size_input = uicontrol('Style', 'edit', 'Position', [130 370 100 20], 'String', '30');
    
    % Function to generate phantom and display k-space
    function generateAndDisplayPhantom
        % Retrieve user inputs for matrix size and structure size
        matrix_size = str2double(get(matrix_size_input, 'String'));
        structure_size = str2double(get(structure_size_input, 'String'));
        
        % Generate phantom with custom parameters
        [phantom_image, k_space] = generateCustomPhantom(matrix_size, structure_size);
        
        % Display k-space
        imshow(log(abs(k_space) + 1), [], 'Parent', hAx1);
        title(hAx1, 'k-Space');
    end

    % Function to reconstruct image based on selected scheme
    function reconstructImageCallback
        % Ensure that k_space has been generated
        if isempty(k_space)
            errordlg('Please generate the phantom first.');
            return;
        end

        % Select acquisition scheme
        scheme = scheme_menu.Value;
        
        % Sample k-space according to the selected scheme
        switch scheme
            case 1
                k_space_sampled = cartesianSampling(k_space, size(k_space, 1));
            case 2
                k_space_sampled = radialSampling(k_space, 180);
            case 3
                k_space_sampled = spiralSampling(k_space, 180);
        end
        
        % Reconstruct and display the image
        reconstructed_image = reconstructImage(k_space_sampled);
        imshow(reconstructed_image, [], 'Parent', hAx2);
        title(hAx2, 'Reconstructed Image');
    end
end

% Function to generate custom phantom and k-space data
function [phantom_image, k_space] = generateCustomPhantom(matrix_size, structure_size)
    % Create an empty matrix of the specified size
    phantom_image = zeros(matrix_size, matrix_size);
    
    % Place a circular "structure" in the center of the phantom
    centerX = matrix_size / 2;
    centerY = matrix_size / 2;
    for x = 1:matrix_size
        for y = 1:matrix_size
            if sqrt((x - centerX)^2 + (y - centerY)^2) <= structure_size
                phantom_image(x, y) = 1;
            end
        end
    end

    % Generate k-space
    k_space = fftshift(fft2(ifftshift(phantom_image)));
end

% Additional helper functions for sampling and reconstruction

function k_space_cartesian = cartesianSampling(k_space, size)
    % Extract Cartesian lines in k-space (simple row-by-row sampling)
    k_space_cartesian = k_space(:, 1:size);
end

function k_space_radial = radialSampling(k_space, num_projections)
    [X, Y] = meshgrid(linspace(-0.5, 0.5, size(k_space,1)), linspace(-0.5, 0.5, size(k_space,2)));
    k_space_radial = zeros(size(k_space));
    
    for theta = linspace(0, pi, num_projections)
        x_line = cos(theta) * X + sin(theta) * Y;
        y_line = -sin(theta) * X + cos(theta) * Y;
        k_space_radial = k_space_radial + interp2(X, Y, k_space, x_line, y_line, 'linear', 0);
    end
end

function k_space_spiral = spiralSampling(k_space, num_spirals)
    t = linspace(0, 1, num_spirals);
    x_spiral = t .* cos(2 * pi * t);
    y_spiral = t .* sin(2 * pi * t);
    
    k_space_spiral = zeros(size(k_space));
    
    for i = 1:length(t)
        x_idx = round(size(k_space,1) / 2 + x_spiral(i) * size(k_space,1) / 2);
        y_idx = round(size(k_space,2) / 2 + y_spiral(i) * size(k_space,2) / 2);
        if x_idx > 0 && x_idx <= size(k_space,1) && y_idx > 0 && y_idx <= size(k_space,2)
            k_space_spiral(x_idx, y_idx) = k_space(x_idx, y_idx);
        end
    end
end

function reconstructed_image = reconstructImage(k_space_sampled)
    % Perform inverse FFT on sampled k-space to obtain the reconstructed image
    reconstructed_image = abs(ifftshift(ifft2(fftshift(k_space_sampled))));
end

% Initialize the GUI
MRI_GUI;