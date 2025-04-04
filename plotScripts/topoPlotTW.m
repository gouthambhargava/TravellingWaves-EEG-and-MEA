function topoPlotTW(chanlocs, values, make_contour, plot_quiver, direction,interpPoints,axes)

% From: https://in.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot?#functions_tab

% Function 'plot_topography' plots a topographical map of the head over the
% desired points given by 'ch_list' and their assigned 'values'.          
%                                                                         
%   Input parameters:                                                     
%       - ch_list:      Channel list is a structure which must contain the following fields - 
%                       theta and radius - available in
%                       .../Programs/ProgramsMAP/Montages/actiCap64.mat
%       - values:       Phase values that are to be plotted for each channel                                  
%       - make_contour: (Optional, default: false) Boolean that controls if
%                       the contour lines should be plotted.              
%       - plot_quiver:  (Optional, default: false) Boolean that controls if
%                       the direction of propagation of waves is to be 
%                       plotted.                 
%       - direction:    Vector of direction values for all channels
%       - plot_clabels: (Optional, default: false) Boolean that controls if
%                       the text labels of each electrode should be
%                       plotted. - currently not enabled
%       - interpPoints:(Optional, default: 1000) No. of interpolation    
%                       points. The lower N, the lower resolution and     
%                       faster computation.    



    % if nargin < 6, plot_clabels = false;
    % end
    if nargin < 6 
        interpPoints = 1000;
    end

    if nargin < 7 % create axis
        axes = getPlotHandles(1,1,[0.1 0.1 0.8 0.8],0.01,0.1,0);
    end
    
    % Finding the desired electrodes
    % ch_list = upper(ch_list);
    elecInfo = [chanlocs(1:64).theta;chanlocs(1:64).radius];
    theta = elecInfo(1,:)';
    radius = elecInfo(2,:)';
    locations = table(theta,radius);
    idx = 1:size(locations,1);
    values = values(:);
    
    % Global parameters
    %   Note: Head radius should be set as 0.6388 if the 10-20 system is used.
    %   This number was calculated taking into account that the distance from Fpz
    %   to Oz is d=2*0.511. Thus, if the circle head must cross the nasion and
    %   the inion, it should be set at 5d/8 = 0.6388.
    %   Note2: When the number of interpolation points rises, the plots become
    %   smoother and more accurate, however, computational time also rises.
    headRadius     = 5*2*0.511/9;  % 1/2  of the nasion-inion distance
    headExtra      = 1*2*0.511/8;  % 1/10 of the nasion-inion distance
    k = 4;                          % Number of nearest neighbors for interpolation
    
    % Interpolating input data
    % Creating the rectangle grid (-1,1)
    [ch_x, ch_y] = pol2cart((pi/180).*((-1).*locations.theta(idx)+90),locations.radius(idx)); % X, Y channel coords
    % Points out of the head to reach more natural interpolation
    r_ext_points = 1.2;
    [add_x, add_y] = pol2cart(0:pi/4:7*pi/4,r_ext_points*ones(1,8));
    linear_grid = linspace(-r_ext_points,r_ext_points,interpPoints);         % Linear grid (-1,1)
    [interp_x, interp_y] = meshgrid(linear_grid, linear_grid);
        
    % Interpolate and create the mask
    outer_rho = max(locations.radius(idx));
    if outer_rho > headRadius 
        mask_radius = outer_rho + headExtra;
    else                      
        mask_radius = headRadius;
    end
    mask = (sqrt(interp_x.^2 + interp_y.^2) <= mask_radius); 
    add_values = compute_nearest_values([add_x(:), add_y(:)], [ch_x(:), ch_y(:)], values(:), k);
    interp_z = griddata([ch_x(:); add_x(:)], [ch_y(:); add_y(:)], [values; add_values(:)], interp_x, interp_y, 'natural');
    interp_z(mask == 0) = NaN;

    % Plotting the final interpolation
    pcolor(axes,interp_x, interp_y, interp_z);
    shading(axes,'interp');
    hold(axes,'on');
        
    % Contour
    if make_contour==1
       [~, hfigc] = contour(interp_x, interp_y, interp_z,axes); 
       set(hfigc, 'LineWidth',0.75, 'Color', [0.2 0.2 0.2]); 
       hold on;
    end

    % Plotting the head limits as a circle         
    head_rho    = headRadius;                      % Head radius
    head_theta  = linspace(0,2*pi,interpPoints);   % From 0 to 360Âº
    head_x      = head_rho.*cos(head_theta);        % Cartesian X of the head
    head_y      = head_rho.*sin(head_theta);        % Cartesian Y of the head
    plot(axes,head_x, head_y, 'Color', 'k', 'LineWidth',4);
    hold(axes,'on');

    % Plotting the nose
    nt = 0.15;      % Half-nose width (in percentage of pi/2)
    nr = 0.22;      % Nose length (in radius units)
    nose_rho   = [head_rho, head_rho+head_rho*nr, head_rho];
    nose_theta = [(pi/2)+(nt*pi/2), pi/2, (pi/2)-(nt*pi/2)];
    nose_x     = nose_rho.*cos(nose_theta);
    nose_y     = nose_rho.*sin(nose_theta);
    plot(axes,nose_x, nose_y, 'Color', 'k', 'LineWidth',4);
    hold(axes,'on');

    % Plotting the ears as ellipses
    ellipse_a = 0.08;                               % Horizontal exentricity
    ellipse_b = 0.16;                               % Vertical exentricity
    ear_angle = 0.9*pi/8;                           % Mask angle
    offset    = 0.05*headRadius;                   % Ear offset
    ear_rho   = @(ear_theta) 1./(sqrt(((cos(ear_theta).^2)./(ellipse_a^2)) ...
        +((sin(ear_theta).^2)./(ellipse_b^2))));    % Ellipse formula in polar coords
    ear_theta_right = linspace(-pi/2-ear_angle,pi/2+ear_angle,interpPoints);
    ear_theta_left  = linspace(pi/2-ear_angle,3*pi/2+ear_angle,interpPoints);
    ear_x_right = ear_rho(ear_theta_right).*cos(ear_theta_right);          
    ear_y_right = ear_rho(ear_theta_right).*sin(ear_theta_right); 
    ear_x_left  = ear_rho(ear_theta_left).*cos(ear_theta_left);         
    ear_y_left  = ear_rho(ear_theta_left).*sin(ear_theta_left); 
    plot(axes,ear_x_right+head_rho+offset, ear_y_right, 'Color', 'k', 'LineWidth',4);
    plot(axes,ear_x_left-head_rho-offset, ear_y_left, 'Color', 'k', 'LineWidth',4);
   
    % plot the channels
    [ch_x, ch_y] = pol2cart((pi/180).*(locations.theta(idx)+90), locations.radius(idx));
    scatter(axes,ch_x, ch_y, 8,'white','filled', 'LineWidth',1.5);
    % Plotting the quiver plot
    if plot_quiver==1 
        quiver(axes,ch_x, ch_y, cos(direction),sin(direction),'color','white','AutoScaleFactor',0.5,'LineWidth',1.5) 
    end

    
    % Last considerations
    max_height = max([max(nose_y), mask_radius]);
    min_height = -mask_radius;
    max_width  = max([max(ear_x_right+head_rho+offset), mask_radius]);
    min_width  = -max_width;
    L = max([min_height, max_height, min_width, max_width]);
    xlim(axes,[-L, L]);
    ylim(axes,[-L, L]);  
    % colorbar;   % Feel free to modify caxis after calling the function
    axis(axes,'square');
    axis(axes,'off');
    hold(axes,'off');
    % h = gcf;
end

% This function compute the mean values of the k-nearest neighbors
%   - coor_add:     XY coordinates of the virtual electrodes
%   - coor_neigh:   XY coordinates of the real electrodes
%   - val_neigh:    Values of the real electrodes
%   - k:            Number of neighbors to consider
function add_val = compute_nearest_values(coor_add, coor_neigh, val_neigh, k)
    
    add_val = NaN(size(coor_add,1),1);
    L = length(add_val);
    
    for i = 1:L
        % Distances between the added electrode and the original ones
        target = repmat(coor_add(i,:),size(coor_neigh,1),1);
        d = sqrt(sum((target-coor_neigh).^2,2));
        
        % K-nearest neighbors
        [~, idx] = sort(d,'ascend');
        idx = idx(2:1+k);
        
        % Final value as the mean value of the k-nearest neighbors
        add_val(i) = mean(val_neigh(idx));
    end
    
end