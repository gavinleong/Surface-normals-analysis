
function [lichen_sample] = lichen_distribution(sample_size)

%%% samples from a lichen thickness distribution obtained in a laser scan 
%%% of stone 30 of ramalina siliquosa.

ptCloud = (pcread('rock_layer_stone30_orient.ply')); % import the point cloud 

% assign point coordinates to individual variables
x = double(ptCloud.Location(1:1:end,1));
y = double(ptCloud.Location(1:1:end,2));
z = double(ptCloud.Location(1:1:end,3));

%% fit a surface to the point cloud using polynomial fit parameters
% There were 25 options for polynomial fitting parameters, so a for loop it
% stores the error sum of squares and we then look for the smallest value 
% and use the polynomial fit that generated it. The fitting function only 
% used 1/100 of the points as the fit didn't change very much from no 
% decimation to ~100 decimation. This made the code run quicker. The 
% surface fit is then plotted. Plotting of data points is suppressed.

% Creates a "poly" string array that contains "poly11", "poly12" etc. until 
% "poly55" via triangle number 25, to be used for comparing the best poly
% fit, i.e. the smallest sum of squared error, a goodness of fit parameter
poly = [];
for k = 1:5
    for j = 1:5
        poly = string(cat(1,poly,strcat('poly',num2str(k),num2str(j))));
    end
end

sse_poly = []; %create array of sse values for minimisation
poly = poly'; %fit can only take 1xn arrays of character vectors  

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% may need to adjust the decimation for larger point clouds!
%%%%%%%%%%%%%%%%%%%%%%%%%%%

decimate = 1;

for i = poly
    [fitobject_0,gof_0,output_0] = fit([x(1:decimate:end,1),y(1:decimate:end,1)],z(1:decimate:end,1),char(i));
    sse_poly = cat(1,sse_poly,gof_0.sse);
end

% minimise the sse_poly to find best poly value for use as actual fit
% method. The min function handily has an index output, so we can use the
% min_index later on to index poly
[minimum,min_index] = min(sse_poly);

% use the min_index poly for actual fitting
[fitobject,gof,output] = fit([x(1:decimate:end,1),y(1:decimate:end,1)],z(1:decimate:end,1),char(poly(min_index)));

% store and tell us the formula for the polynomial fitting as a string
formula_poly = formula(fitobject);

% store and tell us the coefficient values "p00", "p01" etc. 
coeffvals = coeffvalues(fitobject);

%% import the lichen layer for projection onto a polynomial surface fit.
% for artifical lichen creation:

% offset the point cloud by the lichen thickness of interest. If we only
% add the noise to the z-axis, we won't get the right offset when we
% compare the normals. Thus we need to offset the point cloud with
% noise in a direction normal to the local general rock surface. We have 
% that information in the matrix 'normals' already. Since the normals are 
% normalised, we know that applying the normals to the vectors as it is by 
% default will give you a modulus of 1, that is, the vector will always 
% move a total distance of 1 (metre). To scale the normal offset properly 
% we only move the offset in the direction of the normal by 1 mm (scale 
% factor 0.001) or any other lichen_thickness. Then we add it on to the 
% x y z to get x_2 y_2 z_2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT THE LICHEN LAYER HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the second point cloud, the lichen layer
lichenCloud = (pcread('lichen_layer_stone30_orient.ply')); % import the point cloud 

l_x = double(lichenCloud.Location(1:1:end,1));
l_y = double(lichenCloud.Location(1:1:end,2));
l_z = double(lichenCloud.Location(1:1:end,3));
% store lichen points in single array
xyz_lichen = [l_x l_y l_z]; 

% convert the lichen xyz points into a point cloud
xyz_lichen_ptcld = pointCloud(xyz_lichen);

%% project the lichen onto the surface.
% Firstly, store the offset cloud (lichen) and the cloud from which the 
% surface was generated (stone) in a single matrix.
% Then create a coefficient matrix * the x and y values function to 
% calculate the z value of the closest point on the surface relative to the
% offset cloud point. The minimisation MATLAB function will try to solve
% this.
% Then create a distance-calculating function for getting the shortest
% distance between two 3D points. Again the minimisation MATLAB function
% will try to minimise this distance, so that the x and y coordinates of
% where the closest point on the surface is, is solved.
% The loop for minimising the distance and getting the closest point to
% surface uses fminunc. It has many outputs, but we're only interested in
% the solution (x,y) and perhaps the distance. The function to be
% minimised involves both of the previous functions.
% The minimisation function needs x and y points as suggestions, so that
% the minimisation can start from somewhere. It is not assured that the
% global minimum is the same as the local minimum unless solution is found

% store lichen x and y points in a single vector
xy_surface_lichen = [double(l_x) double(l_y)].';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for calculating the x and y values of the polynomial fit. The
% coefficients don't need to be parameters in function as they're constant.
% m(1) are the x coords of a point cloud
% m(2) are the y coords of a point cloud
% coeffvals are the coefficients calculated from a fitting using a
% polynomial method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff_xy_func = @(m)sum(coeffvals.*[m(1).^0.*m(2).^0 m(1).^1.*m(2).^0 m(1).^0.*m(2).^1 m(1).^2.*m(2).^0 m(1).^1.*m(2).^1 m(1).^0.*m(2).^2 m(1).^3.*m(2).^0 m(1).^2.*m(2).^1 m(1).^1.*m(2).^2 m(1).^0.*m(2).^3 m(1).^4.*m(2).^0 m(1).^3.*m(2).^1 m(1).^2.*m(2).^2 m(1).^1.*m(2).^3 m(1).^0.*m(2).^4 m(1).^5.*m(2).^0 m(1).^4.*m(2).^1 m(1).^3.*m(2).^2 m(1).^2.*m(2).^3 m(1).^1.*m(2).^4 m(1).^0.*m(2).^5],2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for loop to solve minimisation to get closest point from l_x, l_y, z_lichen to
% % the surface. This is for the lichen surface only!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for qk = 1:length(xy_surface_lichen) % run loop over index of lichen surface points
    jj = xy_surface_lichen(:,qk); % store a single lichen surface point location
    k1 = l_x(qk,1); % store x, y and z of single lichen surface location
    k2 = l_y(qk,1);
    k3 = l_z(qk,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function to return the shortest distance between a point on the lichen
    % point cloud and the function/surface generated from the polynomial fit to
    % the cleaned surface. It also outputs the xy coordinates of the point 
    % where the shortest distance touches the function/surface.
    % The z coordinate of the shortest point will still need calculating after
    % this, which requires plugging it in to the coeff_xy_func formula.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % values for m(1) and m(2) for which the distance is the minimum will
    % be found using fminunc on distance function. 'coeff_xy_func' function
    % is needed to find the z value 
    distance_function = @(m)sqrt((m(1)-k1).^2 + (m(2)-k2).^2 + (coeff_xy_func(m)-k3).^2);
    % minimsation function for 'distance_function', outputs 'solution' are 
    % the x,y values on the surface where the minimum distance from point
    % to surface is found, and the 'distance' for each corresponding
    % solution
    % the 'distance' takes input in opposite index form to 'solution'
    options = optimoptions('fminunc','Display','off');
    if k3 < coeff_xy_func([k1 k2])
        [solution(:,qk),distance(qk,:)] = fminunc(distance_function,jj,options);
    end
    if k3 > coeff_xy_func([k1 k2])
        [solution(:,qk),distance(qk,:)] = fminunc(distance_function,jj,options);
        distance(qk,:) = -distance(qk,:);
    end
    [solution_full(:,qk),distance_full(qk,:)] = fminunc(distance_function,jj,options);
end

% comment out the following line to avoid having data written over!
save('lichenlayer.mat','solution');

%fit_layer_lichen = load('lichenlayer.mat');

% since it's not possible for the lichen to go into the rock surface (but
% here the rock surface is actually a mathematical fit), we shift all
% values of distance to non-negative
distance2 = distance - min(distance); 
% we remove the zero value as it doesn't behave well with the fitting (or
% does it? test this)
distance2 = nonzeros(distance2);

% may not be necessary, but ensures that the values are displayed as long
% format: many decimal places
format long

figure()
% fit histogram to lichen thickness values
histfit(distance2, 200, 'kernel') 
% fit kernel probability distribution to the lichen data
pd_final = fitdist(distance2,'kernel'); 
% perform one-samples Kolmogorov-Smirnov test at 5% significance level. If
% h_kernel = 1, we reject the null hypothesis that the data comes from this 
% kernel distribution, if = 0, we fail to reject the null hypothesis.
[h_kernel,p] = kstest(distance2,'CDF',pd_final)
title('kernel kstest')

figure()
lichen_sample = random(pd_final, [sample_size,1]); % randomly populate the kernel (pd_final) probability distribution function
histogram(lichen_sample) % plot the randomly populated kernel as a histogram (no distribution fit)

hold on

x_values = 0:0.0001:0.03; % x values for plotting the kernel distribution fit
y = pdf(pd_final,x_values); % give probability values for each individual x_value using kernel distribution fit
plot(x_values,y,'LineWidth',2) % plot the distribution fit
pd_final = fitdist(distance2,'kernel') % 
shg