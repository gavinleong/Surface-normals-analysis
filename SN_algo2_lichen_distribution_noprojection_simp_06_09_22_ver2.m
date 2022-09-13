clear all; 
%close all % close and clear previous figures and workspace
total_time = tic();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure to change the rock_folder, lichen_folder and lichen_write
% version increment so as not to overwrite previous outputs


%% 1. read in a point cloud in .ply format.
% If using openMVS as the last step convert the .ply into a more standard 
% .ply by opening in MeshLab and exporting it as .ply again. Make sure to
% crop out only the section of the point cloud of interest, as the smaller 
% the point cloud, the quicker the code will run. Also, scale the point
% cloud to obtain units in metres. Plotting suppressed.
% Creates folders to store the rock and lichen visualisation outputs.
% Creates a folder to store the lichen point cloud for visualising the
% depth map approach for comparison
tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%         This is for the rock surface          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptCloud = (pcread('model_dense_mesh_for_paper_scaled_cleaned_poissondsksamp_79256_0,2_0,002228.ply')); % import the point cloud 

rock_folder = char('rock_22_09_06_x1_ver8\');
lichen_folder = char('lichen_22_09_06_x1_ver8\');

mkdir(fullfile('C:\Users\gavin\OneDrive - University College London\Surface normals paper\All_figures_13.5.21\from_matlab\', rock_folder));
mkdir(fullfile('C:\Users\gavin\OneDrive - University College London\Surface normals paper\All_figures_13.5.21\from_matlab\', lichen_folder));

mkdir("lichen_ply");
lichen_write = 'lichen_ply\ADPC_for_depthmap_ver8.ply';

disp("1. importing raw cloud:")
toc()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. work out the normals
% of each point on the point cloud (points_rock) using the 6 nearest
% neighbours. Then calculates a plane fit to the point cloud and orients
% all the points so that its average normal aligns with the z-axis.
% Orients the normals away from a point of origin. Plotting is suppressed
tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%         This is for the rock surface          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gets xyz of point cloud without any decimation and turns the values into
% a double which is required by the fitting step. uvw assigns the 
% normals calculated earlier to single vectors. 
x = double(ptCloud.Location(1:1:end,1));
y = double(ptCloud.Location(1:1:end,2));
z = double(ptCloud.Location(1:1:end,3));

points_rock = [x y z];

% fit a plane to points_rock to determine necessary rotation
% using using Object-oriented tools for fitting conics and quadrics by
% Matt J (planarFit())
pfit=planarFit(points_rock.');
x0=mean(points_rock);
Rotate=pfit.vecrot(pfit.normal,[0,0,1]); % compute 3x3 rotation matrix
points_rock_orient = (points_rock-x0)*Rotate.'+x0; % rotated point cloud

% for use when a point cloud version is needed
points_rock_orient_cloud = pointCloud(points_rock_orient);
% so I don't need to change name of variables further on
points_rock = points_rock_orient; 

% this is to ensure that any use of xyz is from the oriented point cloud,
% such as with the displacing of the lichen in the  x_lichen = double(x -
% (u.*lichen_sample*2)); steps.
x = points_rock(:,1);
y = points_rock(:,2);
z = points_rock(:,3);

% work out the normals of the point cloud, 6 nearest neighbours
normals_rock = pcnormals(points_rock_orient_cloud, 6);
u = normals_rock(1:1:end,1);
v = normals_rock(1:1:end,2);
w = normals_rock(1:1:end,3);

% Makes all the normals point away from the "sensorCenter"
sensorCenter = [0,0,0]; 
for k = 1 : numel(x)
   p1 = sensorCenter - [x(k),y(k),z(k)];
   p2 = [u(k),v(k),w(k)];
   % Flip the normal vector if it is pointing towards the sensor.
   angle = atan2(norm(cross(p1,p2)),p1*p2');
   if angle < pi/2 && angle > -pi/2
       u(k) = -u(k);
       v(k) = -v(k);
       w(k) = -w(k);
   end
end

normals_rock = [u v w];

disp("2. normals:")
toc()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. import the lichen layer for projection onto a polynomial surface fit.
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

% version 2 of lichen creation: using the lichen_distribution function,
% obtain a distribution of lichen growth thicknesses using a laser scan of
% a surface of Stonehenge with lichen. Lichen_sample then represents a
% random sampling from this distribution, which is used to create the
% lichen point cloud. This point cloud is processed by denoising using 
% pcdenoise function. This filter can be tuned by adjusting the
% NumNeighbours and Threshold varibles. This cleaned point cloud is then
% saved as .ply. The rock and lichen point clouds are plotted. The normals
% of the lichen point cloud are calculated and normals oriented away from
% origin. If the lichen point cloud has already been created, then reuse
% it.

tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%         This is for the LICHEN surface        %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT THE LICHEN LAYER HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the second point cloud, the lichen layer

sample_size = length(x);
% lichen_sample, samples from a distribution of distances from the rock
% surface that represents ramalina siliquosa thicknesses. 

% if lichen layer has already been created, then reuse the .ply file
if isfile(lichen_write)
     xyz_lichen_ptcld = (pcread(string(lichen_write)));
     x_lichen = double(xyz_lichen_ptcld.Location(1:1:end,1));
     y_lichen = double(xyz_lichen_ptcld.Location(1:1:end,2));
     z_lichen = double(xyz_lichen_ptcld.Location(1:1:end,3));
     points_lichen = [x_lichen y_lichen z_lichen];
% if lichen layer has not been created yet, generate it from the
% lichen_distribution function, making sure we have the right lichen type
% that we want (fruticose, foliose, crutose?)
else
    [lichen_sample] = lichen_distribution(sample_size);
    % with the scalar of distance obtained from the lichen_distribution
    % function, we multiply with the normals of the rock surface and add it to
    % the rock surface points to get the lichen cloud.
    x_lichen = double(x - (u.*lichen_sample*2));
    y_lichen = double(y - (v.*lichen_sample*2));
    z_lichen = double(z - (w.*lichen_sample*2));

    % denoising the lichen point cloud to improve visualisation, this will
    % reduce the number of lichen points in the final visualisation but may
    % improve the ability to outine carvings
    xyz_lichen_ptcld = pointCloud([x_lichen y_lichen z_lichen]);
    [ptCloudOut,inlierIndices,outlierIndices] = pcdenoise(xyz_lichen_ptcld, NumNeighbors=4, Threshold = 0.3);
    points_lichen = [x_lichen(inlierIndices) y_lichen(inlierIndices) z_lichen(inlierIndices)]; 
    
    % redefining the x, y and z_lichen to be the denoised points
    x_lichen = points_lichen(:,1);
    y_lichen = points_lichen(:,2);
    z_lichen = points_lichen(:,3);

    pcwrite(pointCloud([x_lichen,y_lichen,z_lichen]),lichen_write)
    
    % store lichen points in single array
    
end

% plot the rock and lichen clouds
set(0,'DefaultFigureVisible','on');

figure(5000)
%quiver3(x.*100,y.*100,z.*100,u,v,w,0.5,'b');
%hold on
plot3(points_rock_orient(:,1).*100,points_rock_orient(:,2).*100, points_rock_orient(:,3).*100,'k.', 'MarkerSize', 8);
hold on
plot3(x_lichen.*100,y_lichen.*100,z_lichen.*100,'g.','MarkerSize', 8);
axis equal


% convert the lichen xyz points into a point cloud
xyz_lichen_ptcld = pointCloud(points_lichen);
% work out the normals of the point cloud, 6 nearest neighbours
% this is done before any projection is applied. So the neighbours are
% obtained from 3D space
normals_lichen = pcnormals(xyz_lichen_ptcld, 6);

% separate the normals into xyz vector directions to reorient them
u_lichen = normals_lichen(1:1:end,1);
v_lichen = normals_lichen(1:1:end,2);
w_lichen = normals_lichen(1:1:end,3);

% set a general direction of the normals so they all point in the same 
% direction 
sensorCenter = [0,0,0]; 
for k = 1 : numel(x_lichen)
   p1 = sensorCenter - [x_lichen(k),y_lichen(k),z_lichen(k)];
   p2 = [u_lichen(k),v_lichen(k),w_lichen(k)];
   % Flip the normal vector if it is pointing towards the sensor.
   angle = atan2(norm(cross(p1,p2)),p1*p2');
   if angle < pi/2 && angle > -pi/2
       u_lichen(k) = -u_lichen(k);
       v_lichen(k) = -v_lichen(k);
       w_lichen(k) = -w_lichen(k);
   end
end

normals_lichen = [u_lichen v_lichen w_lichen];

disp("3. lichen layer creation:")

toc()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
set(0,'DefaultFigureVisible','off');
nn_var = 4;
tic()

% indices_counter = zeros(length(x),1);
% for nearest_neighbour = nn_var % try changing it around for a better distribution
%     angle_diff_lichen = Algo2only_func_fixed(nearest_neighbour, points_rock, normals_rock, points_rock, normals_rock, indices_counter);
%     for i = 1:30
%         points_plot_rock4plotalgo2AD = angle_diff_lichen > i;
%         plot_algo2_AD_noaxis(points_plot_rock4plotalgo2AD,points_rock);
%         print('-dtiff','-r600','C:\Users\gavin\OneDrive - University College London\Surface normals paper\All_figures_13.5.21\from_matlab\' + string(rock_folder) + sprintf("%d",i) + '_a_' + sprintf("%d",nearest_neighbour) + '_nn_rock.tif')
%     end 
% end

%nn_var = 4;
tic()
indices_counter = zeros(length(inlierIndices),1);
for nearest_neighbour = nn_var % try changing it around for a better distribution
    [angle_diff_lichen, indices_counter] = Algo2only_func_fixed_xy(nearest_neighbour, points_rock, normals_rock, points_lichen, normals_lichen, indices_counter);
    for i = 1:50
        points_plot_rock4plotalgo2AD = angle_diff_lichen > i;
        plot_algo2_AD_noaxis(points_plot_rock4plotalgo2AD,points_rock);
        % export the lichen point cloud with selected points
        %pcwrite(pointCloud([x(points_plot_rock4plotalgo2AD),y(points_plot_rock4plotalgo2AD),z(points_plot_rock4plotalgo2AD)]),append('C:\Users\gavin\OneDrive - University College London\Surface normals paper\All_figures_13.5.21\from_matlab\', lichen_folder, string(i), ".ply"));
        % export the rock point cloud with selected points
        %pcwrite(pointCloud([x_lichen,y_lichen,z_lichen]),lichen_write)

        print('-dtiff','-r600','C:\Users\gavin\OneDrive - University College London\Surface normals paper\All_figures_13.5.21\from_matlab\' + string(lichen_folder) + sprintf("%d",i) + '_a_' + sprintf("%d",nearest_neighbour) + '_nn_lichen.tif')
    end 
end
toc()
toc()
disp("4. plot_figures")

toc()