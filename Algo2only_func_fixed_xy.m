function [diff_angle, indices_counter] = Algo2only_func_fixed_xy(nearest_neighbour, points_rock, normals_rock, points_lichen, normals_lichen, indices_counter)

%%% compares the normal of one point on the rock point cloud with a
%%% "no_nearest_n" number of its nearest neighbours. Repeats for all the
%%% points in the rock point cloud.
%%% It uses an algorithm that averages all the angles of the nearest
%%% neighbours. It then compares this averaged angle with the normal of the
%%% rock point cloud. This comparison result is a diff_angle array.
%%% It DOES NOT then compare a minimum angle with all the values in the
%%% diff_angle array.
%%% Each value in the diff_angle array represents the difference in angle
%%% between the normal at the rock point cloud point and its nearest
%%% neighbours on the lichen point cloud or the rock point cloud itself.

% initialise array for cumulatively adding angles from the nearest
% neighbours from the lichen point cloud
angle_cumul_lichen = zeros(1,3);

% array for storing the angle difference between the averaged nearest 
% neighbour normals in the lichen cloud and the single normal from the rock
diff_angle = zeros(1,length(points_rock));

% turn the projected points into a point cloud
%point_cld_lichen = pointCloud(points_lichen);

% for all points in the rock point cloud...
%for ak = 1:length(points_rock)
% find the nearest neighbours in the lichen point cloud. Put the
% indices of all the nearest neighbours in the lichen point cloud in an
% element of the indices array. First input parameter is the lichen 
% point cloud. Second is the query point "vector". Third is number of 
% nearest neighbours. Distance output is unused.
indices = knnsearch(points_lichen(:,1:2),points_rock(:,1:2),'K', nearest_neighbour);
%indices = indices'; % needs to be inverted
for ak = 1:length(indices(:,1))
    for ik = 1: length(indices(1,:))
        % cumulatively add all the nearest neighbour lichen surface normals
        angle_cumul_lichen(1,:) = angle_cumul_lichen(1,:) + normals_lichen(ik,:);
        indices_counter(ik) = indices_counter(ik) + 1;
    end
    % normalise angle_cumul_rock by dividing it by number of nearest
    % neighbours
    total_nn_angle = angle_cumul_lichen./nearest_neighbour;
    % calculate the difference in angle in degrees between the rock point 
    % cloud normal and the average nearest neighbour normals from the 
    % lichen point cloud
    diff_angle(ak) = abs(atan2d(norm(cross(total_nn_angle,normals_rock(ak,:))),dot(total_nn_angle,normals_rock(ak,:))));
end

%end
