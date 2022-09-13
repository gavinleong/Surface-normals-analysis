function plot_algo2_AD_noaxis(rock_trial,new_rock)

rock_trial_index = find(rock_trial);

figure

rock_plot_avg = plot3(new_rock(rock_trial_index,1).*100,new_rock(rock_trial_index,2).*100,new_rock(rock_trial_index,3).*100,'r.', 'MarkerSize', 3);

hold on

camroll(180)
ax = gca;
ax.XAxisLocation = 'top';
ax.YAxisLocation = 'right';
view(0,-90)
% direction = [0 1 0];
%plot([-49; -49], [34; 44], '-k',  [-49; -39], [44; 44], '-k', 'LineWidth', 2)
%text(-50,39, '10 cm', 'HorizontalAlignment','right')
%text(-44, 46, '10 cm', 'HorizontalAlignment','center')
axis tight
%ylim([-5 35])
%xlim([-15 10])
set(gca,'xtick',[],'ytick',[])
set(gca,'XColor','none','YColor','none','TickDir','out')