clc; clear all; close all

bag_list = ["circle", "hover", "sine", "square"];
% bag_list = ["hover"];
% freq_list = [1 2 5 10 20 60 120];
% freq_list = [20 60 120];

for bag_idx = 1:length(bag_list)
% for freq_idx = 1:length(freq_list)
clc; close all; clearvars -except bag_idx bag_list

shape = bag_list(bag_idx);
bag = rosbag(strcat(shape, '.bag'));
shape

statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
marker_locs = [];
time = [];
for i = 1:length(statesMsgs)
    current_state = statesMsgs(i);
    current_state = current_state{1}.Data;
    if current_state(3) < 0 && isempty(states)
        
    else
        states = [states current_state(1:12)];
        target_states = [target_states current_state(13:24)];
        time = [time current_state(end)];
    end
    if i > length(statesMsgs)/1
        break
    end
end

states(4:6, 1) = [0; 0; 0];
states(10:12, 1) = [0; 0; 0];

%% Prepping measurement data
cam1loc = [10, 10, 10]';
cam2loc = [-10, 10, 10]';
cam3loc = [-10, -10, 10]';
cam4loc = [10, -10, 10]';
y1 = vecnorm(states(1:3,:)-cam1loc);
y2 = vecnorm(states(1:3,:)-cam2loc);
y3 = vecnorm(states(1:3,:)-cam3loc);
y4 = vecnorm(states(1:3,:)-cam4loc);

%% Perturb data
perturbed_y = [];
for i = 1:length(y1)
    perturbed_y = [perturbed_y [y1(i); y2(i); y3(i); y4(i); states(7,i); states(8,i); states(9,i)]];
end

writematrix(perturbed_y, strcat(shape, '_measurements.csv'));
writematrix([states; target_states; time], strcat(shape, '_actual_states.csv'));

end