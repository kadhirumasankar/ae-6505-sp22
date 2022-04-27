clc; clear all; close all;

%% Prepping state data

m = 1.545;
Ixx = 0.029125;
Iyy = 0.029125;
Izz = 0.055225;
g = 9.81;

bag = rosbag('imu.bag');

statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
marker_locs = [];
time = [];
for i = 1:length(statesMsgs)
    if mod(i,1)==0
        current_state = statesMsgs(i);
        current_state = current_state{1}.Data;
        states = [states current_state(1:12)];
        target_states = [target_states current_state(13:24)];
%         marker_locs = [marker_locs get_marker_locs(current_state(1:12))];
        time = [time current_state(end)];
%         figure(1)
%         scatter3(current_state(1), current_state(2), current_state(3), 10, 'm', 'filled')
%         hold on
%         scatter3(marker_locs(1, end), marker_locs(2, end), marker_locs(3, end), 5, 'black', 'filled')
%         scatter3(marker_locs(4, end), marker_locs(5, end), marker_locs(6, end), 5, 'black', 'filled')
%         scatter3(marker_locs(7, end), marker_locs(8, end), marker_locs(9, end), 5, 'black', 'filled')
%         scatter3(marker_locs(10, end), marker_locs(11, end), marker_locs(12, end), 5, 'black', 'filled')
%         scatter3(marker_locs(13, end), marker_locs(14, end), marker_locs(15, end), 5, 'black', 'filled')
%         xlim([current_state(1)-.2 current_state(1)+.2])
%         ylim([current_state(2)-.2 current_state(2)+.2])
%         zlim([current_state(3)-.2 current_state(3)+.2])
%         hold off
    end
    if i > length(statesMsgs)/1
        break
    end
end
time = (time)./1e9;

std(states(7,:))^2
std(states(8,:))^2
std(states(9,:))^2
disp("Everything in SI, angles in radians")