clc; clear all; close all;
bag = rosbag('2022-04-20-15-32-08.bag');

%%
clc; clearvars -except bag; close all;
statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
time = [];
for i = 1:length(statesMsgs)
    current_state = statesMsgs(i);
    current_state = current_state{1}.Data;
    states = [states current_state(1:12)];
    target_states = [target_states current_state(13:24)];
    time = [time current_state(end)];
end
time = (time)./1e9;
plot(time, states(1,:))
hold on
plot(time, target_states(1,:))

target_states = [];
time = [];
statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/target_states'),'DataFormat','struct');
for i = 1:length(statesMsgs)
    current_state = statesMsgs(i);
    current_state = current_state{1};
%     states = [states current_state(1:12)];
    target_states = [target_states [current_state.Pose.Pose.Position.X; current_state.Pose.Pose.Position.Y; current_state.Pose.Pose.Position.Z]];
    time = [time current_state.Header.Stamp.Sec + current_state.Header.Stamp.Nsec/1e9 ];
end
time = (time);
plot(time, target_states(1,:))