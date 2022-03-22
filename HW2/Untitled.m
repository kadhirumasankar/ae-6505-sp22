x = ["0:00", "1:53", "3:18", "4:40", "5:55", "16:58", "21:58", "25:58", "32:23"]

for i=2:length(x)
    convert_to_sec(x(i))
end
convert_back(convert_to_sec(x(2)))

function out = convert_to_sec(time)
    splut = split(time(1), ":");
    out = str2num(splut(1))*60 + str2num(splut(2))
end

function out = convert_back(time)
    out = strcat(num2str(floor(time/60)), ":", num2str(rem(time, 60)))
end