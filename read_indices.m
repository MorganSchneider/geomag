fp = fopen('./EEJ_Data/indices.txt');
indices_cell = textscan(fp, '%f %f %f', 'HeaderLines', 3, 'delimiter', ' ');
fclose(fp);

timestamp = indices_cell{1};
kp = indices_cell{2};
rc = indices_cell{3};

startTime = min([swarm(1).time(1), swarm(2).time(1), swarm(3).time(1)]);
z = find(timestamp <= startTime);
startInd = z(end);

endTime = max([swarm(1).time(end), swarm(2).time(end), swarm(3).time(end)]);
z = find(timestamp <= endTime);
endInd = z(end);

timestamp = rot90(timestamp(startInd:endInd));
kp = rot90(kp(startInd:endInd));
rc = rot90(rc(startInd:endInd));









