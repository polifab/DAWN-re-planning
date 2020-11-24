%% decide

% park_E1
% park_E0

%% take into account only last n samples
n = 60;
pos1 = park_E1(1:n,:);
pos0 = park_E0(end-n+1:end,:);

figure(1)
clf
plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'b')
hold on
plot3(pos0(:,1), pos0(:,2), pos0(:,3), 'r')

%% scan
mindiff_true = 1e10;
index_cross = 0;

norm_diff = zeros(n,1);
for i = 1:n
	temp = pos0(i,:);
	
	norm_diff =  zeros(n,1);
	for j = 1:n
		norm_diff(j) =  norm(pos1(j,:) - temp,2);
	end
	
	mindiff = mink(norm_diff,1);
	
	if mindiff < mindiff_true
		mindiff_true = mindiff;
		index_cross = i;
		
		if i ~= 1
			delete(last_cross)
		end

		figure(1)
		last_cross = plot3(pos0(index_cross,1), pos0(index_cross,2), pos0(index_cross,3), '*k');
		drawnow
	end

	
end
%%
clc

figure(2)
clf
plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'b')
hold on
plot3(pos0(:,1), pos0(:,2), pos0(:,3), 'r')

index_crosspoint_0 = find(row_norm2(park_E0) == row_norm2(pos0(index_cross,:)),1,'last');
crosspoint_0 = park_E0(index_crosspoint_0,:);
plot3(crosspoint_0(1,1), crosspoint_0(1,2), crosspoint_0(1,3), 'dr')

index_crosspoint_1 = find(row_norm2(park_E1 - pos0(index_cross,:)) < 100 ,1,'first');
crosspoint_1 = park_E1(index_crosspoint_1,:);
plot3(crosspoint_1(1,1), crosspoint_1(1,2), crosspoint_1(1,3), 'db')

%plot3(pos0(index_cross,1), pos0(index_cross,2), pos0(index_cross,3), 'db')



%%
%{
%%
figure(1)
clf
plot3(park_E1(:,1), park_E1(:,2), park_E1(:,3), 'b')
hold on
plot3(park_E0(:,1), park_E0(:,2), park_E0(:,3), 'r')

grid on
view([30 30]); axis equal

%%
temp = zeros(9000,3);
temp_norm = zeros(9000,1);
for i = 1:9000
	temp(i,:) = park_E0(end-i,:) - park_E1(i,:);
	temp_norm(i) = norm(temp(i,:),2);
end

figure(2)
clf
plot(temp_norm)

figure(3)
clf
plot3(temp(:,1), temp(:,2), temp(:,3), 'b')

%%
last_E0 = park_E0(end-1999:end,:);
first_E1 = park_E1(1:2000,:);


figure(4)
clf
plot3(last_E0(:,1), last_E0(:,2), last_E0(:,3), 'b')
hold on
plot3(first_E1(:,1), first_E1(:,2), first_E1(:,3), 'r')

figure(5)
clf
plot(last_E0-first_E1, 'b')

find(last_E0==first_E1)

%% 

		  
% % CUT E0 TO GET THE PROPER END OF THE PARK ORBIT!
% % change park orbit adjustment
% temporary = zeros(size(park_V2,1),1);
% for i = 1:size(park_V2,1)
% 	temporary(i) = norm(park_V2(i,:) - orb_change_park(1,:) ,2);
% end
% index_last = find(temporary-mink(temporary,1) < 1e-5, 1, 'last');
% %removing the samples after the start of changing orbit
% park_V2(index_last+1:end,:) = [];
%}