function output = ymd_gen(date1, date2)
%YMD_GEN generates a vector nx3 containing a row with 3 columns per each n
%day like:
%	[year, month, day]

% INPUT
% 	date1 = [year1, month1, day1];
% 	date2 = [year2, month2, day2];

% OUTPUT
%	out = [ year1,	month1,	day1;
% 			...,	...,	...;
%			...,	...,	...;
% 			year2,	month2,	day2]
% 

% Example:
% timevector = ymd_gen([2008, 2, 1], [2010, 3, 1])
% computes a time vector cointaining all days from 1 Feb 2008 to 1 March
% 2010

%% intro

%mese = ["Gennaio","Febbraio","Marzo","Aprile","Maggio","Giugno","Luglio", ...
%		"Agosto","Settembre","Ottobre","Novembre","Dicembre"];

months		= 1:12;
days		= [31,28,31,30,31,30,31,31,30,31,30,31];	
dayslp		= [31,29,31,30,31,30,31,31,30,31,30,31];	% days in months during a leap year
years		= date1(1,1):date2(1,1);

%% computation

n = datenum(date2)- datenum(date1) + 1;

%init
vector	= zeros(n,3);
year	= date1(1,1);
month	= date1(1,2);
day		= date1(1,3);

%iteration
for d = 1:n	
	
	vector(d,:) = [year, month, day];
	
	% leap year check
	if mod(year, 4) == 0
		leap = 1;
	else 
		leap = 0;
	end
	
	% year update
	if month == months(end) && day == days(end)
		year = years(find(years == year) + 1);
	end
	% month and day update
	if leap
		if day == dayslp(find(months == month))
			day = 1;
			if month == months(end)
				month = 1;
			else
				month = months(find(months == month) + 1);
			end
		else
			day = day + 1;
		end
	else
		if day == days(find(months == month))
			day = 1;
			if month == months(end)
				month = 1;
			else
				month = months(find(months == month) + 1);
			end
		else
			day = day + 1;
		end
		
	end
	
	% reset month
	if month == 13
		month = 1;
	end
		
	
end

output = vector;
end

