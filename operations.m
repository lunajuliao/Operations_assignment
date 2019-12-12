%% GETTING DATA FROM EXCEL FILE


clear all
close all
clc

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64');
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab');

data = readtable('Small_data.xlsx','ReadRowNames',true,'ReadVariableNames',true);
BayComplianceData = readtable('BayComplianceMatrix.xlsx','ReadRowNames',true,'ReadVariableNames',true);

PN = height(data);

for i = 1:PN
    plane(i).AT = str2num(erase(string(data.Arrival(i)),':'));
    plane(i).DT = str2num(erase(string(data.Departure(i)),':'));
    plane(i).P = data.Passengers(i);
    if plane(i).AT > plane(i).DT
        plane(i).DT = 2350;
    elseif ( plane(i).DT - plane(i).AT < 15)
        plane(i).DT = 15+plane(i).AT;
    end
    
    % Assigning numbers to different aircraft types
    if     strcmp(data.Type(i),'B747')
        plane(i).Type = 1;
    elseif strcmp(data.Type(i),'B777')
        plane(i).Type = 2;
    elseif strcmp(data.Type(i),'B787')
        plane(i).Type = 3;
    elseif strcmp(data.Type(i),'B737')
        plane(i).Type = 4;
    elseif strcmp(data.Type(i),'ATR72')
        plane(i).Type = 5;
    else
        plane(i).Type = -1; % Invalid aircraft type
    end
end

%% ORDERING DATA
p=[];
for i = 1:PN
    plane(PN+1).AT = 2500;

    for j = i:PN
        if(plane(j).AT < plane(PN+1).AT)
            plane(PN+1) = plane(j);
            index = j;
        end
        plane(j)
    end
    plane(index) = plane(i);
    plane(i) = plane(PN+1);
    p =[p,plane(i).P];
end

%% OVERLAPPING MATRIX OV
% We use only the lower triangular matrix.
OV_initial = zeros(PN);
for i = 1:PN
    for k = i:PN
        if plane(i).DT >= plane(k).AT
            OV_initial(i,k) = 1;
        end
    end
end
OV = OV_initial';

    
            
%%

% for i = 1:length(ts)
%     for j = 1 : PN
%         if (plane(j).AT <= ts(i) && plane(j).DT >= ts(i))
%             plane(j).G_NG(i) = 1;
%         else
%             plane(j).G_NG(i) = 2;
%     end
%     end
% end

% F = readtable ('flight data.xlsx');
% D = readtable ('distance.xlsx');
% 


%% Creating matrix for the constraints
OVtest= [1 0 0
    1 1 0
    0 1 1];
Nflights=PN;
Bays=2;

Aeq=[];
    for j=1:Bays
         Aeq=[Aeq,diag(diag(ones(PN)))];
    end



for i=1:Nflights
    NULL(i,i)=0
end
Aineq=[];
for i=1:Bays
    L=[];
    for j=1:Bays
        if j==i
            L=[L,OV];
        else
            L=[L,NULL];
        end
    end
    Aineq=[Aineq;L];
end 

for i=1:Nflights*Bays
        rightside_ineq(i,1) =1;
end

for i=1:Nflights
        rightside_eq(i,1) =1;
end

for i=1:Nflights*Bays
        ub(i,1) =1;
end

for i=1:Nflights*Bays
        lb(i,1) =0;
end

for i=1:Nflights*Bays
    ctype(1,i)='B';
end

d=[7 4];
f=[];
for i=1:Bays
    for j=1:Nflights
        f=[f,d(i)*p(j)*2];
    end
end

 x=cplexmilp(f, Aineq, rightside_ineq, Aeq, rightside_eq, [],[],[],lb, ub, ctype);
