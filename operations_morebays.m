%% GETTING DATA FROM EXCEL FILE


clear all
close all
clc

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64');
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab');

data = readtable('SampleDataWithAircraftType.xlsx','ReadRowNames',true,'ReadVariableNames',true);
BayComplianceData = readtable('BayComplianceMatrix.xlsx','ReadRowNames',true,'ReadVariableNames',true);

% define how many flights (aircrafts) we have scheduled for this airport
PN = height(data);

for i = 1:PN
    %Get the data from the time from the excel, with some corrections it
    %might happen in the schedule
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
    
    %assign a number to each terminal
    if     strcmp(data.Terminal(i),'A')
        plane(i).terminal = 1;
    elseif strcmp(data.Terminal(i),'B')
        plane(i).terminal = 2;
    elseif strcmp(data.Terminal(i),'C')
        plane(i).terminal = 3;
    elseif strcmp(data.Terminal(i),'D')
        plane(i).terminal = 4;
    end
    
    %assign a minimum and maximum number of possible passengers to each
    %type of aircraft
    r=rand;
    switch plane(i).Type
        case 1
           plane(i).Passenger_min =416;
           plane(i).Passenger_max =524;
        case 2
           plane(i).Passenger_min =301;
           plane(i).Passenger_max =368;
        case 3
           plane(i).Passenger_min =242;
           plane(i).Passenger_max =330;
        case 4
           plane(i).Passenger_min =85;
           plane(i).Passenger_max =215;
        case 5
           plane(i).Passenger_min =68;
           plane(i).Passenger_max =78;
        case -1
           plane(i).Passenger_min =0;
           plane(i).Passenger_max =0;
    end
    %assign a random number of passengers to each plane considering the limits
    plane(i).Passenger= floor((plane(i).Passenger_max - plane(i).Passenger_min)*r) + plane(i).Passenger_min;
    
end

%incorporate distance matrix that tells us the distance from a fixed
%terminal, pre-assigned before to every plane.
d=readmatrix('distance.xlsx');
d(:,1) = [];
d(1,:) = [];
%define how many bays the airport has, from the distance data
NBays =size(d,1);

%% ORDERING DATA - order the flights by time of arrival and set them into a matrix
p=[];
for i = 1:PN
    plane(PN+1).AT = 2500;

    for j = i:PN
        if(plane(j).AT < plane(PN+1).AT)
            plane(PN+1) = plane(j);
            index = j;
        end
        plane(j);
    end
    plane(index) = plane(i);
    plane(i) = plane(PN+1);
    p =[p,plane(i).P];
end

%% INCORPORATE THE BUFFER TIME FOR EACH PLANE 
BT=10; %in minutes
for i=1:PN
    plane(i).AT=plane(i).AT-BT;
    plane(i).DT=plane(i).DT+BT;
    if(mod(plane(i).AT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=mod(plane(i).AT,100)-60;
        plane(i).AT=floor(plane(i).AT/100)*100+100+minutes
    else
    end
    
    if(mod(plane(i).DT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=mod(plane(i).DT,100)-60;
        plane(i).DT=floor(plane(i).DT/100)*100+100+minutes
    else
    end
end


%% OVERLAPPING MATRIX OV
% From the time data, we compute a matrix that shows which planes
% overlap and since they are ordered, we use only the lower triangular(+diagonal) part of the matrix.
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
%In cplex, the constraint matrix has to be separted in equality constraint
%matrix and inequality constraint matrix

%Computation of equality constraint matrix: 
%1-plane constraint: there can not be more than one bay with the some value
Aeq=[];
    for j=1:NBays
         Aeq=[Aeq,diag(diag(ones(PN)))];
    end


%Computation of inequality constraint matrix:
%1-bay and plane constraint: whenever there is an extra plane, we rewrite a
%constraint concerning the ground existing planes which have to be assigned to the
%bays (move or stay in the same place)
for i=1:PN
    NULL(i,i)=0;
end
Aineq=[];
for i=1:NBays
    L=[];
    for j=1:NBays
        if j==i
            L=[L,OV];
        else
            L=[L,NULL];
        end
    end
    Aineq=[Aineq;L];
end 

%set the vector of the right side of the inequality constraints
for i=1:PN*NBays
        rightside_ineq(i,1) =1;
end

%set the vector of the right side of the equality constraints
for i=1:PN
        rightside_eq(i,1) =1;
end

%set the vector of the upper bound of the decision variables
%set the vector of the lower bound of the decision variables
%set the vector that states the decision variables are binary or integer
for i=1:PN*NBays
        ub(i,1) =1;
        lb(i,1) =0;
        ctype(1,i)='B';
end
f=[];
% set the vector of the coefficients of the objective function for distance
% between bays and gates
for i=1:NBays
    for j=1:PN
       f=[f,d(i,(plane(j).terminal))*p(j)*2];
    end
end

%cplex implementation 
 x=cplexmilp(f, Aineq, rightside_ineq, Aeq, rightside_eq, [],[],[],lb, ub, ctype);
