%% GETTING DATA FROM EXCEL FILE


clear all
close all
clc

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64');
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab');

AllData = readtable('SampleDataWithAircraftType.xlsx','ReadRowNames',true,'ReadVariableNames',true);
AllBayComplianceData = importdata('BayComplianceMatrix.xlsx');

% Define number of planes and number of bays used
% PN = height(data); % (The old way.)
PN = 5;
NBays = 4;
data = AllData(1:PN,:);
BayComplianceData = AllBayComplianceData(1:NBays,:);
% Distances to bays
% Please indicate which version of MATLAB you are using.
MATLAB_VERSION = 1; % 0 for R2018b, 1 for R2019b.
if MATLAB_VERSION
    d = readmatrix('distance.xlsx');
    d(:,1) = [];
    d(1,:) = [];

else
    temp = importdata('distance_R2018b.xlsx'); % Temporary variable
    d = zeros(size(temp,1),size(temp,2));
    d = str2double(temp);
end

%%

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

%% INCORPORATE THE TOWING TIME FOR EACH PLANE 
TT=30; %in minutes
for i=1:PN
    plane(i).ATT=plane(i).AT+TT;
    plane(i).DTT=plane(i).DT-TT;
    if(mod(plane(i).ATT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=mod(plane(i).ATT,100)-60;
        plane(i).ATT=floor(plane(i).ATT/100)*100+100+minutes;
    else
    end
    
    if(mod(plane(i).DTT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=100-mod(plane(i).DTT,100);
        plane(i).DTT=floor(plane(i).DTT/100)*100+60-minutes;
    else
    end
end

%% INCORPORATE THE BUFFER TIME FOR EACH PLANE 
BT=10; %in minutes
for i=1:PN
    plane(i).AT=plane(i).AT-BT;
    plane(i).DT=plane(i).DT+BT;
    if(mod(plane(i).AT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=100-mod(plane(i).AT,100);
        plane(i).AT=floor(plane(i).AT/100)*100+60-minutes;
    else
    end
    
    if(mod(plane(i).DT,100) >=60)%correction of the time, to be presented as hours:minutes
        minutes=mod(plane(i).DT,100)-60;
        plane(i).DT=floor(plane(i).DT/100)*100+100+minutes;
    else
    end
end


%% OVERLAPPING MATRIX OV (overlap)
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

%% OVERLAPPING MATRIX OVT (overlap towing)
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

%% MATRIX RELATING PLANE TYPE AND BAY COMPATIBILITY VECTOR
PlaneCompatibilityVectors = [];
for i = 1:PN
    PlaneCompatibilityVectors = [PlaneCompatibilityVectors BayComplianceData(:,plane(i).Type)];
end
RHS_comp =[];
for i = 1:NBays
    RHS_comp = [RHS_comp;PlaneCompatibilityVectors(i,:)'];
end
% This right hand side has to be an inequality constraint
% We have to expand the inequality constraint coefficients with a unity
% matrix with the dimensions of the number of decision variables.
% We add RHS_comp as the right hand side.


%% Creating matrix for the constraints
%NBays=4; (Moved to top)
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
%incorporate distance matrix that tells us the distance from a fixed
%terminal, pre-assigned before to every plane.

% d=readmatrix('distance.xlsx'); (moved to top)
% d(:,1) = [];

f=[];
% set the vector of the coefficients of the objective function for distance
% between bays and gates
for i=1:NBays
    for j=1:PN
       f=[f,d(i,(plane(j).terminal))*p(j)*2];
    end
end

%% Adding Bay Compliance constraints
% Number of Decision variables = Bays*PN
Aineq = [Aineq;diag(diag(ones(NBays*PN)))];
rightside_ineq = [rightside_ineq; RHS_comp];

%cplex implementation 
 x=cplexmilp(f, Aineq, rightside_ineq, Aeq, rightside_eq, [],[],[],lb, ub, ctype);

%% DATA DISPLAYING
% x is the bay number

%building x_output - a matrix which has the decision variables varying by
%plane on the columns and by bays on the rows
x_output= [];

for i = 1:NBays
 x_output = [x_output; x([(i-1)*PN+1: (i-1)*PN+5])'];
%  [(i-1)*PN+1: (i-1)*PN+5]
end
% x_output = x_output';
% x_output
a1=[];
for i = 1 : PN
    a1(i, 1:3) = [(plane(i).AT-mod(plane(i).AT, 100))/100, mod(plane(i).AT, 100), 0];
     a1(i+PN, 1:3) = [(plane(i).DT-mod(plane(i).DT, 100))/100, mod(plane(i).DT, 100), 0];
    plane(i).Color = [i*40, i*20, 50*i]./255;
end


for k = 1:NBays

% b is a vector with the number of the bay, for each plane i have to give a 
%     plot command with specify the color of the line

bayn = k;
b = [bayn, bayn];


for i = 1 : PN
    
    if x_output(k,i)==1
        a = [ a1(i, 1:3);  a1(i+PN, 1:3)];
        c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
        d=strcat(c(:,1),':',c(:,2),':',c(:,3));
        plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
        hold on;
        datetick('x','HH:MM:SS')  

        
     end
    
    
end
 
end
grid on
ylim([0, NBays+1]);
legend
