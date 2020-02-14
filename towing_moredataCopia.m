%% GETTING DATA FROM EXCEL FILE


clear all
close all
clc

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64');
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab');

%Import flight schedule data from excel file
AllData = readtable('Data.xlsx','ReadRowNames',true,'ReadVariableNames',true);
AllBayComplianceData = importdata('BayComplianceMatrix.xlsx');

% Import distance bay-terminal data from excel file
d = readmatrix('distance.xlsx');
d(:,1) = [];
d(1,:) = [];

%Import color data, for the display. from excel file
color = string(importdata('Colors.xlsx')); %Picking up some colors

% Define number of planes and number of bays used
AN = height(AllData);
data = AllData(1:AN,:);
NBays = size(d,1);
%Define matrix for bay compliance 
BayComplianceData = AllBayComplianceData(1:NBays,:);
%% Organize data for the plane structure - some possible needed corrrections are also included

for i = 1:AN
    %Get the data from the time from the excel, with some corrections it
    %might happen in the schedule
    plane(i).AT = str2num(erase(string(data.Arrival(i)),':'));
    plane(i).DT = str2num(erase(string(data.Departure(i)),':'));
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
    %assign a number of passengers to each plane considering the limits
    plane(i).Passenger= floor((plane(i).Passenger_max - plane(i).Passenger_min)/2) + plane(i).Passenger_min;
    
end

%% ORDERING DATA - order the flights by time of arrival and set them into a matrix

for i = 1:AN
    plane(AN+1).AT = 2500;

    for j = i:AN
        if(plane(j).AT < plane(AN+1).AT)
            plane(AN+1) = plane(j);
            index = j;
        end
        plane(j);
    end
    plane(index) = plane(i);
    plane(i) = plane(AN+1);
end


%% Incorporate the towing time for each plane 
TT=30; %in minutes
for i=1:AN
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

%% Incorporate the buffer time for each plane 
BT=10; %in minutes
for i=1:AN
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
% OV_initial = zeros(AN);
% for i = 1:AN
%     for k = i:AN
%         if plane(i).DT >= plane(k).AT
%             OV_initial(i,k) = 1;
%         end
%     end
% end
% OVc = OV_initial';

%% Creating matrix for the equality constraints
%In cplex, the constraint matrix has to be separated in equality constraint
%matrix and inequality constraint matrix

%Computation of equality constraint matrix: 
%1-plane constraint: there can not be more than one bay for the same plane
Aeq=[];
    for j=1:NBays
         Aeq=[Aeq,eye(AN)];
    end
%2-towing constraint: if a plane is towed, then is has to have a leaving
%bay, meaning X_i+NP,k
T=[];
Aux=[];
for i=1:NBays
      Aux=[Aux,eye(AN)];
    end
for i=1:NBays
    T=[T;Aux];
end
%final matrix
Aeq = [Aeq, zeros(size(Aeq)),zeros(size(Aeq));zeros(AN, AN*NBays),Aux,-Aux];


%%  Creating matrix for inequality the constraints

OVab = zeros(AN, 3*AN*NBays);
OVdb = zeros(AN, 3*AN*NBays);

% filling the overlapping matrix for arriving planes, comparing the
% arrival and departure times with towing (ATT, DTT) and without, but with buffer
% (AT,DT)
for i=1:AN
    for j=1:AN
        if (plane(i).AT>plane(j).AT && plane(i).AT<plane(j).DT)
            OVab (i,j) = 1;
            if (plane(i).AT > plane(j).ATT)
                OVab(i,j+2*AN*NBays) = -1;
                if (plane(i).AT > plane(j).DTT)
                    OVab (i,j+AN*NBays) = 1;
                end
            end
        end
        if i==j
            OVab (i,j) = 1;
        end
    end
end

% Filling the overlapping matrix for departing planes, comparing the
% arrival and departure times with towing (ATT, DTT) and without, but with buffer
% (AT,DT)

for i=1:AN
    for j=1:AN
        if(plane(i).DTT > plane(j).AT && plane(i).DT < plane(j).DT)
            OVdb(i,j) = 1;
            if (plane(i).DTT > plane(j).ATT)
            OVdb(i, j+2*(AN*NBays)) = -1;
            end
            if (plane(j).DT>plane(i).DTT)
            OVdb(i, j+AN*NBays) = 1;
            end
        end
        if (i==j)
            OVdb (i,j+AN*NBays) = 1;
        end
    end
end

OVd = OVdb;
OVa = OVab;

%creating the Aineq matrix - shifting the base (OVab and OVdb) matrices of
%AN position <-> apllying the same constraint for every bay
for i =1:NBays-1
%     AN*(3*NBays-i)
    OVd = [OVd; OVdb(:,AN*(3*NBays-i)+1:end), OVdb(:,1:AN*(3*NBays-i))];
    OVa = [OVa; OVab(:,AN*(3*NBays-i)+1:end), OVab(:,1:AN*(3*NBays-i))];
end

%The mentioned before, plus two new last constraints
%1- similar to the equality, a plane that tows can only tow from the bay to
%which arrived and has to have any but one leaving bay
%2- One plane can only have a leaving bay maximum
OV = [OVa;OVd;...
    eye(AN*NBays), T, -eye(AN*NBays);...
    zeros(AN,AN*NBays), Aux,zeros(AN,AN*NBays)]; %departing bays for each plane <=1
Aineq=OV;

%% MATRIX RELATING PLANE TYPE AND BAY COMPATIBILITY VECTOR
PlaneComplianceMatrix = [];
for i = 1:AN
     PlaneComplianceMatrix = [PlaneComplianceMatrix, BayComplianceData(:,plane(i).Type)];
end
PlaneComplianceVector =[];
for i = 1:NBays
    PlaneComplianceVector = [PlaneComplianceVector,PlaneComplianceMatrix(i,:)];
end
% This right hand side has to be an inequality constraint
% We have to expand the inequality constraint coefficients with a unity
% matrix with the dimensions of the number of decision variables.
% We add RHS_comp as the right hand side.

PlaneComplianceVector = [PlaneComplianceVector, PlaneComplianceVector, ones(1,length(PlaneComplianceVector))];

for i =1:size(Aeq, 1)
   Aeq(i,:) = Aeq(i,:).*PlaneComplianceVector;
end

for i =1:size(Aineq, 1)
   Aineq(i,:) = Aineq(i,:).*PlaneComplianceVector;
end



%% Build the other components to run the cplex (upper and lower bound of decision variables, rhs)

%set the vector of the right side of the inequality constraints
for i=1:size(Aineq, 1)
        rightside_ineq(i,1) = 1;
end

%set the vector of the right side of the equality constraints
rightside_eq = [ones(AN,1);zeros(AN,1)];
%set the vector of the upper bound of the decision variables
%set the vector of the lower bound of the decision variables
%set the vector that states the decision variables are binary or integer
for i=1:AN*NBays*3
        ub(i,1) =1;
        lb(i,1) =0;
        ctype(1,i)='B';
end

%define the coefficients of the decision variables for the objective
%function for distance between bays and gates
f=[];
c=1;
for i=1:NBays
    for j=1:AN
       f(c)=2* d(i,(plane(j).terminal))*plane(j).Passenger;
       f(c+AN*NBays) = f(c)/2;
       f(c+2*AN*NBays) = - f(c)/2 + 100;
       c=c+1;
    end
end

%% Adding Bay Compliance constraints
% Number of Decision variables = Bays*AN
% Aineq = [Aineq;diag(diag(ones(NBays*AN)))];
% rightside_ineq = [rightside_ineq; RHS_comp];

%cplex implementation 
 x=cplexmilp(f, Aineq, rightside_ineq, Aeq, rightside_eq, [],[],[],lb, ub, ctype);

 for i = 1 : length(x)
     if (x(i)  > 0)
         x(i) = 1;
     end
 end
 
 %% Display the result in matrix
arriving= [];
towings = [];
leaving = [];
x=x';
for i = 1:NBays
arriving = [arriving; x(1+AN*(i-1):AN*i)];
towings = [towings; x(1+AN*(i-1)+AN*NBays*2:AN*i+AN*NBays*2)];
leaving = [leaving; x(1+AN*(i-1)+AN*NBays:AN*i+AN*NBays)];
end


 fprintf('Every column corresponds to a plane, every row corresponds to a bay\n\n');
 fprintf('X_i,k matrix \n\n');
 disp(arriving);
 fprintf('V_i,k matrix \n\n');
 disp(towings);
 fprintf('X_i+AN,k matrix \n');
 disp(leaving);

 
%% DATA DISPLAYING
% x is the bay number

%building x_output - a matrix which has the decision variables varying by
%plane on the columns and by bays on the rows
x_output= [];

for i = 1:NBays
 x_output = [x_output; x([(i-1)*AN+1: (i-1)*AN+5])'];
end
% x_output = x_output';
% x_output
a1=[];
for i = 1 : AN
    plane(i).at = [(plane(i).AT-mod(plane(i).AT, 100))/100, mod(plane(i).AT, 100), 0];
    plane(i).att = [(plane(i).ATT-mod(plane(i).ATT, 100))/100, mod(plane(i).ATT, 100), 0];
    plane(i).dt = [(plane(i).DT-mod(plane(i).DT, 100))/100, mod(plane(i).DT, 100), 0];
    plane(i).dtt = [(plane(i).DTT-mod(plane(i).DTT, 100))/100, mod(plane(i).DTT, 100), 0];
    plane(i).Color = color(i);
end


for k = 1:NBays

% b is a vector with the number of the bay, for each plane i have to give a 
% plot command with specify the color of the line

b = [k,k];

for i = 1 : AN
    
    
        if (towings(k,i) == 0 && arriving(k,i)==1)
            a = [ plane(i).at;  plane(i).dt];
            c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
            d=strcat(c(:,1),':',c(:,2),':',c(:,3));
            plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
            hold on;
            datetick('x','HH:MM:SS')  
            
        elseif (towings(k,i) == 1 && arriving(k,i)==1)
            a = [ plane(i).at;  plane(i).att];
            c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
            d=strcat(c(:,1),':',c(:,2),':',c(:,3));
            plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
            hold on;
            datetick('x','HH:MM:SS')  
            
        end
        
        if (leaving(k,i) == 1)
            a = [ plane(i).dtt;  plane(i).dt];
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