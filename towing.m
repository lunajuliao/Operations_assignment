%% GETTING DATA FROM EXCEL FILE


clear all
close all
clc

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64');
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab');

AllData = readtable('SampleDataWithAircraftType.xlsx','ReadRowNames',true,'ReadVariableNames',true);
AllBayComplianceData = importdata('BayComplianceMatrix.xlsx');

% Define number of planes and number of bays used

PN = height(AllData); % (The old way.)
data = AllData(1:PN,:);
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

NBays = size(d,1);
BayComplianceData = AllBayComplianceData(1:NBays,:);

color = string(importdata('Colors.xlsx')); %Picking up some colors
%%

for i = 1:PN
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
    plane(i).Passenger= floor((plane(i).Passenger_max - plane(i).Passenger_min)/2) + plane(i).Passenger_min;
    
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
%     p =[p,plane(i).P];
end


%% INCORPORATE THE TOWING TIME FOR EACH PLANE 
TT=30; %in minutes
h = floor(TT/60);
min = mod(TT,60);

for i=1:PN
    hap = floor(plane(i).AT / 100);
    minap = mod(plane(i).AT , 100);
    hdp = floor(plane(i).DT / 100);
    mindp = mod(plane(i).DT , 100);
    plane(i).ATT = 100*(h+hap+floor((minap+min)/60)) + mod (min+minap ,60);
    plane(i).DTT = 100*(hdp-h-double(mindp<min)) + double(mindp>min)*(mindp-min) + double(mindp<min)*(60-min+mindp)  ;
%     plane(i).ATT=plane(i).AT+TT;
%     plane(i).DTT=plane(i).DT-TT;
%     if(mod(plane(i).ATT,100) >=60)%correction of the time, to be presented as hours:minutes
%         minutes=mod(plane(i).ATT,100)-60;
%         plane(i).ATT=floor(plane(i).ATT/100)*100+100+minutes;
%     else
%     end
%     
%     if(mod(plane(i).DTT,100) >=60)%correction of the time, to be presented as hours:minutes
%         minutes=100-mod(plane(i).DTT,100);
%         plane(i).DTT=floor(plane(i).DTT/100)*100+60-minutes;
%     else
%     end

    
end
%% INCORPORATE THE BUFFER TIME FOR EACH PLANE 
BT=0; %in minutes
h = floor(BT/60);
min = mod(BT,60);

for i=1:PN
    
    hap = floor(plane(i).AT / 100);
    minap = mod(plane(i).AT , 100);
    hdp = floor(plane(i).DT / 100);
    mindp = mod(plane(i).DT , 100);
    plane(i).AT = 100*(hap-h-double(minap<min)) + double(minap>min)*(minap-min) + double(minap<min)*(60-min+minap)  ;
    plane(i).DT = 100*(h+hdp+floor((mindp+min)/60)) + mod (min+mindp ,60);
%     if(mod(plane(i).AT,100) >=60)%correction of the time, to be presented as hours:minutes
%         minutes=100-mod(plane(i).AT,100);
%         plane(i).AT=floor(plane(i).AT/100)*100+60-minutes;
%     else
%     end
%     
%     if(mod(plane(i).DT,100) >=60)%correction of the time, to be presented as hours:minutes
%         minutes=mod(plane(i).DT,100)-60;
%         plane(i).DT=floor(plane(i).DT/100)*100+100+minutes;
%     else
%     end

    if 2359 < plane(i).DT
        plane(i).DT = 2359;
    end
    if plane(i).AT < 0
        plane(i).AT = 0;
    end
    
end



%% OVERLAPPING MATRIX OV (overlap)
% % From the time data, we compute a matrix that shows which planes
% % overlap and since they are ordered, we use only the lower triangular(+diagonal) part of the matrix.
OV_initial = zeros(PN);
for i = 1:PN
    for k = i:PN
        if plane(i).DT >= plane(k).AT
            OV_initial(i,k) = 1;
        end
    end
end
OVc = OV_initial';

%% Creating matrix for equality the constraints
%NBays=4; (Moved to top)
%In cplex, the constraint matrix has to be separted in equality constraint
%matrix and inequality constraint matrix

%Computation of equality constraint matrix: 
%1-plane constraint: there can not be more than one bay with the some value
Aeq=[];
    for j=1:NBays
         Aeq=[Aeq,eye(PN)];
    end
    
    T=[];
Aux=[];
for i=1:NBays
      Aux=[Aux,eye(PN)];
    end
for i=1:NBays
    T=[T;Aux];
end
% Aeq = [Aeq, zeros(size(Aeq)),zeros(size(Aeq));zeros(PN*NBays),T,-T];
Aeq = [Aeq, zeros(size(Aeq)),zeros(size(Aeq));zeros(PN, PN*NBays),Aux,-Aux];


%%  Creating matrix for inequality the constraints

OVab = zeros(PN, 3*PN*NBays);
OVdb = zeros(PN, 3*PN*NBays);

% filling the overlapping matrix for arriving planes
for i=1:PN
    for j=1:PN
        if (plane(i).AT>plane(j).AT && plane(i).AT<plane(j).DT)
            OVab (i,j) = 1;
            if (plane(i).AT > plane(j).ATT)
                OVab(i,j+2*PN*NBays) = -1;
                if (plane(i).AT > plane(j).DTT)
                    OVab (i,j+PN*NBays) = 1;
                end
            end
        end
        if i==j
            OVab (i,j) = 1;
        end
    end
end




% filling the overlapping matrix for departing planes

for i=1:PN
    for j=1:PN
        if(plane(i).DTT > plane(j).AT && plane(i).DT < plane(j).DTT)
            OVdb(i,j) = 1;
            if (plane(i).DTT > plane(j).ATT)
            OVdb(i, j+2*(PN*NBays)) = -1;
            end
        elseif (plane(j).DTT<plane(i).DTT && plane(j).DT>plane(i).DTT)
            OVdb(i, j+PN*NBays) = 1;
        end
        if (i==j)
            OVdb (i,j+PN*NBays) = 1;
        end
    end
end

OVd = OVdb;
OVa = OVab;

%creating the Aineq matrix - shifting the base (OVab and OVdb) matrices of
%PN position <-> apllying the same constraint for every bay
for i =1:NBays-1
%     PN*(3*NBays-i)
    OVd = [OVd; OVdb(:,PN*(3*NBays-i)+1:end), OVdb(:,1:PN*(3*NBays-i))];
    OVa = [OVa; OVab(:,PN*(3*NBays-i)+1:end), OVab(:,1:PN*(3*NBays-i))];
end



OV = [OVa;OVd;...
    eye(PN*NBays), T, -eye(PN*NBays);...
    zeros(PN,PN*NBays), Aux,zeros(PN,PN*NBays)]; %departing bays for each plane <=1
Aineq=OV;

%% MATRIX RELATING PLANE TYPE AND BAY COMPATIBILITY VECTOR
PlaneComplianceMatrix = [];
for i = 1:PN
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



%%
%Computation of inequality constraint matrix:
%1-bay and plane constraint: whenever there is an extra plane, we rewrite a
%constraint concerning the ground existing planes which have to be assigned to the
%bays (move or stay in the same place)


%set the vector of the right side of the inequality constraints
for i=1:size(Aineq, 1)
        rightside_ineq(i,1) = 1;
end

%set the vector of the right side of the equality constraints
% for i=1:size(Aeq, 1)
%         rightside_eq(i,1) = 1;
% end
rightside_eq = [ones(PN,1);zeros(PN,1)];
%set the vector of the upper bound of the decision variables
%set the vector of the lower bound of the decision variables
%set the vector that states the decision variables are binary or integer
for i=1:PN*NBays*3
        ub(i,1) =1;
        lb(i,1) =0;
        ctype(1,i)='B';
end
%incorporate distance matrix that tells us the distance from a fixed
%terminal, pre-assigned before to every plane.

% d=readmatrix('distance.xlsx'); (moved to top)
% d(:,1) = [];

f=[];
c=1;
% set the vector of the coefficients of the objective function for distance
% between bays and gates
for i=1:NBays
    for j=1:PN
       f(c)=2* d(i,(plane(j).terminal))*plane(j).Passenger;
       f(c+PN*NBays) = f(c)/2;
       f(c+2*PN*NBays) = - f(c)/2 + 100;
       c=c+1;
    end
end

%% Adding Bay Compliance constraints
% Number of Decision variables = Bays*PN
% Aineq = [Aineq;diag(diag(ones(NBays*PN)))];
% rightside_ineq = [rightside_ineq; RHS_comp];

%cplex implementation 
 x=cplexmilp(f, Aineq, rightside_ineq, Aeq, rightside_eq, [],[],[],lb, ub, ctype);
%%
arriving=[];
towings = [];
leaving = [];
x=x';
for i = 1:NBays
arriving = [arriving; x(1+PN*(i-1):PN*i)];
towings = [towings; x(1+PN*(i-1)+PN*NBays*2:PN*i+PN*NBays*2)];
leaving = [leaving; x(1+PN*(i-1)+PN*NBays:PN*i+PN*NBays)];


end


 fprintf('Every column corresponds to a plane, every row corresponds to a bay\n\n');
 fprintf('X_i,k matrix \n\n');
 disp(arriving);
 fprintf('V_i,k matrix \n\n');
 disp(towings);
 fprintf('X_i+PN,k matrix \n');
 disp(leaving);

 
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
    plane(i).at = [(plane(i).AT-mod(plane(i).AT, 100))/100, mod(plane(i).AT, 100), 0];
    plane(i).att = [(plane(i).ATT-mod(plane(i).ATT, 100))/100, mod(plane(i).ATT, 100), 0];
    plane(i).dt = [(plane(i).DT-mod(plane(i).DT, 100))/100, mod(plane(i).DT, 100), 0];
    plane(i).dtt = [(plane(i).DTT-mod(plane(i).DTT, 100))/100, mod(plane(i).DTT, 100), 0];
    plane(i).Color = color(i);
end


for k = 1:NBays

% b is a vector with the number of the bay, for each plane i have to give a 
%     plot command with specify the color of the line

b = [k,k];

for i = 1 : PN
    
    if arriving(k,i)==1
        if (towings(k,i) == 0)
            a = [ plane(i).at;  plane(i).dt];
            c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
            d=strcat(c(:,1),':',c(:,2),':',c(:,3));
            plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
            hold on;
            datetick('x','HH:MM:SS')  
            
        elseif (towings(k,i) == 1)
            a = [ plane(i).at;  plane(i).att];
            c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
            d=strcat(c(:,1),':',c(:,2),':',c(:,3));
            plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
            hold on;
            datetick('x','HH:MM:SS')  
            
        else
            disp(towings(k,i));
            error('towings (i,k) value is not binary');
        end
        
        if (leaving(k,i) == 1)
            a = [ plane(i).dtt;  plane(i).dt];
            c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
            d=strcat(c(:,1),':',c(:,2),':',c(:,3));
            plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
            hold on;
            datetick('x','HH:MM:SS')  
            
        end
        
        
%         c=cellfun(@(x) num2str(x,'%02d'),num2cell(a),'UniformOutput',false);
%         d=strcat(c(:,1),':',c(:,2),':',c(:,3));
%         plot(datenum(d,'HH:MM:SS'),b,'Linewidth', 6,'Color',plane(i).Color);            
%         hold on;
%         datetick('x','HH:MM:SS')  

        
     end
    
    
end
 
end
grid on
ylim([0, NBays+1]);