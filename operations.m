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
end

%% ORDERING DATA
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
Nflights=PN;
Bays=2;
Gates=2;

%% Creating matrix C for the constraints

for i=1:Nflights
    for j=1:Bays*Gates*Nflights
        if ((1+floor((j-1)/(Bays*Gates)))==i)
            C(i,j)=1;
        else
            C(i,j)=0;
        end
    end
end

for i=Nflights+1:(Nflights+Gates)
    aux=0;
    for j=1:Bays*Gates*Nflights
       if (j==i-Nflights) || (j==(i-Nflights)+(Gates*aux))
            C(i,j)=1;
            aux=aux+1;
       else
            C(i,j)=0;
       end
    end
end

aux1=0;
for i=Nflights+Gates+1:(Nflights+Gates+Bays)
    aux2=0;
    for j=1:Bays*Gates*Nflights-2
       if (j==i-Nflights-Gates) || (j==(i-Nflights-Gates)+(Bays*Gates*aux2))
            C(i,[j+(aux1*(Gates - 1)) j+(aux1*(Gates - 1))+1])=[1 1];
            j=j+1;
            aux2=aux2+1;
       else
            C(i,j+(aux1*(Gates - 1))+1)=0;
       end
    end
    aux1=aux1+1;
end

for i=1:Nflights+Gates+Bays
        rightside(i,1) =1;
end

for i=1:Nflights*Gates*Bays
        ub(i,1) =1;
end

for i=1:Nflights*Gates*Bays
        lb(i,1) =0;
end

d=[0 3 5 6];
f=[d d];

for i=1:Nflights*Gates*Bays
    ctype(1,i)='B';
end


Aeq=C([1:2],:);
Aineq=C([3:6],:);
bineq=rightside([3:6]);
beq=rightside([1:2]);

x=cplexmilp(f, Aineq, bineq, Aeq, beq, [],[],[],lb, ub, ctype);
