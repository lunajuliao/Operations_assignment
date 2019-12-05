clc
clear all

F = readtable ('flight data.xlsx');
D = readtable ('distance.xlsx');

Nflights=2;
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