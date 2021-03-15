clear,clc,close all
rng(1);
%problem generation
agN=3;
agM=4;
probdata=progen(agN,agM);
% min{z1...zN} sum{1,N}(ci'*zi)
% subj: 	sum{1,N} (Hi*zi)=b == sum{1,N}(Hi*zi-bi)=0
%       	zi € Pi, i€{1...N}

% Hi*zi-bi==gi(zi)

%we still need to implement the Pi containing the limits on Dizi<=di
%assignment


b=probdata.b;c=probdata.c;d=probdata.d;
D=probdata.D;H=probdata.H;
LB=probdata.LB;UB=probdata.UB;

%call dual subgradient

optcost=dualsub(b,c,d,D,H,LB,UB);

%results