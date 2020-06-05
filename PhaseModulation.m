folder='Z:\RujiaChen\Results\';
cueTy='endo';
monkey='Mikey';
Unit='VisualResp';   %_unSorted

load([folder monkey 'ArrayLFP_' cueTy '_' Unit '.mat']);
% load([folder  monkey  'ArraySpk_' cueTy '_' Unit '.mat']);
load([folder  monkey 'CueLFP_' cueTy '_' Unit '.mat']);
% load([folder  monkey 'CueSpk_' cueTy '_' Unit '.mat']);









%%
%3-10Hz=2cycles; 11-14 Hz=3cycles; 15-20Hz=4cycles; >20Hz = 5cycles;
LB=-5;
UB=5;
N=100;
FB=2;
FC=6;

for FC=3:80
    if FC<=10
        LB=-1000/FC;
        UB=1000/FC;
    elseif FC<=14
        LB=-1000/FC*1.5;
        UB=1000/FC*1.5;
    elseif FC<=20
        LB=-1000/FC*2;
        UB=1000/FC*2;
    else
        LB=-1000/FC*2.5;
        UB=1000/FC*2.5;
    end
    [PSI,X] = cmorwavf(LB,UB,N,FB,FC);
    
    
end
        
[PSI,X] = cmorwavf(LB,UB,N,FB,FC);
