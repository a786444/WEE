function [maxf,incrementf,u,diffheight,gap_height,hubheight,Umeanwind]=ReadInput(filename)
input=textscan(fopen(filename),'%s');
totaltime=str2double(input{1}{find((strcmp(input{1,1},'totaltime:'))==true)+1});
timestep=str2double(input{1}{find((strcmp(input{1,1},'timestep:'))==true)+1});
maxf=1/timestep/2;
incrementf=1/totaltime;
%% meanwind of different height from z=10
Umeanwind=str2double(input{1}{find((strcmp(input{1,1},'Umeanwind:'))==true)+1});
hubheight=str2double(input{1}{find((strcmp(input{1,1},'hubheight:'))==true)+1});
alpha=str2double(input{1}{find((strcmp(input{1,1},'alpha:'))==true)+1});
diffheight=str2double(input{1}{find((strcmp(input{1,1},'diffheight:'))==true)+1});
gap_height=(hubheight/(diffheight-1));
u=zeros(diffheight,1);
for i=1:diffheight
    u(i)=((i-1)*gap_height/hubheight)^(alpha)*Umeanwind;
    %z(i)=((i-1)*(hubheight/(diffheight-1)));
end
u(1)=[];
%plot(u,z) %plot different meanwind
