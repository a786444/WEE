clear 
clc
filename='input.txt';
[maxf,incrementf,u,diffheight,gap_height,hubheight,Umeanwind]=ReadInput(filename);
%% construct coh
coh=cell(diffheight-1);
f_array=(0:incrementf:maxf);
for i=1:diffheight-1
    for j=1:diffheight-1
        if (i+j)*gap_height/2<60
            Lk=0.7*(i+j)*gap_height/2;
        else
            Lk=42;
        end
        U_hub=0.5*(u(i)+u(j));
        r=abs(i-j)*gap_height;
        coh{i,j}=exp(-12*sqrt((r.*f_array/U_hub).^2+(0.12*r/(Lk))^2));
    end
end
%plot(f_array,coh{1,9})
%% PSD
PSD=cell(diffheight-1,1);
for i=1:diffheight-1
    sigma_u=0.14*(0.75*u(i)+5.6);
    if (i)*gap_height<60
        Lk=0.7*(i)*gap_height;
    else
        Lk=42;
    end
    PSD{i}=(4*(sigma_u^2)*Lk*(u(i)^(2/3)))./(u(i)+6.*f_array*Lk).^(5/3);
end
loglog(f_array,PSD{1,:})
hold on
PSD=cat(3,PSD{:,1});
SM=cell(diffheight-1);
for i=1:diffheight-1
    for j=1:diffheight-1
        for k=1:length(f_array)
            SM{i,j}(k)=sqrt(PSD(1,k,i)*PSD(1,k,j))*coh{i,j}(k);
        end
    end
end
SM=cat(3,SM{1:end,1:end}); 
H=cell(maxf/incrementf+1,1);
for i=1:length(f_array)
    H{i}=chol(reshape(SM(1,i,:),[diffheight-1,diffheight-1]),'lower');
end
%% Stochastic process
X=zeros(diffheight-1,length(f_array));
for i=1:length(f_array)
    phase=exp(1i*( 2*pi.*rand(diffheight-1,1)));
    X(:,i)=H{i}*diag(phase)*ones(diffheight-1,1)*sqrt(2*incrementf);
end
%% inverse
steptime=0.05;
X=[X flip(X,2)];
X(:,length(f_array))=[];
random_u=real(ifft(X,length(X),2))*length(f_array)*2^0.5;
U=u+random_u;
% plot((0:steptime:600),U(1,:),'b')
[p_value,p_fre]=pwelch(U(1,:),2^10,[],2^10,20);
loglog(p_fre,p_value)