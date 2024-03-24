clear;
clc;
close all;
[data_entropy,empty] = xlsread('LUAD_entropy_matrix.xlsx');

pnum=[269,50,71,73,11,26];
data_entropy_size=size(data_entropy);

n=1;
for t=1:length(pnum)
    count=0;
    for c=n:sum(pnum(1:t))   %psize(3)
        count=count+1;
        tmp_entropy=sort(data_entropy(:,c),'descend');
        aver_entropy(t,count)=mean(tmp_entropy(1:floor(data_entropy_size(1)*0.05)));
    end
    case_result(t)=mean(aver_entropy(t,1:pnum(t)));
    n=1+sum(pnum(1:t));
end


% result(1)=mean(case_result(1:2));
% result(2:6)=case_result(3:7);
plot(1:6,case_result,'r-*','LineWidth',3);
ylim([2.8 3.2])






