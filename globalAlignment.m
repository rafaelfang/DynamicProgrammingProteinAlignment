close all
clear
clc


load seq1;
load seq2;

% seq1='ATAGAAT';
% seq2='AGATCAGAAATG';
load blosum62;
load blosum62Letters;

col=size(seq1,2);
row=size(seq2,2);
P=seq1;
Q=seq2;
M=zeros(row+1,col+1);

M(1,2:end)=-5:-5:(-5*col);
M(2:end,1)=-5:-5:(-5*row);
graph=cell(row+1,col+1);
for i=2:(row+1)
    for j=2:(col+1)
        A=M(i-1,j-1)+findSValue( Q(i-1),P(j-1),blosum62Letters,blosum62 );
        B=M(i,j-1)-5;
        C=M(i-1,j)-5;
        M(i,j)=findMaxOutofThree(A,B,C);
        if(M(i,j)==A)
            graph{i,j}='\';
        elseif(M(i,j)==B)
            graph{i,j}='|';
        else
            graph{i,j}='-';
        end
    end
end
