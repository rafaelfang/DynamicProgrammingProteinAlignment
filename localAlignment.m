close all
clear
clc
%written by Chao Fang

% load seq1;
% load seq2;

seq1='GAATTCAGTTA';
seq2='GGATCGA';
load blosum62;
load blosum62Letters;

M=size(seq1,2);
N=size(seq2,2);


% A simple scoring scheme is assumed where
% Si,j = 1 if the residue at position i of sequence #1 
% is the same as the residue at position j of sequence
% #2 (match score); otherwise
% Si,j = 0 (mismatch score)
% w (gap penalty)
w = -5;

% 1. Initialization Step
matrix=zeros(N+1,M+1);


matrix(1,2:end)=w:w:(w*M);
matrix(2:end,1)=w:w:(w*N);

% 2. Matrix Fill Step

for i=2:(N+1)
    for j=2:(M+1)
        A=matrix(i-1,j-1)+findSValue( seq2(i-1),seq1(j-1),blosum62Letters,blosum62 );
        B=matrix(i,j-1)+w;
        C=matrix(i-1,j)+w;
        D=0;
        matrix(i,j)=getMax(A,B,C,D);
      
    end
end


% 3.Trace back

[ i,j ] = findMaxValuePositionFrom2DMatrix( matrix );
stack=[];
while(1)
    termA=matrix(i-1,j-1)+findSValue( seq2(i-1),seq1(j-1),blosum62Letters,blosum62 );
    termB=matrix(i,j-1)+w;
    termC=matrix(i-1,j)+w;
    if(termA==matrix(i,j))
        i=i-1;
        j=j-1;
        stack=[strcat(seq2(i),':',seq1(j));stack];
    elseif(termB==matrix(i,j))
        j=j-1;
        stack=[strcat('-',':',seq1(j));stack];
    else
        i=i-1;
        stack=[strcat(seq2(i),':','-');stack];
    end
        
    if(matrix(i,j)==0)
        break;
    end
end






