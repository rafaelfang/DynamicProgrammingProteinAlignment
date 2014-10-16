function [ x,y ] = findMaxValuePositionFrom2DMatrix( Mat )
%find the max value position from 2D matrix
%written by Chao Fang

row=size(Mat,1);
col=size(Mat,2);
maxVal=-9999;
x=1;
y=1;
for i=1:row
    for j=1:col
        if Mat(i,j)>maxVal
            maxVal=Mat(i,j);
            x=i;
            y=j;
        end
    end
end

end

