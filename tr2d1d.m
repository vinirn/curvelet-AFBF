function x1d=tr2d1d(x2d)

%ZIGZAG
x1d=[];
for lin=1:size(x2d,1),
    if mod(lin,2)==1,
        x1d=[x1d x2d(lin,:)];
    else
        x1d=[x1d x2d(lin,end:-1:1)];
    end
end