function  [ P, Q,  phi] = locate(A,I)
%   Locte a set of points on a grid
%   
%   Input arguments:
%   A: a set of points, could be a vector or an array up to three
%   dimensions.
%   I: a strictly increasing vector of values
%
%   Output values:
%   P: positions of the set of points, A, on I


% 
% A=0:0.3:15;
% A=A'
% I=1:10
% I=I'
% p=size(I,1);


[dim1,dim2]=size(A);
L=length(I); % number of points
AA=repmat(A,[1,1,L]);
II=zeros(1,1,L);
II(1,1,:)=I;
III=repmat(II,dim1,dim2);  
P=sum((AA>=III)==1,3);
if P==0
    P=1;
end
Q=length(I)+1-sum((AA<=III)==1,3);
if Q>L
    Q=P;
end

if P<Q
phi=(I(Q)-A)/(I(Q)-I(P));
else
phi=1;
end

end

