function BC=getBCIndices(BC_nbGhost,A)

BC.topRow=(0:size(A,2)-1)*size(A,1)+1;
BC.topRow_2=(0:size(A,2)-1)*size(A,1)+2;
BC.botRow=(1:size(A,2))*size(A,1);
BC.botRow_2=(1:size(A,2))*size(A,1)-1;
BC.leftRow=1:size(A,1);
BC.leftRow_2=size(A,1)+(1:size(A,1));
BC.rightRow=(size(A,2)-1)*size(A,1)+(1:size(A,1));
BC.rightRow_2=(size(A,2)-2)*size(A,1)+(1:size(A,1));

if BC_nbGhost==2
    BC.topRow_3=(0:size(A,2)-1)*size(A,1)+3;
    BC.botRow_3=(1:size(A,2))*size(A,1)-2;
    BC.leftRow_3=2*size(A,1)+(1:size(A,1));
    BC.rightRow_3=(size(A,2)-3)*size(A,1)+(1:size(A,1));
end