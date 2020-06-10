function varargout = randpoi(S,lambda)



% How many expected
N = poissrnd(lambda * info(S,'totallength'));
varargout{:} = randlocs(S,N);



