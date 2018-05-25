function [LR_MET LR_TT] = Johansen_Test(dataset, nlags)
    
    
    % Comment. With T' I mean T properly adjusted for nlags. (T'= T-nlags-1)
    
    T               = size(dataset,1); % time periods
    nvar            = size(dataset,2); % number of variables
    diff_dataset    = diff(dataset); % taking the first difference of the whole dataset
    
    % Taking the lag on the dataset
    for ilag = 1:nlags+1
        eval(['data_', num2str(ilag), ' = diff_dataset(2+nlags-ilag:end+1-ilag,:);'])
    end
    
    % Building the matrix of regressor (which are the lags of the left hand side)
    reg = zeros(T-nlags-1,nlags*nvar);
    for ilag = 2:nlags+1
        for ivar = 1:nvar
            eval(['reg(:,ivar+size(diff_dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
        end
    end
    
    % Defining objects - For now without constant. Notations is standard
    
    % Run the following regression: dY = Pi*Y + Gam*X + U
    % where dY is (nvar,T'), ? is (nvar,nvar), Y is (nvar,T'),
    % Gam is (nvar,nvar*nlags), X is (nvar*nlags,T') and U is (nvar,T')
    Ylag       = dataset(nlags+1:end-1,:)'; %Y(t-1) level of the right-hand side.
    % I am removing the last obs from Y since it is Y(t-1) which implies that
    % Y(T) will never be used.
    X       = reg'; %dT(t-j) first diff of the right-hand side
    dY      = diff_dataset(nlags+1:end,:)'; %left-hand side - Should be equal to data_1
    % dY should be equal to data_1. Otherwise there is something wrong
    M       = eye(T-nlags-1,T-nlags-1) - X'*(X*X')^(-1)*X; %(T',T')
    S00     = (dY*M*dY')/T; %(nvar,navr)
    S01     = (dY*M*Ylag')/T;  %(nvar,navr)
    S11     = (Ylag*M*Ylag')/T;   %(nvar,navr)
    
    % e = eig(A,B) returns a column vector containing the ...
    % generalized eigenvalues of square matrices A and B.
    % det(lambda*B - A)
    A = S01'*(S00^(-1))*S01;
    B = S11;
    eigv = eig(A,B);
    sorted_eig = sort(eigv,'descend');
    
    %Maximum Eigenvalue Test
    %The null hypothesis is that rank(Pi) = i_eigen - 1
    %The alternative hypothesis is that rank(?) = i_eigen
    LR_MET = zeros(nvar,1);
    for i_eigen = 1:nvar %loop over different H0
        LR_MET(i_eigen) = - T * log(1 - sorted_eig(i_eigen));
    end
    
    %Trace Test
    %The null hypothesis is that rank(Pi) = i_eigen - 1
    %The alternative hypothesis is that i_eigen - 1 < rank(Pi) =< nvar,
    %where nvar is also the maximum number of possible cointegrating vectors
    vec_test = zeros(nvar,nvar);
    for i_eigen = 1:nvar %loop over different H0
        for i_test = i_eigen:nvar
            vec_test(i_eigen,i_test) = log(1 - sorted_eig(i_test));
        end
    end
    LR_TT = - T * sum(vec_test,2);
    
    %Under the hypothesis that there are r cointegrating vectors
    %B is a (nvar - r)-dimensional Brownian motion with covariance matrix I.
    
    %This is table1 for Johansen (1988)
    %Raws represnt the dimension of the brownian motion
    %Columns  2.5%   5%   10%   50%   90%   95%   97.5%
    table1 = [0.0   0.0   0.0   0.6   2.9   4.2   5.3
              1.6   1.9   2.5   5.4   10.3  12.0  13.9
              7.0   7.8   8.8   14.0  21.2  23.8  26.1
              16.0  17.4  19.2  26.3  35.6  38.6  41.2
              28.3  30.4  32.8  42.1  53.6  57.2  60.3]; 
    
    i_eigen = 2
    while table1(i_eigen-1,2) < LR_TT(i_eigen)
        i_eigen = i_eigen + 1
    end
    r = nvar - i_eigen + 1;
    statement = ['Trace Test suggests r = ', num2str(r)]
    disp(statement)
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
