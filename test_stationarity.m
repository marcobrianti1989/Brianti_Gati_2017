function [flag] = test_stationarity(beta)
    %the opnly argument is the result of the vector autoregression with the
    %constant. Inside the function we directly remove the constant from beta
    %Dimension of beta is (nvar , nvar*nlags+1)
    % Beware: need to input beta' if dim(beta) = (nvar*nlags + 1, nvar)
    
    %Removing the constant
    beta = beta(:,2:end);
    
    %Extracting nvar and nlags from beta
    nvar     = size(beta,1);
    nlags    = size(beta,2)/nvar;
    
    %Creating the companion matrix
    comp_matrix                               = zeros(nvar*nlags,nvar*nlags);
    comp_matrix(1:nvar,1:nvar*nlags)          = beta;
    comp_matrix(nvar+1:end,1:nvar*(nlags-1))  = eye(nvar*(nlags-1));
    
    %Evaluating the eigenvalues
    eigens = eig(comp_matrix);
    
    %Test if it exists one eigenvalues outside the unit circle
    %zplane(eigens) %hahahahah this is cool!
    distance = abs(eigens);
    loc = distance >= 1 ;
    test = sum(loc);
    flag = 0;
    if test >= 1
        warning('The VAR is not stationary')
        flag = 1;
    end
    
end