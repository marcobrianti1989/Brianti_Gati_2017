%Forni&Gambetti Orthogonality Test
filename_PC       = 'Dataset_test_PC';
sheet_PC          = 'Quarterly';
range_PC          = 'B2:DC287';
first_n_PCs       = 10;
[pvalue_news_shock, pvalue_IT_shock] = ...
      Forni_Gambetti_orthogonality_test(filename_PC,...
      sheet_PC,range_PC,first_n_PCs,A,gam_opt,res)
