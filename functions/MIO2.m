%valuto Palpha ad occhi aperti su due elettrodi ed eseguo un t-test
O1 = 19 ;
C3 = 9;
FP1 = 1;
ic=1;
Palpha_m_O1_1 = Palpha_m(O1,ic); %occhi aperti
Palpha_m_C3_1 = Palpha_m(C3,ic); %occhi aperti
Palpha_m_FP1_1 = Palpha_m(FP1,ic); %occhi aperti



[h1,p1]=ttest(Palpha_m_O1_1,Palpha_m_C3_1); % O1 vs. C3
[h2,p2]=ttest(Palpha_m_O1_1,Palpha_m_FP1_1); % O1 vs. FP1