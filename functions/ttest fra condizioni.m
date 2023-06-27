%% estraggo Ptot_m(i_el,ic) per O1 nelle due condizioni
i_el=18 % O1
Ptot_m_O1_1 = Ptot_m(i_el,1) %occhi aperti
Ptot_m_O1_2 = Ptot_m(i_el,2) %occhi chiusi

disp(['aperti: Ptot_m_O1_1=' num2str(Ptot_m_O1_1)])
disp(['chiusi: Ptot_m_O1_2=' num2str(Ptot_m_O1_2)])
%% estraggo Palpha_m(i_el,ic) per O1 nelle due condizioni
Palpha_m_O1_1 = Palpha_m(i_el,1) %occhi aperti
Palpha_m_O1_2 = Palpha_m(i_el,2) %occhi chiusi

disp(['aperti: Palpha_m_O1_1=' num2str(Palpha_m_O1_1)])
disp(['chiusi: Palpha_m_O1_2=' num2str(Palpha_m_O1_2)])



