%% estraggo Ptot_m(i_el,ic) per O1 nelle due condizioni
i_el=18 % O1
Ptot_m_O1_1 = Ptot_m(i_el,1); %occhi aperti
Ptot_m_O1_2 = Ptot_m(i_el,2); %occhi chiusi


%% estraggo Palpha_m(i_el,ic) per O1 nelle due condizioni
Palpha_m_O1_1 = Palpha_m(i_el,1); %occhi aperti
Palpha_m_O1_2 = Palpha_m(i_el,2); %occhi chiusi


%% estraggo Palpha_n_m(i_el,ic) per O1 nelle due condizioni 
Palpha_n_m_O1_1 = Palpha_n_m(i_el,1)
Palpha_n_m_O1_2 = Palpha_n_m(i_el,2)

disp(['elettrodo O1'])
disp(['occhi aperti'])
disp(['Ptot_m=' num2str(Ptot_m_O1_1)])
disp(['Palpha_m=' num2str(Palpha_m_O1_1)])
disp(['Palpha_n_m=' num2str(Palpha_n_m_O1_1)])
disp(['occhi chiusi'])
disp(['Ptot_m=' num2str(Ptot_m_O1_2)])
disp(['Palpha_m=' num2str(Palpha_m_O1_2)])
disp(['Palpha_n_m=' num2str(Palpha_n_m_O1_2)])