clear
close all

%%%%%experiment controls
sim_ctrl.nsim = 200;
sel_ctrl.methods = {'AIC', 'BIC', 'EBIC', 'SF_PIC_fractional','SCV_none2', 'SCV_none10', 'SCV_none5', 'SCV_plugin', 'SCV_fractional'};
% sel_ctrl.methods = {'AIC', 'AIC_MA', 'BIC', 'BIC_MA', 'jackknife_MA', 'EBIC', 'SF_PIC_fractional','SCV_none2', 'SCV_none10', 'SCV_none5', 'SCV_plugin', 'SCV_fractional'};
type_ctrl = struct;
verbose = 1;
%%%%%

filename = 'paper_sim_res_noMA.xls';
res_range = 'A4';

%%%%%%%large n
%%%information for writing to spreadsheet
largen_n = 100;
largen_p = 60;
largen_J = 30;
largen_m = 15;
largen_r = 5;
largen_b = [.1 .5];
largen_b1 = largen_b(1);
largen_b2 = largen_b(2);
largen_nsim = sim_ctrl.nsim;

largen_meta1 = table(largen_n, largen_p, largen_J, largen_m, largen_r, largen_b1, largen_nsim);
largen_meta2 = table(largen_n, largen_p, largen_J, largen_m, largen_r, largen_b2, largen_nsim);
%%%

size_ctrl.n_true = largen_n;
size_ctrl.p_true = largen_p;
size_ctrl.J_true = largen_J;
size_ctrl.m_true = largen_m;
size_ctrl.r_true = largen_r;
size_ctrl.b_vals = largen_b;

conv_ctrl.card_upper_bnd = 50;
conv_ctrl.rank_upper_bnd = 10;

%%%mild corr
size_ctrl.dataName = 'mynorm-mildcorr';

bign_mild_test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);

writetable(largen_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largen_smallb_mildcorr', 'Range', 'A1')
writetable(bign_mild_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largen_smallb_mildcorr', 'Range', res_range)
writetable(largen_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largen_largeb_mildcorr', 'Range', 'A1')
writetable(bign_mild_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largen_largeb_mildcorr', 'Range', res_range)

%%%med corr
size_ctrl.dataName = 'mynorm-medcorr';

bign_med_test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);

writetable(largen_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largen_smallb_medcorr', 'Range', 'A1')
writetable(bign_med_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largen_smallb_medcorr', 'Range', res_range)
writetable(largen_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largen_largeb_medcorr', 'Range', 'A1')
writetable(bign_med_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largen_largeb_medcorr', 'Range', res_range)
%%%%%%%


%%%%%%%large p
%%%information for writing to spreadsheet
largep_n = 30;
largep_p = 100;
largep_J = 15;
largep_m = 10;
largep_r = 2;
largep_b = [.2 1];
largep_b1 = largep_b(1);
largep_b2 = largep_b(2);
largep_nsim = sim_ctrl.nsim;

largep_meta1 = table(largep_n, largep_p, largep_J, largep_m, largep_r, largep_b1, largep_nsim);
largep_meta2 = table(largep_n, largep_p, largep_J, largep_m, largep_r, largep_b2, largep_nsim);
%%%

size_ctrl.n_true = largep_n;
size_ctrl.p_true = largep_p;
size_ctrl.J_true = largep_J;
size_ctrl.m_true = largep_m;
size_ctrl.r_true = largep_r;
size_ctrl.b_vals = largep_b;

conv_ctrl.card_upper_bnd = 25;
conv_ctrl.rank_upper_bnd = 4;

%%%mild corr
size_ctrl.dataName = 'mynorm-mildcorr';

bigp_mild_test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);

writetable(largep_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_smallb_mildcorr', 'Range', 'A1')
writetable(bigp_mild_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_smallb_mildcorr', 'Range', res_range)
writetable(largep_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_largeb_mildcorr', 'Range', 'A1')
writetable(bigp_mild_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_largeb_mildcorr', 'Range', res_range)
%%%med corr
size_ctrl.dataName = 'mynorm-medcorr';

bigp_med_test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);

writetable(largep_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_smallb_medcorr', 'Range', 'A1')
writetable(bigp_med_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_smallb_medcorr', 'Range', res_range)
writetable(largep_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_largeb_medcorr', 'Range', 'A1')
writetable(bigp_med_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_largeb_medcorr', 'Range', res_range)
%%%%%%%
