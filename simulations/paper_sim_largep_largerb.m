clear
close all

%%%%%experiment controls
sim_ctrl.nsim = 200;
sel_ctrl.methods = {'AIC', 'AIC_MA', 'BIC', 'BIC_MA', 'EBIC', 'SF_PIC_fractional','SCV_none2', 'SCV_none10', 'SCV_none5', 'SCV_plugin', 'SCV_fractional'};
type_ctrl = struct;
verbose = 1;
%%%%%

%%%%%%%large p, larger b
%%%information for writing to spreadsheet
filename = 'paper_sim_largep_largerb_res.xls';
res_range = 'A4';

largep_n = 30;
largep_p = 100;
largep_J = 15;
largep_m = 10;
largep_r = 2;
largep_b = [2 3 4];
largep_b1 = largep_b(1);
largep_b2 = largep_b(2);
largep_b3 = largep_b(3);
largep_nsim = sim_ctrl.nsim;

largep_meta1 = table(largep_n, largep_p, largep_J, largep_m, largep_r, largep_b1, largep_nsim);
largep_meta2 = table(largep_n, largep_p, largep_J, largep_m, largep_r, largep_b2, largep_nsim);
largep_meta3 = table(largep_n, largep_p, largep_J, largep_m, largep_r, largep_b3, largep_nsim);
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

writetable(largep_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b2_mildcorr', 'Range', 'A1')
writetable(bigp_mild_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b2_mildcorr', 'Range', res_range)
writetable(largep_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b3_mildcorr', 'Range', 'A1')
writetable(bigp_mild_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b3_mildcorr', 'Range', res_range)
writetable(largep_meta3, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b4_mildcorr', 'Range', 'A1')
writetable(bigp_mild_test_errs{3}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b4_mildcorr', 'Range', res_range)
%%%med corr
size_ctrl.dataName = 'mynorm-medcorr';

bigp_med_test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);

writetable(largep_meta1, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b2_medcorr', 'Range', 'A1')
writetable(bigp_med_test_errs{1}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b2_medcorr', 'Range', res_range)
writetable(largep_meta2, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b3_medcorr', 'Range', 'A1')
writetable(bigp_med_test_errs{2}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b3_medcorr', 'Range', res_range)
writetable(largep_meta3, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'Sheet', 'largep_b4_medcorr', 'Range', 'A1')
writetable(bigp_med_test_errs{3}, filename, 'FileType', 'spreadsheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', 'largep_b4_medcorr', 'Range', res_range)
%%%%%%%