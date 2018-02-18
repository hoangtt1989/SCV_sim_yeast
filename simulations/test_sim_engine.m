clear
close all


size_ctrl.n_true = 100;
size_ctrl.p_true = 15;
size_ctrl.J_true = 8;
size_ctrl.m_true = 5;
size_ctrl.r_true = 2;
size_ctrl.b_strength = 1;


sim_ctrl.nsim = 3;
sel_ctrl.nfolds = 5;

sel_ctrl.methods = {'AIC', 'BIC', 'SF_PIC_fractional', 'SCV_plugin', 'SCV_none2'};

type_ctrl = struct;

conv_ctrl.card_upper_bnd = 20;
conv_ctrl.rank_upper_bnd = 5;

% [errs, ranks, cards, Ms, FAs] = SRR_sim_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, 0);


size_ctrl.b_vals = [.5 .7];

test_errs = SRR_b_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, 0);

test_errs