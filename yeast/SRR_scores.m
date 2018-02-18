function [ C1, C2 ] = SRR_scores( perf_ctrl )
%%%after getting the output from SRR_tune_perf or SRR_record_tuning,
%%%compute the scores for non-zero predictors

[U, S, V] = svd(perf_ctrl.BOpt);
U2 = U(logical(perf_ctrl.est_patt), :);
C1 = U2 * S(:, diag(S) >= 1e-4);
C2 = V(:, diag(S) >= 1e-4);

end

