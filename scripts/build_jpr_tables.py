# JPR tablolari — tum sayilar dogrulanmis koşu ciktilarindan
T=[]
T.append(r"""\begin{table}[htbp]
\caption{Site-level quantification deposited for each substrate load, and the sites used here. Union across loads: 44 lysines; 20 pass the confidence threshold in all four loads.}
\label{tbl:coverage}
\begin{tabular}{lrrrl}
\toprule
Substrate load (ng) & Rows deposited & Confidence $>90$ & Used & SHA-256 (first 16) \\
\midrule
250  & 35 & 22 & 22 & \texttt{6d9fe8b86b77d205} \\
500  & 40 & 33 & 33 & \texttt{9b8c7d8c5b78c202} \\
1000 & 44 & 37 & 37 & \texttt{0df9e146fd41e4b5} \\
2000 & 45 & 38 & 38 & \texttt{8224464fdb5606f5} \\
\bottomrule
\end{tabular}
\end{table}""")
T.append(r"""\begin{table}[htbp]
\caption{Label diagnostics. Weights are inverse-variance, median-normalised within a load. Kish effective sample size (ESS) is $(\sum w)^2/\sum w^2$. Truncation counts responses falling outside the unit interval before truncation.}
\label{tbl:labels}
\begin{tabular}{lrrrrrrr}
\toprule
Load (ng) & $n$ & $w_{\min}$ & $w_{\max}$ & Kish ESS & ESS/$n$ & Truncated & Response $=0$ / $=1$ \\
\midrule
250  & 22 & 0.0002 & 28.4 & 5.0  & 0.23 & 2 & 3 / 0 \\
500  & 33 & 0.0048 & 18.7 & 9.7  & 0.30 & 1 & 1 / 1 \\
1000 & 37 & 0.0016 & 18.6 & 11.8 & 0.32 & 5 & 5 / 0 \\
2000 & 38 & 0.0171 & 38.0 & 11.4 & 0.30 & 3 & 4 / 0 \\
\bottomrule
\end{tabular}
\end{table}""")
T.append(r"""\begin{table}[htbp]
\caption{Condition-level results under each weighting scheme. Baseline is the training-fold mean predictor; $\Delta$MAE is the advantage over the site-structure permutation control; $P$ is the empirical permutation $P$ value from 200 permutations, for which 0.005 is the resolution limit. Response dispersion is given for the same site sets.}
\label{tbl:replication}
\small
\begin{tabular}{lrrrrrrrrrr}
\toprule
 & & & \multicolumn{3}{c}{Unweighted (primary)} & \multicolumn{3}{c}{Inverse-variance} & \\
\cmidrule(lr){4-6}\cmidrule(lr){7-9}
Load (ng) & $n$ & Baseline & MAE & $\rho$ & $P$ & MAE & $\rho$ & $P$ & SD (IQR) \\
\midrule
250  & 22 & 0.337 & 0.274 & 0.414 & 0.010 & 0.304 & 0.141 & 0.109 & 0.357 (0.574) \\
500  & 33 & 0.275 & 0.187 & 0.616 & 0.005 & 0.184 & 0.681 & 0.005 & 0.310 (0.509) \\
1000 & 37 & 0.268 & 0.221 & 0.399 & 0.005 & 0.233 & 0.298 & 0.010 & 0.304 (0.533) \\
2000 & 38 & 0.270 & 0.188 & 0.581 & 0.005 & 0.209 & 0.574 & 0.005 & 0.298 (0.492) \\
\bottomrule
\end{tabular}
\par\smallskip
\footnotesize Baseline is shown for the unweighted analysis; the inverse-variance baselines are 0.396, 0.300, 0.302 and 0.290. $\Delta$MAE values are 0.090, 0.102, 0.064 and 0.091 unweighted, and 0.048, 0.100, 0.047 and 0.079 weighted.
\end{table}""")

T.append(r"""\begin{table}[htbp]
\caption{Matched twenty-site core, unweighted, with the random-subset reference. The reference distribution is 200 random twenty-site subsets of the 500 ng condition: mean $\rho = 0.400$ (SD 0.265), median 0.461, 5th and 95th percentiles $-0.135$ and 0.708, with 74\% above 0.3 and a mean baseline improvement of 0.062. The matched core reached $\rho = -0.224$ and a baseline improvement of 0.003, the 3rd percentile of that distribution.}
\label{tbl:core}
\begin{tabular}{lrrrrr}
\toprule
Load (ng) & $n$ & MAE & Baseline MAE & $\rho$ & $P$ \\
\midrule
250  & 20 & 0.326 & 0.347 & \phantom{-}0.164 & 0.080 \\
500  & 20 & 0.257 & 0.260 & $-0.224$ & 0.194 \\
1000 & 20 & 0.226 & 0.221 & \phantom{-}0.134 & 0.264 \\
2000 & 20 & 0.231 & 0.248 & \phantom{-}0.405 & 0.129 \\
\bottomrule
\end{tabular}
\end{table}""")
open("tables_jpr.tex","w").write("\n\n".join(T)+"\n")
print("4 tablo yazildi")
