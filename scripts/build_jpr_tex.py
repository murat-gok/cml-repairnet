import re, pathlib
U="/mnt/user-data/uploads/RepairNet/docs/"
ms=pathlib.Path(U+"submission_JPR/MANUSCRIPT_JPR_TR.md").read_text(encoding="utf-8")
me=pathlib.Path(U+"revision_A1/A3_METHODS_PROTEOMICS_TR.md").read_text(encoding="utf-8")

def esc(t):
    t=t.replace("\\","\\textbackslash{}")
    for a,b in [("&","\\&"),("%","\\%"),("$","\\$"),("#","\\#"),("_","\\_")]:
        t=t.replace(a,b)
    t=t.replace("~","\\textasciitilde{}")
    # kesme isaretleri ve tireler
    t=t.replace("“","``").replace("”","''").replace("’","'")
    return t

MATH={
 "N-epsilon-carboxymethyllysine": r"$N^{\\varepsilon}$-carboxymethyllysine",
 "N-epsilon-carboxymethyl": r"$N^{\\varepsilon}$-carboxymethyl",
 "var(A)/C^2 + A^2 var(C)/C^4": r"$\\mathrm{var}(A)/C^{2} + A^{2}\\,\\mathrm{var}(C)/C^{4}$",
}
def inline(t):
    holes=[]
    for k,v in MATH.items():
        if k in t:
            t=t.replace(k, "\x00%d\x00"%len(holes)); holes.append(v)
    # backtick kod parcalari -> \texttt, kirilabilir
    codes=[]
    def _c(m):
        codes.append(m.group(1)); return "\x01%d\x01"%(len(codes)-1)
    t=re.sub(r"`([^`]+)`", _c, t)
    t=esc(t)
    t=t.replace("^","\\textasciicircum{}")
    t=re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
    for i,v in enumerate(holes): t=t.replace("\x00%d\x00"%i, v)
    for i,c in enumerate(codes):
        cc=esc(c).replace("\\_", "\\_\\allowbreak ").replace("{", "\\{").replace("}", "\\}")
        t=t.replace("\x01%d\x01"%i, r"\texttt{%s}"%cc)
    return t

def paras(block):
    out=[]
    for p in re.split(r"\n\s*\n", block.strip()):
        p=" ".join(l.strip() for l in p.splitlines() if l.strip())
        if not p: continue
        if p.startswith("[") and p.endswith("]"): continue   # TR yorum satiri
        out.append(inline(p))
    return "\n\n".join(out)

def grab(txt, start, end):
    i=txt.index(start)+len(start); j=txt.index(end, i)
    return txt[i:j]

# --- baslik ---
title_blk=grab(ms,"## Baslik","## Ozet")
title=re.search(r"\*\*(.+?)\*\*", title_blk, re.S).group(1)
title=" ".join(title.split())

abstract=paras(grab(ms,"## Ozet","## Giris"))
intro   =paras(grab(ms,"## Giris","## Sonuclar (Results)"))

# --- results: alt bolumler ---
res_raw=grab(ms,"## Sonuclar (Results)","## Tartisma (Discussion)")
res=[]
for m in re.finditer(r"### (.+?)\n(.*?)(?=\n### |\Z)", res_raw, re.S):
    head=" ".join(m.group(1).split()); body=paras(m.group(2))
    res.append((inline(head), body))

disc=paras(grab(ms,"## Tartisma (Discussion)","---\n\n## TABLO PLANI"))

# --- methods: JPR sirasi (agirliklama one) ---
me_body=grab(me,"# TASLAK METIN","---\n\n## DEGISMESI GEREKEN YERLER")
msec={}
for m in re.finditer(r"## (2\.\d+) (.+?)\n(.*?)(?=\n## 2\.|\Z)", me_body, re.S):
    msec[m.group(1)]=(" ".join(m.group(2).split()), m.group(3))
ORDER=["2.1","2.2","2.3","2.4","2.5","2.6","2.7","2.8","2.9","2.10"]
JPR_ORDER=["2.1","2.2","2.3","2.4","2.5","2.6","2.7","2.8","2.9","2.10"]
meth=[]
for k in JPR_ORDER:
    if k not in msec: continue
    h,b=msec[k]
    b=b.replace("(two,\none, five and two)","(two, one, five and three)")
    b=b.replace("two,\none, five and two","two, one, five and three")
    meth.append((inline(h), paras(b)))

FIGS=[("figure1_reuse_workflow","Reuse workflow. Site-level quantification deposited for four substrate loads is used as deposited; the condition reported by the published figure is identified, labels are constructed and the weighting scheme evaluated, the structure--response relationship is tested separately in each condition, and site composition is examined with a matched core and random-subset resampling."),
      ("figure2_label_weighting","Label weighting governs whether the relationship is detected. (A) Inverse-variance weights, median-normalised, for each substrate load; the dotted line marks unit weight. (B) Kish effective sample size as a fraction of the nominal site count, with the raw values above each bar. (C) Spearman rank correlation under each weighting scheme. (D) Empirical permutation $P$ values on a logarithmic scale; the dashed line marks $\\alpha = 0.05$."),
      ("figure3_replication_composition","Condition replication and site composition. (A) Permutation null distributions of mean absolute error for each load (grey) with the observed value marked. (B) Out-of-fold predictions against observed response, by load; the dotted line is the identity. (C) Spearman correlations of 200 random twenty-site subsets of the 500 ng condition, with the median and the matched twenty-site core marked. (D) Spearman correlation on the matched core in each condition, with empirical $P$ values.")]

L=[]
A=L.append
A(r"\documentclass[journal=jprobs,manuscript=article,layout=traditional]{achemso}")
A(r"\usepackage[T1]{fontenc}")
A(r"\usepackage[utf8]{inputenc}")
A(r"\usepackage{graphicx}")
A(r"\usepackage{amsmath}")
A(r"\usepackage{booktabs}")
A(r"\usepackage{microtype}")
A(r"\SectionNumbersOn")
A(r"\author{Murat G\"ok}")
A(r"\affiliation{Department of Computer Engineering, Faculty of Engineering, Yalova University, 77200 Yalova, T\"urkiye}")
A(r"\email{murat.gok@yalova.edu.tr}")
A(r"\title{%s}" % inline(title))
A(r"\begin{document}")
A(r"\begin{abstract}")
A(abstract)
A(r"\end{abstract}")
A(r"\section{Introduction}")
A(intro)
A(r"\section{Experimental Section}")
for h,b in meth:
    A(r"\subsection{%s}" % h); A(b)
A(r"\section{Results and Discussion}")
for i,(h,b) in enumerate(res):
    A(r"\subsection{%s}" % h); A(b)
    if i==0:
        A(r"\begin{figure}[htbp]\centering\includegraphics[width=\textwidth]{figures/%s.pdf}" % FIGS[0][0])
        A(r"\caption{%s}\label{fig:workflow}\end{figure}" % FIGS[0][1])
    if h.lower().startswith("label weighting"):
        A(r"\begin{figure}[htbp]\centering\includegraphics[width=\textwidth]{figures/%s.pdf}" % FIGS[1][0])
        A(r"\caption{%s}\label{fig:weighting}\end{figure}" % FIGS[1][1])
    if h.lower().startswith("site composition"):
        A(r"\begin{figure}[htbp]\centering\includegraphics[width=\textwidth]{figures/%s.pdf}" % FIGS[2][0])
        A(r"\caption{%s}\label{fig:composition}\end{figure}" % FIGS[2][1])
A(r"\subsection{Discussion}")
A(disc)
A(r"\begin{acknowledgement}")
A(r"The author thanks the investigators who generated and openly shared the carboxymethyllysine proteomics datasets reused here.")
A(r"\end{acknowledgement}")
A(r"\end{document}")
out="\n\n".join(L)+"\n"
pathlib.Path("jpr_tex/manuscript_jpr.tex").write_text(out,encoding="utf-8")
print("tex yazildi:",len(out.splitlines()),"satir")
print("baslik:",title[:90])
print("results alt bolum:",len(res),"| methods alt bolum:",len(meth))
