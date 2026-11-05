import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---- ACS-uyumlu parametreler (tek yerden degistirilebilir) ----
CM=1/2.54
W1=8.3*CM          # tek sutun
W2=17.8*CM         # cift sutun
DPI_TIFF=600
FONT="Liberation Sans"   # Arial ile metrik uyumlu
plt.rcParams.update({
 "font.family":FONT, "font.size":7, "axes.labelsize":7, "axes.titlesize":7.5,
 "xtick.labelsize":6.5, "ytick.labelsize":6.5, "legend.fontsize":6.5,
 "axes.linewidth":0.6, "xtick.major.width":0.6, "ytick.major.width":0.6,
 "xtick.major.size":2.5, "ytick.major.size":2.5, "lines.linewidth":1.0,
 "pdf.fonttype":42, "ps.fonttype":42, "savefig.bbox":"tight",
 "savefig.pad_inches":0.02, "axes.spines.top":False, "axes.spines.right":False,
})
# renk korlugune dayanikli palet (Okabe-Ito)
C={"250":"#0072B2","500":"#D55E00","1000":"#009E73","2000":"#CC79A7",
   "unw":"#0072B2","wt":"#999999","core":"#D55E00"}

d=np.load("figdata.npz")
CONDS=["0p25","0p5","1p0","2p0"]; LAB={"0p25":"250","0p5":"500","1p0":"1000","2p0":"2000"}
# raporlanan sayilar
UNW={"250":(0.414,0.0100),"500":(0.616,0.0050),"1000":(0.399,0.0050),"2000":(0.581,0.0050)}
WTD={"250":(0.141,0.1095),"500":(0.681,0.0050),"1000":(0.298,0.0100),"2000":(0.574,0.0050)}
ESS={"250":(22,5.0),"500":(33,9.7),"1000":(37,11.8),"2000":(38,11.4)}
CORE={"250":(0.164,0.0796),"500":(-0.224,0.1940),"1000":(0.134,0.2637),"2000":(0.405,0.1294)}

def save(fig,name):
    fig.savefig(name+".pdf")
    try:
        fig.savefig(name+".tif",dpi=DPI_TIFF,pil_kwargs={"compression":"tiff_lzw"})
    except Exception as e:
        fig.savefig(name+".png",dpi=DPI_TIFF); print(f"  (tif yerine png: {e})")
    plt.close(fig); print("  yazildi:",name)

def panel(ax,letter):
    ax.text(-0.22,1.06,letter,transform=ax.transAxes,fontsize=8.5,fontweight="bold",va="top")

# ================= FIGURE 1: workflow =================
fig,ax=plt.subplots(figsize=(W2,W2*0.34)); ax.axis("off")
ax.set_xlim(0,10); ax.set_ylim(0,3)
boxes=[(0.15,"PXD070124\n4 substrate loads\n250/500/1000/2000 ng"),
       (2.15,"Condition identity\nreplicate match,\nsite count, >90% sites"),
       (4.15,"Label construction\nratio response,\nweighting evaluated"),
       (6.15,"Per-condition test\nblocked CV +\nsite-structure permutation"),
       (8.15,"Site composition\nmatched core,\nsubset resampling")]
for x,t in boxes:
    ax.add_patch(FancyBboxPatch((x,0.9),1.7,1.2,boxstyle="round,pad=0.06",
        linewidth=0.7,edgecolor="#333333",facecolor="#F2F2F2"))
    ax.text(x+0.85,1.5,t,ha="center",va="center",fontsize=6.3,linespacing=1.35)
for x in [1.9,3.9,5.9,7.9]:
    ax.add_patch(FancyArrowPatch((x,1.5),(x+0.22,1.5),arrowstyle="-|>",
        mutation_scale=7,linewidth=0.7,color="#333333"))
ax.text(5.0,0.45,"deposited processed quantification reused as deposited; no raw spectra reprocessed",
        ha="center",fontsize=6.0,style="italic",color="#555555")
save(fig,"figure1_reuse_workflow")

# ================= FIGURE 2: weighting =================
fig,axes=plt.subplots(2,2,figsize=(W2,W2*0.74)); fig.subplots_adjust(hspace=0.42,wspace=0.30)
# A weight distributions
ax=axes[0,0]; panel(ax,"A")
for i,c in enumerate(CONDS):
    w=d[f"w_{c}"]; k=LAB[c]
    y=np.full(len(w),i)+np.random.default_rng(i).uniform(-0.17,0.17,len(w))
    ax.scatter(w,y,s=3.5,alpha=0.75,color=C[k],linewidths=0)
ax.set_xscale("log"); ax.set_yticks(range(4)); ax.set_yticklabels([LAB[c] for c in CONDS])
ax.set_xlabel("inverse-variance weight (median-normalised)"); ax.set_ylabel("substrate load (ng)")
ax.axvline(1.0,color="#666666",linestyle=":",linewidth=0.7)
ax.set_title("Weights span up to five orders of magnitude",fontsize=7)
# B ESS/n
ax=axes[0,1]; panel(ax,"B")
ks=[LAB[c] for c in CONDS]; ratio=[ESS[k][1]/ESS[k][0] for k in ks]
ax.bar(range(4),ratio,color=[C[k] for k in ks],width=0.62,edgecolor="none")
ax.axhline(1.0,color="#333333",linewidth=0.7)
for i,k in enumerate(ks):
    ax.text(i,ratio[i]+0.03,f"{ESS[k][1]:.1f}/{ESS[k][0]}",ha="center",fontsize=6.0)
ax.set_xticks(range(4)); ax.set_xticklabels(ks); ax.set_ylim(0,1.12)
ax.set_xlabel("substrate load (ng)"); ax.set_ylabel("Kish ESS / nominal n")
ax.set_title("Effective sample size under weighting",fontsize=7)
# C Spearman unweighted vs weighted
ax=axes[1,0]; panel(ax,"C")
x=np.arange(4); wd=0.36
ax.bar(x-wd/2,[UNW[k][0] for k in ks],wd,label="unweighted",color=C["unw"],edgecolor="none")
ax.bar(x+wd/2,[WTD[k][0] for k in ks],wd,label="inverse-variance",color=C["wt"],edgecolor="none")
ax.axhline(0,color="#333333",linewidth=0.6)
ax.set_xticks(x); ax.set_xticklabels(ks); ax.set_xlabel("substrate load (ng)")
ax.set_ylabel("Spearman $\\rho$"); ax.set_ylim(0,0.84)
ax.legend(frameon=False,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.02),handlelength=1.2,columnspacing=1.0)
ax.set_title("Rank correlation by weighting scheme",fontsize=7)
# D permutation P
ax=axes[1,1]; panel(ax,"D")
ax.bar(x-wd/2,[UNW[k][1] for k in ks],wd,label="unweighted",color=C["unw"],edgecolor="none")
ax.bar(x+wd/2,[WTD[k][1] for k in ks],wd,label="inverse-variance",color=C["wt"],edgecolor="none")
ax.axhline(0.05,color="#D55E00",linestyle="--",linewidth=0.8)
ax.text(3.45,0.056,"$\\alpha$ = 0.05",fontsize=6.0,color="#D55E00",ha="right")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(ks)
ax.set_xlabel("substrate load (ng)"); ax.set_ylabel("permutation $P$")
ax.set_title("4/4 unweighted vs 3/4 weighted",fontsize=7)
ax.legend(frameon=False,loc="upper left")
save(fig,"figure2_label_weighting")

# ================= FIGURE 3: replication + composition =================
fig,axes=plt.subplots(2,2,figsize=(W2,W2*0.74)); fig.subplots_adjust(hspace=0.42,wspace=0.30)
# A permutation nulls
ax=axes[0,0]; panel(ax,"A")
for i,c in enumerate(CONDS):
    nl=d[f"null_{c}"]; k=LAB[c]
    parts=ax.violinplot([nl],positions=[i],widths=0.68,showextrema=False)
    for b in parts["bodies"]: b.set_facecolor("#CCCCCC"); b.set_alpha(0.85); b.set_edgecolor("none")
    ax.plot([i-0.3,i+0.3],[float(d[f"mae_{c}"])]*2,color=C[k],linewidth=1.6)
ax.set_xticks(range(4)); ax.set_xticklabels([LAB[c] for c in CONDS])
ax.set_xlabel("substrate load (ng)"); ax.set_ylabel("mean absolute error")
ax.set_title("Observed MAE vs permutation null",fontsize=7)
# B observed vs predicted
ax=axes[0,1]; panel(ax,"B")
for c in CONDS:
    k=LAB[c]; ax.scatter(d[f"y_{c}"],d[f"oof_{c}"],s=7,alpha=0.8,color=C[k],linewidths=0,label=k)
ax.plot([0,1],[0,1],color="#666666",linestyle=":",linewidth=0.7)
ax.set_xlim(-0.03,1.03); ax.set_ylim(-0.03,1.03)
ax.set_xlabel("observed response"); ax.set_ylabel("out-of-fold prediction")
ax.legend(frameon=False,title="load (ng)",title_fontsize=6.3,loc="upper left",ncol=2)
ax.set_title("Out-of-fold predictions by condition",fontsize=7)
# C subset histogram
ax=axes[1,0]; panel(ax,"C")
sp=d["subset_rho"]; cr=float(d["core_rho"])
ax.hist(sp,bins=22,color="#CCCCCC",edgecolor="white",linewidth=0.4)
ax.axvline(cr,color=C["core"],linewidth=1.4)
pct=(sp<cr).mean()*100
suf = "th" if 10<=int(round(pct))%100<=20 else {1:"st",2:"nd",3:"rd"}.get(int(round(pct))%10,"th")
ax.annotate(f"matched core\n$\\rho$ = {cr:.3f}\n{pct:.0f}{suf} percentile",
    xy=(cr,ax.get_ylim()[1]*0.55),xycoords="data",
    xytext=(0.04,0.93),textcoords="axes fraction",
    fontsize=6.2,color=C["core"],ha="left",va="top",linespacing=1.35,
    arrowprops=dict(arrowstyle="-|>",color=C["core"],linewidth=0.7,mutation_scale=6,
                    shrinkA=2,shrinkB=2))
ax.axvline(np.median(sp),color="#333333",linestyle="--",linewidth=0.7)
ax.text(np.median(sp)+0.02,ax.get_ylim()[1]*0.95,f"median {np.median(sp):.3f}",fontsize=6.0)
ax.set_xlabel("Spearman $\\rho$, random 20-site subsets (n = 200)")
ax.set_ylabel("count"); ax.set_title("Power is not the explanation",fontsize=7)
# D matched core per condition
ax=axes[1,1]; panel(ax,"D")
vals=[CORE[k][0] for k in ks]; ps=[CORE[k][1] for k in ks]
bars=ax.bar(range(4),vals,color=[C[k] for k in ks],width=0.62,edgecolor="none")
ax.axhline(0,color="#333333",linewidth=0.6)
for i,(v,pv) in enumerate(zip(vals,ps)):
    off=0.035 if v>=0 else -0.075
    ax.text(i,v+off,f"P = {pv:.2f}",ha="center",fontsize=6.0)
ax.set_xticks(range(4)); ax.set_xticklabels(ks); ax.set_ylim(-0.33,0.56)
ax.set_xlabel("substrate load (ng)"); ax.set_ylabel("Spearman $\\rho$, matched 20-site core")
ax.set_title("No condition retains the signal",fontsize=7)
save(fig,"figure3_replication_composition")
print("\nparametreler: tek sutun %.2f in | cift sutun %.2f in | TIFF %d dpi | font %s"%(W1,W2,DPI_TIFF,FONT))
