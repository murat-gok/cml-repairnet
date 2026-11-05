import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
CM=1/2.54; W2=17.8*CM
plt.rcParams.update({"font.family":"Liberation Sans","font.size":7,"axes.labelsize":7,
 "axes.titlesize":7.5,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"legend.fontsize":6.5,
 "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "xtick.major.size":2.5,"ytick.major.size":2.5,"pdf.fonttype":42,
 "axes.spines.top":False,"axes.spines.right":False,"savefig.bbox":"tight"})
C={"casein":"#0072B2","hemoglobin":"#D55E00","collagen":"#009E73"}
fig,ax=plt.subplots(1,3,figsize=(W2,W2*0.30)); fig.subplots_adjust(wspace=0.34)
def P(a,l): a.text(-0.20,1.08,l,transform=a.transAxes,fontsize=8.5,fontweight="bold",va="top")

# A gozlenen bulk yanit, seviye bazli
obs={"casein":[1.34,0.97,0.92],"hemoglobin":[0.89,0.87,0.82],"collagen":[0.80,0.65,0.60]}
sd ={"casein":[0.61,0.06,0.01],"hemoglobin":[0.07,0.06,0.02],"collagen":[0.06,0.05,0.05]}
x=np.arange(3); w=0.26; P(ax[0],"A")
for i,(k,v) in enumerate(obs.items()):
    ax[0].bar(x+(i-1)*w,v,w,yerr=sd[k],color=C[k],label=k,edgecolor="none",
              error_kw=dict(lw=0.6,capsize=1.6,capthick=0.6))
ax[0].axhline(1.0,color="#666666",linestyle=":",linewidth=0.7)
ax[0].set_xticks(x); ax[0].set_xticklabels(["low","mid","high"])
ax[0].set_xlabel("modification level"); ax[0].set_ylabel("observed bulk repair fraction")
ax[0].set_ylim(0,2.35)
ax[0].legend(frameon=False,loc="upper center",ncol=3,handlelength=0.9,
              columnspacing=0.8,handletextpad=0.4,fontsize=5.8,bbox_to_anchor=(0.5,1.02))
ax[0].set_title("Ordering holds at every level",fontsize=7)

# B ongorulen siralama dagilimi
orders=["Hb > Cas > Col","Col > Cas > Hb","Cas > Col > Hb",
        "Cas > Hb > Col","Col > Hb > Cas"]
cnt=[39,35,12,7,3]; cols=["#CCCCCC"]*5; cols[3]="#D55E00"
P(ax[1],"B")
ax[1].barh(range(5)[::-1],cnt,color=cols,edgecolor="none",height=0.62)
ax[1].set_yticks(range(5)[::-1]); ax[1].set_yticklabels(orders,fontsize=6.0)
ax[1].axvline(96/6,color="#333333",linestyle="--",linewidth=0.7)
ax[1].text(96/6+1.0,0.15,"chance\n(1 in 6)",fontsize=5.8,color="#333333",linespacing=1.2)
ax[1].set_xlabel("pipeline configurations (of 96 complete)")
ax[1].set_title("Predicted ordering rarely matches",fontsize=7)

# C lizin uygunlugu
P(ax[2],"C")
sc=["casein","hemoglobin","collagen"]; tot=[58,24,69]; strict=[8,24,0]; dom=[18,12,1]
xx=np.arange(3); w2=0.26
ax[2].bar(xx-w2,tot,w2,color="#CCCCCC",label="lysines",edgecolor="none")
ax[2].bar(xx,strict,w2,color="#888888",label="strict quality",edgecolor="none")
ax[2].bar(xx+w2,dom,w2,color="#0072B2",label="in applicability domain",edgecolor="none")
ax[2].set_xticks(xx); ax[2].set_xticklabels(sc,fontsize=6.2)
ax[2].set_ylabel("lysines"); ax[2].set_ylim(0,88)
ax[2].legend(frameon=False,loc="upper center",ncol=1,handlelength=0.9,fontsize=5.8)
ax[2].set_title("Collagen has no eligible site",fontsize=7)
fig.savefig("figureS1_cross_scaffold.pdf")
fig.savefig("figureS1_cross_scaffold.tif",dpi=600,pil_kwargs={"compression":"tiff_lzw"})
print("Figure S1 uretildi")
