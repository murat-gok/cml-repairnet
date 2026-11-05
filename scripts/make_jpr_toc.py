import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
CM=1/2.54
W,H=3.25,1.75          # ACS TOC azami olcusu (inc)
plt.rcParams.update({"font.family":"Liberation Sans","font.size":6.2,
 "axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
 "xtick.major.size":1.8,"ytick.major.size":1.8,"pdf.fonttype":42,
 "axes.spines.top":False,"axes.spines.right":False,
 "xtick.labelsize":5.6,"ytick.labelsize":5.6})
BLUE="#0072B2"; GREY="#999999"; OR="#D55E00"

fig=plt.figure(figsize=(W,H))
# --- sol: arsivdeki dort kosul ---
axl=fig.add_axes([0.015,0.10,0.40,0.78]); axl.axis("off")
axl.set_xlim(0,1); axl.set_ylim(0,1)
axl.text(0.5,0.97,"one archive, four conditions",ha="center",va="top",
         fontsize=6.4,fontweight="bold")
loads=[("2000 ng",0.66),("1000 ng",0.48),("500 ng",0.30),("250 ng",0.12)]
for lab,y in loads:
    pub = lab.startswith("500")
    axl.add_patch(FancyBboxPatch((0.10,y),0.62,0.135,boxstyle="round,pad=0.012",
        linewidth=0.9 if pub else 0.5,
        edgecolor=OR if pub else "#555555",
        facecolor="#FDE7D9" if pub else "#F2F2F2"))
    axl.text(0.41,y+0.068,lab,ha="center",va="center",fontsize=6.0,
             color=OR if pub else "#222222", fontweight="bold" if pub else "normal")
axl.annotate("published\nfigure",xy=(0.74,0.368),xytext=(0.80,0.18),
    fontsize=5.5,color=OR,ha="left",va="center",linespacing=1.2,
    arrowprops=dict(arrowstyle="-|>",color=OR,linewidth=0.6,mutation_scale=5))
# ok
axm=fig.add_axes([0.405,0.10,0.09,0.78]); axm.axis("off")
axm.set_xlim(0,1); axm.set_ylim(0,1)
axm.add_patch(FancyArrowPatch((0.15,0.46),(0.88,0.46),arrowstyle="-|>",
    mutation_scale=7,linewidth=0.8,color="#333333"))
# --- sag: agirliklama sonucu ---
axr=fig.add_axes([0.575,0.235,0.405,0.56])
x=np.arange(4); wd=0.38
unw=[0.010,0.005,0.005,0.005]; wtd=[0.109,0.005,0.010,0.005]
axr.bar(x-wd/2,unw,wd,color=BLUE,edgecolor="none",label="unweighted")
axr.bar(x+wd/2,wtd,wd,color=GREY,edgecolor="none",label="inv-variance")
axr.axhline(0.05,color=OR,linestyle="--",linewidth=0.7)
axr.text(3.52,0.057,"0.05",fontsize=5.0,color=OR,ha="right",va="bottom")
axr.set_yscale("log"); axr.set_xticks(x); axr.set_xticklabels(["250","500","1k","2k"])
axr.set_ylabel("permutation $P$",fontsize=5.8,labelpad=1.5)
axr.set_xlabel("substrate load (ng)",fontsize=5.8,labelpad=1.5)
axr.tick_params(pad=1.2)
axr.legend(frameon=False,fontsize=5.2,loc="upper center",ncol=2,
           bbox_to_anchor=(0.5,1.30),handlelength=0.9,columnspacing=0.7,handletextpad=0.4)
fig.text(0.78,0.055,"label weighting decides the answer",ha="center",
         fontsize=6.0,fontweight="bold")
fig.savefig("toc_graphic.pdf")
fig.savefig("toc_graphic.tif",dpi=600,pil_kwargs={"compression":"tiff_lzw"})
print("TOC uretildi: %.2f x %.2f inc"%(W,H))
