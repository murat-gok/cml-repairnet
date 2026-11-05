import json, numpy as np, pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr

SP="/tmp/claude-0/-home-claude/2bc12b7f-2a5c-59c2-826e-1fb99cbd7646/scratchpad/"
UP="/mnt/user-data/uploads/RepairNet/"
FEATS=["residue_sasa","relative_sasa_lysine","distance_to_protein_centroid",
 "neighbor_hydrophobic_fraction_8a","neighbor_cysteine_count_8a","neighbor_aromatic_count_8a",
 "ca_contact_count_8a","ca_contact_count_10a","ca_contact_count_12a"]
CONDS=["0p25","0p5","1p0","2p0"]; RATIO={"0p25":0.25,"0p5":0.5,"1p0":1.0,"2p0":2.0}
VFLOOR=1e-4

raw=json.load(open(SP+"pxd070124_parsed.json"))
struct=pd.read_csv(SP+"structure_features_union44.csv")

# ---- build long label table (same formula as prepare_labels.py) ----
rows=[]
for c in CONDS:
    for r in raw[c]:
        res,conf=int(r[0]),float(r[1]); a=np.array(r[2:5],float); k=np.array(r[5:8],float)
        am,km=a.mean(),k.mean(); asd,ksd=a.std(ddof=1),k.std(ddof=1)
        if km==0: continue
        raw_resp=1.0-am/km
        var=(asd**2/3)/(km**2)+(am**2)*(ksd**2/3)/(km**4)
        rows.append(dict(cond=c,ratio=RATIO[c],site_uniprot=res,conf=conf,
            active_mean=am,control_mean=km,resp_raw=raw_resp,
            resp=float(np.clip(raw_resp,0,1)),var=var))
lab=pd.DataFrame(rows)
lab["var"]=lab["var"].fillna(1.0)
lab["w_raw"]=1.0/(lab["var"]+VFLOOR)
lab["weight"]=lab.groupby("cond")["w_raw"].transform(lambda s:s/s.median())
hi=lab[lab.conf>90].copy()

def blocked(df,seed=0,nperm=200,ngroup=4,wmode='none'):
    X=df[FEATS].to_numpy(float); y=df["resp"].to_numpy(float); w=df["weight"].to_numpy(float)
    if wmode=="none": w=np.ones(len(w))
    XYZ=df[["ca_x","ca_y","ca_z"]].to_numpy(float)
    g=AgglomerativeClustering(n_clusters=ngroup,linkage="ward").fit_predict(XYZ)
    def run(Xm):
        oof=np.full(len(y),np.nan)
        for gg in np.unique(g):
            te=g==gg; tr=~te
            if tr.sum()<4 or te.sum()<1: continue
            sc=StandardScaler().fit(Xm[tr])
            m=ExtraTreesRegressor(n_estimators=500,max_depth=4,min_samples_leaf=2,
                                  random_state=seed,n_jobs=1)
            m.fit(sc.transform(Xm[tr]),y[tr],sample_weight=w[tr])
            oof[te]=m.predict(sc.transform(Xm[te]))
        ok=~np.isnan(oof)
        return mean_absolute_error(y[ok],oof[ok]), spearmanr(y[ok],oof[ok]).statistic, oof,ok
    mae,rho,oof,ok=run(X)
    # null: mean of training folds
    nb=[]
    for gg in np.unique(g):
        te=g==gg; tr=~te
        if tr.sum()<4 or te.sum()<1: continue
        nb.append(np.abs(y[te]-np.average(y[tr],weights=w[tr])))
    base=float(np.concatenate(nb).mean())
    rng=np.random.default_rng(seed)
    pm=[]
    for _ in range(nperm):
        idx=rng.permutation(len(y))
        pm.append(run(X[idx])[0])
    pm=np.array(pm)
    p=(np.sum(pm<=mae)+1)/(nperm+1)
    return dict(n=len(y),groups=int(len(np.unique(g))),mae=mae,spearman=float(rho),
                baseline=base,rel_impr=(base-mae)/base,perm_mae_mean=float(pm.mean()),
                dMAE=float(pm.mean()-mae),p_perm=float(p))


NG={"0p25":250,"0p5":500,"1p0":1000,"2p0":2000}
avail={}
for c in CONDS:
    h=hi[hi.cond==c]
    avail[c]=h.merge(struct[["site_uniprot"]+FEATS+["ca_x","ca_y","ca_z"]],on="site_uniprot")

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr

avail={}
for c in CONDS:
    h=hi[hi.cond==c]
    avail[c]=h.merge(struct[["site_uniprot"]+FEATS+["ca_x","ca_y","ca_z"]],on="site_uniprot")
NG={"0p25":250,"0p5":500,"1p0":1000,"2p0":2000}

def core_run(df,seed=0,ngroup=4,unw=True,nperm=0,ret_null=False):
    X=df[FEATS].to_numpy(float); y=df["resp"].to_numpy(float)
    w=np.ones(len(y)) if unw else df["weight"].to_numpy(float)
    XYZ=df[["ca_x","ca_y","ca_z"]].to_numpy(float)
    g=AgglomerativeClustering(n_clusters=ngroup,linkage="ward").fit_predict(XYZ)
    def run(Xm):
        oof=np.full(len(y),np.nan)
        for gg in np.unique(g):
            te=g==gg; tr=~te
            if tr.sum()<4 or te.sum()<1: continue
            sc=StandardScaler().fit(Xm[tr])
            m=ExtraTreesRegressor(n_estimators=500,max_depth=4,min_samples_leaf=2,random_state=seed,n_jobs=1)
            m.fit(sc.transform(Xm[tr]),y[tr],sample_weight=w[tr])
            oof[te]=m.predict(sc.transform(Xm[te]))
        ok=~np.isnan(oof)
        return mean_absolute_error(y[ok],oof[ok]), spearmanr(y[ok],oof[ok]).statistic, oof, ok
    mae,rho,oof,ok=run(X)
    out=dict(mae=mae,rho=float(rho),oof=oof,y=y)
    if nperm:
        rng=np.random.default_rng(seed); pm=[]
        for _ in range(nperm):
            pm.append(run(X[rng.permutation(len(y))])[0])
        out["null"]=np.array(pm)
    return out

save={}
# A) permutasyon bos dagilimlari (agirliksiz, 4 kosul)
for c in CONDS:
    r=core_run(avail[c],nperm=200)
    save[f"null_{c}"]=r["null"]; save[f"mae_{c}"]=r["mae"]
    save[f"oof_{c}"]=r["oof"]; save[f"y_{c}"]=r["y"]
    print(f"null {NG[c]}: real={r['mae']:.3f} nullmean={r['null'].mean():.3f}",flush=True)
# B) 200 rastgele 20-altkume Spearman (agirliksiz, 0p5)
d5=avail["0p5"].reset_index(drop=True)
core=sorted(set.intersection(*[set(avail[c].site_uniprot) for c in CONDS]))
rng=np.random.default_rng(7); sp=[]
for i in range(200):
    idx=rng.choice(len(d5),20,replace=False)
    sp.append(core_run(d5.iloc[idx].reset_index(drop=True),ngroup=3)["rho"])
save["subset_rho"]=np.array(sp)
save["core_rho"]=core_run(d5[d5.site_uniprot.isin(core)].reset_index(drop=True),ngroup=3)["rho"]
print(f"subset: mean={np.mean(sp):.3f} core={save['core_rho']:.3f}",flush=True)
# C) agirliklar, guven-yanit, cekirdek/disi
for c in CONDS:
    save[f"w_{c}"]=avail[c]["weight"].to_numpy(float)
    save[f"conf_{c}"]=avail[c]["conf"].to_numpy(float)
    save[f"resp_{c}"]=avail[c]["resp"].to_numpy(float)
save["core_sites"]=np.array(core)
save["d5_sites"]=d5.site_uniprot.to_numpy()
save["d5_resp"]=d5["resp"].to_numpy(float)
np.savez("figdata.npz",**save)
print("figdata.npz yazildi")
