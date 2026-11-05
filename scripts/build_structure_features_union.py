import json, numpy as np, pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1
HYD=set("AILMFWVY"); ARO=set("FWY")
L=lambda n: protein_letters_3to1.get(n.upper(),"X")

PDB="/mnt/user-data/uploads/RepairNet/data/raw/4F5S.pdb"
OLD="/mnt/user-data/uploads/RepairNet/data/processed/bsa_structure_features.csv"
raw=json.load(open("pxd070124_parsed.json"))
CONDS=["0p25","0p5","1p0","2p0"]
union=sorted({int(r[0]) for c in CONDS for r in raw[c] if float(r[1])>90})
print(f"conf>90 birlesim: {len(union)} bolge")

st=PDBParser(QUIET=True).get_structure("p",PDB)
model=next(st.get_models())
chain=model["A"]
ShrakeRupley().compute(model, level="R")
residues=[r for r in chain if r.id[0]==" " and "CA" in r]
by={int(r.id[1]):r for r in residues}
allca=np.vstack([r["CA"].coord for r in residues])
centroid=allca.mean(axis=0)
print(f"chain A standart rezidu: {len(residues)} | model zincir sayisi: {len(list(model))}")

rows=[]
for su in union:
    sp=su-24
    res=by.get(sp)
    if res is None:
        rows.append(dict(site_uniprot=su,site_pdb=sp,structure_missing=1)); continue
    if L(res.resname)!="K":
        rows.append(dict(site_uniprot=su,site_pdb=sp,structure_missing=1,note=f"not K: {res.resname}")); continue
    ca=np.asarray(res["CA"].coord,float)
    d=np.linalg.norm(allca-ca,axis=1)
    nb=[o for o in residues if o is not res and np.linalg.norm(o["CA"].coord-ca)<=8.0]
    nl=[L(o.resname) for o in nb]
    sasa=float(getattr(res,"sasa",np.nan))
    r=dict(site_uniprot=su,site_pdb=sp,structure_missing=0,
        ca_x=ca[0],ca_y=ca[1],ca_z=ca[2],bfactor_ca=float(res["CA"].bfactor),
        residue_sasa=sasa, relative_sasa_lysine=min(sasa/205.0,1.5),
        distance_to_protein_centroid=float(np.linalg.norm(ca-centroid)),
        neighbor_hydrophobic_fraction_8a=float(np.mean([a in HYD for a in nl])) if nl else 0.0,
        neighbor_cysteine_count_8a=nl.count("C"),
        neighbor_aromatic_count_8a=sum(a in ARO for a in nl))
    for rad in (8.0,10.0,12.0):
        r[f"ca_contact_count_{rad:g}a"]=int(np.sum((d>0)&(d<=rad)))
    rows.append(r)
new=pd.DataFrame(rows)
print(f"cikarilan: {(new.structure_missing==0).sum()} / {len(new)}")
miss=new[new.structure_missing==1]
if len(miss): print("cikarilamayan:", miss.to_dict("records"))

# DOGRULAMA
old=pd.read_csv(OLD)
FE=["residue_sasa","relative_sasa_lysine","distance_to_protein_centroid",
    "neighbor_hydrophobic_fraction_8a","neighbor_cysteine_count_8a","neighbor_aromatic_count_8a",
    "ca_contact_count_8a","ca_contact_count_10a","ca_contact_count_12a","ca_x","ca_y","ca_z","bfactor_ca"]
m=old[["site_uniprot"]+FE].merge(new[["site_uniprot"]+FE],on="site_uniprot",suffixes=("_old","_new"))
print(f"\n=== DOGRULAMA: ortak {len(m)} bolge ===")
ok=True
for f in FE:
    a=m[f+"_old"].to_numpy(float); b=m[f+"_new"].to_numpy(float)
    mx=np.nanmax(np.abs(a-b)); same=mx<1e-6
    ok&=same
    print(f"  {f:>34} max|fark|={mx:.2e} {'OK' if same else '<<< FARKLI'}")
print(f"\nTUM OZNITELIKLER AYNI: {ok}")
new.to_csv("structure_features_union44.csv",index=False)
newonly=sorted(set(union)-set(old.site_uniprot))
print(f"\nYENI bolgeler ({len(newonly)}): {newonly}")
print(new[new.site_uniprot.isin(newonly)][["site_uniprot","site_pdb","structure_missing"]].to_string(index=False))
