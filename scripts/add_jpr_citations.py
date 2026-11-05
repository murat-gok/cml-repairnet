import re, pathlib, sys
T=pathlib.Path("jpr_tex/manuscript_jpr.tex"); s=T.read_text(encoding="utf-8")
# (arama dizgisi, cite anahtarlari) — ILK gecise eklenir
M=[("ProteomeXchange PXD070124","trabosh2026a,trabosh2026b"),
   ("an enzyme that reverses","trabosh2026a"),
   ("Shrake--Rupley algorithm","shrake1973"),
   ("Shrake-Rupley algorithm","shrake1973"),
   ("maximum accessibility scale for lysine","tien2013"),
   ("extremely randomised trees ensemble","geurts2006"),
   ("PDB 4F5S","bujacz2012"),
   ("scikit-learn","pedregosa2011"),
   ("SciPy","virtanen2020"),
   ("random site-level cross-validation would place","roberts2017"),
   ("glyoxalase system and fructosamine-3-kinase","schalkwijk2020"),
   ("remains contested","richarme2015"),
   ("Kish effective sample size","kish1965"),
   ("recent transformer-based models","johansen2006,reddy2018,yu2019,liu2022,jia2024"),
   ("abundant in long-lived proteins","chaudhuri2018"),
   ("Site-resolved proteomics increasingly supplies","disanzo2021a"),
   ("reconciled against UniProt","uniprot2025"),
   ("UniProt numbering","uniprot2025"),
   ("Biopython","cock2009"),
   ("deposited in ProteomeXchange","perezriverol2025"),("files ","perezriverol2025")]
hit=[];miss=[]
for pat,keys in M:
    i=s.find(pat)
    if i<0: miss.append(pat); continue
    j=i+len(pat)
    s=s[:j]+("~\\cite{%s}"%keys)+s[j:]
    hit.append((pat,keys))
# bibliyografya
s=s.replace(r"\end{document}", "\\bibliography{references_jpr}\n\n\\end{document}")
T.write_text(s,encoding="utf-8")
print("EKLENEN (%d):"%len(hit))
for p,k in hit: print(f"   {k:<48} <- \"{p[:46]}\"")
print("\nTUTMAYAN (%d):"%len(miss))
for p in miss: print("   ",p)
