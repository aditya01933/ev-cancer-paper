"""
cross_species_te_density.py  v2
Use NCBI FTP _rm.out.gz (RepeatMasker output) per species.
Parse TE intervals on Chr1, compute fraction of 5kb windows >50% TE covered.
Correlate with Zone3% from HANDOFF.
"""
import os, gzip, json, subprocess, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats

WORK = os.path.expanduser("~/ai-project/AncientKeyGen1/imp-research")
WIN, TE_THRESH = 5000, 0.50

# NCBI FTP rm.out.gz URLs + local filename + Chr1 RefSeq accession + chr1 size
SPECIES = {
    "Rice":         ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rm.out.gz",
                     "rice_rm.out.gz",    "NC_029256.1", 43270923),
    "Maize":        ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rm.out.gz",
                     "maize_rm.out.gz",   "NC_050096.1", 308452471),
    "Sorghum":      ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_rm.out.gz",
                     "sorghum_rm.out.gz", "NC_012870.2", 80884392),
    "Brachypodium": ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/505/GCF_000005505.3_Brachypodium_distachyon_v3.0/GCF_000005505.3_Brachypodium_distachyon_v3.0_rm.out.gz",
                     "brachy_rm.out.gz",  "NC_016131.3", 75070568),
    "Arabidopsis":  ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_rm.out.gz",
                     "ath_rm.out.gz",     "NC_003070.9", 30427671),
    "Soybean":      ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v2.1/GCF_000004515.6_Glycine_max_v2.1_rm.out.gz",
                     "soybean_rm.out.gz", "NC_016088.3", 59478403),
}

ZONE3 = {"Arabidopsis":8.4,"Rice":25.2,"Maize":35.1,"Sorghum":32.4,
         "Brachypodium":37.7,"Soybean":12.4}

def download(url, dest):
    if os.path.exists(dest):
        print(f"  cached: {dest}"); return True
    print(f"  downloading {os.path.basename(dest)}...", flush=True)
    r = subprocess.run(["curl","-fsSL","--retry","3","-o",dest,url])
    ok = r.returncode == 0 and os.path.getsize(dest) > 1000
    if not ok:
        print(f"  FAILED"); os.remove(dest) if os.path.exists(dest) else None
    else:
        print(f"  ok ({os.path.getsize(dest)//1024}KB)")
    return ok

def parse_rm_out(path, chrom):
    """Return sorted (start,end) 0-based intervals on chrom from rm.out(.gz)."""
    ivs = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path,"rt",errors="replace") as fh:
        for i,line in enumerate(fh):
            if i < 3: continue
            f = line.split()
            if len(f) < 7: continue
            if f[4] != chrom: continue
            try: ivs.append((int(f[5])-1, int(f[6])))
            except ValueError: continue
    return sorted(ivs)

def dense_fraction(ivs, chr_size):
    n_total = chr_size // WIN
    if not ivs or n_total == 0: return 0.0
    n_dense, iv_idx = 0, 0
    for w in range(n_total):
        ws, we = w*WIN, (w+1)*WIN
        cov = 0
        for s,e in ivs[iv_idx:]:
            if s >= we: break
            cov += min(e,we) - max(s,ws)
        if cov >= TE_THRESH * WIN: n_dense += 1
        while iv_idx < len(ivs) and ivs[iv_idx][1] <= we: iv_idx += 1
    return n_dense / n_total

def main():
    os.chdir(WORK)
    results = []
    for sp,(url,fname,chrom,sz) in SPECIES.items():
        path = os.path.join(WORK,fname)
        if not download(url,path): continue
        print(f"  [{sp}] parsing {chrom}...", flush=True)
        ivs = parse_rm_out(path, chrom)
        print(f"    intervals={len(ivs)}", flush=True)
        if not ivs:
            print(f"    WARNING: 0 intervals — check chrom name"); continue
        frac = dense_fraction(ivs, sz)
        z3 = ZONE3[sp]
        print(f"    dense_frac={frac*100:.1f}%  Zone3={z3}%")
        results.append({"species":sp,"te_dense_frac":frac,"zone3_pct":z3})

    if len(results) < 4:
        print("Too few species — check downloads/chrom names"); return

    xs = np.array([r["te_dense_frac"]*100 for r in results])
    ys = np.array([r["zone3_pct"]          for r in results])
    labels = [r["species"] for r in results]

    r,p = stats.pearsonr(xs,ys)
    rho,p_rho = stats.spearmanr(xs,ys)
    print(f"\nPearson r={r:.3f} p={p:.4f}  Spearman rho={rho:.3f} p={p_rho:.4f}")

    fig,ax = plt.subplots(figsize=(7,5.5))
    sl,ic,*_ = stats.linregress(xs,ys)
    xr = np.linspace(xs.min()-1,xs.max()+1,200)
    ax.plot(xr,sl*xr+ic,color="#555",lw=1.8,alpha=0.7,label=f"r={r:.3f}  p={p:.3f}")
    for x,y,sp in zip(xs,ys,labels):
        ax.scatter(x,y,s=90,color="#2e8b57",edgecolors="#2e8b57",zorder=3)
        ax.annotate(sp,(x,y),xytext=(x+0.3,y+0.4),fontsize=9,color="#2e8b57",
                    path_effects=[pe.withStroke(linewidth=2,foreground="white")])
    ax.set_xlabel("Windows with >50% TE coverage (%)",fontsize=10)
    ax.set_ylabel("Zone3 % (Ev scan)",fontsize=10)
    ax.set_title("Cross-species: direct TE density vs Zone3%\n(NCBI RepeatMasker, Chr1)",fontsize=11)
    ax.legend(fontsize=9); ax.grid(color="white"); ax.set_facecolor("#f8f9fa")
    plt.tight_layout()
    for ext in ("png","pdf"):
        plt.savefig(f"{WORK}/cross_species_te_density_plot.{ext}",
                    dpi=180 if ext=="png" else 150,bbox_inches="tight")
    plt.close()
    print("Saved plots.")
    json.dump({"results":results,"pearson_r":round(r,4),"pearson_p":round(p,6),
               "spearman_rho":round(rho,4),"spearman_p":round(p_rho,6)},
              open(f"{WORK}/te_density_results.json","w"),indent=2)

if __name__=="__main__":
    main()
