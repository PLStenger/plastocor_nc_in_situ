#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob
BASE=/nvme/bio/data_fungi/plastocor_nc_in_situ
OLD="$BASE/00_NEW_RES"
RAW="$BASE/01_raw_data/Original_reads"
OUT="${OUT:-$OLD/11_all_samples_identity_comparison}"
source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir -p "$OUT"/{metadata,trimmed,qiime2,logs,tmp}
export TMPDIR="$OUT/tmp"
exec 9>"$OUT/.lock"
flock -n 9 || { echo 'Execution deja en cours' >&2; exit 1; }
exec > >(tee -a "$OUT/logs/run_$(date +%Y%m%d_%H%M%S).log") 2>&1
conda activate qiime2-amplicon-2025.7
# Reprendre la correspondance explicitement confirmee lors du premier pipeline.
CONTROL_STEM="$(python - "$OLD/00_metadata/identity_lock.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['control'])
PY
)"
WRONG_STEM=FILTRE-STERILE_S60
[[ -f "$RAW/${CONTROL_STEM}_R1_001.fastq.gz" ]]
[[ -f "$RAW/${WRONG_STEM}_R1_001.fastq.gz" ]]
# Tous les fichiers originaux, SED et Temoin-pool compris.
python - "$RAW" "$OUT" "$CONTROL_STEM" "$WRONG_STEM" <<'PY'
from pathlib import Path
import sys,csv,json
raw,out,control,wrong=sys.argv[1:]; raw=Path(raw); out=Path(out)
rows=[]
control_label=control.rsplit('_S',1)[0]
for r1 in sorted(raw.glob('*_R1_001.fastq.gz')):
 sid=r1.name.removesuffix('_R1_001.fastq.gz'); r2=raw/f'{sid}_R2_001.fastq.gz'
 if not r2.exists(): raise SystemExit(f'R2 absent: {r2}')
 original=sid.rsplit('_S',1)[0]
 corrected='FILTRE-STERILE' if sid==control else control_label if sid==wrong else original
 def info(label):
  if label=='FILTRE-STERILE': return 'sterile_control','not_applicable','filter'
  tokens=label.split('-'); site='Nenouire' if 'NC' in tokens else 'Tontouta' if 'C' in tokens else 'not_specified'
  typ='pooled_control' if 'Temoin-pool' in label else 'sample'
  return typ,site,tokens[0]
 rows.append([sid,original,corrected,*info(original),*info(corrected),str(r1),str(r2)])
if not rows: raise SystemExit('Aucun FASTQ')
config={'samples':rows,'control':control,'wrong':wrong}
lock=out/'metadata/input_lock.json'
if lock.exists() and json.loads(lock.read_text())!=config:
 raise SystemExit('Entrees modifiees: utiliser un autre OUT')
lock.write_text(json.dumps(config,indent=2))
header=['sample-id','original_name','corrected_name','original_type','original_site','original_material','corrected_type','corrected_site','corrected_material','r1','r2']
with (out/'metadata/all_samples.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(header); w.writerows(rows)
print('Nombre total de paires de fichiers:',len(rows))
print('Vrai controle selon correspondance existante:',control)
PY
# Reprendre exactement les adaptateurs du premier pipeline.
ADAPTERS="$OLD/00_metadata/library-adapters.fasta"
[[ -s "$ADAPTERS" ]] || { echo "Adaptateurs absents: $ADAPTERS" >&2; exit 1; }
if [[ ! -f "$OUT/trimmed.done" ]]; then
 conda activate trimmomatic
 while IFS=$'\t' read -r sid original corrected ot os om ct cs cm r1 r2; do
  [[ "$sid" == sample-id ]] && continue
  trimmomatic -Xmx4G PE -threads 1 -phred33 "$r1" "$r2" \
   "$OUT/trimmed/${sid}_R1_P.fastq.gz" "$OUT/trimmed/${sid}_R1_U.fastq.gz" \
   "$OUT/trimmed/${sid}_R2_P.fastq.gz" "$OUT/trimmed/${sid}_R2_U.fastq.gz" \
   "ILLUMINACLIP:$ADAPTERS:2:30:10" MINLEN:100 \
   2> "$OUT/logs/trimmomatic_${sid}.log"
 done < "$OUT/metadata/all_samples.tsv"
 touch "$OUT/trimmed.done"
fi
conda activate qiime2-amplicon-2025.7
Q="$OUT/qiime2"
if [[ ! -f "$OUT/dada2.done" ]]; then
 # Apres interruption, archiver/supprimer les sorties QIIME partielles avant reprise.
 python - "$OUT" <<'PY'
from pathlib import Path
import csv,sys
p=Path(sys.argv[1])
with (p/'metadata/manifest.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 w.writerow(['sample-id','forward-absolute-filepath','reverse-absolute-filepath'])
 for r in csv.DictReader((p/'metadata/all_samples.tsv').open(),delimiter='\t'):
  sid=r['sample-id']; w.writerow([sid,*[str((p/'trimmed'/f'{sid}_R{m}_P.fastq.gz').resolve()) for m in (1,2)]])
PY
 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$OUT/metadata/manifest.tsv" --input-format PairedEndFastqManifestPhred33V2 --output-path "$Q/imported.qza"
 qiime cutadapt trim-paired --i-demultiplexed-sequences "$Q/imported.qza" \
  --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0.1 --p-overlap 15 --p-discard-untrimmed --p-cores 8 \
  --o-trimmed-sequences "$Q/no-primers.qza" --verbose 2>&1 | tee "$OUT/logs/cutadapt.log"
 qiime demux summarize --i-data "$Q/no-primers.qza" --o-visualization "$Q/no-primers.qzv"
 qiime dada2 denoise-paired --i-demultiplexed-seqs "$Q/no-primers.qza" \
  --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 \
  --p-max-ee-f 2 --p-max-ee-r 2 --p-n-threads 8 --p-chimera-method consensus \
  --o-table "$Q/table-all.qza" --o-representative-sequences "$Q/rep-seqs-all.qza" --o-denoising-stats "$Q/stats.qza"
 qiime metadata tabulate --m-input-file "$Q/stats.qza" --o-visualization "$Q/stats.qzv"
 qiime feature-table summarize --i-table "$Q/table-all.qza" \
  --m-sample-metadata-file "$OUT/metadata/all_samples.tsv" --o-visualization "$Q/table-all.qzv"
 touch "$OUT/dada2.done"
fi
python - "$OUT" <<'PY'
from pathlib import Path
import sys,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage,dendrogram,leaves_list
from qiime2 import Artifact
p=Path(sys.argv[1]); meta=pd.read_csv(p/'metadata/all_samples.tsv',sep='\t',index_col=0)
x=Artifact.load(str(p/'qiime2/table-all.qza')).view(pd.DataFrame)
x=x.reindex(meta.index,fill_value=0)
x.to_csv(p/'ASV_counts_original_file_ids.tsv',sep='\t',index_label='sample-id')
depths=x.sum(axis=1)
meta.assign(reads_after_DADA2=depths).to_csv(p/'sample_inventory.tsv',sep='\t')
valid=depths.gt(0)
meta.loc[~valid].to_csv(p/'samples_without_ASVs.tsv',sep='\t')
x=x.loc[valid]; x=x.loc[:,x.sum(axis=0).gt(0)]
if len(x)<2: raise SystemExit('Moins de deux echantillons non vides')
rel=x.div(x.sum(axis=1),axis=0)
d=pdist(rel.to_numpy(),metric='braycurtis'); dm=squareform(d)
z=linkage(d,method='average'); order=leaves_list(z)
# Calcul unique, memes coordonnees et ordre dans les deux versions.
n=len(x); h=np.eye(n)-np.ones((n,n))/n
b=-0.5*h@(dm**2)@h
values,vectors=np.linalg.eigh(b); idx=np.argsort(values)[::-1]
values=values[idx]; vectors=vectors[:,idx]
positive=values>1e-12
coords=np.zeros((n,2)); pos=np.flatnonzero(positive)[:2]
coords[:,:len(pos)]=vectors[:,pos]*np.sqrt(values[pos])
variance=np.zeros(2)
if positive.any(): variance[:len(pos)]=100*values[pos]/values[positive].sum()
pd.DataFrame({'eigenvalue':values}).to_csv(p/'PCoA_eigenvalues.tsv',sep='\t',index=False)
for version,col in [('01_original_names','original_name'),('02_corrected_names','corrected_name')]:
 out=p/version; out.mkdir(exist_ok=True)
 labels=[f'{meta.loc[s,col]} [{s}]' for s in x.index]
 pd.DataFrame(dm,index=labels,columns=labels).to_csv(out/'braycurtis_distances.tsv',sep='\t')
 counts=x.copy(); counts.index=labels; counts.to_csv(out/'ASV_counts.tsv',sep='\t',index_label='sample')
 fig,ax=plt.subplots(figsize=(15,max(9,len(x)*0.4)))
 dendrogram(z,labels=labels,orientation='right',ax=ax)
 ax.set_xlabel('Bray-Curtis distance (relative abundances)'); fig.tight_layout()
 fig.savefig(out/'dendrogram.png',dpi=200); plt.close(fig)
 fig,ax=plt.subplots(figsize=(15,13))
 im=ax.imshow(dm[np.ix_(order,order)],vmin=0,vmax=1,cmap='viridis')
 ax.set_xticks(range(n)); ax.set_yticks(range(n))
 ax.set_xticklabels([labels[i] for i in order],rotation=90,fontsize=7)
 ax.set_yticklabels([labels[i] for i in order],fontsize=7)
 fig.colorbar(im,ax=ax,label='Bray-Curtis distance'); fig.tight_layout()
 fig.savefig(out/'heatmap.png',dpi=200); plt.close(fig)
 fig,ax=plt.subplots(figsize=(14,11)); ax.scatter(coords[:,0],coords[:,1],s=35)
 for i,label in enumerate(labels): ax.annotate(label,coords[i],xytext=(4,3),textcoords='offset points',fontsize=7)
 ax.set_xlabel(f'PCoA 1 ({variance[0]:.1f}% of positive inertia)')
 ax.set_ylabel(f'PCoA 2 ({variance[1]:.1f}% of positive inertia)')
 fig.tight_layout(); fig.savefig(out/'PCoA.png',dpi=200); plt.close(fig)
 pd.DataFrame(coords,index=labels,columns=['PCoA1','PCoA2']).to_csv(out/'PCoA_coordinates.tsv',sep='\t')
 neighbors=[]
 for i,sid in enumerate(x.index):
  others=[j for j in np.argsort(dm[i]) if j!=i][:3]
  for rank,j in enumerate(others,1): neighbors.append([labels[i],rank,labels[j],dm[i,j]])
 pd.DataFrame(neighbors,columns=['sample','rank','neighbor','distance']).to_csv(out/'nearest_neighbors.tsv',sep='\t',index=False)
print('Echantillons enumeres:',len(meta),'avec ASV:',len(x))
print('Aucune suppression des ASV du controle. Aucune rarefaction.')
print('Valeurs propres PCoA negatives:',int((values < -1e-10).sum()))
print('Deux presentations des MEMES profils, distances et coordonnees:',p)
PY
