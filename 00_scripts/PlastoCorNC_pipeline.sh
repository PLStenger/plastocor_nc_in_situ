#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob
umask 0022
PROJECT_NAME=PlastoCorNC
RAW_ROOT_DIR=/nvme/bio/data_fungi/plastocor_nc_in_situ/01_raw_data/Original_reads
RESULTS_DIR="${RESULTS_DIR:-/nvme/bio/data_fungi/plastocor_nc_in_situ/00_NEW_RES}"
STAGE="${1:-qc}" # qc | dada2 | taxonomy | diversity | all
THREADS=1
QIIME_THREADS=8
TRIMMOMATIC_HEAP=4G
FASTQC_ENV=fastqc
MULTIQC_ENV=multiqc
TRIMMOMATIC_ENV=trimmomatic
QIIME2_ENV=qiime2-amplicon-2025.7
# Obligatoire : confirmer la correspondance entre numeros terrain et fichiers.
SWAP_CONFIRMED="${SWAP_CONFIRMED:-YES}"
ACTUAL_CONTROL_STEM="${ACTUAL_CONTROL_STEM:-SED-TANNE-C_S55}"
MISLABELLED_CONTROL_STEM=FILTRE-STERILE_S60

# Ex.: ACTUAL_CONTROL_STEM=SED-TANNE-C_S55 UNIQUEMENT apres verification.
# Adaptateurs de preparation de librairie, pas les amorces 16S.
ADAPTERS_FA="${ADAPTERS_FA:-/nvme/bio/data_fungi/plastocor_nc_in_situ/99_softwares/adapters_sequences.fasta}"
PRIMERS_PRESENT="${PRIMERS_PRESENT:-YES}" 
FWD=CCTACGGGNGGCWGCAG
REV=GACTACHVGGGTATCTAATCC

# Valeurs a fixer APRES inspection du QZV post-amorces.
TRUNC_F="${TRUNC_F:-}"
TRUNC_R="${TRUNC_R:-}"
DEPTH="${DEPTH:-}"
SILVA_VERSION=144
SILVA_URL="${SILVA_URL:-https://ftp.arb-silva.de/release_144/Exports/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz}"
die(){ echo "ERREUR: $*" >&2; exit 1; }
log(){ echo "[$(date -Is)] $*"; }
CONDA_BASE="$(conda info --base 2>/dev/null || true)"
[[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]] || die 'Conda introuvable'
source "$CONDA_BASE/etc/profile.d/conda.sh"
[[ -d "$RAW_ROOT_DIR" ]] || die "Raw root absent : $RAW_ROOT_DIR"
[[ "$SWAP_CONFIRMED" == YES && -n "$ACTUAL_CONTROL_STEM" ]] || die 'Confirmer inversion : SWAP_CONFIRMED=YES et ACTUAL_CONTROL_STEM=<prefixe exact sans _R1_001.fastq.gz>'
[[ "$ACTUAL_CONTROL_STEM" == SED* ]] || die 'Le vrai controle doit ici etre le fichier initialement etiquete SED concerne par inversion'
[[ "$PRIMERS_PRESENT" == YES || "$PRIMERS_PRESENT" == NO ]] || die 'PRIMERS_PRESENT doit etre YES ou NO'
case "$STAGE" in qc|dada2|taxonomy|diversity|all) ;; *) die 'Etape inconnue';; esac
mkdir -p "$RESULTS_DIR"/{00_metadata,01_raw_links,02_fastqc_raw,03_trimmomatic,04_fastqc_clean,05_qiime2,06_controls,07_taxonomy,08_diversity,09_exports,logs,tmp,state}
LOG_DIR="$RESULTS_DIR/logs"; TMPDIR_BASE="$RESULTS_DIR/tmp"; export TMPDIR="$TMPDIR_BASE"
exec > >(tee -a "$LOG_DIR/${STAGE}_$(date +%Y%m%d_%H%M%S).log") 2>&1
exec 9>"$RESULTS_DIR/.lock"
flock -n 9 || die 'Une autre execution utilise ce repertoire'
activate(){ conda activate "$1"; }
Q="$RESULTS_DIR/05_qiime2"; C="$RESULTS_DIR/06_controls"; T="$RESULTS_DIR/07_taxonomy"; D="$RESULTS_DIR/08_diversity"
M="$RESULTS_DIR/00_metadata"
# Aucune modification des fichiers bruts. Les exclusions suivent la correction.
activate "$QIIME2_ENV"
python - "$RAW_ROOT_DIR" "$RESULTS_DIR" "$ACTUAL_CONTROL_STEM" "$MISLABELLED_CONTROL_STEM" <<'PY'
import sys, pathlib, csv, json
raw,out,control,wrong=sys.argv[1:]; raw=pathlib.Path(raw); out=pathlib.Path(out)
config={'raw':str(raw),'control':control,'wrong_control':wrong}
lock=out/'00_metadata'/'identity_lock.json'
if lock.exists() and json.loads(lock.read_text())!=config:
    raise SystemExit('Correspondances modifiees : utiliser un autre RESULTS_DIR')
for stem in (control,wrong):
    for mate in (1,2):
        assert (raw/f'{stem}_R{mate}_001.fastq.gz').is_file(), (stem,mate)
lock.write_text(json.dumps(config,indent=2))
rows=[]; audit=[]; seen=set()
for r1 in sorted(raw.glob('*_R1_001.fastq.gz')):
    stem=r1.name.removesuffix('_R1_001.fastq.gz'); r2=raw/f'{stem}_R2_001.fastq.gz'
    assert r2.is_file(), f'Mate absent: {r2}'
    label=stem.rsplit('_S',1)[0]
    if stem==control:
        sid='CTRL-STERILE'; typ='negative_control'; site='not_applicable'; mat='filter'; status='not_applicable'; keep=True
        reason='Inversion corrigee: vrai filtre sterile'
    elif stem==wrong:
        sid='SED-TANNE-C-corrected'; keep=False; reason='Inversion corrigee: vrai sediment, exclu'
    elif label.startswith('SED') or 'Temoin-pool' in label:
        sid=label; keep=False; reason='SED ou Temoin-pool exclu'
    elif label.startswith(('PE-','PET-','EAU-')):
        sid=label; typ='sample'; mat=label.split('-')[0]
        tokens=label.split('-'); status='NC' if 'NC' in tokens else 'C' if 'C' in tokens else None
        assert status is not None, f'Statut inconnu: {stem}'
        site='Nenouire' if status=='NC' else 'Tontouta'; keep=True; reason='Echantillon biologique'
    else: raise SystemExit(f'Echantillon non reconnu: {stem}')
    audit.append([stem,sid,'yes' if keep else 'no',reason])
    if not keep: continue
    assert sid not in seen; seen.add(sid)
    paths=[]
    for mate,src in ((1,r1),(2,r2)):
        dest=out/'01_raw_links'/f'{sid}_R{mate}_001.fastq.gz'
        if dest.is_symlink(): assert dest.resolve()==src.resolve()
        elif dest.exists(): raise SystemExit(f'Fichier preexistant non lien: {dest}')
        else: dest.symlink_to(src.resolve())
        paths.append(str(dest))
    rows.append([sid,typ,site,mat,status,stem,*paths])
def write(name,header,data):
    with (out/'00_metadata'/name).open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(header); w.writerows(data)
write('sample_audit.tsv',['original_stem','corrected_id','included','reason'],audit)
write('samples.tsv',['sample-id','sample_type','site','material','contamination','original_stem'],[r[:6] for r in rows])
write('biological.tsv',['sample-id','sample_type','site','material','contamination','original_stem'],[r[:6] for r in rows if r[1]=='sample'])
write('selected_reads.tsv',['sample-id','r1','r2'],[[r[0],r[6],r[7]] for r in rows])
print('Retenus:',len(rows),'dont controles:',sum(r[1]=='negative_control' for r in rows))
PY
qiime info > "$LOG_DIR/qiime-info.txt"
# Chaque etape est marquee seulement apres succes. Un repertoire partiel n'est
# jamais considere comme termine; le supprimer avant de relancer son etape.
qc(){
 [[ ! -e "$RESULTS_DIR/state/qc.done" ]] || { log 'QC deja terminee'; return; }
 [[ -f "$ADAPTERS_FA" ]] || die 'Fournir ADAPTERS_FA : FASTA des adaptateurs de la librairie utilisee'
 cp "$ADAPTERS_FA" "$M/library-adapters.fasta"
 activate "$FASTQC_ENV"
 while IFS=$'\t' read -r sid r1 r2; do
  [[ "$sid" == sample-id ]] && continue
  fastqc --threads "$THREADS" --outdir "$RESULTS_DIR/02_fastqc_raw" "$r1" "$r2"
 done < "$M/selected_reads.tsv"
 activate "$MULTIQC_ENV"
 multiqc "$RESULTS_DIR/02_fastqc_raw" -o "$RESULTS_DIR/02_fastqc_raw/multiqc" --force
 activate "$TRIMMOMATIC_ENV"
 # DADA2 gere la qualite : pas de SLIDINGWINDOW agressif ici.
 while IFS=$'\t' read -r sid r1 r2; do
  [[ "$sid" == sample-id ]] && continue
  trimmomatic -Xmx"$TRIMMOMATIC_HEAP" PE -threads "$THREADS" -phred33 \
   "$r1" "$r2" \
   "$RESULTS_DIR/03_trimmomatic/${sid}_R1_P.fastq.gz" "$RESULTS_DIR/03_trimmomatic/${sid}_R1_U.fastq.gz" \
   "$RESULTS_DIR/03_trimmomatic/${sid}_R2_P.fastq.gz" "$RESULTS_DIR/03_trimmomatic/${sid}_R2_U.fastq.gz" \
   "ILLUMINACLIP:$M/library-adapters.fasta:2:30:10" MINLEN:100 \
   2> "$LOG_DIR/trimmomatic_${sid}.log"
 done < "$M/selected_reads.tsv"
 activate "$FASTQC_ENV"
 fastqc --threads "$THREADS" --outdir "$RESULTS_DIR/04_fastqc_clean" "$RESULTS_DIR"/03_trimmomatic/*_P.fastq.gz
 activate "$MULTIQC_ENV"
 multiqc "$RESULTS_DIR/04_fastqc_clean" "$LOG_DIR" -o "$RESULTS_DIR/04_fastqc_clean/multiqc" --force
 activate "$QIIME2_ENV"
 python - "$RESULTS_DIR" <<'PY'
import sys,pathlib,csv,zipfile
p=pathlib.Path(sys.argv[1])
with (p/'09_exports/fastqc_modules.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t'); w.writerow(['stage','file','module','status'])
 for stage,folder in [('raw','02_fastqc_raw'),('clean','04_fastqc_clean')]:
  for z in sorted((p/folder).glob('*_fastqc.zip')):
   with zipfile.ZipFile(z) as a:
    key=next(x for x in a.namelist() if x.endswith('/summary.txt'))
    for line in a.read(key).decode().splitlines():
     status,module,name=line.split('\t'); w.writerow([stage,name,module,status])
with (p/'00_metadata/manifest_clean.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['sample-id','forward-absolute-filepath','reverse-absolute-filepath'])
 for row in csv.DictReader((p/'00_metadata/selected_reads.tsv').open(),delimiter='\t'):
  sid=row['sample-id']; w.writerow([sid,*[str(p/'03_trimmomatic'/f'{sid}_R{m}_P.fastq.gz') for m in (1,2)]])
PY
 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$M/manifest_clean.tsv" --input-format PairedEndFastqManifestPhred33V2 --output-path "$Q/imported.qza"
 if [[ "$PRIMERS_PRESENT" == YES ]]; then
  qiime cutadapt trim-paired --i-demultiplexed-sequences "$Q/imported.qza" \
--p-front-f "$FWD" --p-front-r "$REV" \
--p-error-rate 0.1 --p-overlap 15 \
--p-discard-untrimmed --p-cores "$QIIME_THREADS" \
   --o-trimmed-sequences "$Q/no-primers.qza" 2>&1 | tee "$LOG_DIR/cutadapt.log"
 else
  cp "$Q/imported.qza" "$Q/no-primers.qza"
 fi
 qiime demux summarize --i-data "$Q/no-primers.qza" --o-visualization "$Q/no-primers.qzv"
 printf '%s\n' "$PRIMERS_PRESENT" > "$M/primers_present.txt"
 touch "$RESULTS_DIR/state/qc.done"
 log 'Inspecter no-primers.qzv puis definir TRUNC_F/TRUNC_R. Les amorces sont supposees en 5 prime sans spacer.'
}
denoise(){
 [[ -e "$RESULTS_DIR/state/qc.done" ]] || die 'Lancer qc avant dada2'
 [[ ! -e "$RESULTS_DIR/state/dada2.done" ]] || { log 'DADA2 deja terminee'; return; }
 [[ "$TRUNC_F" =~ ^[0-9]+$ && "$TRUNC_R" =~ ^[0-9]+$ ]] || die 'Definir TRUNC_F et TRUNC_R apres inspection (0 = pas de troncature fixe)'
 if (( TRUNC_F>0 && TRUNC_R>0 && TRUNC_F+TRUNC_R<445 )); then
  die 'Somme des longueurs <445 : recouvrement V3-V4 probablement insuffisant; revoir la qualite / strategie'
 fi
 activate "$QIIME2_ENV"
 qiime dada2 denoise-paired --i-demultiplexed-seqs "$Q/no-primers.qza" \
  --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f "$TRUNC_F" --p-trunc-len-r "$TRUNC_R" \
  --p-max-ee-f 2 --p-max-ee-r 2 --p-n-threads "$QIIME_THREADS" --p-chimera-method consensus \
  --o-table "$Q/table.qza" --o-representative-sequences "$Q/rep-seqs.qza" --o-denoising-stats "$Q/denoising-stats.qza"
 qiime metadata tabulate --m-input-file "$Q/denoising-stats.qza" --o-visualization "$Q/denoising-stats.qzv"
 qiime feature-table summarize --i-table "$Q/table.qza" --m-sample-metadata-file "$M/samples.tsv" --o-visualization "$Q/table.qzv"
 python - "$Q/table.qza" "$C" <<'PY'
import sys,pathlib,pandas as pd
from qiime2 import Artifact
x=Artifact.load(sys.argv[1]).view(pd.DataFrame); out=pathlib.Path(sys.argv[2])
if 'CTRL-STERILE' not in x.index:
 raise SystemExit('Controle absent apres DADA2 : verifier les stats, ne pas assimiler une perte du controle a une absence de contamination.')
ctrl=x.loc['CTRL-STERILE']; bad=ctrl[ctrl>0]
bad.rename('control_reads').to_csv(out/'control_ASVs_counts.tsv',sep='\t',index_label='feature-id')
x.sum(axis=1).rename('reads_before').to_csv(out/'reads_before.tsv',sep='\t',index_label='sample-id')
x.loc[:,~x.columns.isin(bad.index)].sum(axis=1).rename('reads_after_ASV_removal').to_csv(out/'reads_after.tsv',sep='\t',index_label='sample-id')
if len(bad)==x.shape[1]: raise SystemExit('Tous les ASV sont presents dans le controle: analyse arretee.')
print('ASV a enlever:',len(bad),'/',x.shape[1])
PY
 # Pas de liste vide transmise au filtre QIIME.
 if [[ $(wc -l < "$C/control_ASVs_counts.tsv") -gt 1 ]]; then
  qiime feature-table filter-features --i-table "$Q/table.qza" --m-metadata-file "$C/control_ASVs_counts.tsv" \
   --p-exclude-ids --o-filtered-table "$C/table-no-control-ASVs.qza"
 else
  cp "$Q/table.qza" "$C/table-no-control-ASVs.qza"
 fi
 qiime feature-table filter-samples --i-table "$C/table-no-control-ASVs.qza" --m-metadata-file "$M/biological.tsv" \
  --p-min-frequency 1 --o-filtered-table "$C/table-bio.qza"
 qiime feature-table filter-seqs --i-data "$Q/rep-seqs.qza" --i-table "$C/table-bio.qza" --o-filtered-data "$C/rep-seqs-bio.qza"
 qiime feature-table summarize --i-table "$C/table-bio.qza" --m-sample-metadata-file "$M/biological.tsv" --o-visualization "$C/table-bio.qzv"
 printf '%s\t%s\n' "$TRUNC_F" "$TRUNC_R" > "$M/dada2_truncation.tsv"
 touch "$RESULTS_DIR/state/dada2.done"
}
taxonomy(){
 [[ -e "$RESULTS_DIR/state/dada2.done" ]] || die 'DADA2 requise'
 [[ ! -e "$RESULTS_DIR/state/taxonomy.done" ]] || { log 'Taxonomie deja terminee'; return; }
 activate "$QIIME2_ENV"
 mkdir -p "$T/reference"
 local ref="$T/reference"
 curl --fail --location --retry 4 "$SILVA_URL" -o "$ref/silva.fasta.gz.part"
 gzip -t "$ref/silva.fasta.gz.part"
 mv "$ref/silva.fasta.gz.part" "$ref/silva.fasta.gz"
 sha256sum "$ref/silva.fasta.gz" > "$ref/sha256.txt"
 printf '%s\n' "$SILVA_URL" > "$ref/source-url.txt"
 # Import direct : compatible avec l'ancien environnement QIIME2, sans exiger
 # que RESCRIPt 2025.7 connaisse le numero de version 144.
 # Taxonomie native SILVA complete issue des en-tetes; pas de normalisation
 # artificielle en sept rangs. U->T; sequences ambigues ou courtes eliminees.
 python - "$ref" <<'PY'
import gzip,sys,pathlib
from Bio import SeqIO
p=pathlib.Path(sys.argv[1]); kept=0; skipped=0
with gzip.open(p/'silva.fasta.gz','rt') as src,(p/'dna.fasta').open('w') as fa,(p/'taxonomy.tsv').open('w') as tax:
 for r in SeqIO.parse(src,'fasta'):
  s=str(r.seq).upper().replace('U','T').replace('.','').replace('-','')
  lineage=r.description.split(None,1)
  if len(s)<900 or set(s)-set('ACGT') or len(lineage)!=2:
   skipped+=1; continue
  fa.write(f'>{r.id}\n{s}\n'); tax.write(f'{r.id}\t{lineage[1]}\n'); kept+=1
if not kept: raise SystemExit('Aucune sequence reference importable')
(p/'curation_counts.tsv').write_text(f'kept\tskipped\n{kept}\t{skipped}\n')
PY
 qiime tools import --type 'FeatureData[Sequence]' --input-path "$ref/dna.fasta" --output-path "$ref/dna.qza"
 qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat \
  --input-path "$ref/taxonomy.tsv" --output-path "$ref/taxonomy.qza"
 qiime feature-classifier extract-reads --i-sequences "$ref/dna.qza" --p-f-primer "$FWD" --p-r-primer "$REV" \
  --p-min-length 350 --p-max-length 500 --p-n-jobs "$QIIME_THREADS" --o-reads "$ref/v34.qza"
 qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads "$ref/v34.qza" \
  --i-reference-taxonomy "$ref/taxonomy.qza" --o-classifier "$ref/classifier-v34.qza"
 qiime feature-classifier classify-sklearn --i-classifier "$ref/classifier-v34.qza" --i-reads "$C/rep-seqs-bio.qza" \
  --p-n-jobs "$QIIME_THREADS" --p-read-orientation same --o-classification "$T/taxonomy.qza"
 qiime metadata tabulate --m-input-file "$T/taxonomy.qza" --o-visualization "$T/taxonomy.qzv"
 # Analyse bacteries + archees, hors organites et eucaryotes.
 qiime taxa filter-table --i-table "$C/table-bio.qza" --i-taxonomy "$T/taxonomy.qza" \
  --p-include 'Bacteria,Archaea' --p-exclude 'mitochondria,chloroplast' --o-filtered-table "$T/table-prokaryotes.qza"
 qiime feature-table filter-seqs --i-data "$C/rep-seqs-bio.qza" --i-table "$T/table-prokaryotes.qza" --o-filtered-data "$T/rep-seqs.qza"
 qiime taxa barplot --i-table "$T/table-prokaryotes.qza" --i-taxonomy "$T/taxonomy.qza" \
  --m-metadata-file "$M/biological.tsv" --o-visualization "$T/taxa-barplot.qzv"
 qiime feature-table summarize --i-table "$T/table-prokaryotes.qza" --m-sample-metadata-file "$M/biological.tsv" --o-visualization "$T/table-prokaryotes.qzv"
 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$T/rep-seqs.qza" --p-n-threads "$QIIME_THREADS" \
  --o-alignment "$T/aligned.qza" --o-masked-alignment "$T/masked.qza" --o-tree "$T/unrooted.qza" --o-rooted-tree "$T/rooted.qza"
 for spec in "table-prokaryotes:table" "rep-seqs:sequences" "taxonomy:taxonomy"; do
  qiime tools export --input-path "$T/${spec%%:*}.qza" --output-path "$RESULTS_DIR/09_exports/${spec##*:}"
 done
 python - "$T/table-prokaryotes.qza" "$RESULTS_DIR/09_exports" <<'PY'
import sys,pathlib,pandas as pd
from qiime2 import Artifact
x=Artifact.load(sys.argv[1]).view(pd.DataFrame); p=pathlib.Path(sys.argv[2])
x.to_csv(p/'ASV_counts_samples_by_features.tsv',sep='\t',index_label='sample-id')
x.sum(axis=1).sort_values().rename('reads').to_csv(p/'depths_before_rarefaction.tsv',sep='\t',index_label='sample-id')
PY
 touch "$RESULTS_DIR/state/taxonomy.done"
}
diversity(){
 [[ -e "$RESULTS_DIR/state/taxonomy.done" ]] || die 'Taxonomie requise'
 [[ ! -e "$RESULTS_DIR/state/diversity.done" ]] || { log 'Diversite deja terminee'; return; }
 [[ "$DEPTH" =~ ^[1-9][0-9]*$ ]] || die 'DEPTH doit etre un entier positif choisi apres inspection des effectifs'
 activate "$QIIME2_ENV"
 python - "$T/table-prokaryotes.qza" "$DEPTH" "$D" <<'PY'
import sys,pathlib,pandas as pd
from qiime2 import Artifact
x=Artifact.load(sys.argv[1]).view(pd.DataFrame).sum(axis=1); depth=int(sys.argv[2]); out=pathlib.Path(sys.argv[3])
pd.DataFrame({'reads':x,'retained':x>=depth}).to_csv(out/'rarefaction_sample_retention.tsv',sep='\t',index_label='sample-id')
if (x>=depth).sum()<2: raise SystemExit('Moins de deux echantillons conserves a cette profondeur')
PY
 qiime diversity alpha-rarefaction --i-table "$T/table-prokaryotes.qza" --i-phylogeny "$T/rooted.qza" \
  --p-max-depth "$DEPTH" --p-steps 20 --p-iterations 10 --m-metadata-file "$M/biological.tsv" --o-visualization "$D/alpha-rarefaction.qzv"
 qiime diversity core-metrics-phylogenetic --i-table "$T/table-prokaryotes.qza" --i-phylogeny "$T/rooted.qza" \
  --p-sampling-depth "$DEPTH" --m-metadata-file "$M/biological.tsv" --output-dir "$D/core-metrics"
 qiime tools export --input-path "$D/core-metrics/rarefied_table.qza" --output-path "$RESULTS_DIR/09_exports/rarefied_table"
 printf '%s\n' "$DEPTH" > "$M/rarefaction_depth.txt"
 touch "$RESULTS_DIR/state/diversity.done"
}
log "Demarrage $PROJECT_NAME, etape $STAGE; correction des identites AVANT exclusions."
case "$STAGE" in
 qc) qc;; dada2) denoise;; taxonomy) taxonomy;; diversity) diversity;;
 all) qc; denoise; taxonomy; diversity;;
esac
log "Etape $STAGE terminee. Resultats : $RESULTS_DIR"
