# Epimachinery blast 

Blast queries from: https://github.com/hputnam/Epi_Machine/tree/main. I used git clone to download this repo to my computer and then uploaded the query fastas to URI HPC Unity.

Download CDS file:

```
wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.cds.fna.gz
gunzip Pocillopora_acuta_HIv2.genes.cds.fna.gz
```

Combine query fastas:

```
cd /project/pi_hputnam_uri_edu/estrand/Ploidy_HoloInt/epimachinery_queries
cat * > combined.fasta

## check duplicate headers
grep "^>" combined.fasta | sort | uniq -d

>HIRA-201 peptide: ENSP00000263208 pep:protein_coding
>MST1-202 peptide: ENSP00000414287 pep:protein_coding
>PRKAA2-201 peptide: ENSP00000360290 pep:protein_coding
>STK4-202 peptide: ENSP00000361892 pep:protein_coding

## remove the duplicates (the sequences are true duplicates)
awk '/^>/ {if(seen[$0]++) skip=1; else skip=0} !skip' combined.fasta > combined_no_dups.fasta
grep "^>" combined_no_dups.fasta | sort | uniq -d

## count number of queries: 499
grep -c "^>" combined.fasta
```

Create meta file for group of proteins and protein name

```
echo "fasta_file,header" > ../fasta_headers.csv

for f in *; do
    [ -e "$f" ] || continue
    [[ "$f" == "combined.fasta" || "$f" == "combined_no_dups.fasta" ]] && continue
    awk -v file="$f" '/^>/ {print file "," substr($0,2)}' "$f"
done >> ../fasta_headers.csv
```



Run Blast on combined query fasta against Pacuta CDS file:

```
#!/bin/bash
#SBATCH --partition=uri-cpu
#SBATCH --job-name=tblastn_cds_orthologs
#SBATCH --output=logs/tblastn_cds_orthologs_%j.out
#SBATCH --error=logs/tblastn_cds_orthologs_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

############################
# User settings
############################
QUERY_PROT="/project/pi_hputnam_uri_edu/estrand/Ploidy_HoloInt/epimachinery_queries/combined_no_dups.fasta"
TARGET_CDS="/project/pi_hputnam_uri_edu/estrand/Ploidy_HoloInt/Pocillopora_acuta_HIv2.genes.cds.fna"
SPECIES_TAG="Pocillopora_acuta"
EVALUE="1e-5"

############################
# Load software
############################
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/blast

############################
# Output names
############################
DB_NAME="${SPECIES_TAG}_cds_db"
ARCHIVE="${SPECIES_TAG}.tblastn.asn"
TABLE_ALL="${SPECIES_TAG}.tblastn_all_hits.tsv"
TABLE_TOP="${SPECIES_TAG}.tblastn_top_hits.tsv"
FASTA_HIT_REGIONS="${SPECIES_TAG}.top_hit_regions.fna"
FASTA_FULL_CDS="${SPECIES_TAG}.top_hit_full_cds.fna"

############################
# 1. Build BLAST database
############################
makeblastdb \
    -in "${TARGET_CDS}" \
    -dbtype nucl \
    -parse_seqids \
    -out "${DB_NAME}"

############################
# 2. Search query proteins vs CDS database
############################
tblastn \
    -query "${QUERY_PROT}" \
    -db "${DB_NAME}" \
    -evalue "${EVALUE}" \
    -num_threads "${SLURM_CPUS_PER_TASK}" \
    -outfmt 11 \
    -out "${ARCHIVE}"

############################
# 3. Convert archive to tabular format
############################
blast_formatter \
    -archive "${ARCHIVE}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
    -out "${TABLE_ALL}"

############################
# 4. Keep best hit per query
############################
sort -k1,1 -k12,12gr -k11,11g -k4,4gr "${TABLE_ALL}" | \
awk '!seen[$1]++' > "${TABLE_TOP}"

############################
# 5. Extract aligned nucleotide region
############################
: > "${FASTA_HIT_REGIONS}"

while IFS=$'\t' read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe
do
    if [[ "${sstart}" -le "${send}" ]]; then
        RANGE="${sstart}-${send}"
        STRAND="plus"
    else
        RANGE="${send}-${sstart}"
        STRAND="minus"
    fi

    blastdbcmd \
        -db "${DB_NAME}" \
        -entry "${sseqid}" \
        -range "${RANGE}" \
        -strand "${STRAND}" \
        -outfmt "%f" | \
    awk -v q="${qseqid}" -v s="${sseqid}" -v r="${RANGE}" '
        /^>/ {print ">" q "|" s "|region:" r; next}
        {print}
    ' >> "${FASTA_HIT_REGIONS}"

done < "${TABLE_TOP}"

############################
# 6. Extract full CDS for each top hit
############################
: > "${FASTA_FULL_CDS}"

while IFS=$'\t' read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe
do
    blastdbcmd \
        -db "${DB_NAME}" \
        -entry "${sseqid}" \
        -outfmt "%f" | \
    awk -v q="${qseqid}" -v s="${sseqid}" '
        /^>/ {print ">" q "|" s; next}
        {print}
    ' >> "${FASTA_FULL_CDS}"

done < "${TABLE_TOP}"

echo "Finished."
echo "Archive file: ${ARCHIVE}"
echo "All hits table: ${TABLE_ALL}"
echo "Top hits table: ${TABLE_TOP}"
echo "Hit regions fasta: ${FASTA_HIT_REGIONS}"
echo "Full CDS fasta: ${FASTA_FULL_CDS}"
```

The best hit is chosen by:
- k1,1 → query ID (qseqid)  
- k12,12gr → within each query, sort by bitscore descending  
- k11,11g → then by evalue ascending  
- k4,4gr → then by alignment length descending  

