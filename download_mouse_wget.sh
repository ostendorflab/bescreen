OUTPUT_FASTA="bescreen/resources/Ensembl/release-112/fasta/mus_musculus/dna/"
OUTPUT_GTF="bescreen/resources/Ensembl/release-112/gtf/mus_musculus/"
OUTPUT_VEP="bescreen/resources/Ensembl/release-112/variation/indexed_vep_cache/"

URL_FASTA="ftp://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
URL_GTF="ftp://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz"
URL_VEP="ftp://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/mus_musculus_vep_112_GRCm39.tar.gz"

mkdir -p "$OUTPUT_FASTA"
mkdir -p "$OUTPUT_GTF"
mkdir -p "$OUTPUT_VEP"

wget "$URL_FASTA" -P "$OUTPUT_FASTA"
gunzip "${OUTPUT_FASTA}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
wget "$URL_GTF" -P "$OUTPUT_GTF"
wget "$URL_VEP" -P "$OUTPUT_VEP"
