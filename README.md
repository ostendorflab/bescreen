# How to use BEscreen
The easiest way is to use the public webserver at [bescreen.ostendorflab.org](https://bescreen.ostendorflab.org). If you want to use BEscreen locally on your computer follow this guide to use its command-line interface.

## Local installation
To run BEscreen, you need to run the script `bescreen.py` in a Python environment containing several modules. If you want to use the annotations of Ensembl's Variant Effect Predictor (VEP) and NCBI's BLAST you also need the command-line version of these tools ([VEP](https://www.ensembl.org/info/docs/tools/vep/index.html); [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)). The easiest way to generate such an environment is using `conda`. We provide multiple environment files for the above mentioned combinations:
- `environment_basic.yaml (bescreen_basic)`: The basic environment for the command line version.
- `environment_full.yaml (bescreen_full)`: The environment for the command line version with Ensembl's VEP and NCBI's BLAST.
- `environment_full_vepfixed.yaml (bescreen_full_vepfixed)`: The environment for the command line version with Ensembl's VEP with fixed version and NCBI's BLAST (use this if the one without version fixing doesn't work).

To install one of these environments run:
```
# install a conda distribution like mini-forge
git clone git@github.com:bepeople/bescreen.git
cd bescreen
# Install with conda:
conda env create -f envs/environment.yaml # replace environment.yaml with one of the above mentioned files
```
Before you use BEscreen you also need to activate this environment:
```
# Activate conda environment:
conda activate bescreen # replace bescreen with one of the above mentioned environment names in the parentheses
```

### Reference files
You also need to download (and install) some reference files from Ensembl's and NCBI's FTP servers. BEscreen was tested with the Ensembl release 112, but others should also work. If you are working with macOS or Linux, you can use from of the download scripts (`download_(species)_(wget/curl).sh`) to directly download the necessary files to the default locations. Most Linux distributions provide `wget` per default and macOS provides `curl`. So, if you don't have the other tool installed, simply use the `wget` scripts on Linux and the `curl` scripts on macOS. You need to execute the scripts in the root folder of the repository:
```
git clone git@github.com:bepeople/bescreen.git # if not already done
cd bescreen # if not already done
./download_(species)_(wget/curl).sh
```

The necessary files are:

#### The reference genome
- human: ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
- mouse: ftp://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
- Download the one you need and unzip it.
- For best compatibility store at `(bescreen/)bescreen/resources/Ensembl/release-112/fasta/(homo_sapiens|mus_musculus)/dna/`

#### The corresponding GTF annotation file
- human: ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
- mouse: ftp://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
- Download the one you need.
- For best compatibility store at `(bescreen/)bescreen/resources/Ensembl/release-112/gtf/(homo_sapiens|mus_musculus)/`
- To use the annotation file in IGV for viewing your guide BAM you need to sort the GTF (you need tabix for this):
```
bgzip -d gtf_file.gtf.gz > gtf_file.gtf
sort -k1,1 -k4,4n -s gtf_file.gtf > gtf_file.sorted.gtf
bgzip -c gtf_file.sorted.gtf > gtf_file.sorted.gtf.gz
tabix gtf_file.sorted.gtf.gz
```

#### The corresponding VEP cache
- human: ftp://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz
- mouse: ftp://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/mus_musculus_vep_112_GRCm39.tar.gz
- Download the one you need.
- For best compatibility store at `(bescreen/)bescreen/resources/Ensembl/variation/indexed_vep_cache/`
- Install it by executing this command in the download folder using your environment (the human cache takes around 27 GB):
```
# # general command
# vep_install --AUTO c \
# --SPECIES {species} \
# --ASSEMBLY {build} \
# --CACHE_VERSION {release} \
# --CACHEURL . \
# --CACHEDIR cache \
# --NO_UPDATE

# human
vep_install --AUTO c \
--SPECIES homo_sapiens \
--ASSEMBLY GRCh38 \
--CACHE_VERSION 112 \
--CACHEURL . \
--CACHEDIR cache \
--NO_UPDATE

# mouse
vep_install --AUTO c \
--SPECIES mus_musculus \
--ASSEMBLY GRCm39 \
--CACHE_VERSION 112 \
--CACHEURL . \
--CACHEDIR cache \
--NO_UPDATE
```

#### The dbSNP VCF
To use rsIDs from dbSNP as input you need to build a local dbSNP database. To do this first download the dbSNP VCF from here: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz to here `(bescreen/)bescreen/resources/dbSNP/snp/organisms/human_9606/VCF/`.

If you have done so use the the script `(bescreen/)bescreen/dbsnp_sqlite3.py` to create the local database by executing this command in the download folder using your environment:
```
git clone git@github.com:bepeople/bescreen.git # if not already done
cd bescreen # if not already done
python dbsnp_sqlite3.py --dbsnp-vcf resources/dbSNP/snp/organisms/human_9606/VCF/00-All.vcf.gz --rsid-db resources/dbSNP/snp/organisms/human_9606/VCF/rsid.db # adjust the paths according to your needs
```
This might take several hours and needs a lot of memory. If someone you know already generated this file, it might be best to ask them for the file. The generated file will be around 20 GB in size. For best compatibility store the generated `rsid.db` at `(bescreen/)bescreen/resources/dbSNP/snp/organisms/human_9606/VCF/`.

Now you should be good to go!

## Using the command line
To run BEscreen from the command line interface cd to `(bescreen/)bescreen/` and execute `python bescreen.py -h` for help (output of this see below):

BEscreen takes either variants (in form of genomic positions, rsIDs or protein amino acid changes) or whole genes as input. The input can be provided either directly or within a file. If you use a variant as input the output will contain all base editing guides for this variant. If you use a gene as input the output will be all editing and non-editing guides of this gene's CDS.

### Variants as input
Variant inputs need to be provided either as:
- genomic positions:<br>
  [chromosome]\_[genomic_position]\_[alternative_base] or [chromosome]\_[genomic_position]\_[reference_base]\_[alternative_base]<br>
  e.g.: 12_6537866_C_T or 12_6537866_T

- rsID:<br>
  e.g.: rs1062436

- protein amino acid changes<br>
  [gene_symbol]-[transcript_number]-[WTAApositionMUTAA] or [gene_symbol]-[transcript_number]-[WTAApositionMUTAA]<br>
  e.g.: GAPDH-L270F or GAPDH-201-L270F<br>
  If you don't provide a transcript number the MANE select transcript will be used.

These can be provided either directly using the `--variant (-v)` flag. If you provide multiple variants, they need to be separated by commas without any spaces. The variants can also be provided using a CSV file containing the column `variant` with the variants in one or a mixture of the above mentioned formats. Alternatively for file you can also use the genomic coordinates directly using the columns `chrom`, `pos`, and `alt` (optionally also using the column `ref`) with the genomic coordinates of the variants as used in a VCF file (just in minor letters for columns names). You cannot mix the columns `variant` and `chrom`, `pos`, (`ref`) and `alt`.

CSV files containing variants need to be formatted as either of the following examples:
- columns `variant`
    ```
    variant
    12_6537866_C_T
    12_6537866_T
    rs1062436
    GAPDH-L270F
    GAPDH-201-L270F
    ```

- columns `chrom`, `pos`, `ref` and `alt`:
    ```
    chrom,pos,ref,alt
    12,6537866,C,T
    ```

The reference base can be omitted in any case (in the file you can leave the column out completely), but then BEscreen will not check, if your reference base is the one it also finds. If these differ this hints to the usage of different reference genomes.

### Genes as input
Gene inputs need to be provided as HGNC gene symbol (e.g. GAPDH). These can be provided either directly using the `--gene-symbols (-g)` flag. If you provide multiple genes, they need to be separated by commas without any spaces. The genes can also be provided using a CSV file containing the column `symbol`.

CSV files containing genes need to be formatted as either of the following examples:
- columns `variant`:
    ```
    symbol
    GAPDH
    ```

### Regions as input
If you want to input regions directly, you need to use the form [CHROM]:[START]-[END] (e.g. like this: `12:6536490-6537490`) with the `--regions (-z)` fladg. If you choose to input from file you can use the column `region` as shown in the following example:
- columns `region`:
    ```
    region
    12:6536490-6537490
    ```


### Minimal command
The minimal command for variants is this:
```
python bescreen.py --ref-genome path/to/reference_genome.fa --output path/to/output_file_with_output_prefix --annotation-file path/to/annotation_file.gtf.gz --variant variant_one,variant_two,variant_three
```

The minimal command for genes is this:
```
python bescreen.py --ref-genome path/to/reference_genome.fa --output path/to/output_file_with_output_prefix --annotation-file path/to/annotation_file.gtf.gz --gene-symbols gene_one,gene_two,gene_three
```

The minimal command for regions is this:
```
python bescreen.py --ref-genome path/to/reference_genome.fa --output path/to/output_file_with_output_prefix --annotation-file path/to/annotation_file.gtf.gz --region region_one,region_two,region_three
```

### Help
More options can be found in the help:
```
usage: python bescreen.py [-h] -r <PATH> -o <PATH> [-t <PATH>] (-i <PATH> | -v VARIANT | -g GENE_SYMBOLS | -z REGIONS) [-p {NGG,NG}] [-s WINDOW_START] [-e WINDOW_END] [-l GUIDE_LENGTH] [-b {both,ABE,CBE}] [-q WINDOW_START_PLUS] [-w WINDOW_END_PLUS] [-d {collapsed,exploded}] [-x]
                          [--be-preset {default,A3A-BE3,ABE7.10,ABE7.10*,ABE7.9,BE-PLUS,BE1,BE2,BE3,BE3-R33A,BE4,BE4-Gam,BE4max,CP-ABEmax variants,eBE-S3,EE-BE3,evoAPOBEC1-BE4max,evoFERNY-BE4max,hA3A-eBE-Y130F,hA3A-eBE-Y132D,hA3A-eBE3,HF-BE3,NG-ABEmax,Target-AID,Target-AID-NG,xABE,xBE3,YE1-BE3,YE2-BE3,YEE-BE3}] [--filter-synonymous] [--filter-splice-site] [--filter-specific] [--filter-missense] [--write-bam]
                          [-n IGNORE_STRING] [-a] [-f {variant,vcf,symbol,region}] [-y] [--add-vep-annotations] [--vep-species VEP_SPECIES] [--vep-assembly VEP_ASSEMBLY] [--vep-dir_cache VEP_DIR_CACHE] [--vep-dir_plugins VEP_DIR_PLUGINS] [--vep-cache_version VEP_CACHE_VERSION] [--vep-flags VEP_FLAGS] [--dbsnp-db DBSNP_DB] [--add-blast-counts] [--blast-main-chromosomes-only]

bescreen

options:
  -h, --help            show this help message and exit
  -i <PATH>, --input <PATH>
                        Path to the input file
  -v VARIANT, --variant VARIANT
                        One or more variants you want to design guides for using variant strings ([chromosome]_[genomic_position]_[alternative_base] or [chromosome]_[genomic_position]_[reference_base]_[alternative_base]) or rsIDs or protein amino acid change ([gene_symbol]-[transcript_number]-[WTpositionMUT] or [gene_symbol]-[WTpositionMUT])
  -g GENE_SYMBOLS, --gene-symbols GENE_SYMBOLS
                        One or more HGNC gene symbols separated by commas
  -z REGIONS, --regions REGIONS
                        One or more regions of the form [CHROM]:[START]-[END] separated by commas
  -p {NGG,NG}, --pam-site {NGG,NG}
                        Sequence of the PAM site
  -s WINDOW_START, --window-start WINDOW_START
                        Starting position of editing window
  -e WINDOW_END, --window-end WINDOW_END
                        End position of editing window
  -l GUIDE_LENGTH, --guide-length GUIDE_LENGTH
                        Length of guide
  -b {both,ABE,CBE}, --base-editor {both,ABE,CBE}
                        Base editor to design guides for
  -q WINDOW_START_PLUS, --window-start-plus WINDOW_START_PLUS
                        Flanking region before the editing window for broader specificity
  -w WINDOW_END_PLUS, --window-end-plus WINDOW_END_PLUS
                        Flanking region after the editing window for broader specificity
  -d {collapsed,exploded}, --aspect {collapsed,exploded}
                        Collapsed (default) or exploded aspect of the output
  -x, --write-parquet   Force writing a new parquet file
  --be-preset {default,A3A-BE3,ABE7.10,ABE7.10*,ABE7.9,BE-PLUS,BE1,BE2,BE3,BE3-R33A,BE4,BE4-Gam,BE4max,CP-ABEmax variants,eBE-S3,EE-BE3,evoAPOBEC1-BE4max,evoFERNY-BE4max,hA3A-eBE-Y130F,hA3A-eBE-Y132D,hA3A-eBE3,HF-BE3,NG-ABEmax,Target-AID,Target-AID-NG,xABE,xBE3,YE1-BE3,YE2-BE3,YEE-BE3}
                        Use a known base editor to set the guide options
  --filter-synonymous   Pre-filter the output for synonymous guides
  --filter-splice-site  Pre-filter the output for splice_site guides
  --filter-specific     Pre-filter the output for specific guides
  --filter-missense     Pre-filter the output for missense guides
  --write-bam           Write a BAM file to view your guides in IGV (DO NOT USE THE SEQUENCES FOR YOUR LIBRARY SINCE REVERSE COMPLEMENT GUIDES WILL BE WRONG)
  -n IGNORE_STRING, --ignore-string IGNORE_STRING
                        Substring to be ignored in the variant string
  -a, --all-possible    Design for all ALTs
  -f {variant,vcf,symbol,region}, --input-format {variant,vcf,symbol,region}
                        Use column(s) "variant" (variant) or "chr", "pos", ("ref",) and "alt" (vcf) or "symbol" (symbol) in your input
  -y, --no-splice-sites
                        Exclude splice sites
  --add-vep-annotations, --vep
                        Add Ensembl's VEP annotations
  --vep-species VEP_SPECIES
                        Species for Ensembl's VEP annotations
  --vep-assembly VEP_ASSEMBLY
                        Assembly for Ensembl's VEP annotations
  --vep-dir_cache VEP_DIR_CACHE
                        Path to cache directory for Ensembl's VEP annotations
  --vep-dir_plugins VEP_DIR_PLUGINS
                        Path to plugin directory for Ensembl's VEP annotations
  --vep-cache_version VEP_CACHE_VERSION
                        Species for Ensembl's VEP annotations
  --vep-flags VEP_FLAGS
                        Modify the output of --add-vep-annotations (not all flags might work). Please provide the flags as string within quotes (e.g. --vep-flags '--mane --pubmed'). If you only want to use one flag you have to end the string with a space (e.g. --vep-flags '--everything '). The following flags should work: --sift [p|s|b], --polyphen [p|s|b], --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory,
                        --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna
  --dbsnp-db DBSNP_DB   Path to dbSNP database for rsID input
  --add-blast-counts, --blast
                        Add counts of genome wide (main and mitochondrial chromosomes only) perfect matches using NCBI's BLAST
  --blast-main-chromosomes-only, --no-contigs
                        BLAST only the main and mitochondrial chromosomes of the reference genome

required arguments:
  -r <PATH>, --ref-genome <PATH>
                        Path to the reference genome
  -o <PATH>, --output <PATH>
                        Path to the output file with output prefix
  -t <PATH>, --annotation-file <PATH>
                        Path to the GTF annotation file
```
