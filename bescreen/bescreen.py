import os
import sys
import argparse
import polars as pl
import bedesigner
import besaturate
import beregion
import shared


def arguments():

    bescreendir = os.path.dirname(__file__)

    parser = argparse.ArgumentParser(description='bescreen',
                                     prog='python bescreen.py')

    # required
    required = parser.add_argument_group("required arguments")

    # for both
    required.add_argument('-r', '--ref-genome', help='Path to the reference genome',
                          metavar='<PATH>',
                        #   default='resources/Ensembl/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
                          required=True)
    required.add_argument('-o', '--output', help='Path to the output file with output prefix',
                          metavar='<PATH>', required=True)
    required.add_argument('-t', '--annotation-file', help='Path to the GTF annotation file',
                          metavar='<PATH>',
                        #   default='resources/Ensembl/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.sorted.gtf.gz',
                          required=False)


    # mutually exclusive
    input_group = parser.add_mutually_exclusive_group(required=True)
    # for both
    input_group.add_argument('-i', '--input', help='Path to the input file',
                             metavar='<PATH>')

    # for bedesigner
    input_group.add_argument('-v', '--variant', help='One or more variants you want to design guides for using variant strings ([chromosome]_[genomic_position]_[alternative_base] or [chromosome]_[genomic_position]_[reference_base]_[alternative_base]) or rsIDs or protein amino acid change ([gene_symbol]-[transcript_number]-[WTpositionMUT] or [gene_symbol]-[WTpositionMUT])',
                             type=str)

    # for besaturate
    input_group.add_argument('-g', '--gene-symbols', help='One or more HGNC gene symbols separated by commas',
                             type=str)

    # for beregion
    input_group.add_argument('-z', '--regions', help='One or more regions of the form [CHROM]:[START]-[END] separated by commas',
                             type=str)

    # optional
    # for both
    parser.add_argument('-p', '--pam-site', help='Sequence of the PAM site',
                        choices=['NGG', 'NG'], default='NG', type=str)
    parser.add_argument('-s', '--window-start', help='Starting position of editing window',
                        default=4, type=int)
    parser.add_argument('-e', '--window-end', help='End position of editing window',
                        default=8, type=int)
    parser.add_argument('-l', '--guide-length', help='Length of guide',
                        default=20, type=int)
    parser.add_argument('-b', '--base-editor', help='Base editor to design guides for',
                        choices=['both', 'ABE', 'CBE'], default='both', type=str)
    parser.add_argument('-q', '--window-start-plus', help='Flanking region before the editing window for broader specificity',
                        default=0, type=int)
    parser.add_argument('-w', '--window-end-plus', help='Flanking region after the editing window for broader specificity',
                        default=0, type=int)
    parser.add_argument('-d', '--aspect', help='Collapsed (default) or exploded aspect of the output',
                        choices=['collapsed', 'exploded'], default='collapsed', type=str)
    parser.add_argument('-x', '--write-parquet', help='Force writing a new parquet file',
                        action='store_true')
    base_editor_presets = shared.get_be_presets_dict(os.path.join(bescreendir, 'base_editors/be_presets.tsv'))
    parser.add_argument('--be-preset', help='Use a known base editor to set the guide options',
                        choices=list(base_editor_presets.keys()), type=str)
    parser.add_argument('--filter-synonymous', help='Pre-filter the output for synonymous guides',
                        action='store_true')
    parser.add_argument('--filter-splice-site', help='Pre-filter the output for splice_site guides',
                        action='store_true')
    parser.add_argument('--filter-specific', help='Pre-filter the output for specific guides',
                        action='store_true')
    parser.add_argument('--filter-missense', help='Pre-filter the output for missense guides',
                        action='store_true')
    parser.add_argument('--filter-nonsense', help='Pre-filter the output for nonsense (stop gain) guides',
                        action='store_true')
    parser.add_argument('--filter-stoplost', help='Pre-filter the output for stop lost guides',
                        action='store_true')
    parser.add_argument('--filter-startlost', help='Pre-filter the output for start lost guides',
                        action='store_true')
    parser.add_argument('--write-bam', help='Write a BAM file to view your guides in IGV (DO NOT USE THE SEQUENCES FOR YOUR LIBRARY SINCE REVERSE COMPLEMENT GUIDES WILL BE WRONG)',
                        action='store_true')

    # for bedesigner
    parser.add_argument('-n', '--ignore-string', help='Substring to be ignored in the variant string',
                        type=str)
    parser.add_argument('-a', '--all-possible', help='Design for all ALTs',
                        action='store_true')
    parser.add_argument('-f', '--input-format', help='Use column(s) "variant" (variant) or "chr", "pos", ("ref",) and "alt" (vcf) or "symbol" (symbol) in your input',
                        choices=['variant', 'vcf', 'symbol', 'region'], type=str)

    # for besaturate
    parser.add_argument('-y', '--no-splice-sites', help='Exclude splice sites',
                        action='store_false') # change this to store_true, but then adjust whole script

    # vep
    parser.add_argument('--add-vep-annotations', '--vep', help="Add Ensembl's VEP annotations",
                        action='store_true')
    parser.add_argument('--vep-species', help="Species for Ensembl's VEP annotations",
                        default='homo_sapiens', type=str)
    parser.add_argument('--vep-assembly', help="Assembly for Ensembl's VEP annotations",
                        default='GRCh38', type=str)
    parser.add_argument('--vep-dir_cache', help="Path to cache directory for Ensembl's VEP annotations",
                        default='resources/Ensembl/release-112/variation/indexed_vep_cache/cache', type=str)
    parser.add_argument('--vep-dir_plugins', help="Path to plugin directory for Ensembl's VEP annotations",
                        default='resources/Ensembl/release-112/variation/indexed_vep_cache/plugins', type=str) # currently not in use
    parser.add_argument('--vep-cache_version', help="Species for Ensembl's VEP annotations",
                        default='112', type=str)
    parser.add_argument('--vep-flags', help='Modify the output of --add-vep-annotations (not all flags might work). Please provide the flags as string within quotes (e.g. --vep-flags \'--mane --pubmed\'). If you only want to use one flag you have to end the string with a space (e.g. --vep-flags \'--everything \').\nThe following flags should work: --sift [p|s|b], --polyphen [p|s|b], --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna ', default='', type=str)

    # dbSNP
    parser.add_argument('--dbsnp-db', help="Path to dbSNP database for rsID input",
                        default='resources/dbSNP/snp/organisms/human_9606/VCF/rsid.db', type=str)

    # blast
    parser.add_argument('--add-blast-counts', '--blast', help="Add counts of genome wide (main and mitochondrial chromosomes only) perfect matches using NCBI's BLAST",
                        action='store_true')
    parser.add_argument('--blast-main-chromosomes-only', '--no-contigs', help="BLAST only the main and mitochondrial chromosomes of the reference genome",
                        action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if args.input and not os.path.isfile(args.input):
        sys.exit('Input file not found!')

    if not os.path.isfile(args.ref_genome):
        sys.exit('Reference genome not found!')

    if args.annotation_file and not os.path.isfile(args.annotation_file): # for besaturate (should be set as requried now)
        sys.exit('GTF annotation file not found!') # required for besaturate


    # input

    # optional
    # both
    pamsite = args.pam_site
    edit_window_start = args.window_start
    edit_window_end = args.window_end
    guidelength = args.guide_length
    baseeditor = args.base_editor
    aspect = args.aspect
    edit_window_start_plus = args.window_start_plus
    edit_window_end_plus = args.window_end_plus
    write_parquet = args.write_parquet
    filter_synonymous = args.filter_synonymous
    filter_splice_site = args.filter_splice_site
    filter_specific = args.filter_specific
    filter_missense = args.filter_missense
    filter_nonsense = args.filter_nonsense
    filter_stoplost = args.filter_stoplost
    filter_startlost = args.filter_startlost
    write_bam = args.write_bam
    # bedesigner
    ignorestring = args.ignore_string
    allpossible = args.all_possible
    input_format = args.input_format
    # besaturate
    splice_sites = args.no_splice_sites # False, if flag is used; True, if flag is not used

    # vep
    vep = args.add_vep_annotations
    vep_species = args.vep_species
    vep_assembly = args.vep_assembly
    vep_dir_cache = args.vep_dir_cache
    vep_dir_plugins = args.vep_dir_plugins
    vep_cache_version = args.vep_cache_version
    vep_flags = args.vep_flags

    # dbSNP
    dbsnp_db = args.dbsnp_db

    # blast
    blast = args.add_blast_counts
    no_contigs = args.blast_main_chromosomes_only

    # required
    # both
    ref_genome = args.ref_genome # named refgenome in bedesigner
    output_file = args.output + "_p" + pamsite + "_e" + \
                  str(edit_window_start) + "-" + str(edit_window_end) + \
                  "_l" + str(guidelength)
    # besaturate
    annotation_file = args.annotation_file

    # mutually exclusive
    # both
    input_file = args.input
    # bedesigner
    input_variant = args.variant
    # besaturate
    gene_symbols = args.gene_symbols
    # beregion
    regions = args.regions

    be_preset = args.be_preset
    if be_preset:
        sys_argvs = sys.argv
        if (not '-p' in sys_argvs) and (not '--pam-site' in sys_argvs):
            pamsite = base_editor_presets[be_preset]['pam']
        if (not '-s' in sys_argvs) and (not '--window-start' in sys_argvs):
            edit_window_start = base_editor_presets[be_preset]['window_start']
        if (not '-e' in sys_argvs) and (not '--window-end' in sys_argvs):
            edit_window_end = base_editor_presets[be_preset]['window_end']
        if (not '-l' in sys_argvs) and (not '--guide-length' in sys_argvs):
            guidelength = base_editor_presets[be_preset]['guide_length']
        if (not '-b' in sys_argvs) and (not '--base-editor' in sys_argvs):
            baseeditor = base_editor_presets[be_preset]['class']
        if (not '-q' in sys_argvs) and (not '--window-start-plus' in sys_argvs):
            edit_window_start_plus = base_editor_presets[be_preset]['plus_start']
        if (not '-w' in sys_argvs) and (not '--window-end-plus' in sys_argvs):
            edit_window_end_plus = base_editor_presets[be_preset]['plus_end']

    return (ref_genome,
            output_file,
            annotation_file,
            input_file,
            input_variant,
            gene_symbols,
            regions,
            pamsite,
            edit_window_start,
            edit_window_end,
            guidelength,
            baseeditor,
            vep,
            ignorestring, # deprecate!
            allpossible,
            input_format, # deprecate or evolve!
            edit_window_start_plus,
            edit_window_end_plus,
            splice_sites,
            write_parquet,
            aspect,
            vep_flags,
            vep_species,
            vep_assembly,
            vep_dir_cache,
            vep_dir_plugins,
            vep_cache_version,
            dbsnp_db,
            blast,
            no_contigs,
            filter_synonymous,
            filter_splice_site,
            filter_specific,
            filter_missense,
            filter_nonsense,
            filter_stoplost,
            filter_startlost,
            write_bam)


def check_mode(input_file, input_variant, gene_symbols, regions, input_format):

    mode = 'not_set'

    if input_file:
        if not input_format:

            input_columns = pl.read_csv(input_file, separator=',').columns

            if all(cols in input_columns for cols in ['variant', 'chr', 'pos', 'alt', 'symbol', 'region']): # 1
                sys.exit('Input file contains columns "variant", "chr", "pos", "alt", "symbol" and "region"!\nPlease use either "variant" or "chr", "pos"(, "ref") and "alt", "symbol" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['chr', 'pos', 'alt', 'symbol', 'region']): # 2
                sys.exit('Input file contains columns "chr", "pos", "alt", "symbol"  and "region"!\nPlease use either "chr", "pos"(, "ref") and "alt", "symbol" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'symbol', 'region']): # 3
                sys.exit('Input file contains columns "variant", "symbol" and "region"!\nPlease use either "variant", "symbol" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'chr', 'pos', 'alt', 'region']): # 4
                sys.exit('Input file contains columns "variant", "chr", "pos", "alt" and "region"!\nPlease use either "variant", "chr", "pos"(, "ref") and "alt" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'chr', 'pos', 'alt', 'symbol']): # 5
                sys.exit('Input file contains columns "variant", "chr", "pos", "alt" and "symbol"!\nPlease use either "variant" or "chr", "pos"(, "ref") and "alt" or "symbol" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['region', 'symbol']): # 6
                sys.exit('Input file contains columns "region" and "symbol"!\nPlease use either "region" or "symbol" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['chr', 'pos', 'alt', 'region']): # 7
                sys.exit('Input file contains columns "chr", "pos", "alt" and "region"!\nPlease use either "chr", "pos"(, "ref") and "alt" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['chr', 'pos', 'alt', 'symbol']): # 8
                sys.exit('Input file contains columns "chr", "pos", "alt" and "symbol"!\nPlease use either "chr", "pos"(, "ref") and "alt" or "symbol" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'region']): # 9
                sys.exit('Input file contains columns "variant" and "region"!\nPlease use either "variant" or "region" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'symbol']): # 10
                sys.exit('Input file contains columns "variant" and "symbol"!\nPlease use either "variant" or "symbol" or manually define --input-format!')

            elif all(cols in input_columns for cols in ['variant', 'chr', 'pos', 'alt']): # 11
                sys.exit('Input file contains columns "variant", "chr", "pos" and "alt"!\nPlease use either "variant" or "chr", "pos"(, "ref") and "alt" or manually define --input-format!')

            elif 'region' in input_columns: # 12
                mode = 'beregion'

            elif 'symbol' in input_columns: # 13
                mode = 'besaturate'

            elif all(cols in input_columns for cols in ['chr', 'pos', 'alt']) or 'variant' in input_columns: # 14 and 15
                mode = 'bedesigner'

        elif input_format in ['variant', 'vcf']:
            mode = 'bedesigner'
        elif input_format == 'symbol':
            mode = 'besaturate'
        elif input_format == 'region':
            mode = 'beregion'

    elif input_variant:
        mode = 'bedesigner'
    elif gene_symbols:
        mode = 'besaturate'
    elif regions:
        mode = 'beregion'

    return mode


if __name__ == "__main__":

    ref_genome_arg, \
    output_file_arg, \
    annotation_file_arg, \
    input_file_arg, \
    input_variant_arg, \
    gene_symbols_arg, \
    regions_arg, \
    pamsite_arg, \
    edit_window_start_arg, \
    edit_window_end_arg, \
    guidelength_arg, \
    baseeditor_arg, \
    vep_arg, \
    ignorestring_arg, \
    allpossible_arg, \
    input_format_arg, \
    edit_window_start_plus_arg, \
    edit_window_end_plus_arg, \
    splice_sites_arg, \
    write_parquet_arg, \
    aspect_arg, \
    vep_flags_arg, \
    vep_species_arg, \
    vep_assembly_arg, \
    vep_dir_cache_arg, \
    vep_dir_plugins_arg, \
    vep_cache_version_arg, \
    dbsnp_db_arg, \
    blast_arg, \
    no_contigs_arg, \
    filter_synonymous_arg, \
    filter_splice_site_arg, \
    filter_specific_arg, \
    filter_missense_arg, \
    filter_nonsense_arg, \
    filter_stoplost_arg, \
    filter_startlost_arg, \
    write_bam_arg = arguments()

    mode_arg = check_mode(input_file_arg, input_variant_arg, gene_symbols_arg, regions_arg, input_format_arg)

    if mode_arg == 'bedesigner':
        # try:
            guides, sam_df = bedesigner.design_bes(annotation_file_arg,
                                                   ref_genome_arg,
                                                   input_variant_arg,
                                                   input_file_arg,
                                                   pamsite_arg,
                                                   edit_window_start_arg,
                                                   edit_window_end_arg,
                                                   guidelength_arg,
                                                   ignorestring_arg,
                                                   baseeditor_arg,
                                                   edit_window_start_plus_arg,
                                                   edit_window_end_plus_arg,
                                                   allpossible_arg,
                                                   input_format_arg,
                                                   write_parquet_arg,
                                                   vep_arg,
                                                   aspect_arg,
                                                   vep_flags_arg,
                                                   vep_species_arg,
                                                   vep_assembly_arg,
                                                   vep_dir_cache_arg,
                                                   vep_dir_plugins_arg,
                                                   vep_cache_version_arg,
                                                   dbsnp_db_arg,
                                                   blast_arg,
                                                   no_contigs_arg,
                                                   filter_synonymous_arg,
                                                   filter_splice_site_arg,
                                                   filter_specific_arg,
                                                   filter_missense_arg,
                                                   filter_nonsense_arg,
                                                   filter_stoplost_arg,
                                                   filter_startlost_arg)

            bedesigner.output_sgrnas(guides, output_file_arg)

            if write_bam_arg:
                bedesigner.output_guides_sam(sam_df, output_file_arg, ref_genome_arg)
        # except Exception as error:
        #     sys.exit(str(error))

    elif mode_arg == 'besaturate':
        # try:
            guides, guides_ne, \
            sam_df, sam_ne_df, \
            not_found_genes = besaturate.saturate_bes(annotation_file_arg,
                                                      ref_genome_arg,
                                                      gene_symbols_arg,
                                                      input_file_arg,
                                                      pamsite_arg,
                                                      edit_window_start_arg,
                                                      edit_window_end_arg,
                                                      guidelength_arg,
                                                      baseeditor_arg,
                                                      edit_window_start_plus_arg,
                                                      edit_window_end_plus_arg,
                                                      splice_sites_arg,
                                                      write_parquet_arg,
                                                      vep_arg,
                                                      aspect_arg,
                                                      vep_flags_arg,
                                                      vep_species_arg,
                                                      vep_assembly_arg,
                                                      vep_dir_cache_arg,
                                                      vep_dir_plugins_arg,
                                                      vep_cache_version_arg,
                                                      blast_arg,
                                                      no_contigs_arg,
                                                      filter_synonymous_arg,
                                                      filter_splice_site_arg,
                                                      filter_specific_arg,
                                                      filter_missense_arg,
                                                      filter_nonsense_arg,
                                                      filter_stoplost_arg,
                                                      filter_startlost_arg)

            if not_found_genes: # redundant
                if len(not_found_genes) == 1:
                    print(f'WARNING: The gene {', '.join(not_found_genes)} was not found in the annotation file. It will not be mentioned in the output tables.')
                elif len(not_found_genes) > 1:
                    print(f'WARNING: The {len(not_found_genes)} genes {', '.join(not_found_genes)} were not found in the annotation file. They will not be mentioned in the output tables.')

            besaturate.output_sgrnas(guides, guides_ne, output_file_arg)

            if write_bam_arg:
                besaturate.output_guides_sam(sam_df, sam_ne_df, output_file_arg, ref_genome_arg)
        # except Exception as error:
        #     sys.exit(str(error))

    if mode_arg == 'beregion':
        # try:
            guides, guides_ne, \
            sam_df, sam_ne_df, \
            not_found_sequences = beregion.saturate_region(ref_genome_arg,
                                                           regions_arg,
                                                           input_file_arg,
                                                           pamsite_arg,
                                                           edit_window_start_arg,
                                                           edit_window_end_arg,
                                                           guidelength_arg,
                                                           baseeditor_arg,
                                                           edit_window_start_plus_arg,
                                                           edit_window_end_plus_arg,
                                                           vep_arg,
                                                           aspect_arg,
                                                           vep_flags_arg,
                                                           vep_species_arg,
                                                           vep_assembly_arg,
                                                           vep_dir_cache_arg,
                                                           vep_dir_plugins_arg,
                                                           vep_cache_version_arg,
                                                           blast_arg,
                                                           no_contigs_arg)

            if not_found_sequences: # redundant
                if len(not_found_sequences) == 1:
                    print(f'WARNING: The region {', '.join(not_found_sequences)} was not found in the reference genome. It will not be mentioned in the output tables.')
                elif len(not_found_sequences) > 1:
                    print(f'WARNING: The {len(not_found_sequences)} regions {', '.join(not_found_sequences)} were not found in the reference genome. They will not be mentioned in the output tables.')

            beregion.output_sgrnas(guides, guides_ne, output_file_arg)

            if write_bam_arg:
                beregion.output_guides_sam(sam_df, sam_ne_df, output_file_arg, ref_genome_arg)
        # except Exception as error:
        #     sys.exit(str(error))
