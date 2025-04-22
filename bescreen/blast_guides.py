import os
import subprocess
import polars as pl
import io
import shlex


def check_blastdb(reference_genome, remake_blastdb):

    extensions = ['.njs',
                  '.ntf',
                  '.nto',
                  '.not',
                  '.ndb',
                  '.nos',
                  '.nog',
                  '.nsq',
                  '.nhr',
                  '.nin']
    
    blastfb_files = [reference_genome + extension for extension in extensions]
    for blastfb_file in blastfb_files:
        if not os.path.isfile(blastfb_file):
            remake_blastdb = True
    
    if remake_blastdb:
        make_blastdb(reference_genome)


def make_blastdb(reference_genome):

    cmd = ['makeblastdb',
           '-in', reference_genome,
           '-dbtype', 'nucl',
           '-parse_seqids'] # this will fail, if there are any spaces in the path; use symlink without spaces in that case

    process_blastdb = subprocess.Popen(cmd)
    process_blastdb.wait()


def guide_blast(guides,
                guidelength,
                reference_genome,
                mode,
                main_chroms_only): # exclude contigs

    # can the if-elif blocks be simplified?
    if mode == 'variants':
        guides = guides.rename({'index': 'indexvar'})
        guides = guides.with_columns(indexguide = pl.int_range(pl.len()).over('indexvar'))
        guides = guides.with_columns(pl.col('indexvar').cast(pl.String))
        guides = guides.with_columns(pl.col('indexguide').cast(pl.String))
        guides = guides.with_columns(indexboth = pl.concat_str([pl.col('indexvar'), pl.col('indexguide')], separator='_'))
        guides = guides.with_columns(indexboth_blast = pl.concat_str('>' + pl.col('indexboth')))
        guides_valid = guides.filter(~pl.col('guide').is_in(['be_not_usable', 'no_guides_found']))
        if guides_valid.is_empty(): # no guides found at all
            blast_results = pl.DataFrame({'blastcount': []})
            return (blast_results)
        guides_toblast = '\n'.join(guides_valid.with_columns(pl.concat_str([pl.col('indexboth_blast'), pl.col('guide')], separator='\n').alias("guides_toblast"))["guides_toblast"].to_list())

    elif mode == 'genes':
        guides = guides.with_columns(pl.col('index').cast(pl.String))
        guides = guides.with_columns(index_blast = pl.concat_str('>' + pl.col('index'))) # change to (index = pl.concat_str(['>' + pl.col('index')], separator='_')) # maybe with []?
        guides_toblast = '\n'.join(guides.with_columns(pl.concat_str([pl.col('index_blast'), pl.col('guide')], separator='\n').alias("guides_toblast"))["guides_toblast"].to_list())

    cmd_blast = ['blastn',
                 '-db', reference_genome,
                 '-task', 'blastn',
                 '-dust', 'no',
                 '-word_size', str(guidelength),
                 '-outfmt', '6']

    process_blast = subprocess.Popen(cmd_blast, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    blast_stdout = process_blast.communicate(input=guides_toblast.encode('utf-8'))[0].decode('utf-8').splitlines()
    blast_csv = io.StringIO()

    header_outfmt6 = ['qseqid',
                      'sseqid',
                      'pident',
                      'length',
                      'mismatch',
                      'gapopen',
                      'qstart',
                      'qend',
                      'sstart',
                      'send',
                      'evalue',
                      'bitscore'] # this is fixed; better use a tuple? could also use dtypes_outfmt6.keys(), but would that be sorted?

    dtypes_outfmt6 = {'qseqid': pl.String,
                      'sseqid': pl.String,
                      'pident': pl.Float64,
                      'length': pl.Int64,
                      'mismatch': pl.Int64,
                      'gapopen': pl.Int64,
                      'qstart': pl.Int64,
                      'qend': pl.Int64,
                      'sstart': pl.Int64,
                      'send': pl.Int64,
                      'evalue': pl.Float64,
                      'bitscore': pl.Float64}

    for line in blast_stdout:
        blast_csv.write(line + '\n')

    blast_csv.seek(0) # what does that do exactly? Does it set the cursor to line 0?
    blast_results = pl.read_csv(blast_csv,
                                separator='\t',
                                has_header=False,
                                new_columns=header_outfmt6,
                                schema_overrides=dtypes_outfmt6)
    blast_csv.close()

    # blast_results.to_csv('BLASTOUT.tsv', sep='\t') # for debugging

    if main_chroms_only: # exclude contigs
        blast_results = blast_results.filter(pl.col('sseqid').is_in(main_chroms))

    blast_results = blast_results['qseqid'].value_counts()

    # blast_results.to_csv('BLASTCOUNTOUT.tsv', sep='\t') # for debugging

    blast_results = blast_results.with_columns(pl.col('count').cast(pl.String)) # easier handling of strings then ints in mainscript

    if blast_results.is_empty():
        blast_results = pl.DataFrame({'blastcount': []})
        return (blast_results)

    # can the if-elif blocks be simplified?
    if mode == 'variants':
        blast_results = blast_results.rename({'qseqid': 'indexboth'})
        blast_results = guides.select('indexboth').join(blast_results, on='indexboth', how='left').with_columns(pl.col('count').fill_null(0))
        blast_results = blast_results.with_columns(pl.col('indexboth').str.split_exact('_', 1).struct.rename_fields(['indexvar', 'indexguide']).alias('fields')).unnest('fields')
        blast_results = blast_results.with_columns(pl.col('indexvar').cast(pl.UInt32))
        blast_results = blast_results.with_columns(pl.col('indexguide').cast(pl.UInt32))
        blast_results = blast_results.sort(['indexvar', 'indexguide'])
        blast_results = blast_results.drop('indexboth')
        blast_results = blast_results.group_by('indexvar', maintain_order=True).agg(pl.all())
        blast_results = blast_results.drop('indexguide') # can I do this earlier?
        blast_results = blast_results.rename({'count': 'blastcount'})

    elif mode == 'genes':
        blast_results = blast_results.rename({'qseqid': 'index'})
        if main_chroms_only: # if a guide if found on a contig only, but contigs are excluded for blasting
            blast_results = guides.select('index').join(blast_results, on='index', how='left').with_columns(pl.col('count').fill_null(0))
        blast_results = blast_results.with_columns(pl.col('index').cast(pl.UInt32))
        blast_results = blast_results.sort('index')
        blast_results = blast_results.rename({'count': 'blastcount'})

    return (blast_results)


main_chroms = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
               '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'MT',
               'X', 'Y', # Ensembl
               'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21',
               'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chrM', 'chrX', 'chrY', # UCSC
               'NC_000001.11', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12',
               'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10',
               'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000002.12',
               'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000003.12',
               'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14',
               'NC_000008.11', 'NC_000009.12', 'NC_012920.1', 'NC_000023.11',
               'NC_000024.10', # RefSeq; how to ignore behind dot?
               'CM000663.2', 'CM000672.2', 'CM000673.2', 'CM000674.2',
               'CM000675.2', 'CM000676.2', 'CM000677.2', 'CM000678.2',
               'CM000679.2', 'CM000680.2', 'CM000681.2', 'CM000664.2',
               'CM000682.2', 'CM000683.2', 'CM000684.2', 'CM000665.2',
               'CM000666.2', 'CM000667.2', 'CM000668.2', 'CM000669.2',
               'CM000670.2', 'CM000671.2', 'J01415.2', 'CM000685.2',
               'CM000686.2'] # GenBak; how to ignore behind dot?
# taken from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt


if __name__ == "__main__":
    pass
