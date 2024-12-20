import subprocess
import polars as pl
import io
import shlex

def get_vep_annotation(variants,
                       species,
                       assembly,
                       dir_cache,
                    #    dir_plugins, # currently not in use
                       cache_version,
                       flags=''):

    variants_transformed = []
    for variant in variants:
        variant_split = variant.split('_')
        if len(variant_split) == 4:
            variant_transformed = variant_split[0] + '\t' + variant_split[1] + '\t.\t' + variant_split[2] + '\t' + variant_split[3]
        if len(variant_split) == 3:
            variant_transformed = variant_split[0] + '\t' + variant_split[1] + '\t.\t.\t' + variant_split[2]
        variants_transformed.append(variant_transformed)
    variants_vep = '\t.\t.\t.\n'.join(variants_transformed) + '\t.\t.\t.'

    cmd_vep = ['vep',
               '--species', species,
               '--assembly', assembly,
               '--dir_cache', dir_cache,
            #    '--dir_plugins', dir_plugin, # currently not in use
               '--cache',
               '--cache_version', cache_version,
               '--force',
               '--format', 'vcf',
               '--no_stats',
               '--offline',
               '--output_file', 'STDOUT',
               '--quiet',
               '--vcf'] + shlex.split(flags)


    process_vep = subprocess.Popen(cmd_vep, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    vep_stdout = process_vep.communicate(input=variants_vep.encode('utf-8'))[0].decode('utf-8').splitlines()
    vep_csv = io.StringIO()

    dtypes_vep = {'#CHROM': pl.String,
                  'POS': pl.Int64,
                  'ID': pl.String,
                  'REF': pl.String,
                  'ALT': pl.String,
                  'QUAL': pl.String,
                  'FILTER': pl.String,
                  'INFO': pl.String}

    for line in vep_stdout:
        if line.startswith('##INFO'):
            info = line
        elif not line.startswith('##'):
            vep_csv.write(line + '\n')

    vep_csv.seek(0) # what does that do exactly? Does it set the cursor to line 0?
    vep_annotation = pl.read_csv(vep_csv,
                                 separator='\t',
                                 schema_overrides=dtypes_vep) # .with_row_index()
    vep_csv.close()

    return (info, vep_annotation)


if __name__ == "__main__":
    pass
