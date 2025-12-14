import sys
import polars as pl
import shared
import copy

def get_variant_from_protein(transcript,
                             mutation,
                             bescreen_annotation,
                             bescreen_ref_genome,
                             snvs_only):

    codon_sun_one_letter = copy.deepcopy(shared.codon_sun_one_letter)

    # generate genome cds transformational information
    cds = '' # full cds for translation
    chromosomes = [] # chromosome of transcript (== chromosome of gene)
    genome_ranges = [] # base position on genome
    strands = [] # strand of transcript (== strand of gene)
    position_in_cds = 0 # position of first base
    exon_ranges = [] # base positions in cds

    if transcript == '':
        return ['no_input_gene_given']

    if transcript == '-':
        return ['no_input_transcript_given']

    if mutation == '':
        return ['no_input_mutation_given']

    if len(transcript.split('-')) == 2:
        transcript_exons = bescreen_annotation.filter(pl.col('transcript_name') == transcript).sort('exon_number') # maybe they are already sorted, but make sure
        if transcript_exons.is_empty():
            return ['input_transcript_not_found']

    elif len(transcript.split('-')) == 1:
        transcript_exons = bescreen_annotation.filter((pl.col('gene_name') == transcript) & pl.col('MANE_Select')).sort('exon_number') # maybe they are already sorted, but make sure
        if transcript_exons.is_empty():
            return ['input_gene_not_found']

    for row in transcript_exons.iter_rows(named=True):
        chromosomes.append(row['Chromosome']) # chromosome name
        genome_ranges.append(range(row['Start'], row['End'])) # range of this cds exon in genome
        strand = row['Strand']
        strands.append(strand)
        position_start = position_in_cds
        if strand == '+':
            cds_exon = str(bescreen_ref_genome[row['Chromosome']][row['Start']:row['End']])
        elif strand == '-':
            cds_exon = shared.revcom(str(bescreen_ref_genome[row['Chromosome']][row['Start']:row['End']]))
        cds += cds_exon
        position_in_cds += len(cds_exon) # end of this cds exon and start of next cds exon
        exon_ranges.append(range(position_start, position_in_cds)) # range of this cds exon in cds

    # translate full cds to protein
    aaseq = '' # full aa sequence

    for i in range(3, len(cds)+1, 3):
        aaseq += codon_sun_one_letter[cds[i-3:i]]

    # get genomic position of mutation
    # format needs to be: [one original aa][position in protein][one mutated aa]

    aa_ref = mutation[0]
    if not aa_ref in codon_sun_one_letter.values():
        return ['reference_not_amino_acid']
    aa_mut = mutation[-1]
    if not aa_mut in codon_sun_one_letter.values():
        return ['mutation_not_amino_acid']
    position_mut = mutation[1:-1]
    if position_mut.isnumeric():
        position_mut = int(position_mut)
    else:
        return ['input_position_not_numeric']
    alt_codons = [k for k,v in codon_sun_one_letter.items() if v == aa_mut]

    exon_found = False

    if aaseq[position_mut-1] == aa_ref: # check, reference aa is in transcript; cds and aaseq only necessary for this; simplify?
        start_cds = (position_mut-1) * 3 # transform bases to aa
        end_cds = ((position_mut-1) * 3) + 3 # transform bases to aa

        while not exon_found:

            for index, exon_range in enumerate(exon_ranges): # find exon with overlap; does this also find partial overlaps?

                if set(range(start_cds, end_cds)).issubset(set(exon_range)):
                    index_found = index
                    exon_range_found = exon_range
                    exon_found = True
                    break

        if strands[index_found] == '+':
            genome_start = genome_ranges[index_found][0] + (start_cds-exon_range_found[0])
            genome_end = genome_ranges[index_found][0] + (end_cds-exon_range_found[0])
            codon_genome = str(bescreen_ref_genome[chromosomes[index_found]][genome_start:genome_end])

        elif strands[index_found] == '-':
            genome_start = (genome_ranges[index_found][-1]+1) - (end_cds-exon_range_found[0])
            genome_end = (genome_ranges[index_found][-1]+1) - (start_cds-exon_range_found[0])
            codon_genome = shared.revcom(str(bescreen_ref_genome[chromosomes[index_found]][genome_start:genome_end]))

        generated_variants = []

        for alt_codon in alt_codons:
            edited_codon = []

            for index_base, base in enumerate(alt_codon):
                if codon_genome[index_base] != base:
                    if strands[index_found] == '+':
                        edited_codon.append(f'{chromosomes[index_found]}_{genome_start+index_base+1}_{codon_genome[index_base]}_{base}')
                    elif strands[index_found] == '-':
                        edited_codon.append(f'{chromosomes[index_found]}_{genome_end-index_base}_{shared.revcom(codon_genome[index_base])}_{shared.revcom(base)}')


            generated_variants.append(edited_codon)

        if snvs_only:
            return [''.join(generated_variant) for generated_variant in generated_variants if len(generated_variant) == 1] # list of snvs leading to desired aa change

        else:
            return generated_variants # list of variants leading to desired aa change

    else:
        return ['wrong_reference_amino_acid']


if __name__ == "__main__":
    pass
