import os
import sys
import polars as pl
import numpy as np
import pyfaidx
import pysam
import get_vep
import shared
import dbsnp_sqlite3
import protein_variant
import blast_guides
import itertools


def design_bes(annotation_file,
               refgenome,
               input_variant,
               input_file,
               pamsite,
               fiveprimepam,
               edit_window_start,
               edit_window_end,
               guidelength,
               ignorestring,
               baseeditor,
               edit_window_start_plus,
               edit_window_end_plus,
               allpossible,
               input_format,
               mane_select_only,
               write_parquet,
               vep,
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
               filter_startlost):

    if guidelength < 17:
        raise ValueError('Please set the guide length to at least 17 bp!')

    bes = shared.bes

    if fiveprimepam:
         # swich fwd and rev for edits
        bes = shared.fiveprimepam_bes(bes)

         # switch to reverse PAM
        pamsite = shared.revcom(pamsite)

        # switch window start and end and count from end
        edit_window_start_new = guidelength - edit_window_end + 1
        edit_window_end = guidelength - edit_window_start + 1
        edit_window_start = edit_window_start_new

        # switch window start plus and end plus
        edit_window_start_plus_new = edit_window_end_plus
        edit_window_end_plus = edit_window_start_plus
        edit_window_start_plus = edit_window_start_plus_new

    # select base editors
    if not baseeditor:
        raise ValueError('Please select at least one base editor!')
    elif not 'all' in baseeditor:
        bes = {key: bes[key] for key in baseeditor}

    ref_genome_pyfaidx = pyfaidx.Fasta(refgenome)

    parquet_file = shared.check_parquet(annotation_file, write_parquet)

    cdss = pl.read_parquet(parquet_file)

    if mane_select_only:
        cdss = cdss.filter(pl.col('MANE_Select'))

    if input_file:

        input_df = pl.read_csv(input_file, separator=',')
        input_columns = input_df.columns

        if not input_format:

            if all(cols in input_columns for cols in ['variant', 'chr', 'pos', 'alt']):
                raise ValueError('Input file contains columns "variant", "chr", "pos" and "alt"!\nPlease use either "variant" or "chr", "pos"(, "ref") and "alt" or manually define --input-format!')

            elif 'variant' in input_columns:
                detected_mode = 'variant'

            elif all(cols in input_columns for cols in ['chr', 'pos', 'alt']):
                detected_mode = 'vcf'

            else:
                raise ValueError('Input file does not contain the column "variant" or the columns "chr", "pos"(, "ref") and "alt"!\nPlease provide either "variant" or "chr", "pos"(, "ref") and "alt"!')

            input_format = detected_mode

        if input_format == "variant":
            if 'variant' in input_columns:
                variant_list = input_df["variant"].to_list()
            else:
                raise ValueError('Input file does not contain the column "variant"!\nPlease provide it if selecting --input-format variant!')

        elif input_format == "vcf":
            if all(cols in input_columns for cols in ['chr', 'pos', 'alt']):
                if "ref" in input_columns:
                    variant_list = input_df.with_columns(variant = pl.concat_str([pl.col("chr"), pl.col("pos"), pl.col("ref"), pl.col("alt")], separator='_'))['variant'].to_list()
                else:
                    variant_list = input_df.with_columns(variant = pl.concat_str([pl.col("chr"), pl.col("pos"), pl.col("alt")], separator='_'))['variant'].to_list()
            else:
                raise ValueError('Input file does not contain the columns "chr", "pos"(, "ref") and "alt"!\nPlease provide them if selecting --input-format vcf!')

    elif input_variant:
        variant_list = input_variant.replace(' ', '').split(",")

    # final lists to fill after looping
    all_variant = []
    all_variant_real = []
    all_editable = []
    all_be_strings = []
    all_original_alt = []
    all_target_seq_ref = []
    all_target_seq_ref_match = []
    all_target_seq = []
    all_target_base_ref = []
    all_target_base = []
    all_possible_guides_with_pam = []
    all_edit_strings = []
    all_edit_pos_strings = []
    all_possible_guides = []
    all_possible_pams = []
    all_rev_com = []
    all_guide_starts = []
    all_guide_ends = []
    all_chroms = []

    # final lists for shared.analyze_guide()
    all_edit_window = []
    all_num_edits = []
    all_specific = []
    all_edit_window_plus = []
    all_num_edits_plus = []
    all_specific_plus = []
    all_safety_region = []
    all_num_edits_safety = []
    all_additional_in_safety = []
    all_specificity = []
    all_distance_median_variant = []
    # all_quality_scores_variant = []
    all_distance_median_all = []
    # all_quality_scores_all = []

    # final lists for manual annotation
    all_gene_symbols = []
    all_transcript_symbols = []
    all_exon_numbers = []
    all_first_transcript_exons = []
    all_last_transcript_exons = []
    all_codonss = []
    all_codonss_edited = []
    all_aass = []
    all_aass_edited = []
    all_aa_positionss = []
    all_splice_sites_included = []
    all_synonymouss = []
    all_consequences = []

    # rsIDs to transform
    rsids = []

    # identify rsIDs in variants
    for variant in variant_list:
        if len(variant.split("_")) == 1 and variant.startswith('rs'): # rsID
            rsids.append(variant)

    for rsid in rsids:
        # if rsid in variant_list: # this is not really necessary
        #     variant_list.remove(rsid)
        variant_list.remove(rsid) # remove rsIDs from variant_list

    if rsids:
        rsidvars = dbsnp_sqlite3.transform_locations(dbsnp_sqlite3.query_rsids(rsids,
                                                                               dbsnp_db))
        rsids_df = pl.DataFrame({'rsid': rsids})
        rsidvars = rsidvars.join(rsids_df, on='rsid', how='right')
        rsidvars = rsidvars.with_columns(pl.col("variant").fill_null('non_existent_input_rsID'))
        rsidvars = rsidvars.with_columns(rsid_variant = pl.concat_str([pl.col('rsid'), pl.col('variant')], separator=':'))['rsid_variant'].to_list() # better in the function?

        variant_list += rsidvars # readd the transformed rsIDs

    # transcript mutations to transform
    tmuts = []

    # identify transcript mutations in variants
    for variant in variant_list:
        if len(variant.split("_")) == 1 and len(variant.split("-")) in [2, 3]: # tmuts
            tmuts.append(variant)

    for tmut in tmuts:
        # if tmut in variant_list: # this is not really necessary
        #     variant_list.remove(tmut)
        variant_list.remove(tmut) # remove transcript mutations from variant_list

    if tmuts:
        tmutvars = []
        for tmut in tmuts:
            tr_mut_tmutvar = tmut.split("-")
            tr_tmutvar = '-'.join(tr_mut_tmutvar[0:-1])
            mut_tmutvar = tr_mut_tmutvar[-1]

            tmutvar = protein_variant.get_variant_from_protein(tr_tmutvar,
                                                               mut_tmutvar,
                                                               cdss,
                                                               ref_genome_pyfaidx,
                                                               True)

            for tmutsnp in tmutvar:
                tmutvars.append(f'{tmut}:{tmutsnp}')

        variant_list += tmutvars # readd the transformed transcript mutations

    # precalculation of quality scores for positions in editing window
    # distance_median_dict, quality_scores_dict = shared.qc_precalc(edit_window_start, edit_window_end)
    distance_median_dict = shared.qc_precalc(edit_window_start, edit_window_end)

    # only need to do this once
    non_stop_aas = list(set([x for x in shared.codon_sun_one_letter.values() if x != 'Stop']))
    any_filter = any([filter_synonymous,
                      filter_splice_site,
                      filter_specific,
                      filter_missense,
                      filter_nonsense,
                      filter_stoplost,
                      filter_startlost])

    for be in bes:

        edits = {
            'A': 'be_not_usable',
            'T': 'be_not_usable',
            'C': 'be_not_usable',
            'G': 'be_not_usable'
        }
        for direction in bes[be]:
            edits[bes[be][direction]['REF']] = bes[be][direction]['ALT']

        for variant in variant_list:
            variant_coords = variant.split("_")
            if len(variant_coords) == 4 and variant_coords[0].startswith('rs'):
                rsidchrom, position, ref, alt = variant_coords
                rsid, chrom = rsidchrom.split(":")
            elif len(variant_coords) == 4 and len(variant_coords[0].split("-")) in [2, 3]:
                tmutchrom, position, ref, alt = variant_coords
                tmut, chrom = tmutchrom.split(":")
            elif len(variant_coords) == 4: # REF given
                chrom, position, ref, alt = variant_coords
                alt = alt.upper()
                ref = ref.upper()
            elif len(variant_coords) == 3: # no REF given
                chrom, position, alt = variant_coords
                alt = alt.upper()
                ref = ""
            else:
                variant += ':variant_is_improperly_formatted'
                chrom, position, ref, alt = 'variant_is_improperly_formatted'.split("_")
            if ignorestring:
                chrom = chrom.replace(ignorestring, "")

            if not (variant.endswith('no_input_gene_given') or
                    variant.endswith('no_input_transcript_given') or
                    variant.endswith('no_input_mutation_given') or
                    variant.endswith('input_transcript_not_found') or
                    variant.endswith('input_gene_not_found') or
                    variant.endswith('reference_not_amino_acid') or
                    variant.endswith('mutation_not_amino_acid') or
                    variant.endswith('input_position_not_numeric') or
                    variant.endswith('wrong_reference_amino_acid') or
                    variant.endswith('non_existent_input_rsID') or
                    variant.endswith('variant_is_improperly_formatted')):
                if position.isnumeric():
                    position = int(position) - 1 # position based on VCF POS, which is 1-based
                else:
                    variant += ':genomic_position_not_numeric'
                    chrom, position, ref, alt = 'genomic_position_not_numeric'.split("_")

            length = guidelength # unnecessary

            position_edit_window_start = edit_window_start # unnecessary
            position_edit_window_end = edit_window_end # unnecessary

            pam = pamsite # unnecessary
            pam_length = len(pam) # 2 or 3

            pam_relevant = pam.lstrip('N') # needs to be rstrip() for 5' PAM

            pamlist = list(pam_relevant)
            pamlist_real = [shared.iupac_nt_code.get(item, item) for item in pamlist]
            pamlist_real_all = list(itertools.product(*pamlist_real))
            pamlist_real_string = [''.join(pam_real_string) for pam_real_string in pamlist_real_all]

            total_length = length + pam_length
            length_edit_window = position_edit_window_end - (position_edit_window_start - 1)

            bases_before_ew = position_edit_window_start - 1
            bases_after_ew = length - position_edit_window_end
            bases_after_ew_with_pam = total_length - position_edit_window_end

            bases_before_variant = bases_before_ew + length_edit_window - 1 # - 1 since variant always need to stay in the edit window
            bases_after_variant_with_pam = bases_after_ew_with_pam + length_edit_window - 1 # - 1 since variant always need to stay in the edit window

            if variant.endswith('variant_is_improperly_formatted'):
                target_base_ref = 'variant_is_improperly_formatted'
            elif variant.endswith('no_input_gene_given'):
                target_base_ref = 'no_input_gene_given'
            elif variant.endswith('no_input_transcript_given'):
                target_base_ref = 'no_input_transcript_given'
            elif variant.endswith('no_input_mutation_given'):
                target_base_ref = 'no_input_mutation_given'
            elif variant.endswith('input_transcript_not_found'):
                target_base_ref = 'input_transcript_not_found'
            elif variant.endswith('input_gene_not_found'):
                target_base_ref = 'input_gene_not_found'
            elif variant.endswith('reference_not_amino_acid'):
                target_base_ref = 'reference_not_amino_acid'
            elif variant.endswith('mutation_not_amino_acid'):
                target_base_ref = 'mutation_not_amino_acid'
            elif variant.endswith('input_position_not_numeric'):
                target_base_ref = 'input_position_not_numeric'
            elif variant.endswith('wrong_reference_amino_acid'):
                target_base_ref = 'wrong_reference_amino_acid'
            elif variant.endswith('non_existent_input_rsID'):
                target_base_ref = 'non_existent_input_rsID'
            elif variant.endswith('genomic_position_not_numeric'):
                target_base_ref = 'genomic_position_not_numeric'
            else:
                try:
                    target_base_ref = str(ref_genome_pyfaidx[chrom][position])
                except:
                    variant += ':genomic_coordinates_not_found'
                    chrom, position, ref, alt = 'genomic_coordinates_not_found'.split("_")
                    target_base_ref = 'genomic_coordinates_not_found'

            if ref:
                if ref == target_base_ref:
                    target_seq_ref_match = "ref_match"
                else:
                    target_seq_ref_match = "REF_MISMATCH"
            else:
                target_seq_ref_match = "no_ref_input"

            if (target_base_ref == bes[be]['fwd']['REF'] and alt == bes[be]['fwd']['ALT']) or \
               (target_base_ref == bes[be]['fwd']['REF'] and allpossible):
                editable = True
                rev_com = False
                be_string = be
                if alt == bes[be]['fwd']['ALT']:
                    original_alt = True
                else:
                    original_alt = False

                target_base = target_base_ref

                target_seq_ref_start = position - bases_before_variant
                target_seq_ref_end = position + bases_after_variant_with_pam + 1
                target_seq_ref = str(ref_genome_pyfaidx[chrom][target_seq_ref_start:target_seq_ref_end]) # + 1 since last position is excluded
                target_seq = target_seq_ref

            elif (target_base_ref == bes[be]['rev']['REF'] and alt == bes[be]['rev']['ALT']) or \
                 (target_base_ref == bes[be]['rev']['REF'] and allpossible):
                editable = True
                rev_com = True
                be_string = be
                if alt == bes[be]['rev']['ALT']:
                    original_alt = True
                else:
                    original_alt = False

                target_base = shared.revcom(target_base_ref)

                target_seq_ref_end = position - bases_after_variant_with_pam
                target_seq_ref_start = position + bases_before_variant + 1
                target_seq_ref = str(ref_genome_pyfaidx[chrom][target_seq_ref_end:target_seq_ref_start]) # + 1 since last position is excluded
                target_seq = shared.revcom(target_seq_ref)

            else:
                editable = False

            alt_edited = 'be_not_usable' if target_base_ref in ['variant_is_improperly_formatted',
                                                                'no_input_gene_given',
                                                                'no_input_transcript_given',
                                                                'no_input_mutation_given',
                                                                'input_transcript_not_found',
                                                                'input_gene_not_found',
                                                                'reference_not_amino_acid',
                                                                'mutation_not_amino_acid',
                                                                'input_position_not_numeric',
                                                                'wrong_reference_amino_acid',
                                                                'non_existent_input_rsID',
                                                                'genomic_position_not_numeric',
                                                                'genomic_coordinates_not_found'] else edits[target_base_ref] if allpossible else alt
            variant_real = variant
            rsidtmutvariant = variant_real.split(":")
            if (variant.endswith('no_input_gene_given') or
                variant.endswith('no_input_transcript_given') or
                variant.endswith('no_input_mutation_given') or
                variant.endswith('input_transcript_not_found') or
                variant.endswith('input_gene_not_found') or
                variant.endswith('reference_not_amino_acid') or
                variant.endswith('mutation_not_amino_acid') or
                variant.endswith('input_position_not_numeric') or
                variant.endswith('wrong_reference_amino_acid') or
                variant.endswith('non_existent_input_rsID') or
                variant.endswith('genomic_position_not_numeric') or
                variant.endswith('genomic_coordinates_not_found') or
                variant.endswith('variant_is_improperly_formatted')):
                variant_real = rsidtmutvariant[1]
            elif len(rsidtmutvariant) == 2 and (variant_real.startswith('rs') or len(rsidtmutvariant[0].split("-")) in [2, 3]):
                variant_real = rsidtmutvariant[1]
            if allpossible and target_base_ref not in ['variant_is_improperly_formatted',
                                                    'no_input_gene_given',
                                                    'no_input_transcript_given',
                                                    'no_input_mutation_given',
                                                    'input_transcript_not_found',
                                                    'input_gene_not_found',
                                                    'reference_not_amino_acid',
                                                    'mutation_not_amino_acid',
                                                    'input_position_not_numeric',
                                                    'wrong_reference_amino_acid',
                                                    'non_existent_input_rsID',
                                                    'genomic_position_not_numeric',
                                                    'genomic_coordinates_not_found']:
                variant_real = variant_real.split('_')
                if edits[target_base_ref] != 'be_not_usable':
                    variant_real[-1] = edits[target_base_ref]
                variant_real = '_'.join(variant_real)


            if editable:

                # lists per variant
                possible_guides_with_pam = []
                possible_guides = []
                possible_pams = []
                edit_strings = []
                edit_pos_strings = []
                possible_starts = []
                possible_ends = []
                possible_chroms = []

                # lists per variant for shared.analyze_guide()
                edit_windows = []
                num_editss = []
                specifics = []
                edit_window_pluss = []
                num_edits_pluss = []
                specific_pluss = []
                safety_regions = []
                num_edits_safetys = []
                additional_in_safetys = []
                specificitys = []
                distance_median_variants = []
                # quality_scores_variants = []
                distance_median_alls = []
                # quality_scores_alls = []

                # lists per variant for manual annotation
                gene_symbolss = []
                transcript_symbolss = []
                exon_numberss = []
                first_transcript_exonss = []
                last_transcript_exonss = []
                codonsss = []
                codonsss_edited = []
                aasss = []
                aasss_edited = []
                aa_positionsss = []
                splice_sitess_included = []
                synonymousss = []
                consequences = []

                count_ingored_N_in_pam = pam_length - len(pam_relevant)
                count_of_other_in_pam = pam_length - count_ingored_N_in_pam

                start_of_pam_search = bases_after_variant_with_pam - bases_after_ew - count_ingored_N_in_pam # the N of NGG/NG is irrelevant for the search here
                end_of_pam_search = count_of_other_in_pam - 1 # the N of NGG/NG is irrelevant for the search here; stop at beginning of last pam
                length_of_pam_search = count_of_other_in_pam # the N of NGG/NG is irrelevant for the search here

                start_before_hit = length + count_ingored_N_in_pam
                end_after_hit_with_pam = count_of_other_in_pam
                end_after_hit_guide_only = count_ingored_N_in_pam

                variant_position = bases_before_variant

                for i in range(len(target_seq) - start_of_pam_search,
                            len(target_seq) - end_of_pam_search):
                    if any(target_seq[i:i + length_of_pam_search] == pam_real_string for pam_real_string in pamlist_real_string):
                        possible_guide_with_pam = target_seq[i - start_before_hit:i + end_after_hit_with_pam]
                        possible_guide = target_seq[i - start_before_hit:i - end_after_hit_guide_only]
                        possible_pam = target_seq[i - count_ingored_N_in_pam:i + count_of_other_in_pam]

                        # shared.analyze_guide() per guide
                        edit_window, \
                        num_edits, \
                        specific, \
                        edit_window_plus, \
                        num_edits_plus, \
                        specific_plus, \
                        safety_region, \
                        num_edits_safety, \
                        additional_in_safety, \
                        edit_string, \
                        edit_pos_string, \
                        specificity, \
                        distance_median_variant, \
                        distance_median_all = shared.analyze_guide(possible_guide,
                                                                position_edit_window_start,
                                                                position_edit_window_end,
                                                                edit_window_start_plus,
                                                                edit_window_end_plus,
                                                                target_base,
                                                                variant_position,
                                                                distance_median_dict,
                                                                fiveprimepam)
                        # the following lines have been removed from the statement above:
                        # quality_scores_variant, \
                        # distance_median_all, \
                        # quality_scores_all = shared.analyze_guide(possible_guide,
                                                                #   distance_median_dict,
                                                                #   quality_scores_dict)

                        edit_string_an = edit_string
                        if rev_com:
                            edit_string_an = edit_string_an[::-1]

                        edit_list = []
                        variant_pos_rel = None
                        for i in range(len(edit_string_an)):
                            if edit_string_an[i] == '*':
                                edit_list.append(i)
                            elif edit_string_an[i] == 'V':
                                edit_list.append(i)
                                variant_pos_rel = i

                        total_poss = [position - (variant_pos_rel - m) for m in edit_list]

                        cdss_filtered = cdss.filter((pl.col('Chromosome') == chrom) &
                                                    (pl.col('Start') - 2 <= position) &
                                                    (position < pl.col('End') + 2)) # find CDSs with splice sites at the position of the variant

                        if not cdss_filtered.is_empty():

                            # lists per guide for manual annotation
                            gene_symbols = []
                            transcript_symbols = []
                            exon_numbers = []
                            first_transcript_exons = []
                            last_transcript_exons = []
                            codonss = []
                            codonss_edited = []
                            aass = []
                            aass_edited = []
                            aa_positionss = []
                            splice_sites_included = []
                            synonymouss = []

                            for row in cdss_filtered.iter_rows(named=True): # could be a lot since overlapping exons (multiple transcripts (common), overlapping CDS of different genes (rare),...)

                                # fixed annotations per overlapping exons for manual annotation
                                gene_symbol = row['gene_name']
                                transcript_symbol = row['transcript_name']
                                exon_number = row['exon_number']
                                first_transcript_exon = row['first_transcript_exon']
                                last_transcript_exon = row['last_transcript_exon']

                                # lists for codons and aas per overlapping exons per guide/editing window
                                codons = []
                                codons_edited = []
                                aas = []
                                aas_edited = []
                                aa_positions = []
                                includes_splice_site = False
                                synonymous = False

                                for edit in total_poss:

                                    if str(row['Strand']) == '+':
                                        index_edit_cds = edit - row['Start']
                                    elif str(row['Strand']) == '-':
                                        index_edit_cds = row['End'] - 1 - edit

                                    if row["Start"] <= edit < row["End"]: # edit is in cds

                                        offset = shared.get_offset(str(row['Strand']), int(row['Frame']), index_edit_cds)

                                        if (row["Start"] <= (edit - offset) < row["End"]) and (row["Start"] <= (edit + 3 - offset) < row["End"]): # whole codon in cds
                                            codon = str(ref_genome_pyfaidx[row["Chromosome"]][edit - offset:edit + 3 - offset])
                                            codon_edited = shared.replace_str_index(codon, offset, alt_edited)

                                            if str(row['Strand']) == '-':
                                                codon = shared.revcom(codon)
                                                codon_edited = shared.revcom(codon_edited)

                                            aa = shared.codon_sun_one_letter[codon]
                                            aa_edited = shared.codon_sun_one_letter[codon_edited]

                                            if aa == shared.codon_sun_one_letter["ATG"] and exon_number == first_transcript_exon: # startcodon
                                                if (row['Strand'] == '+' and (row["Start"] <= edit < (row["Start"] + 3))) or (row['Strand'] == '-' and ((row["End"] - 3) <= edit < row["End"])):
                                                    codon = "Start" + codon
                                                    aa = "Start" + aa

                                        else:
                                            codon = "incomplete_codon"
                                            codon_edited = "incomplete_codon"
                                            aa = "codon_incomplete"
                                            aa_edited = "codon_incomplete"

                                        if len(total_poss) == 1:
                                            if aa == aa_edited:
                                                synonymous = True
                                            else:
                                                synonymous = False

                                        if row['Strand'] == "+":
                                            nt_position = row['transcript_length_before'] + ((edit + 1) - row['Start'])
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                        elif row['Strand'] == "-":
                                            nt_position = row['transcript_length_before'] + (row['End'] - (edit)) # (row['End'] - (edit + 1))?
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                    elif (row["Start"] - 2) <= edit < (row["End"] + 2):
                                        includes_splice_site = True
                                        synonymous = False

                                        if (row["Start"] - 2) <= edit < row["Start"]:
                                            if row["Strand"] == '+':
                                                codon = "5prime_splice_site"
                                                codon_edited = "5prime_splice_site"
                                                aa = "5prime_splice_site"
                                                aa_edited = "5prime_splice_site"
                                                if exon_number == first_transcript_exon:
                                                    codon = "5prime_UTR"
                                                    codon_edited = "5prime_UTR"
                                                    aa = "5prime_UTR"
                                                    aa_edited = "5prime_UTR"
                                                    includes_splice_site = False
                                            elif row["Strand"] == '-':
                                                codon = "3prime_splice_site"
                                                codon_edited = "3prime_splice_site"
                                                aa = "3prime_splice_site"
                                                aa_edited = "3prime_splice_site"
                                                if exon_number == last_transcript_exon:
                                                    codon = "3prime_UTR"
                                                    codon_edited = "3prime_UTR"
                                                    aa = "3prime_UTR"
                                                    aa_edited = "3prime_UTR"
                                                    includes_splice_site = False
                                        elif row["End"] <= edit < (row["End"] + 2):
                                            if row["Strand"] == '-':
                                                codon = "5prime_splice_site"
                                                codon_edited = "5prime_splice_site"
                                                aa = "5prime_splice_site"
                                                aa_edited = "5prime_splice_site"
                                                if exon_number == first_transcript_exon:
                                                    codon = "5prime_UTR"
                                                    codon_edited = "5prime_UTR"
                                                    aa = "5prime_UTR"
                                                    aa_edited = "5prime_UTR"
                                                    includes_splice_site = False
                                            elif row["Strand"] == '+':
                                                codon = "3prime_splice_site"
                                                codon_edited = "3prime_splice_site"
                                                aa = "3prime_splice_site"
                                                aa_edited = "3prime_splice_site"
                                                if exon_number == last_transcript_exon:
                                                    codon = "3prime_UTR"
                                                    codon_edited = "3prime_UTR"
                                                    aa = "3prime_UTR"
                                                    aa_edited = "3prime_UTR"
                                                    includes_splice_site = False

                                        aa_position = aa

                                    else:
                                        codon = "not_in_CDS"
                                        codon_edited = "not_in_CDS"
                                        aa = "not_in_CDS"
                                        aa_edited = "not_in_CDS"

                                        aa_position = aa

                                    if edit != position:
                                        codon = codon.lower()
                                        codon_edited = codon_edited.lower()
                                        aa = aa.lower()
                                        aa_edited = aa_edited.lower()

                                    # lists for codons and aas per overlapping exons per guide/editing window
                                    codons.append(codon)
                                    codons_edited.append(codon_edited)
                                    aas.append(aa)
                                    aas_edited.append(aa_edited)

                                    aa_positions.append(aa_position)

                                # lists per guide for manual annotation
                                gene_symbols.append(gene_symbol) # mostly only one, but could be multiple, if CDSs overlapp
                                transcript_symbols.append(transcript_symbol) # multiple, categorized by according gene_symbol
                                exon_numbers.append(str(int(exon_number))) # multiple, categorized by according transcript_symbol; consider calling int (and maybe str) earlier
                                first_transcript_exons.append(str(int(first_transcript_exon))) # multiple, categorized by according transcript_symbol; consider calling int (and maybe str) earlier
                                last_transcript_exons.append(str(int(last_transcript_exon))) # multiple, categorized by according transcript_symbol; consider calling int (and maybe str) earlier
                                codonss.append(codons) # multiple, categorized by according frame according to transcript_symbol; consider different join character; tuple() not necessary
                                codonss_edited.append(codons_edited) # multiple, categorized by according frame according to transcript_symbol; consider different join character; tuple() not necessary
                                aass.append(aas) # multiple, categorized by according frame according to transcript_symbol; consider different join character; tuple() not necessary
                                aass_edited.append(aas_edited) # multiple, categorized by according frame according to transcript_symbol; consider different join character; tuple() not necessary
                                aa_positionss.append(aa_positions)
                                splice_sites_included.append(str(includes_splice_site)) # str() is necessary for collapsing later, find a better solution
                                synonymouss.append(str(synonymous)) # str() is necessary for collapsing later, find a better solution

                        else:
                            gene_symbols = ['no_CDS_found']
                            transcript_symbols = ['no_CDS_found']
                            exon_numbers = ['no_CDS_found']
                            first_transcript_exons = ['no_CDS_found']
                            last_transcript_exons = ['no_CDS_found']
                            codonss = [['no_CDS_found']] # [tuple(['no_CDS_found'])] (not necessary anymore with polars)
                            codonss_edited = [['no_CDS_found']] # [tuple(['no_CDS_found'])] (not necessary anymore with polars)
                            aass = [['no_CDS_found']] # [tuple(['no_CDS_found'])] (not necessary anymore with polars)
                            aass_edited = [['no_CDS_found']] # [tuple(['no_CDS_found'])] (not necessary anymore with polars)
                            aa_positionss = [['no_CDS_found']]
                            splice_sites_included = ['no_CDS_found'] # [str(False)]
                            synonymouss = ['no_CDS_found'] # [str(False)]

                        # deprecated version
                        # gene_symbolss.append('|'.join(gene_symbols))
                        # transcript_symbolss.append('|'.join(trasncript_symbols))
                        # exon_numberss.append('|'.join(exon_numbers))
                        # last_transcript_exonss.append('|'.join(last_transcript_exons))
                        # codonsss.append('|'.join(codonss))
                        # codonsss_edited.append('|'.join(codonss_edited))
                        # aasss.append('|'.join(aass))
                        # aasss_edited.append('|'.join(aass_edited))

                        # combine lists per guide for manual annotation in pl.DataFrame for better aggregation
                        annotations = pl.DataFrame({'gene_symbolss': gene_symbols,
                                                    'transcript_symbolss': transcript_symbols,
                                                    'exon_numberss': exon_numbers,
                                                    'first_transcript_exonss': first_transcript_exons,
                                                    'last_transcript_exonss': last_transcript_exons,
                                                    'codonsss': codonss,
                                                    'codonsss_edited': codonss_edited,
                                                    'aasss': aass,
                                                    'aasss_edited':aass_edited,
                                                    'aa_positionsss': aa_positionss,
                                                    'splice_sitess_included': splice_sites_included,
                                                    'synonymousss': synonymouss})

                        # group everything except exon_numberss, transcript_symbolss, first_transcript_exonss, last_transcript_exonss; should or could this be extended
                        annotations = annotations.group_by([col for col in annotations.columns if col not in ['aa_positionsss', 'exon_numberss', 'transcript_symbolss', 'first_transcript_exonss', 'last_transcript_exonss']], maintain_order=True).agg(pl.all())

                        aas_for_filtering = []
                        for aa_list in annotations['aasss'].to_list():
                            aas_for_filtering += aa_list

                        aas_edited_for_filtering = []
                        for aa_edited_list in annotations['aasss_edited'].to_list():
                            aas_edited_for_filtering += aa_edited_list

                        guide_is_missense = False
                        guide_is_synonymous = False
                        guide_is_stopgain = False
                        guide_is_startlost = False
                        guide_is_stoplost = False

                        for i, aa in enumerate(aas_for_filtering):
                            aa_editied = aas_edited_for_filtering[i]
                            if aa in non_stop_aas and aa_editied in non_stop_aas and aa != aa_editied:
                                guide_is_missense = True
                                # consequence.append('missense')
                            elif aa in non_stop_aas and aa_editied in non_stop_aas and aa == aa_editied:
                                guide_is_synonymous = True
                                # consequence.append('synonymous')

                            elif aa != 'Stop' and aa_editied == 'Stop':
                                guide_is_stopgain = True
                                # consequence.append('stopgain')
                            elif aa == 'Stop' and aa_editied != 'Stop':
                                guide_is_stoplost = True
                                # consequence.append('stoplost')
                            elif aa == 'Stop' and aa_editied == 'Stop': # can this happen?
                                guide_is_synonymous = True
                                # consequence.append('synonymous')

                            elif aa == 'StartM' and aa_editied !='M': # StartM not used in aas_edited yet; Any edit in StartM would be a startlost
                                guide_is_startlost = True
                                # consequence.append('startlost')
                            elif aa == 'StartM' and aa_editied =='M': # can't happen
                                guide_is_synonymous = True
                                # consequence.append('synonymous')

                            # else: # should not happen; only necessary, if using consequence list
                            #     consequence.append('complex')

                        consequence = []
                        if "True" in annotations['synonymousss'].to_list(): # list always only has 1 element
                            consequence.append('synonymous')
                        if "True" in annotations['splice_sitess_included'].to_list(): # list always only has 1 element
                            consequence.append('splice_site')
                        if specific:
                            consequence.append('specific')
                        if guide_is_missense:
                            consequence.append('missense')
                        if guide_is_stopgain:
                            consequence.append('stopgain')
                        if guide_is_stoplost:
                            consequence.append('stoplost')
                        if guide_is_startlost:
                            consequence.append('startlost')
                        if not consequence:
                            consequence.append('complex')

                        if not ((filter_synonymous) and ("True" not in annotations['synonymousss'].to_list()) or # list always only has 1 element
                                (filter_splice_site) and ("True" not in annotations['splice_sitess_included'].to_list()) or # list always only has 1 element
                                (filter_specific) and (not specific) or
                                (filter_missense) and (not guide_is_missense) or
                                (filter_nonsense) and (not guide_is_stopgain) or
                                (filter_stoplost) and (not guide_is_stoplost) or
                                (filter_startlost) and (not guide_is_startlost)):

                            if fiveprimepam:
                                possible_guide = shared.revcom(possible_guide)
                                possible_guide_with_pam = shared.revcom(possible_guide_with_pam)
                                possible_pam = shared.revcom(possible_pam)
                                edit_window = shared.revcom(edit_window)
                                edit_window_plus = shared.revcom(edit_window_plus)
                                safety_region = shared.revcom(safety_region)
                                edit_string = edit_string[::-1]
                                edit_pos_string = edit_pos_string[::-1]
                                distance_median_variant = distance_median_variant[::-1]
                                distance_median_all = distance_median_all[::-1]

                            # main fields
                            possible_guides_with_pam.append(possible_guide_with_pam)
                            possible_guides.append(possible_guide)
                            possible_pams.append(possible_pam)
                            possible_starts.append(str(target_seq_ref_start - guidelength + 1 if rev_com else target_seq_ref_start + 1)) # + 1 to switch to 1-based for SAM file
                            possible_ends.append(str(target_seq_ref_start + 1 if rev_com else target_seq_ref_start + guidelength + 1)) # + 1 to switch to 1-based for SAM file
                            possible_chroms.append(str(chrom))

                            # combine annotations, if multiple gene symbol or frames appeared to get annotations per guide
                            gene_symbolss.append(annotations['gene_symbolss'].to_list()) # consider different join character
                            transcript_symbolss.append(annotations['transcript_symbolss'].to_list()) # consider different join character
                            exon_numberss.append(annotations['exon_numberss'].to_list()) # consider different join character
                            first_transcript_exonss.append(annotations['first_transcript_exonss'].to_list()) # consider different join character
                            last_transcript_exonss.append(annotations['last_transcript_exonss'].to_list()) # consider different join character
                            codonsss.append(annotations['codonsss'].to_list()) # consider different join character
                            codonsss_edited.append(annotations['codonsss_edited'].to_list()) # consider different join character
                            aasss.append(annotations['aasss'].to_list()) # consider different join character
                            aasss_edited.append(annotations['aasss_edited'].to_list()) # consider different join character
                            aa_positionsss.append(annotations['aa_positionsss'].to_list())
                            splice_sitess_included.append(annotations['splice_sitess_included'].to_list())
                            synonymousss.append(annotations['synonymousss'].to_list())
                            consequences.append(consequence)

                            # lists per variant for shared.analyze_guide()
                            edit_windows.append(edit_window)
                            num_editss.append(str(num_edits)) # str() is necessary for collapsing later, find a better solution
                            specifics.append(str(specific)) # str() is necessary for collapsing later, find a better solution
                            edit_window_pluss.append(edit_window_plus)
                            num_edits_pluss.append(str(num_edits_plus)) # str() is necessary for collapsing later, find a better solution
                            specific_pluss.append(str(specific_plus)) # str() is necessary for collapsing later, find a better solution
                            safety_regions.append(safety_region)
                            num_edits_safetys.append(str(num_edits_safety)) # str() is necessary for collapsing later, find a better solution
                            additional_in_safetys.append(str(additional_in_safety)) # str() is necessary for collapsing later, find a better solution
                            edit_strings.append(edit_string)
                            edit_pos_strings.append(edit_pos_string)
                            specificitys.append(str(specificity)) # str() is necessary for collapsing later, find a better solution
                            distance_median_variants.append(distance_median_variant)
                            # quality_scores_variants.append(quality_scores_variant)
                            distance_median_alls.append(distance_median_all)
                            # quality_scores_alls.append(quality_scores_all)

                    # decreasing variant position for next loop; is there a more direct way? This seems unintuitive
                    variant_position -= 1
                    if rev_com:
                        target_seq_ref_start -= 1
                    else:
                        target_seq_ref_start += 1

                if possible_guides == []:
                    gene_symbolss.append(["no_guides_found"])
                    possible_guides.append("no_guides_found")
                    possible_guides_with_pam.append("no_guides_found")
                    edit_windows.append("no_guides_found")
                    num_editss.append("no_guides_found")
                    specifics.append("no_guides_found")
                    edit_window_pluss.append("no_guides_found")
                    num_edits_pluss.append("no_guides_found")
                    specific_pluss.append("no_guides_found")
                    safety_regions.append("no_guides_found")
                    num_edits_safetys.append("no_guides_found")
                    additional_in_safetys.append("no_guides_found")
                    synonymousss.append(["no_guides_found"])
                    consequences.append(["no_guides_found"])
                    codonsss.append([["no_guides_found"]])
                    codonsss_edited.append([["no_guides_found"]])
                    aasss.append([["no_guides_found"]])
                    aasss_edited.append([["no_guides_found"]])
                    aa_positionsss.append([[["no_guides_found"]]])
                    splice_sitess_included.append(["no_guides_found"])
                    edit_strings.append("no_guides_found")
                    edit_pos_strings.append("no_guides_found")
                    specificitys.append("no_guides_found")
                    distance_median_variants.append("no_guides_found")
                    # quality_scores_variants.append("no_guides_found")
                    distance_median_alls.append("no_guides_found")
                    # quality_scores_alls.append("no_guides_found")
                    transcript_symbolss.append([["no_guides_found"]])
                    exon_numberss.append([["no_guides_found"]])
                    first_transcript_exonss.append([["no_guides_found"]])
                    last_transcript_exonss.append([["no_guides_found"]])
                    possible_starts.append("no_guides_found")
                    possible_ends.append("no_guides_found")
                    possible_chroms.append("no_guides_found")

                # final lists for editbale variants
                all_variant.append(variant)
                all_variant_real.append(variant_real)
                all_editable.append(editable)
                all_be_strings.append(be_string)
                all_original_alt.append(original_alt)
                all_target_seq_ref.append(target_seq_ref)
                all_target_seq_ref_match.append(target_seq_ref_match)
                all_target_seq.append(target_seq)
                all_target_base_ref.append(target_base_ref)
                all_target_base.append(target_base)
                all_possible_guides_with_pam.append(possible_guides_with_pam) # list
                all_edit_strings.append(edit_strings) # list
                all_edit_pos_strings.append(edit_pos_strings) # list
                all_possible_guides.append(possible_guides) # list
                all_possible_pams.append(possible_pams)
                all_rev_com.append("{}".format("-" if rev_com and not fiveprimepam else
                                            "+" if not rev_com and not fiveprimepam else
                                            '+' if rev_com and fiveprimepam else
                                            '-' if not rev_com and fiveprimepam else
                                            ''))
                all_edit_window.append(edit_windows) # list
                all_num_edits.append(num_editss)
                all_specific.append(specifics) # list
                all_edit_window_plus.append(edit_window_pluss) # list
                all_num_edits_plus.append(num_edits_pluss)
                all_specific_plus.append(specific_pluss) # list
                all_safety_region.append(safety_regions)
                all_num_edits_safety.append(num_edits_safetys)
                all_additional_in_safety.append(additional_in_safetys)
                all_specificity.append(specificitys) # list
                all_distance_median_variant.append(distance_median_variants) # list
                # all_quality_scores_variant.append(quality_scores_variants) # list
                all_distance_median_all.append(distance_median_alls) # list
                # all_quality_scores_all.append(quality_scores_alls) # list
                all_gene_symbols.append(gene_symbolss) # list
                all_transcript_symbols.append(transcript_symbolss) # list
                all_exon_numbers.append(exon_numberss) # list
                all_first_transcript_exons.append(first_transcript_exonss) # list
                all_last_transcript_exons.append(last_transcript_exonss) # list
                all_codonss.append(codonsss) # list
                all_codonss_edited.append(codonsss_edited) # list
                all_aass.append(aasss) # list
                all_aass_edited.append(aasss_edited) # list
                all_aa_positionss.append(aa_positionsss) # list
                all_splice_sites_included.append(splice_sitess_included)
                all_synonymouss.append(synonymousss)
                all_consequences.append(consequences)
                all_guide_starts.append(possible_starts)
                all_guide_ends.append(possible_ends)
                all_chroms.append(possible_chroms)

            else:
                # final lists for non-editbale variants
                all_variant.append(variant)
                all_variant_real.append(variant_real) # not in output
                all_editable.append(editable) # not in output
                all_be_strings.append(be)
                all_original_alt.append("be_not_usable")
                all_target_seq_ref.append("") # not in output
                all_target_seq_ref_match.append(target_seq_ref_match)
                all_target_seq.append("") # not in output
                all_target_base_ref.append(target_base_ref) # not in output
                all_target_base.append("") # not in output
                all_possible_guides_with_pam.append(["be_not_usable"])
                all_edit_strings.append(["be_not_usable"])
                all_edit_pos_strings.append(["be_not_usable"])
                all_possible_guides.append(["be_not_usable"])
                all_possible_pams.append("") # not in output
                all_rev_com.append("{}".format("be_not_usable"))
                all_edit_window.append(["be_not_usable"])
                all_num_edits.append(["be_not_usable"])
                all_specific.append(["be_not_usable"])
                all_edit_window_plus.append(["be_not_usable"])
                all_num_edits_plus.append(["be_not_usable"])
                all_specific_plus.append(["be_not_usable"])
                all_safety_region.append(["be_not_usable"])
                all_num_edits_safety.append(["be_not_usable"])
                all_additional_in_safety.append(["be_not_usable"])
                all_specificity.append(["be_not_usable"])
                all_distance_median_variant.append(["be_not_usable"])
                # all_quality_scores_variant.append(["be_not_usable"])
                all_distance_median_all.append(["be_not_usable"])
                # all_quality_scores_all.append(["be_not_usable"])
                all_gene_symbols.append([["be_not_usable"]])
                all_transcript_symbols.append([[["be_not_usable"]]])
                all_exon_numbers.append([[["be_not_usable"]]])
                all_first_transcript_exons.append([[["be_not_usable"]]])
                all_last_transcript_exons.append([[["be_not_usable"]]])
                all_codonss.append([[["be_not_usable"]]])
                all_codonss_edited.append([[["be_not_usable"]]])
                all_aass.append([[["be_not_usable"]]])
                all_aass_edited.append([[["be_not_usable"]]])
                all_aa_positionss.append([[[["be_not_usable"]]]])
                all_splice_sites_included.append([["be_not_usable"]])
                all_synonymouss.append([["be_not_usable"]])
                all_consequences.append([["be_not_usable"]])
                all_guide_starts.append(["be_not_usable"])
                all_guide_ends.append(["be_not_usable"])
                all_chroms.append(["be_not_usable"])

    # transform lists to pl.DataFrame for better handling
    sgrnas = pl.DataFrame({"variant": all_variant,
                           "base_editor": all_be_strings,
                           "symbol": all_gene_symbols,
                           "guide": all_possible_guides,
                           "guide_chrom": all_chroms,
                           "guide_start": all_guide_starts,
                           "guide_end": all_guide_ends,
                           "guide_with_pam": all_possible_guides_with_pam,
                           "edit_window": all_edit_window,
                           "num_edits": all_num_edits,
                           "specific": all_specific,
                           "edit_window_plus": all_edit_window_plus,
                           "num_edits_plus": all_num_edits_plus,
                           "specific_plus": all_specific_plus,
                           "safety_region": all_safety_region,
                           "num_edits_safety": all_num_edits_safety,
                           "additional_in_safety": all_additional_in_safety,
                           "ne_plus": "NA_for_variants",
                           "synonymous": all_synonymouss,
                           "consequence": all_consequences,
                           "strand": all_rev_com,
                           "codon_ref": all_codonss,
                           "aa_ref": all_aass,
                           "aa_pos": all_aa_positionss,
                           "codon_edit": all_codonss_edited,
                           "aa_edit": all_aass_edited,
                           "splice_site_included": all_splice_sites_included,
                           "originally_intended_ALT": all_original_alt,
                           "ref_match": all_target_seq_ref_match,
                           "off_target_bases": all_edit_strings,
                           "edited_positions": all_edit_pos_strings,
                           "specificity": all_specificity,
                           "distance_median_variant": all_distance_median_variant,
                        #    "efficiency_scores_variant": all_quality_scores_variant,
                           "distance_median_all": all_distance_median_all,
                        #    "efficiency_scores_all": all_quality_scores_all,
                           "transcript": all_transcript_symbols,
                           "exon_number": all_exon_numbers,
                           "first_transcript_exon": all_first_transcript_exons,
                           "last_transcript_exon": all_last_transcript_exons},
                           strict=False) # get rid of this

    if blast: # needs to be adjusted to polars

        blast_guides.check_blastdb(refgenome, False)

        sgrnas = sgrnas.with_row_index('index')

        guides = sgrnas.select('index', 'guide').explode('guide') # should never be deeper then list in list before explode()

        blast_results = blast_guides.guide_blast(guides,
                                                 guidelength,
                                                 refgenome,
                                                 'variants',
                                                 no_contigs)
        if not blast_results.is_empty():
            sgrnas = sgrnas.join(blast_results, left_on='index', right_on='indexvar', how='left')
            sgrnas = sgrnas.with_columns(blastcount = pl.when((pl.col('guide') == ["no_guides_found"]) |
                                                              (pl.col('guide') == ["be_not_usable"]))
                                                        .then(pl.col('guide'))
                                                        .otherwise(pl.col('blastcount')))
        else:
            sgrnas = sgrnas.with_columns(blastcount = pl.col('guide'))
        sgrnas = sgrnas.drop('index')

    # add vep annotations, if wanted
    if vep:

        variants_vep = pl.DataFrame({'variant': all_variant_real}).with_row_index('for_sorting_later')
        variants_non_vep = variants_vep.filter(pl.col('variant').is_in(['variant_is_improperly_formatted',
                                                                        'no_input_gene_given',
                                                                        'no_input_transcript_given',
                                                                        'no_input_mutation_given',
                                                                        'input_transcript_not_found',
                                                                        'input_gene_not_found',
                                                                        'reference_not_amino_acid',
                                                                        'mutation_not_amino_acid',
                                                                        'input_position_not_numeric',
                                                                        'wrong_reference_amino_acid',
                                                                        'non_existent_input_rsID',
                                                                        'genomic_position_not_numeric',
                                                                        'genomic_coordinates_not_found']))
        variants_vep = variants_vep.filter(~pl.col('variant').is_in(['variant_is_improperly_formatted',
                                                                     'no_input_gene_given',
                                                                     'no_input_transcript_given',
                                                                     'no_input_mutation_given',
                                                                     'input_transcript_not_found',
                                                                     'input_gene_not_found',
                                                                     'reference_not_amino_acid',
                                                                     'mutation_not_amino_acid',
                                                                     'input_position_not_numeric',
                                                                     'wrong_reference_amino_acid',
                                                                     'non_existent_input_rsID',
                                                                     'genomic_position_not_numeric',
                                                                     'genomic_coordinates_not_found']))
        variants_vep_sorted = shared.sort_variantsdf(variants_vep) # new vep algorithm allows optimization here
        variants_vep_sorted_input = variants_vep_sorted['variant'].unique(maintain_order=True).to_list()

        vep_annotations = get_vep.get_vep_annotation(variants_vep_sorted_input,
                                                               species=vep_species,
                                                               assembly=vep_assembly,
                                                               dir_cache=vep_dir_cache,
                                                            #    dir_plugins=vep_dir_plugins, # currently not in use
                                                               cache_version=vep_cache_version,
                                                               flags=vep_flags)

        vep_annotations = vep_annotations.rename(lambda column_name: "VEP_" + column_name)
        vep_annotations = vep_annotations.group_by('VEP_#Uploaded_variation').agg(pl.all().str.join(","))
        vep_annotations = vep_annotations.with_columns(
            pl.col('VEP_#Uploaded_variation').str.replace('./', '', literal=True).str.replace('/', '_', literal=True) # this needs to be tested
            )
        variants_vep_sorted = variants_vep_sorted.join(vep_annotations, left_on='variant', right_on='VEP_#Uploaded_variation', how='left')

        variants_vep_resorted = shared.resort_variantsdf(variants_vep_sorted) # new vep algorithm allows optimization here

        for col in [col for col in vep_annotations.columns if col in variants_vep_resorted.columns]: # is there a better way?
            variants_non_vep = variants_non_vep.with_columns(pl.lit("not_suitable_for_VEP").alias(col))

        variants_vep_resorted = pl.concat([variants_vep_resorted, variants_non_vep]).sort('for_sorting_later').drop('for_sorting_later')

        variants_vep_resorted = variants_vep_resorted.select(pl.exclude(["original_index", "variant"]))
        sgrnas = sgrnas.with_columns(variants_vep_resorted)

    symbols_to_contract = ['symbol']

    splice_sites_to_modify = ['splice_site_included',
                              'synonymous']

    consequences_to_modify = ['consequence']

    aa_pos_to_modify_very_first = ['aa_pos']

    transcript_cols_to_modify_first = ['exon_number',
                                       'transcript',
                                       'first_transcript_exon',
                                       'last_transcript_exon']

    variant_cols_to_modify_first = ['codon_ref',
                                    'aa_ref',
                                    'codon_edit',
                                    'aa_edit']

    non_list_columns = ['variant',
                        'symbol',
                        'base_editor',
                        'ne_plus',
                        'strand',
                        'originally_intended_ALT',
                        'ref_match']

    if vep:
        non_list_columns += variants_vep_resorted.columns

    columns_to_modify_last = [col for col in sgrnas.columns if col not in non_list_columns]

    sgrnas = sgrnas.with_columns(
        pl.col(aa_pos_to_modify_very_first).list.eval(pl.element().list.eval(pl.element().list.eval(pl.element().list.join(";")).list.join("~")).list.join("^")),
        pl.col(symbols_to_contract).list.eval(pl.element().list.unique(maintain_order=True).list.join("?")), # make nested symbol list unique and flatten; join with '?', if multiple symbols appear
        pl.col(transcript_cols_to_modify_first).list.eval(pl.element().list.eval(pl.element().list.join("~")).list.join("^")),
        pl.col(variant_cols_to_modify_first).list.eval(pl.element().list.eval(pl.element().list.join(";")).list.join("^")), # this separator should be something different (the first one; old comment?)
        pl.col(splice_sites_to_modify).list.eval(pl.element().list.join("^")), # splice_site_included on same level as transcript groups
        pl.col(consequences_to_modify).list.eval(pl.element().list.join("&")) # consequences on same level as transcript groups
    )

    if allpossible:
        all_variant = sgrnas['variant'].to_list()
        sgrnas = sgrnas.with_columns(pl.Series('variant', [f'{all_variant[i]}({all_variant_real[i]})' if all_variant[i].split(':')[-1] != all_variant_real[i] else all_variant[i] for i in range(len(all_variant))]))

    if aspect == 'exploded': # one line per guide
        sgrnas = sgrnas.explode(symbols_to_contract + columns_to_modify_last) # explode per guide (and the associated columns)

    elif aspect == 'collapsed': # one line per variant
        sgrnas = sgrnas.with_columns(
            pl.col(symbols_to_contract).list.unique(maintain_order=True).list.join("?"), # make nested symbol list unique and flatten; join with '?', if multiple symbols appear
            pl.col(columns_to_modify_last).list.join("•") # join guides (and the associated columns)
        )

    sgrnas = sgrnas.unique().sort(by=sgrnas.columns) # duplicates are unlikely due to all the annotations; consider revision

    sam_df = sgrnas.filter((pl.col("guide").cast(str) != "no_guides_found") & # editbale, but no guides found
                           (pl.col("guide").cast(str) != "be_not_usable")) # not editable

    sam_df = sam_df.select(
        ['variant', 'guide', 'guide_chrom', 'guide_start', 'strand']
        ).with_columns(
            pl.col("guide").str.split("•").alias("guide"),
            pl.col("guide_chrom").str.split("•").alias("guide_chrom"),
            pl.col("guide_start").str.split("•").alias("guide_start")
            ).explode(['guide', 'guide_chrom', 'guide_start']) # force explosion

    sam_df = sam_df.with_columns(
        pl.col('variant').alias('QNAME'),
        pl.when(pl.col("strand") == '+').then(0).otherwise(
            pl.when(pl.col("strand") == '-').then(16)
            ).alias('FLAG'),
        pl.col('guide_chrom').alias('RNAME'),
        pl.col('guide_start').cast(pl.Int64).alias('POS'),
        pl.lit(255).alias('MAPQ'),
        pl.lit(str(guidelength) + 'M').alias('CIGAR'),
        pl.lit('*').alias('RNEXT'),
        pl.lit(0).alias('PNEXT'),
        pl.lit(0).alias('TLEN'),
        pl.when(pl.col("strand") == '+').then(pl.col('guide')).otherwise(
            pl.when(pl.col("strand") == '-').then(pl.col('guide').map_elements(shared.revcom, return_dtype=pl.String))
            ).alias('SEQ'),
        pl.lit('*').alias('QUAL')
    ).select(['QNAME',
              'FLAG',
              'RNAME',
              'POS',
              'MAPQ',
              'CIGAR',
              'RNEXT',
              'PNEXT',
              'TLEN',
              'SEQ',
              'QUAL']).sort(by=[
                  'RNAME',
                  'POS',
                  'QNAME',
                  'FLAG',
                  'MAPQ',
                  'CIGAR',
                  'RNEXT',
                  'PNEXT',
                  'TLEN',
                  'SEQ',
                  'QUAL'
                  ])

    sgrnas = sgrnas.drop(['ne_plus',
                          'off_target_bases',
                          'specificity'])

    if not allpossible:
        sgrnas = sgrnas.drop(['originally_intended_ALT'])

    if edit_window_start_plus == 0 and edit_window_end_plus == 0:
        sgrnas = sgrnas.drop(['edit_window_plus',
                              'num_edits_plus',
                              'specific_plus',
                              'safety_region',
                              'num_edits_safety',
                              'additional_in_safety'])

    return (sgrnas, sam_df)


def output_sgrnas(sgrnas, output_file):

    sgrnas.write_csv(output_file + ".tsv", separator='\t')

    # filter for variant, where guides were found; is thiss still the best way?
    sgrnas_filtered = sgrnas.filter((pl.col("guide").cast(str) != "no_guides_found") & # editbale, but no guides found
                                    (pl.col("guide").cast(str) != "be_not_usable")) # not editable

    sgrnas_filtered.write_csv(output_file + "_filtered.tsv", separator='\t')


def output_guides_sam(sam_df, output_file, refgenome):

    sam_df.write_csv(output_file + "_filtered.sam", separator='\t', include_header=False)
    pysam.view(output_file + "_filtered.sam",
               '-b',
               '-o', output_file + "_filtered_unsorted.bam",
               '-t', refgenome + '.fai',
               catch_stdout=False)
    os.remove(output_file + "_filtered.sam")
    pysam.sort('-o', output_file + "_filtered.bam", output_file + "_filtered_unsorted.bam") # should already be sorted, but this acts as a failsafe
    os.remove(output_file + "_filtered_unsorted.bam")
    pysam.index(output_file + "_filtered.bam")


if __name__ == "__main__":
    pass
