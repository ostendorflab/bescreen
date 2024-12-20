import os
import re
import polars as pl
import pyfaidx
import pysam
import get_vep
import shared
import blast_guides


def saturate_region(ref_genome,
                    regions,
                    input_file,
                    pamsite,
                    edit_window_start,
                    edit_window_end,
                    guidelength,
                    baseeditor,
                    edit_window_start_plus,
                    edit_window_end_plus,
                    vep,
                    aspect,
                    vep_flags,
                    vep_species,
                    vep_assembly,
                    vep_dir_cache,
                    vep_dir_plugins,
                    vep_cache_version,
                    blast,
                    no_contigs):

    bes = {
        'ABE': {'fwd': {'REF': 'A', 'ALT': 'G'},
                'rev': {'REF': 'T', 'ALT': 'C'}},
        'CBE': {'fwd': {'REF': 'C', 'ALT': 'T'},
                'rev': {'REF': 'G', 'ALT': 'A'}}
    }

    # modified input
    if baseeditor == "ABE":
        del bes['CBE']
    elif baseeditor == "CBE":
        del bes['ABE']

    if input_file:

        input_df = pl.read_csv(input_file, separator=',')
        input_columns = input_df.columns

        if 'region' in input_columns:
            regions_list = pl.read_csv(input_file, separator=',')['region'].to_list()
        else:
            raise ValueError('Input file does not contain the column "region"!\nPlease provide it!')

    elif regions:
        regions_list = regions.replace(' ', '').split(",")

    ref_genome_pyfaidx = pyfaidx.Fasta(ref_genome)

    pamsite_relevant = pamsite.replace('N', '')

    # ---[***}************NG <- guide
    # >TAA----------------NG <- sequence
    # CN************{***]--- <- shared.revcom(guide)
    # CN----------------TAC< <- shared.revcom(sequence)
    flanklength_pamsided = guidelength + len(pamsite) - edit_window_start

    # ***{***]------------NG <- guide
    # -------ATG>>>>>>>>>>NG <- sequence
    # CN------------[***}*** <- shared.revcom(guide)
    # CN>>>>>>>>>>ATT------- <- shared.revcom(sequence)
    flanklength_nonpamsided = edit_window_end - 1

    # real code
    # create empty lists to fill for editing guides
    be_output = []
    variant_output = []
    guide_output = []
    guide_with_pam_output = []
    edit_window_output = []
    guide_output_strand = []
    guide_with_pam_output_strand = []
    edit_window_output_strand = []
    direction_output = []
    specific_output = []
    edit_window_plus_output = []
    edit_window_plus_output_strand = []
    specific_plus_output = []
    num_edits_output = []
    num_edits_plus_output = []
    safety_region_output = []
    num_edits_safety_output = []
    additional_in_safety_output = []
    edit_string_output = []
    edit_pos_string_output = []
    distance_median_all_output = []
    quality_scores_all_output = []
    guide_starts_output = []
    guide_ends_output = []
    chroms_output = []

    # create empty lists to fill for non-editing guides
    be_output_ne = []
    guide_output_ne = []
    guide_with_pam_output_ne = []
    edit_window_output_ne = []
    guide_output_strand_ne = []
    guide_with_pam_output_strand_ne = []
    edit_window_output_strand_ne = []
    edit_window_plus_output_ne = []
    edit_window_plus_output_strand_ne = []
    ne_plus_output_ne = []
    direction_output_ne = []
    safety_region_output_ne = []
    guide_starts_output_ne = []
    guide_ends_output_ne = []
    chroms_output_ne = []

    # precalculation of quality scores for positions in editing window
    distance_median_dict, quality_scores_dict = shared.qc_precalc(edit_window_start, edit_window_end)

    sequences_not_found = []

    # iterate through regions
    for region in regions_list:

        try:
            chrom, positions = region.split(':')
            start, end = positions.split('-')
            start = int(start)
            end = int(end)

            # get region sequence
            sequence = ref_genome_pyfaidx[chrom][start:end] # pyranges transforms gtf file into 0-based

            #generate forward search string with:
            # 5' overhang to account for first base is last base of editing window
            flanking_upstream_forward = ref_genome_pyfaidx[chrom][start - flanklength_nonpamsided:start]
            # 3' overhang to account for last base is first base of editing window plus PAM site
            flanking_downstream_forward = ref_genome_pyfaidx[chrom][end:end + flanklength_pamsided]
            # concat substring
            forward_search_string = str(flanking_upstream_forward) + str(sequence) + str(flanking_downstream_forward)

            #generate reverse search string with:
            # 5' overhang to account for first base is first base of editing window plus PAM site
            flanking_upstream_reverse = ref_genome_pyfaidx[chrom][start - flanklength_pamsided:start]
            # 3' overhang to account for last base is last base of editing window
            flanking_downstream_reverse = ref_genome_pyfaidx[chrom][end:end + flanklength_nonpamsided]
            # concat substring
            reverse_search_string = str(flanking_upstream_reverse) + str(sequence) + str(flanking_downstream_reverse)

            sequence_found = True
        except:
            sequences_not_found.append(region)
            sequence_found = False

        if sequence_found:
            for be in bes: # for ABE and CBE
                startpos_fwd = start - flanklength_nonpamsided
                startpos_rev = start - flanklength_pamsided

                for i in range(len(forward_search_string) - (guidelength + len(pamsite)) + 1): # start foward search; + 1 correct?

                    possible_guide_with_pam = forward_search_string[i:i + guidelength + len(pamsite)] # sliding window for guide with PAM
                    possible_guide = forward_search_string[i:i + guidelength] # sliding window for guide without PAM
                    # edit_window_bases = possible_guide_with_pam[edit_window_start - 1:edit_window_end] # sliding window for edit window
                    edit_window_bases = possible_guide[edit_window_start - 1:edit_window_end] # sliding window for edit window
                    edit_window_plus_bases = possible_guide[edit_window_start - 1 - edit_window_start_plus:edit_window_end + edit_window_end_plus] # sliding window for edit window plus
                    safety_region = possible_guide[edit_window_start - 1 - edit_window_start_plus:edit_window_start - 1] + edit_window_bases.lower() + possible_guide[edit_window_end:edit_window_end + edit_window_end_plus]

                    if possible_guide_with_pam.endswith(pamsite_relevant): # PAM found at 5' end

                        if bes[be]['fwd']['REF'] in edit_window_bases: # editable base in edit window

                            # find all indices for possible edits in the edit window in search string
                            indices_edit = [i + (edit_window_start - 1) + m.start() for m in re.finditer(bes[be]['fwd']['REF'], edit_window_bases)]
                            indices_plus_edit = [i + (edit_window_start - 1 - edit_window_start_plus) + m.start() for m in re.finditer(bes[be]['fwd']['REF'], edit_window_plus_bases)] # is this correct?

                            variants = [] # reset variant list for sliding window

                            poss = [start - flanklength_nonpamsided + index_edit for index_edit in indices_edit]
                            if any(start <= pos < end for pos in poss):

                                for index_edit in indices_edit:

                                    pos = start - flanklength_nonpamsided + index_edit # pos for fasta

                                    if bes[be]['fwd']['REF'] == ref_genome_pyfaidx[chrom][pos]: # failsafe; can be removed

                                        variants.append(f'{chrom}_{pos + 1}_{ref_genome_pyfaidx[chrom][pos]}_{bes[be]["fwd"]["ALT"]}') # append variant to variantlist 1-based as VCF; add here already to not miss off-target edits outside CDS

                            if variants != []: # is this redundant?

                                edit_window_bases_unused, \
                                num_edits, \
                                specific, \
                                edit_window_plus_bases_unused, \
                                num_edits_plus, \
                                specific_plus, \
                                safety_region_unused, \
                                num_edits_safety, \
                                additional_in_safety, \
                                edit_string, \
                                edit_pos_string, \
                                specificity_unused, \
                                distance_median_variant_unused, \
                                quality_scores_variant_unused, \
                                distance_median_all, \
                                quality_scores_all = shared.analyze_guide(possible_guide,
                                                                        edit_window_start,
                                                                        edit_window_end,
                                                                        edit_window_start_plus,
                                                                        edit_window_end_plus,
                                                                        bes[be]['fwd']['REF'],
                                                                        None,
                                                                        distance_median_dict,
                                                                        quality_scores_dict) # _unused variables just for testing purposes; can be removed later

                                # if  (edit_window_bases_unused != edit_window_bases or
                                #     edit_window_plus_bases_unused != edit_window_plus_bases): # just for testing purposes; can be removed later
                                #     raise Exception("output of analyze_guide() wrong in forward search")

                                be_output.append(be)
                                variant_output.append(variants) # list
                                guide_output.append(possible_guide)
                                guide_with_pam_output.append(possible_guide_with_pam)
                                edit_window_output.append(edit_window_bases)
                                guide_output_strand.append(possible_guide)
                                guide_with_pam_output_strand.append(possible_guide_with_pam)
                                edit_window_output_strand.append(edit_window_bases)
                                specific_output.append(specific)
                                edit_window_plus_output.append(edit_window_plus_bases)
                                edit_window_plus_output_strand.append(edit_window_plus_bases)
                                specific_plus_output.append(specific_plus)
                                direction_output.append('+') # forward
                                num_edits_output.append(num_edits)
                                num_edits_plus_output.append(num_edits_plus)
                                safety_region_output.append(safety_region)
                                num_edits_safety_output.append(num_edits_safety)
                                additional_in_safety_output.append(additional_in_safety)
                                edit_string_output.append(edit_string)
                                edit_pos_string_output.append(edit_pos_string)
                                distance_median_all_output.append(distance_median_all)
                                quality_scores_all_output.append(quality_scores_all)
                                guide_starts_output.append(str(startpos_fwd + 1))
                                guide_ends_output.append(str(startpos_fwd + guidelength + 1))
                                chroms_output.append(str(chrom))

                        else:
                            if bes[be]['fwd']['REF'] in edit_window_plus_bases:
                                ne_plus = False
                            else:
                                ne_plus = True
                            be_output_ne.append(be)
                            guide_output_ne.append(possible_guide)
                            guide_with_pam_output_ne.append(possible_guide_with_pam)
                            edit_window_output_ne.append(edit_window_bases)
                            guide_output_strand_ne.append(possible_guide)
                            guide_with_pam_output_strand_ne.append(possible_guide_with_pam)
                            edit_window_output_strand_ne.append(edit_window_bases)
                            edit_window_plus_output_ne.append(edit_window_plus_bases)
                            edit_window_plus_output_strand_ne.append(edit_window_plus_bases)
                            ne_plus_output_ne.append(ne_plus)
                            direction_output_ne.append('+') # forward
                            safety_region_output_ne.append(safety_region)
                            guide_starts_output_ne.append(startpos_fwd + 1)
                            guide_ends_output_ne.append(startpos_fwd + guidelength + 1)
                            chroms_output_ne.append(chrom)

                    startpos_fwd += 1

                for i in range(len(reverse_search_string) - (guidelength + len(pamsite)) + 1): # start reverse search; + 1 correct?

                    possible_guide_with_pam = reverse_search_string[i:i + guidelength + len(pamsite)] # sliding window for guide with PAM
                    possible_guide = reverse_search_string[i + len(pamsite):i + guidelength + len(pamsite)] # sliding window for guide without PAM
                    # edit_window_bases = possible_guide_with_pam[guidelength + len(pamsite) - edit_window_end:guidelength + len(pamsite) - edit_window_start + 1] # sliding window for edit window
                    edit_window_bases = possible_guide[guidelength - edit_window_end:guidelength - edit_window_start + 1] # sliding window for edit window
                    edit_window_plus_bases = possible_guide[guidelength - edit_window_end - edit_window_end_plus:guidelength - edit_window_start + 1 + edit_window_start_plus] # sliding window for edit window plus
                    safety_region = possible_guide[guidelength - edit_window_end - edit_window_end_plus:guidelength - edit_window_end] + edit_window_bases.lower() + possible_guide[guidelength - edit_window_start + 1:guidelength - edit_window_start + 1 + edit_window_start_plus]

                    if possible_guide_with_pam.startswith(shared.revcom(pamsite_relevant)): # shared.revcom(PAM) found at 3' end

                        if bes[be]['rev']['REF'] in edit_window_bases: # shared.revcom(editable base) in edit window

                            # find all indices for possible edits in the edit window in search string
                            indices_edit = [i + (guidelength + len(pamsite) - edit_window_end) + m.start() for m in re.finditer(bes[be]['rev']['REF'], edit_window_bases)]
                            indices_plus_edit = [i + (guidelength + len(pamsite) - edit_window_end + edit_window_end_plus) + m.start() for m in re.finditer(bes[be]['rev']['REF'], edit_window_plus_bases)] # is this correct?

                            variants = [] # reset variant list for sliding window
                            codons = [] # reset variant list for sliding window
                            codons_edited = [] # reset variant list for sliding window
                            aas = [] # reset variant list for sliding window
                            aas_edited = [] # reset variant list for sliding window
                            # includes_splice_site = "not_tested" # also check for splice-sites even, if not explicitely added
                            # if splice_sites:
                            #     includes_splice_site = False
                            includes_splice_site = False

                            poss = [start - flanklength_pamsided + index_edit for index_edit in indices_edit]
                            if any(start <= pos < end for pos in poss):

                                for index_edit in indices_edit:

                                    pos = start - flanklength_pamsided + index_edit # pos for fasta

                                    if bes[be]['rev']['REF'] == ref_genome_pyfaidx[chrom][pos]: # failsafe; can be removed

                                        variants.append(f'{chrom}_{pos + 1}_{ref_genome_pyfaidx[chrom][pos]}_{bes[be]["rev"]["ALT"]}') # append variant to variantlist 1-based as VCF; add here already to not miss off-target edits outside CDS

                            if variants != []: # is this redundant?

                                edit_window_bases_unused, \
                                num_edits, \
                                specific, \
                                edit_window_plus_bases_unused, \
                                num_edits_plus, \
                                specific_plus, \
                                safety_region_unused, \
                                num_edits_safety, \
                                additional_in_safety, \
                                edit_string, \
                                edit_pos_string, \
                                specificity_unused, \
                                distance_median_variant_unused, \
                                quality_scores_variant_unused, \
                                distance_median_all, \
                                quality_scores_all = shared.analyze_guide(shared.revcom(possible_guide),
                                                                        edit_window_start,
                                                                        edit_window_end,
                                                                        edit_window_start_plus,
                                                                        edit_window_end_plus,
                                                                        bes[be]['fwd']['REF'],
                                                                        None,
                                                                        distance_median_dict,
                                                                        quality_scores_dict) # _unused variables just for testing purposes; can be removed later

                                # edit_window_bases_unused = shared.revcom(edit_window_bases_unused)
                                # edit_window_plus_bases_unused = shared.revcom(edit_window_plus_bases_unused)

                                # if  (edit_window_bases_unused != edit_window_bases or
                                #     edit_window_plus_bases_unused != edit_window_plus_bases): # just for testing purposes; can be removed later
                                #     raise Exception("output of analyze_guide() wrong in reverse search")

                                be_output.append(be)
                                variant_output.append(variants) # list
                                guide_output.append(shared.revcom(possible_guide))
                                guide_with_pam_output.append(shared.revcom(possible_guide_with_pam))
                                edit_window_output.append(shared.revcom(edit_window_bases))
                                guide_output_strand.append(possible_guide)
                                guide_with_pam_output_strand.append(possible_guide_with_pam)
                                edit_window_output_strand.append(edit_window_bases)
                                specific_output.append(specific)
                                edit_window_plus_output.append(shared.revcom(edit_window_plus_bases))
                                edit_window_plus_output_strand.append(edit_window_plus_bases)
                                specific_plus_output.append(specific_plus)
                                direction_output.append('-') # reverse
                                num_edits_output.append(num_edits)
                                num_edits_plus_output.append(num_edits_plus)
                                safety_region_output.append(shared.revcom(safety_region))
                                num_edits_safety_output.append(num_edits_safety)
                                additional_in_safety_output.append(additional_in_safety)
                                edit_string_output.append(edit_string)
                                edit_pos_string_output.append(edit_pos_string)
                                distance_median_all_output.append(distance_median_all)
                                quality_scores_all_output.append(quality_scores_all)
                                guide_starts_output.append(str(startpos_rev + len(pamsite) + 1))
                                guide_ends_output.append(str(startpos_rev + guidelength + len(pamsite) + 1))
                                chroms_output.append(str(chrom))

                        else:
                            if bes[be]['rev']['REF'] in edit_window_plus_bases:
                                ne_plus = False
                            else:
                                ne_plus = True
                            be_output_ne.append(be)
                            guide_output_ne.append(shared.revcom(possible_guide))
                            guide_with_pam_output_ne.append(shared.revcom(possible_guide_with_pam))
                            edit_window_output_ne.append(shared.revcom(edit_window_bases))
                            guide_output_strand_ne.append(possible_guide)
                            guide_with_pam_output_strand_ne.append(possible_guide_with_pam)
                            edit_window_output_strand_ne.append(edit_window_bases)
                            edit_window_plus_output_ne.append(shared.revcom(edit_window_plus_bases))
                            edit_window_plus_output_strand_ne.append(edit_window_plus_bases)
                            ne_plus_output_ne.append(ne_plus)
                            direction_output_ne.append('-') # reverse
                            safety_region_output_ne.append(shared.revcom(safety_region))
                            guide_starts_output_ne.append(startpos_rev + len(pamsite) + 1)
                            guide_ends_output_ne.append(startpos_rev + guidelength + len(pamsite) + 1)
                            chroms_output_ne.append(chrom)

                    startpos_rev += 1

    if not variant_output:
        if len(regions_list) == 1:
            raise ValueError("Your region was not found.")
        if len(regions_list) > 1:
            raise ValueError("None of your regions were found.")

    sgrnas = pl.DataFrame({"variant": variant_output,
                           "base_editor": be_output,
                           "guide": guide_output,
                           "guide_chrom": chroms_output,
                           "guide_start": guide_starts_output,
                           "guide_end": guide_ends_output,
                           "guide_with_pam": guide_with_pam_output,
                           "edit_window": edit_window_output,
                           "num_edits": num_edits_output,
                           "specific": specific_output,
                           "edit_window_plus": edit_window_plus_output,
                           "num_edits_plus": num_edits_plus_output,
                           "specific_plus": specific_plus_output,
                           "safety_region": safety_region_output,
                           "num_edits_safety": num_edits_safety_output,
                           "additional_in_safety": additional_in_safety_output,
                           "strand": direction_output,
                           "off_target_bases": edit_string_output,
                           "edited_positions": edit_pos_string_output,
                           "distance_median_all": distance_median_all_output,
                           "efficiency_scores_all": quality_scores_all_output})

    sgrnas_ne = pl.DataFrame({"variant":  "NA_for_non_editing_guides",
                              "base_editor": be_output_ne,
                              "guide": guide_output_ne,
                              "guide_chrom": chroms_output_ne,
                              "guide_start": guide_starts_output_ne,
                              "guide_end": guide_ends_output_ne,
                              "guide_with_pam": guide_with_pam_output_ne,
                              "edit_window": edit_window_output_ne,
                              "edit_window_plus": edit_window_plus_output_ne,
                              "safety_region": safety_region_output_ne,
                              "ne_plus": ne_plus_output_ne,
                              "strand": direction_output_ne})

    if blast:

        blast_guides.check_blastdb(ref_genome, False)

        sgrnas = sgrnas.with_row_index('index')
        sgrnas_ne = sgrnas_ne.with_row_index('index')

        guides_blast = sgrnas.select('index', 'guide')
        blast_results = blast_guides.guide_blast(guides_blast,
                                                 guidelength,
                                                 ref_genome,
                                                 'genes',
                                                 no_contigs)
        sgrnas = sgrnas.join(blast_results, on='index', how='left')

        guides_ne_blast = sgrnas_ne.select('index', 'guide')
        blast_results_ne = blast_guides.guide_blast(guides_ne_blast,
                                                    guidelength,
                                                    ref_genome,
                                                    'genes',
                                                    no_contigs)
        sgrnas_ne = sgrnas_ne.join(blast_results_ne, on='index', how='left')

        sgrnas = sgrnas.drop('index')
        sgrnas_ne =sgrnas_ne.drop('index')

    if vep: # better join by index

        variants_vep = sgrnas.select('variant').with_row_index('original_index').explode('variant')
        variants_vep_sorted = shared.sort_variantsdf(variants_vep)
        variants_vep_sorted_input = variants_vep_sorted['variant'].to_list()

        vep_info, vep_annotations = get_vep.get_vep_annotation(variants_vep_sorted_input,
                                                               species=vep_species,
                                                               assembly=vep_assembly,
                                                               dir_cache=vep_dir_cache,
                                                            #    dir_plugins=vep_dir_plugins, # currently not in use
                                                               cache_version=vep_cache_version,
                                                               flags=vep_flags)

        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_chrom = vep_annotations['#CHROM'].cast(str))
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_pos = vep_annotations['POS'].cast(str))
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_id = vep_annotations['ID'])
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_ref = vep_annotations['REF'])
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_alt = vep_annotations['ALT'])
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_qual = vep_annotations['QUAL'])
        # variants_vep_sorted = variants_vep_sorted.with_columns(vep_filter = vep_annotations['FILTER'])
        variants_vep_sorted = variants_vep_sorted.with_columns(vep_annotations['INFO'].alias('vep_info'))

        variants_vep_resorted = shared.resort_variantsdf(variants_vep_sorted)

        variants_vep_resorted = variants_vep_resorted.group_by('original_index', maintain_order=True).agg(pl.all())

        # sgrnas = sgrnas.with_columns(vep_chrom = variants_vep_resorted['vep_chrom'])
        # sgrnas = sgrnas.with_columns(vep_pos = variants_vep_resorted['vep_pos'])
        # sgrnas = sgrnas.with_columns(vep_id = variants_vep_resorted['vep_id'])
        # sgrnas = sgrnas.with_columns(vep_ref = variants_vep_resorted['vep_ref'])
        # sgrnas = sgrnas.with_columns(vep_alt = variants_vep_resorted['vep_alt'])
        # sgrnas = sgrnas.with_columns(vep_qual = variants_vep_resorted['vep_qual'])
        # sgrnas = sgrnas.with_columns(vep_filter = variants_vep_resorted['vep_filter'])
        sgrnas = sgrnas.with_columns(variants_vep_resorted['vep_info'].alias(f'vep_info ({vep_info})'))

        # sgrnas_ne = sgrnas_ne.with_columns(vep_chrom = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_pos = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_id = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_ref = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_alt = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_qual = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(vep_filter = pl.lit("NA_for_non_editing_guides"))
        # sgrnas_ne = sgrnas_ne.with_columns(pl.lit("NA_for_non_editing_guides").alias(f'vep_info ({vep_info})'))

    variant_cols_to_modify_second = ['variant']

    if vep:
        # variant_cols_to_modify_second += ['vep_chrom',
        #                                   'vep_pos',
        #                                   'vep_id',
        #                                   'vep_ref',
        #                                   'vep_alt',
        #                                   'vep_qual',
        #                                   'vep_filter',
        #                                   f'vep_info ({vep_info})']
        variant_cols_to_modify_second += [f'vep_info ({vep_info})']

    if aspect == 'exploded': # maybe also show off target edits in codons in small letters as for bedesigner
        sgrnas = sgrnas.explode(variant_cols_to_modify_second)

    elif aspect == 'collapsed':
        sgrnas = sgrnas.with_columns(
            pl.col(variant_cols_to_modify_second).list.join(";")
        )

    sgrnas = sgrnas.unique().sort(by=sgrnas.columns)
    sgrnas_ne = sgrnas_ne.unique().sort(by=sgrnas_ne.columns)

    sam_df = sgrnas.select(
        ['variant', 'guide', 'guide_chrom', 'guide_start', 'strand']
        ) # will always be exploded, but maybe unify somehow different?

    sam_ne_df = sgrnas_ne.select(
        ['variant', 'guide', 'guide_chrom', 'guide_start', 'strand']
        ) # will always be exploded, but maybe unify somehow different?

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

    sam_ne_df = sam_ne_df.with_columns(
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

    sgrnas = sgrnas.drop(['off_target_bases'])

    sgrnas_ne = sgrnas_ne.drop(['variant'])

    if edit_window_start_plus == 0 and edit_window_end_plus == 0:
        sgrnas = sgrnas.drop(['edit_window_plus',
                              'num_edits_plus',
                              'specific_plus',
                              'safety_region',
                              'num_edits_safety',
                              'additional_in_safety'])

        sgrnas_ne = sgrnas_ne.drop(['edit_window_plus',
                                    'safety_region',
                                    'ne_plus'])

    return (sgrnas, sgrnas_ne, sam_df, sam_ne_df, sequences_not_found)


def output_sgrnas(sgrnas, sgrnas_ne, output_file):
    sgrnas.write_csv(output_file + ".tsv", separator='\t')
    sgrnas_ne.write_csv(output_file + "_ne.tsv", separator='\t')


def output_guides_sam(sam_df, sam_ne_df, output_file, ref_genome):

    sam_df.write_csv(output_file + ".sam", separator='\t', include_header=False)
    pysam.view(output_file + ".sam",
               '-b',
               '-o', output_file + "_unsorted.bam",
               '-t', ref_genome + '.fai',
               catch_stdout=False)
    os.remove(output_file + ".sam")
    pysam.sort('-o', output_file + ".bam", output_file + "_unsorted.bam") # should already be sorted, but this acts as a failsafe
    os.remove(output_file + "_unsorted.bam")
    pysam.index(output_file + ".bam")

    sam_ne_df.write_csv(output_file + "_ne.sam", separator='\t', include_header=False)
    pysam.view(output_file + "_ne.sam",
               '-b',
               '-o', output_file + "_ne_unsorted.bam",
               '-t', ref_genome + '.fai',
               catch_stdout=False)
    os.remove(output_file + "_ne.sam")
    pysam.sort('-o', output_file + "_ne.bam", output_file + "_ne_unsorted.bam") # should already be sorted, but this acts as a failsafe
    os.remove(output_file + "_ne_unsorted.bam")
    pysam.index(output_file + "_ne.bam")


if __name__ == "__main__":
    pass