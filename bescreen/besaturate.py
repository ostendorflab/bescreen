import os
import re
import polars as pl
import pyfaidx
import pysam
import get_vep
import shared
import blast_guides


def saturate_bes(annotation_file,
                 ref_genome,
                 gene_symbols,
                 input_file,
                 pamsite,
                 edit_window_start,
                 edit_window_end,
                 guidelength,
                 baseeditor,
                 edit_window_start_plus,
                 edit_window_end_plus,
                 splice_sites,
                 write_parquet,
                 vep,
                 aspect,
                 vep_flags,
                 vep_species,
                 vep_assembly,
                 vep_dir_cache,
                 vep_dir_plugins,
                 vep_cache_version,
                 blast,
                 no_contigs,
                 filter_synonymous,
                 filter_splice_site,
                 filter_specific,
                 filter_missense):

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

        if 'symbol' in input_columns:
            gene_symbols_list = pl.read_csv(input_file, separator=',')['symbol'].to_list()
        else:
            raise ValueError('Input file does not contain the column "symbol"!\nPlease provide it!')

    elif gene_symbols:
        gene_symbols_list = gene_symbols.replace(' ', '').split(",")

    ref_genome_pyfaidx = pyfaidx.Fasta(ref_genome)

    parquet_file = shared.check_parquet(annotation_file, write_parquet)

    cdss = pl.read_parquet(parquet_file)
    cdss_gene = cdss.filter(pl.col('gene_name').is_in(gene_symbols_list))
    if cdss_gene.is_empty():
        if len(gene_symbols_list) == 1:
            raise ValueError("Your gene was not found.")
        if len(gene_symbols_list) > 1:
            raise ValueError("None of your genes were found.")

    genes_found = cdss_gene['gene_name'].unique().to_list()
    genes_not_found = [gene for gene in gene_symbols_list if gene not in genes_found]

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
    symbol_output = []
    transcript_output = []
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
    synonymous_output = []
    codon_ref = []
    aa_ref = []
    codon_edit = []
    aa_edit = []
    splice_site_included = []
    exon_number_output = []
    first_transcript_exon_output = []
    last_transcript_exon_output = []
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
    symbol_output_ne = []
    transcript_output_ne = []
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
    exon_number_output_ne = []
    first_transcript_exon_output_ne = []
    last_transcript_exon_output_ne = []
    safety_region_output_ne = []
    guide_starts_output_ne = []
    guide_ends_output_ne = []
    chroms_output_ne = []

    # precalculation of quality scores for positions in editing window
    distance_median_dict, quality_scores_dict = shared.qc_precalc(edit_window_start, edit_window_end)

    # iterate through genes and cds-exons
    for row in cdss_gene.iter_rows(named=True):

        gene_symbol = row['gene_name']
        transcript_symbol = row['transcript_name']
        exon_number = row['exon_number']
        first_transcript_exon = row['first_transcript_exon']
        last_transcript_exon = row['last_transcript_exon']
        chrom = row['Chromosome']

         # get cds-exon sequence
        sequence = ref_genome_pyfaidx[row["Chromosome"]][row["Start"]:row["End"]] # pyranges transforms gtf file into 0-based

        #generate forward search string with:
        # 5' overhang to account for first base is last base of editing window
        flanking_upstream_forward = ref_genome_pyfaidx[row["Chromosome"]][row["Start"] - flanklength_nonpamsided:row["Start"]]
        # 3' overhang to account for last base is first base of editing window plus PAM site
        flanking_downstream_forward = ref_genome_pyfaidx[row["Chromosome"]][row["End"]:row["End"] + flanklength_pamsided]
        # concat substring
        forward_search_string = str(flanking_upstream_forward) + str(sequence) + str(flanking_downstream_forward)
        if splice_sites:
            flanking_upstream_forward_splice = ref_genome_pyfaidx[row["Chromosome"]][row["Start"] - flanklength_nonpamsided - 2:row["Start"] - flanklength_nonpamsided]
            flanking_downstream_forward_splice = ref_genome_pyfaidx[row["Chromosome"]][row["End"] + flanklength_pamsided:row["End"] + flanklength_pamsided + 2]
            forward_search_string = str(flanking_upstream_forward_splice) + forward_search_string + str(flanking_downstream_forward_splice)

        #generate reverse search string with:
        # 5' overhang to account for first base is first base of editing window plus PAM site
        flanking_upstream_reverse = ref_genome_pyfaidx[row["Chromosome"]][row["Start"] - flanklength_pamsided:row["Start"]]
        # 3' overhang to account for last base is last base of editing window
        flanking_downstream_reverse = ref_genome_pyfaidx[row["Chromosome"]][row["End"]:row["End"] + flanklength_nonpamsided]
        # concat substring
        reverse_search_string = str(flanking_upstream_reverse) + str(sequence) + str(flanking_downstream_reverse)
        if splice_sites:
            flanking_upstream_reverse_splice = ref_genome_pyfaidx[row["Chromosome"]][row["Start"] - flanklength_pamsided - 2:row["Start"] - flanklength_pamsided]
            flanking_downstream_reverse_splice = ref_genome_pyfaidx[row["Chromosome"]][row["End"] + flanklength_nonpamsided:row["End"] + flanklength_nonpamsided + 2]
            reverse_search_string = str(flanking_upstream_reverse_splice) + reverse_search_string + str(flanking_downstream_reverse_splice)

        for be in bes: # for ABE and CBE
            startpos_fwd = row["Start"] - flanklength_nonpamsided
            startpos_rev = row["Start"] - flanklength_pamsided
            if splice_sites:
                startpos_fwd = row["Start"] - flanklength_nonpamsided - 2
                startpos_rev = row["Start"] - flanklength_pamsided - 2


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
                        codons = [] # reset variant list for sliding window
                        codons_edited = [] # reset variant list for sliding window
                        aas = [] # reset variant list for sliding window
                        aas_edited = [] # reset variant list for sliding window
                        # includes_splice_site = "not_tested" # also check for splice-sites even, if not explicitely added
                        # if splice_sites:
                        #     includes_splice_site = False
                        includes_splice_site = False

                        poss = [row["Start"] - flanklength_nonpamsided + index_edit - 2 if splice_sites else row["Start"] - flanklength_nonpamsided + index_edit for index_edit in indices_edit]
                        start = row["Start"] - 2 if splice_sites else row["Start"]
                        end = row["End"] + 2 if splice_sites else row["End"]
                        if any(start <= pos < end for pos in poss):

                            for index_edit in indices_edit:

                                pos = row["Start"] - flanklength_nonpamsided + index_edit # pos for fasta

                                if splice_sites:
                                    pos -= 2

                                if bes[be]['fwd']['REF'] == ref_genome_pyfaidx[row["Chromosome"]][pos]: # failsafe; can be removed

                                    variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["fwd"]["ALT"]}') # append variant to variantlist 1-based as VCF; add here already to not miss off-target edits outside CDS

                                    if row["Start"] <= pos < row["End"]: # assure pos is in cds-exon; <= or <?: test showed < for row["End"]
                                        # variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["fwd"]["ALT"]}') # append variant to variantlist 1-based as VCF

                                        index_edit_cds = index_edit - flanklength_nonpamsided # index cds-exon sequence
                                        if splice_sites:
                                            index_edit_cds -= 2

                                        if index_edit_cds < 0: # failsafe; can be removed
                                            raise Exception("index_edit_cds < 0")

                                        if str(row['Strand']) == '-':
                                            index_edit_cds = len(str(sequence)) - 1 - index_edit_cds
                                        offset = shared.get_offset(str(row['Strand']), int(row['Frame']), index_edit_cds) # with new offset, if should be independent of search direction

                                        if (row["Start"] <= (pos - offset) < row["End"]) and (row["Start"] <= (pos + 3 - offset) < row["End"]): # <= or <?: test showed < for row["End"]
                                            codon = str(ref_genome_pyfaidx[row["Chromosome"]][pos - offset:pos + 3 - offset]) # with new offset, if should be independent of search direction
                                            codon_edited = shared.replace_str_index(codon, offset, bes[be]['fwd']['ALT'])

                                            if str(row['Strand']) == '-':
                                                codon = shared.revcom(codon)
                                                codon_edited = shared.revcom(codon_edited)
                                            aa = shared.codon_sun_one_letter[codon]
                                            aa_edited = shared.codon_sun_one_letter[codon_edited]

                                            if aa == shared.codon_sun_one_letter["ATG"] and exon_number == first_transcript_exon:
                                                if (row['Strand'] == '+' and (row["Start"] <= pos < (row["Start"] + 3))) or (row['Strand'] == '-' and ((row["End"] - 3) <= pos < row["End"])):
                                                    codon = "Start" + codon
                                                    aa = "Start" + aa

                                        else:
                                            codon = "incomplete_codon"
                                            codon_edited = "incomplete_codon"
                                            aa = "codon_incomplete"
                                            aa_edited = "codon_incomplete"

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        if len(indices_edit) == 1: # only 1 editable base in editing window
                                            # specific = True # will be analyzed by analyze_guide()
                                            # if len(indices_plus_edit) == 1:
                                            #     specific_plus = True
                                            # else:
                                            #     specific_plus = False
                                            if aa == aa_edited:
                                                synonymous = True
                                            else:
                                                synonymous = False
                                        else: # if more than 1 editable base in editing window
                                            # specific = False # will be analyzed by analyze_guide()
                                            synonymous = False
                                            # specific_plus = False # will be analyzed by analyze_guide()

                                    elif (row["Start"] - 2) <= pos < (row["End"] + 2):
                                        # if splice_sites: # also check for splice-sites even, if not explicitely added
                                            # variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["fwd"]["ALT"]}') # append variant to variantlist 1-based as VCF
                                            includes_splice_site = True

                                            if (row["Start"] - 2) <= pos < row["Start"]:
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
                                            elif row["End"] <= pos < (row["End"] + 2):
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

                                            codons.append(codon)
                                            codons_edited.append(codon_edited)
                                            aas.append(aa)
                                            aas_edited.append(aa_edited)

                                            # if len(indices_edit) == 1: # only 1 editable base in editing window # whole block not necessary anymore since using analyze_guide(); replaced with synonymous = False
                                            #     # specific = True # will be analyzed by analyze_guide()
                                            #     synonymous = False
                                            #     # if len(indices_plus_edit) == 1: # will be analyzed by analyze_guide()
                                            #     #     specific_plus = True
                                            #     # else:
                                            #     #     specific_plus = False
                                            # else: # if more than 1 editable base in editing window
                                            #     # specific = False # will be analyzed by analyze_guide()
                                            #     synonymous = False
                                            #     # specific_plus = False # will be analyzed by analyze_guide()
                                            synonymous = False

                                    else:
                                        codons.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        codons_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        aas.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        aas_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?

                                else: # failsafe; can be removed
                                    raise Exception

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

                            if not ((filter_synonymous) and (not synonymous) or
                                    (filter_splice_site) and (not includes_splice_site) or
                                    (filter_specific) and (not specific) or
                                    (filter_missense) and (aas == aas_edited)):

                                symbol_output.append(gene_symbol)
                                transcript_output.append(transcript_symbol)
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
                                synonymous_output.append(synonymous)
                                direction_output.append('+') # forward
                                codon_ref.append(codons) # list
                                aa_ref.append(aas) # list
                                codon_edit.append(codons_edited) # list
                                aa_edit.append(aas_edited) # list
                                splice_site_included.append(includes_splice_site)
                                exon_number_output.append(exon_number)
                                first_transcript_exon_output.append(first_transcript_exon)
                                last_transcript_exon_output.append(last_transcript_exon)
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
                        symbol_output_ne.append(gene_symbol)
                        transcript_output_ne.append(transcript_symbol)
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
                        exon_number_output_ne.append(exon_number)
                        first_transcript_exon_output_ne.append(first_transcript_exon)
                        last_transcript_exon_output_ne.append(last_transcript_exon)
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

                        poss = [row["Start"] - flanklength_pamsided + index_edit - 2 if splice_sites else row["Start"] - flanklength_pamsided + index_edit for index_edit in indices_edit]
                        start = row["Start"] - 2 if splice_sites else row["Start"]
                        end = row["End"] + 2 if splice_sites else row["End"]
                        if any(start <= pos < end for pos in poss):

                            for index_edit in indices_edit:

                                pos = row["Start"] - flanklength_pamsided + index_edit # pos for fasta

                                if splice_sites:
                                    pos -= 2

                                if bes[be]['rev']['REF'] == ref_genome_pyfaidx[row["Chromosome"]][pos]: # failsafe; can be removed

                                    variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["rev"]["ALT"]}') # append variant to variantlist 1-based as VCF; add here already to not miss off-target edits outside CDS

                                    if row["Start"] <= pos < row["End"]: # assure pos is in cds-exon; <= or <?: test showed < for row["End"]
                                        # variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["rev"]["ALT"]}') # append variant to variantlist 1-based as VCF

                                        index_edit_cds = index_edit - flanklength_pamsided # index cds-exon sequence
                                        if splice_sites:
                                            index_edit_cds -= 2

                                        if index_edit_cds < 0: # failsafe; can be removed
                                            raise Exception("index_edit_cds < 0")

                                        if str(row['Strand']) == '-':
                                            index_edit_cds = len(str(sequence)) - 1 - index_edit_cds
                                        offset = shared.get_offset(str(row['Strand']), int(row['Frame']), index_edit_cds) # with new offset, if should be independent of search direction

                                        if (row["Start"] <= (pos - offset) < row["End"]) and (row["Start"] <= (pos + 3 - offset) < row["End"]): # <= or <?: test showed < for row["End"]
                                            codon = str(ref_genome_pyfaidx[row["Chromosome"]][pos - offset:pos + 3 - offset]) # with new offset, if should be independent of search direction
                                            codon_edited = shared.replace_str_index(codon, offset, bes[be]['rev']['ALT'])

                                            if str(row['Strand']) == '-':
                                                codon = shared.revcom(codon)
                                                codon_edited = shared.revcom(codon_edited)
                                            aa = shared.codon_sun_one_letter[codon]
                                            aa_edited = shared.codon_sun_one_letter[codon_edited]

                                            if aa == shared.codon_sun_one_letter["ATG"] and exon_number == first_transcript_exon:
                                                if (row['Strand'] == '+' and (row["Start"] <= pos < (row["Start"] + 3))) or (row['Strand'] == '-' and ((row["End"] - 3) <= pos < row["End"])):
                                                    codon = "Start" + codon
                                                    aa = "Start" + aa

                                        else:
                                            codon = "incomplete_codon"
                                            codon_edited = "incomplete_codon"
                                            aa = "codon_incomplete"
                                            aa_edited = "codon_incomplete"

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        if len(indices_edit) == 1: # only 1 editable base in editing window
                                            # specific = True # will be analyzed by analyze_guide()
                                            # if len(indices_plus_edit) == 1:
                                            #     specific_plus = True
                                            # else:
                                            #     specific_plus = False
                                            if aa == aa_edited:
                                                synonymous = True
                                            else:
                                                synonymous = False
                                        else: # if more than 1 editable base in editing window
                                            # specific = False # will be analyzed by analyze_guide()
                                            synonymous = False
                                            # specific_plus = False # will be analyzed by analyze_guide()

                                    elif (row["Start"] - 2) <= pos < (row["End"] + 2):
                                        # if splice_sites: # also check for splice-sites even, if not explicitely added
                                            # variants.append(f'{row["Chromosome"]}_{pos + 1}_{ref_genome_pyfaidx[row["Chromosome"]][pos]}_{bes[be]["rev"]["ALT"]}') # append variant to variantlist 1-based as VCF
                                            includes_splice_site = True

                                            if (row["Start"] - 2) <= pos < row["Start"]:
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
                                            elif row["End"] <= pos < (row["End"] + 2):
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

                                            codons.append(codon)
                                            codons_edited.append(codon_edited)
                                            aas.append(aa)
                                            aas_edited.append(aa_edited)

                                            # if len(indices_edit) == 1: # only 1 editable base in editing window # whole block not necessary anymore since using analyze_guide(); replaced with synonymous = False
                                            #     # specific = True # will be analyzed by analyze_guide()
                                            #     synonymous = False
                                            #     # if len(indices_plus_edit) == 1: # will be analyzed by analyze_guide()
                                            #     #     specific_plus = True
                                            #     # else:
                                            #     #     specific_plus = False
                                            # else: # if more than 1 editable base in editing window
                                            #     # specific = False # will be analyzed by analyze_guide()
                                            #     synonymous = False
                                            #     # specific_plus = False # will be analyzed by analyze_guide()
                                            synonymous = False

                                    else:
                                        codons.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        codons_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        aas.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        aas_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?

                                else: # failsafe; can be removed
                                    raise Exception

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

                            if not ((filter_synonymous) and (not synonymous) or
                                    (filter_splice_site) and (not includes_splice_site) or
                                    (filter_specific) and (not specific) or
                                    (filter_missense) and (aas == aas_edited)):

                                symbol_output.append(gene_symbol)
                                transcript_output.append(transcript_symbol)
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
                                synonymous_output.append(synonymous)
                                direction_output.append('-') # reverse
                                codon_ref.append(codons) # list
                                aa_ref.append(aas) # list
                                codon_edit.append(codons_edited) # list
                                aa_edit.append(aas_edited) # list
                                splice_site_included.append(includes_splice_site)
                                exon_number_output.append(exon_number)
                                first_transcript_exon_output.append(first_transcript_exon)
                                last_transcript_exon_output.append(last_transcript_exon)
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
                        symbol_output_ne.append(gene_symbol)
                        transcript_output_ne.append(transcript_symbol)
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
                        exon_number_output_ne.append(exon_number)
                        first_transcript_exon_output_ne.append(first_transcript_exon)
                        last_transcript_exon_output_ne.append(last_transcript_exon)
                        safety_region_output_ne.append(shared.revcom(safety_region))
                        guide_starts_output_ne.append(startpos_rev + len(pamsite) + 1)
                        guide_ends_output_ne.append(startpos_rev + guidelength + len(pamsite) + 1)
                        chroms_output_ne.append(chrom)

                startpos_rev += 1

    # cols_do_modify = ["variant",
    #                   "codon_ref",
    #                   "aa_ref",
    #                   "codon_edit",
    #                   "aa_edit"]

    sgrnas = pl.DataFrame({"variant": variant_output,
                           "base_editor": be_output,
                           "symbol": symbol_output,
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
                           "ne_plus": "NA_for_editing_guides",
                           "synonymous": synonymous_output,
                           "strand": direction_output,
                           "codon_ref": codon_ref,
                           "aa_ref": aa_ref,
                           "codon_edit": codon_edit,
                           "aa_edit": aa_edit,
                           "splice_site_included": splice_site_included,
                           "originally_intended_ALT": "NA_for_genes",
                           "ref_match": "NA_for_genes",
                           "off_target_bases": edit_string_output,
                           "edited_positions": edit_pos_string_output,
                           "specificity":  "NA_for_genes",
                           "distance_median_variant":  "NA_for_genes",
                           "efficiency_scores_variant":  "NA_for_genes",
                           "distance_median_all": distance_median_all_output,
                           "efficiency_scores_all": quality_scores_all_output,
                           "transcript": transcript_output,
                           "exon_number": exon_number_output,
                           "first_transcript_exon": first_transcript_exon_output,
                           "last_transcript_exon": last_transcript_exon_output})

    sgrnas_cols = sgrnas.columns

    # for col in cols_do_modify:
    #     sgrnas = sgrnas.with_columns(sgrnas[col].alias(col))

    sgrnas = sgrnas.with_columns(exon_number = sgrnas['exon_number'].cast(int).cast(str))
    sgrnas = sgrnas.with_columns(first_transcript_exon = sgrnas['first_transcript_exon'].cast(int).cast(str))
    sgrnas = sgrnas.with_columns(last_transcript_exon = sgrnas['last_transcript_exon'].cast(int).cast(str))
    sgrnas = sgrnas.with_columns(synonymous = sgrnas['synonymous'].cast(str))
    sgrnas = sgrnas.with_columns(splice_site_included = sgrnas['splice_site_included'].cast(str))
    sgrnas = sgrnas.with_columns(variant = sgrnas['variant'])
    sgrnas = sgrnas.with_columns(codon_ref = sgrnas['codon_ref'])
    sgrnas = sgrnas.with_columns(aa_ref = sgrnas['aa_ref'])
    sgrnas = sgrnas.with_columns(codon_edit = sgrnas['codon_edit'])
    sgrnas = sgrnas.with_columns(aa_edit = sgrnas['aa_edit'])
    # sgrnas = sgrnas.with_columns(variant = sgrnas['variant'].list.join('&'))
    # sgrnas = sgrnas.with_columns(codon_ref = sgrnas['codon_ref'].list.join('&'))
    # sgrnas = sgrnas.with_columns(aa_ref = sgrnas['aa_ref'].list.join('&'))
    # sgrnas = sgrnas.with_columns(codon_edit = sgrnas['codon_edit'].list.join('&'))
    # sgrnas = sgrnas.with_columns(aa_edit = sgrnas['aa_edit'].list.join('&'))
    sgrnas = sgrnas.group_by([col for col in sgrnas.columns if col not in ['exon_number',
                                                                          'transcript',
                                                                          'first_transcript_exon',
                                                                          'last_transcript_exon']],
                              maintain_order=True).agg(pl.all())
    sgrnas = sgrnas.group_by([col for col in sgrnas.columns if col not in ['synonymous',
                                                                          'codon_ref',
                                                                          'aa_ref',
                                                                          'codon_edit',
                                                                          'aa_edit',
                                                                          'splice_site_included',
                                                                          'exon_number',
                                                                          'transcript',
                                                                          'first_transcript_exon',
                                                                          'last_transcript_exon']],
                              maintain_order=True).agg(pl.all())

    sgrnas = sgrnas[sgrnas_cols]

    sgrnas_ne = pl.DataFrame({"variant":  "NA_for_non_editing_guides",
                              "base_editor": be_output_ne,
                              "symbol": symbol_output_ne,
                              "guide": guide_output_ne,
                              "guide_chrom": chroms_output_ne,
                              "guide_start": guide_starts_output_ne,
                              "guide_end": guide_ends_output_ne,
                              "guide_with_pam": guide_with_pam_output_ne,
                              "edit_window": edit_window_output_ne,
                              "num_edits": "NA_for_non_editing_guides",
                              "specific": "NA_for_non_editing_guides",
                              "edit_window_plus": edit_window_plus_output_ne,
                              "num_edits_plus": "NA_for_non_editing_guides",
                              "specific_plus": "NA_for_non_editing_guides",
                              "safety_region": safety_region_output_ne,
                              "num_edits_safety": "NA_for_non_editing_guides",
                              "additional_in_safety": "NA_for_non_editing_guides",
                              "ne_plus": ne_plus_output_ne,
                              "synonymous": "NA_for_non_editing_guides",
                              "strand": direction_output_ne,
                              "codon_ref": "NA_for_non_editing_guides",
                              "aa_ref": "NA_for_non_editing_guides",
                              "codon_edit": "NA_for_non_editing_guides",
                              "aa_edit": "NA_for_non_editing_guides",
                              "splice_site_included": "NA_for_non_editing_guides",
                              "originally_intended_ALT": "NA_for_genes",
                              "ref_match": "NA_for_genes",
                              "off_target_bases": "NA_for_non_editing_guides",
                              "edited_positions": "NA_for_non_editing_guides",
                              "specificity": "NA_for_genes",
                              "distance_median_variant": "NA_for_genes",
                              "efficiency_scores_variant": "NA_for_genes",
                              "distance_median_all": "NA_for_non_editing_guides",
                              "efficiency_scores_all": "NA_for_non_editing_guides",
                              "transcript": transcript_output_ne,
                              "exon_number": exon_number_output_ne,
                              "first_transcript_exon": first_transcript_exon_output_ne,
                              "last_transcript_exon": last_transcript_exon_output_ne})

    sgrnas_ne_cols = sgrnas_ne.columns

    # for col in cols_do_modify:
    #     sgrnas_ne = sgrnas_ne.with_columns(sgrnas_ne[col].alias(col))

    sgrnas_ne = sgrnas_ne.with_columns(exon_number = sgrnas_ne['exon_number'].cast(int).cast(str))
    sgrnas_ne = sgrnas_ne.with_columns(first_transcript_exon = sgrnas_ne['first_transcript_exon'].cast(int).cast(str))
    sgrnas_ne = sgrnas_ne.with_columns(last_transcript_exon = sgrnas_ne['last_transcript_exon'].cast(int).cast(str))
    sgrnas_ne = sgrnas_ne.with_columns(synonymous = sgrnas_ne['synonymous'].cast(str))
    sgrnas_ne = sgrnas_ne.with_columns(splice_site_included = sgrnas_ne['splice_site_included'].cast(str))
    sgrnas_ne = sgrnas_ne.group_by([col for col in sgrnas_ne.columns if col not in ['exon_number',
                                                                                   'transcript',
                                                                                   'first_transcript_exon',
                                                                                   'last_transcript_exon']],
                                    maintain_order=True).agg(pl.all().str.join("~"))
    sgrnas_ne = sgrnas_ne.group_by([col for col in sgrnas_ne.columns if col not in ['synonymous',
                                                                                   'codon_ref',
                                                                                   'aa_ref',
                                                                                   'codon_edit',
                                                                                   'aa_edit',
                                                                                   'splice_site_included',
                                                                                   'exon_number',
                                                                                   'transcript',
                                                                                   'first_transcript_exon',
                                                                                   'last_transcript_exon']],
                                    maintain_order=True).agg(pl.all().str.join("^")) # should not happen since no annotations

    sgrnas_ne = sgrnas_ne[sgrnas_ne_cols]

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

    transcript_cols_to_modify_first = ['exon_number',
                                       'transcript',
                                       'first_transcript_exon',
                                       'last_transcript_exon']

    variant_cols_to_modify_first = ['codon_ref',
                                    'aa_ref',
                                    'codon_edit',
                                    'aa_edit']

    variant_cols_to_modify_second = ['variant']

    transcript_cols_to_modify_second = ['synonymous',
                                        'codon_ref',
                                        'aa_ref',
                                        'codon_edit',
                                        'aa_edit',
                                        'splice_site_included',
                                        'exon_number',
                                        'transcript',
                                        'first_transcript_exon',
                                        'last_transcript_exon']

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
        sgrnas = sgrnas.with_columns(
            pl.col(transcript_cols_to_modify_first).list.eval(pl.element().list.join("~")).list.join("^"),
            pl.col(variant_cols_to_modify_first).map_batches(lambda s:
                s.to_frame()
                .with_row_index("row")
                .explode(pl.last())
                .with_row_index()
                .explode(pl.last())
                .with_columns(pl.int_ranges(pl.col("index").rle().struct.field("len")).flatten().alias("index"))
                .group_by("index", "row", maintain_order=True)
                .agg(pl.last())
                .group_by("row", maintain_order=True)
                .agg(pl.last())
                .select(pl.last())
                .to_series() # this transposes the columns in variant_cols_to_modify_first
                .list.eval(pl.element().list.join("^"))),
            pl.col([col for col in transcript_cols_to_modify_second if col not in transcript_cols_to_modify_first + variant_cols_to_modify_first]).list.join("^")
        ).explode(variant_cols_to_modify_second + variant_cols_to_modify_first)

    elif aspect == 'collapsed':
        sgrnas = sgrnas.with_columns(
            pl.col(transcript_cols_to_modify_first).list.eval(pl.element().list.join("~")).list.join("^"),
            pl.col(variant_cols_to_modify_first).list.eval(pl.element().list.join(";")).list.join("^"),
            pl.col(variant_cols_to_modify_second).list.join(";"),
            pl.col([col for col in transcript_cols_to_modify_second if col not in transcript_cols_to_modify_first + variant_cols_to_modify_first]).list.join("^")
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

    sgrnas = sgrnas.drop(['ne_plus',
                          'originally_intended_ALT',
                          'ref_match',
                          'off_target_bases',
                          'specificity',
                          'distance_median_variant',
                          'efficiency_scores_variant'])

    sgrnas_ne = sgrnas_ne.drop(['variant',
                                'num_edits',
                                'specific',
                                'num_edits_plus',
                                'specific_plus',
                                'num_edits_safety',
                                'additional_in_safety',
                                'synonymous',
                                'codon_ref',
                                'aa_ref',
                                'codon_edit',
                                'aa_edit',
                                'splice_site_included',
                                'originally_intended_ALT',
                                'ref_match',
                                'off_target_bases',
                                'edited_positions',
                                'specificity',
                                'distance_median_variant',
                                'efficiency_scores_variant',
                                'distance_median_all',
                                'efficiency_scores_all'])

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

    return (sgrnas, sgrnas_ne, sam_df, sam_ne_df, genes_not_found)


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
