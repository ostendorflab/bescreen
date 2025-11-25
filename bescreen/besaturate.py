import os
import re
import polars as pl
import pyfaidx
import pysam
import get_vep
import shared
import blast_guides
import itertools


def saturate_bes(annotation_file,
                 ref_genome,
                 gene_symbols,
                 input_file,
                 pamsite,
                 fiveprimepam,
                 edit_window_start,
                 edit_window_end,
                 guidelength,
                 basechange,
                 edit_window_start_plus,
                 edit_window_end_plus,
                 splice_sites,
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
    if not basechange:
        raise ValueError('Please select at least one base editor!')
    elif not 'all' in basechange:
        bes = {key: bes[key] for key in basechange}

    if input_file:

        input_df = pl.read_csv(input_file, separator=',')
        input_columns = input_df.columns

        if 'symbol' in input_columns:
            gene_symbols_list = pl.read_csv(input_file, separator=',')['symbol'].to_list()
        else:
            raise ValueError('Input file does not contain the column "symbol"!\nPlease provide it!')

    elif gene_symbols:
        gene_symbols_list = gene_symbols.replace(' ', '').split(",")

    # real_gene_symbols_list = [gene for gene in gene_symbols_list if len(gene.split('-')) == 1] # some genes have hyphen in their name, so this solution is wrong
    # transcript_symbols_list = [transcript for transcript in gene_symbols_list if len(transcript.split('-')) == 2] # some genes have hyphen in their name, so this solution is wrong
    # genes_malformed = [gene for gene in gene_symbols_list if len(gene.split('-')) not in [1, 2]] # some genes have hyphen in their name, so this solution is wrong

    ref_genome_pyfaidx = pyfaidx.Fasta(ref_genome)

    parquet_file = shared.check_parquet(annotation_file, write_parquet)

    cdss = pl.read_parquet(parquet_file)

    all_genes = cdss.get_column('gene_name').unique().to_list()
    all_transcripts = cdss.get_column('transcript_name').unique().to_list()
    all_genes_transcripts = all_genes + all_transcripts
    all_genes = set(all_genes)
    all_transcripts = set(all_transcripts)
    all_genes_transcripts = set(all_genes_transcripts)

    # real_gene_symbols_list = [gene for gene in gene_symbols_list if gene in all_genes] # this works and keeps duplicates; not sure what's better
    # transcript_symbols_list = [transcript for transcript in gene_symbols_list if transcript in all_transcripts] # this works and keeps duplicates; not sure what's better
    # genes_malformed = [malformed for malformed in gene_symbols_list if malformed not in all_genes_transcripts] # this works and keeps duplicates; not sure what's better

    gene_symbols_set = set(gene_symbols_list) # this works and removes duplicates; not sure what's better

    real_gene_symbols_list = list(gene_symbols_set & all_genes) # this works and removes duplicates; not sure what's better
    transcript_symbols_list = list(gene_symbols_set & all_transcripts) # this works and removes duplicates; not sure what's better
    genes_malformed = list(gene_symbols_set - all_genes_transcripts) # this works and removes duplicates; not sure what's better

    if mane_select_only:
        cdss = cdss.filter(pl.col('MANE_Select'))

    cdss_real_gene = cdss.filter(pl.col('gene_name').is_in(real_gene_symbols_list))
    cdss_transcript = cdss.filter(pl.col('transcript_name').is_in(transcript_symbols_list))
    cdss_gene = pl.concat([cdss_real_gene, cdss_transcript]).unique()

    if cdss_gene.is_empty():
        if len(gene_symbols_list) == 1:
            raise ValueError(f"The gene {', '.join(gene_symbols_list)} was not found in the annotation file.")
        if len(gene_symbols_list) > 1:
            raise ValueError(f"None of the {len(gene_symbols_list)} genes {', '.join(gene_symbols_list)} were found in the annotation file.")

    real_genes_found = cdss_real_gene['gene_name'].unique().to_list()
    real_genes_not_found = [gene for gene in real_gene_symbols_list if gene not in real_genes_found]
    transcripts_found = cdss_transcript['transcript_name'].unique().to_list()
    transcripts_not_found = [transcript for transcript in transcript_symbols_list if transcript not in transcripts_found]
    genes_not_found = real_genes_not_found + transcripts_not_found + genes_malformed

    pamsite_relevant = pamsite.lstrip('N') # needs to be rstrip() for 5' PAM

    pamlist = list(pamsite_relevant)
    pamlist_real = [shared.iupac_nt_code.get(item, item) for item in pamlist]
    pamlist_real_all = list(itertools.product(*pamlist_real))
    pamlist_real_string = [''.join(pam_real_string) for pam_real_string in pamlist_real_all]

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
    consequence_output = []
    codon_ref = []
    aa_ref = []
    aa_pos = []
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
    distance_to_center_output = []
    # quality_scores_all_output = []
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

                if any(possible_guide_with_pam.endswith(pam_real_string) for pam_real_string in pamlist_real_string): # PAM found at 5' end

                    if bes[be]['fwd']['REF'] in edit_window_bases: # editable base in edit window

                        # find all indices for possible edits in the edit window in search string
                        indices_edit = [i + (edit_window_start - 1) + m.start() for m in re.finditer(bes[be]['fwd']['REF'], edit_window_bases)]
                        indices_plus_edit = [i + (edit_window_start - 1 - edit_window_start_plus) + m.start() for m in re.finditer(bes[be]['fwd']['REF'], edit_window_plus_bases)] # is this correct?

                        variants = [] # reset variant list for sliding window
                        codons = [] # reset variant list for sliding window
                        codons_edited = [] # reset variant list for sliding window
                        aas = [] # reset variant list for sliding window
                        aas_edited = [] # reset variant list for sliding window
                        consequence = [] # reset variant list for sliding window
                        aa_positions = [] # reset variant list for sliding window
                        # includes_splice_site = "not_tested" # also check for splice-sites even, if not explicitely added
                        # if splice_sites:
                        #     includes_splice_site = False
                        includes_splice_site = False
                        synonymous = False

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
                                            try:
                                                codon_list = []
                                                for ntpos in range(pos - offset, pos + 3 - offset):
                                                    if ntpos < row["Start"]:
                                                        new_exon = cdss_gene.filter(
                                                            (pl.col('gene_name') == gene_symbol) &
                                                            (pl.col('transcript_name') == transcript_symbol) &
                                                            (pl.col('exon_number') == exon_number - 1) &
                                                            (pl.col('first_transcript_exon') == first_transcript_exon) &
                                                            (pl.col('last_transcript_exon') == last_transcript_exon) &
                                                            (pl.col('Chromosome') == chrom)
                                                            )
                                                        if not len(new_exon) == 0:
                                                            raise AssertionError('More then one previous exon!')
                                                        new_exon = new_exon.to_dict(as_series=False)
                                                        new_ntpos = new_exon['End'][0] - (row["Start"] - ntpos)
                                                        codon_list.append(str(ref_genome_pyfaidx[new_exon['Chromosome'][0]][new_ntpos]))
                                                    elif row["Start"] <= ntpos < row["End"]:
                                                        codon_list.append(str(ref_genome_pyfaidx[row["Chromosome"]][ntpos]))
                                                    elif row["End"] <= ntpos:
                                                        new_exon = cdss_gene.filter(
                                                            (pl.col('gene_name') == gene_symbol) &
                                                            (pl.col('transcript_name') == transcript_symbol) &
                                                            (pl.col('exon_number') == exon_number + 1) &
                                                            (pl.col('first_transcript_exon') == first_transcript_exon) &
                                                            (pl.col('last_transcript_exon') == last_transcript_exon) &
                                                            (pl.col('Chromosome') == chrom)
                                                            )
                                                        if not len(new_exon) == 0:
                                                            raise AssertionError('More then one next exon!')
                                                        new_exon = new_exon.to_dict(as_series=False)
                                                        new_ntpos = new_exon['Start'][0] + (ntpos - row["End"])
                                                        codon_list.append(str(ref_genome_pyfaidx[new_exon['Chromosome'][0]][new_ntpos]))
                                                codon = ''.join(codon_list)
                                                codon_edited = shared.replace_str_index(codon, offset, bes[be]['fwd']['ALT']) # specific for search direction

                                                if str(row['Strand']) == '-':
                                                    codon = shared.revcom(codon)
                                                    codon_edited = shared.revcom(codon_edited)
                                                aa = shared.codon_sun_one_letter[codon]
                                                aa_edited = shared.codon_sun_one_letter[codon_edited]

                                                if aa == shared.codon_sun_one_letter["ATG"] and exon_number == first_transcript_exon:
                                                    if (row['Strand'] == '+' and (row["Start"] <= pos < (row["Start"] + 3))) or (row['Strand'] == '-' and ((row["End"] - 3) <= pos < row["End"])):
                                                        codon = "Start" + codon
                                                        aa = "Start" + aa
                                            except:
                                                codon = "incomplete_codon"
                                                codon_edited = "incomplete_codon"
                                                aa = "codon_incomplete"
                                                aa_edited = "codon_incomplete"

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        if aa in non_stop_aas and aa_edited in non_stop_aas and aa != aa_edited:
                                            consequenc = 'missense'
                                        elif aa in non_stop_aas and aa_edited in non_stop_aas and aa == aa_edited:
                                            consequenc = 'synonymous'

                                        elif aa != 'Stop' and aa_edited == 'Stop':
                                            consequenc = 'stopgain'
                                        elif aa == 'Stop' and aa_edited != 'Stop':
                                            consequenc = 'stoplost'
                                        elif aa == 'Stop' and aa_edited == 'Stop': # can this happen?
                                            consequenc = 'synonymous'

                                        elif aa == 'StartM' and aa_edited !='M': # StartM not used in aas_edited yet; Any edit in StartM would be a startlost
                                            consequenc = 'startlost'
                                        elif aa == 'StartM' and aa_edited =='M': # can't happen
                                            consequenc = 'synonymous'

                                        elif aa == "codon_incomplete" and aa_edited == "codon_incomplete":
                                            consequenc = 'indefinite'

                                        else: # can this happen?
                                            consequenc = 'complex'

                                        consequence.append(consequenc)

                                        if row['Strand'] == "+":
                                            nt_position = row['transcript_length_before'] + ((pos + 1) - row['Start'])
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                        elif row['Strand'] == "-":
                                            nt_position = row['transcript_length_before'] + (row['End'] - (pos)) # (row['End'] - (pos + 1))?
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                        aa_positions.append(aa_position)

                                        if len(indices_edit) == 1: # only 1 editable base in editing window
                                            # specific = True # will be analyzed by analyze_guide()
                                            # if len(indices_plus_edit) == 1:
                                            #     specific_plus = True
                                            # else:
                                            #     specific_plus = False
                                            if aa == aa_edited and aa != "codon_incomplete":
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
                                            consequenc = 'splice_site'

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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'

                                            codons.append(codon)
                                            codons_edited.append(codon_edited)
                                            aas.append(aa)
                                            aas_edited.append(aa_edited)

                                            consequence.append(consequenc)

                                            aa_positions.append(aa)

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
                                        # codons.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # codons_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # aas.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # aas_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # this happens due to off-target effects in introns or the UTRs
                                        codon = "intron"
                                        codon_edited = "intron"
                                        aa = "intron"
                                        aa_edited = "intron"
                                        consequenc = "intron"

                                        if pos < (row["Start"] - 2):
                                            if row["Strand"] == '+' and exon_number == first_transcript_exon:
                                                codon = "5prime_UTR"
                                                codon_edited = "5prime_UTR"
                                                aa = "5prime_UTR"
                                                aa_edited = "5prime_UTR"
                                                consequenc = 'utr'
                                            elif row["Strand"] == '-' and exon_number == last_transcript_exon:
                                                codon = "3prime_UTR"
                                                codon_edited = "3prime_UTR"
                                                aa = "3prime_UTR"
                                                aa_edited = "3prime_UTR"
                                                consequenc = 'utr'

                                        elif (row["End"] + 2) <= pos:
                                            if row["Strand"] == '-' and exon_number == first_transcript_exon:
                                                codon = "5prime_UTR"
                                                codon_edited = "5prime_UTR"
                                                aa = "5prime_UTR"
                                                aa_edited = "5prime_UTR"
                                                consequenc = 'utr'
                                            elif row["Strand"] == '+' and exon_number == last_transcript_exon:
                                                codon = "3prime_UTR"
                                                codon_edited = "3prime_UTR"
                                                aa = "3prime_UTR"
                                                aa_edited = "3prime_UTR"
                                                consequenc = 'utr'

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        consequence.append(consequenc)

                                        aa_positions.append(aa)
                                        synonymous = False

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
                            distance_to_center_variant_unused, \
                            distance_to_center = shared.analyze_guide(possible_guide,
                                                                      edit_window_start,
                                                                      edit_window_end,
                                                                      edit_window_start_plus,
                                                                      edit_window_end_plus,
                                                                      bes[be]['fwd']['REF'],
                                                                      None,
                                                                      distance_median_dict,
                                                                      fiveprimepam)
                            # the following lines have been removed from the statement above:
                            # quality_scores_variant_unused, \
                            # distance_to_center, \
                            # quality_scores_all = shared.analyze_guide(possible_guide,
                                                                    #   distance_median_dict,
                                                                    #   quality_scores_dict) # _unused variables just for testing purposes; can be removed later

                            # if  (edit_window_bases_unused != edit_window_bases or
                            #     edit_window_plus_bases_unused != edit_window_plus_bases): # just for testing purposes; can be removed later
                            #     raise Exception("output of analyze_guide() wrong in forward search")

                            # guide_is_missense = False
                            # guide_is_synonymous = False
                            # guide_is_stopgain = False
                            # guide_is_startlost = False
                            # guide_is_stoplost = False

                            # for i, aa in enumerate(aas):
                            #     aa_edited = aas_edited[i]
                            #     if aa in non_stop_aas and aa_edited in non_stop_aas and aa != aa_edited:
                            #         guide_is_missense = True
                            #         # consequence.append('missense')
                            #     elif aa in non_stop_aas and aa_edited in non_stop_aas and aa == aa_edited:
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     elif aa != 'Stop' and aa_edited == 'Stop':
                            #         guide_is_stopgain = True
                            #         # consequence.append('stopgain')
                            #     elif aa == 'Stop' and aa_edited != 'Stop':
                            #         guide_is_stoplost = True
                            #         # consequence.append('stoplost')
                            #     elif aa == 'Stop' and aa_edited == 'Stop': # can this happen?
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     elif aa == 'StartM' and aa_edited !='M': # StartM not used in aas_edited yet; Any edit in StartM would be a startlost
                            #         guide_is_startlost = True
                            #         # consequence.append('startlost')
                            #     elif aa == 'StartM' and aa_edited =='M': # can't happen
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     # else: # should not happen; only necessary, if using consequence list
                            #     #     consequence.append('complex')

                            # consequence = []
                            # if synonymous:
                            #     consequence.append('synonymous')
                            # if includes_splice_site:
                            #     consequence.append('splice_site')
                            # if specific:
                            #     consequence.append('specific')
                            # if guide_is_missense:
                            #     consequence.append('missense')
                            # if guide_is_stopgain:
                            #     consequence.append('stopgain')
                            # if guide_is_stoplost:
                            #     consequence.append('stoplost')
                            # if guide_is_startlost:
                            #     consequence.append('startlost')
                            # if not consequence:
                            #     consequence.append('complex')

                            # if not ((filter_synonymous) and (not synonymous) or
                            #         (filter_splice_site) and (not includes_splice_site) or
                            #         (filter_specific) and (not specific) or
                            #         (filter_missense) and (not guide_is_missense) or
                            #         (filter_nonsense) and (not guide_is_stopgain) or
                            #         (filter_stoplost) and (not guide_is_stoplost) or
                            #         (filter_startlost) and (not guide_is_startlost)):

                            if fiveprimepam:
                                possible_guide = shared.revcom(possible_guide)
                                possible_guide_with_pam = shared.revcom(possible_guide_with_pam)
                                edit_window_bases = shared.revcom(edit_window_bases)
                                edit_window_plus_bases = shared.revcom(edit_window_plus_bases)
                                safety_region = shared.revcom(safety_region)
                                edit_string = edit_string[::-1]
                                edit_pos_string = edit_pos_string[::-1]
                                distance_to_center = distance_to_center[::-1]

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
                            consequence_output.append(consequence)
                            direction_output.append('+' if not fiveprimepam else '-') # forward
                            codon_ref.append(codons) # list
                            aa_ref.append(aas) # list
                            aa_pos.append(aa_positions)
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
                            distance_to_center_output.append(distance_to_center)
                            # quality_scores_all_output.append(quality_scores_all)
                            guide_starts_output.append(str(startpos_fwd + 1))
                            guide_ends_output.append(str(startpos_fwd + guidelength + 1))
                            chroms_output.append(str(chrom))

                    else:
                        if bes[be]['fwd']['REF'] in edit_window_plus_bases:
                            ne_plus = False
                        else:
                            ne_plus = True

                        if fiveprimepam:
                            possible_guide = shared.revcom(possible_guide)
                            possible_guide_with_pam = shared.revcom(possible_guide_with_pam)
                            edit_window_bases = shared.revcom(edit_window_bases)
                            edit_window_plus_bases = shared.revcom(edit_window_plus_bases)
                            safety_region = shared.revcom(safety_region)

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
                        direction_output_ne.append('+' if not fiveprimepam else '-') # forward
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

                if any(possible_guide_with_pam.startswith(shared.revcom(pam_real_string)) for pam_real_string in pamlist_real_string): # shared.revcom(PAM) found at 3' end

                    if bes[be]['rev']['REF'] in edit_window_bases: # shared.revcom(editable base) in edit window

                        # find all indices for possible edits in the edit window in search string
                        indices_edit = [i + (guidelength + len(pamsite) - edit_window_end) + m.start() for m in re.finditer(bes[be]['rev']['REF'], edit_window_bases)]
                        indices_plus_edit = [i + (guidelength + len(pamsite) - edit_window_end + edit_window_end_plus) + m.start() for m in re.finditer(bes[be]['rev']['REF'], edit_window_plus_bases)] # is this correct?

                        variants = [] # reset variant list for sliding window
                        codons = [] # reset variant list for sliding window
                        codons_edited = [] # reset variant list for sliding window
                        aas = [] # reset variant list for sliding window
                        aas_edited = [] # reset variant list for sliding window
                        consequence = [] # reset variant list for sliding window
                        aa_positions = [] # reset variant list for sliding window
                        # includes_splice_site = "not_tested" # also check for splice-sites even, if not explicitely added
                        # if splice_sites:
                        #     includes_splice_site = False
                        includes_splice_site = False
                        synonymous = False

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
                                            try:
                                                codon_list = []
                                                for ntpos in range(pos - offset, pos + 3 - offset):
                                                    if ntpos < row["Start"]:
                                                        new_exon = cdss_gene.filter(
                                                            (pl.col('gene_name') == gene_symbol) &
                                                            (pl.col('transcript_name') == transcript_symbol) &
                                                            (pl.col('exon_number') == exon_number - 1) &
                                                            (pl.col('first_transcript_exon') == first_transcript_exon) &
                                                            (pl.col('last_transcript_exon') == last_transcript_exon) &
                                                            (pl.col('Chromosome') == chrom)
                                                            )
                                                        if not len(new_exon) == 0:
                                                            raise AssertionError('More then one previous exon!')
                                                        new_exon = new_exon.to_dict(as_series=False)
                                                        new_ntpos = new_exon['End'][0] - (row["Start"] - ntpos)
                                                        codon_list.append(str(ref_genome_pyfaidx[new_exon['Chromosome'][0]][new_ntpos]))
                                                    elif row["Start"] <= ntpos < row["End"]:
                                                        codon_list.append(str(ref_genome_pyfaidx[row["Chromosome"]][ntpos]))
                                                    elif row["End"] <= ntpos:
                                                        new_exon = cdss_gene.filter(
                                                            (pl.col('gene_name') == gene_symbol) &
                                                            (pl.col('transcript_name') == transcript_symbol) &
                                                            (pl.col('exon_number') == exon_number + 1) &
                                                            (pl.col('first_transcript_exon') == first_transcript_exon) &
                                                            (pl.col('last_transcript_exon') == last_transcript_exon) &
                                                            (pl.col('Chromosome') == chrom)
                                                            )
                                                        if not len(new_exon) == 0:
                                                            raise AssertionError('More then one next exon!')
                                                        new_exon = new_exon.to_dict(as_series=False)
                                                        new_ntpos = new_exon['Start'][0] + (ntpos - row["End"])
                                                        codon_list.append(str(ref_genome_pyfaidx[new_exon['Chromosome'][0]][new_ntpos]))
                                                codon = ''.join(codon_list)
                                                codon_edited = shared.replace_str_index(codon, offset, bes[be]['rev']['ALT']) # specific for search direction

                                                if str(row['Strand']) == '-':
                                                    codon = shared.revcom(codon)
                                                    codon_edited = shared.revcom(codon_edited)
                                                aa = shared.codon_sun_one_letter[codon]
                                                aa_edited = shared.codon_sun_one_letter[codon_edited]

                                                if aa == shared.codon_sun_one_letter["ATG"] and exon_number == first_transcript_exon:
                                                    if (row['Strand'] == '+' and (row["Start"] <= pos < (row["Start"] + 3))) or (row['Strand'] == '-' and ((row["End"] - 3) <= pos < row["End"])):
                                                        codon = "Start" + codon
                                                        aa = "Start" + aa
                                            except:
                                                codon = "incomplete_codon"
                                                codon_edited = "incomplete_codon"
                                                aa = "codon_incomplete"
                                                aa_edited = "codon_incomplete"

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        if aa in non_stop_aas and aa_edited in non_stop_aas and aa != aa_edited:
                                            consequenc = 'missense'
                                        elif aa in non_stop_aas and aa_edited in non_stop_aas and aa == aa_edited:
                                            consequenc = 'synonymous'

                                        elif aa != 'Stop' and aa_edited == 'Stop':
                                            consequenc = 'stopgain'
                                        elif aa == 'Stop' and aa_edited != 'Stop':
                                            consequenc = 'stoplost'
                                        elif aa == 'Stop' and aa_edited == 'Stop': # can this happen?
                                            consequenc = 'synonymous'

                                        elif aa == 'StartM' and aa_edited !='M': # StartM not used in aas_edited yet; Any edit in StartM would be a startlost
                                            consequenc = 'startlost'
                                        elif aa == 'StartM' and aa_edited =='M': # can't happen
                                            consequenc = 'synonymous'

                                        elif aa == "codon_incomplete" and aa_edited == "codon_incomplete":
                                            consequenc = 'indefinite'

                                        else: # can this happen?
                                            consequenc = 'complex'

                                        consequence.append(consequenc)

                                        if row['Strand'] == "+":
                                            nt_position = row['transcript_length_before'] + ((pos + 1) - row['Start'])
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                        elif row['Strand'] == "-":
                                            nt_position = row['transcript_length_before'] + (row['End'] - (pos)) # (row['End'] - (pos + 1))?
                                            aa_position = -(-nt_position // 3)
                                            aa_position = str(aa_position)

                                        aa_positions.append(aa_position)

                                        if len(indices_edit) == 1: # only 1 editable base in editing window
                                            # specific = True # will be analyzed by analyze_guide()
                                            # if len(indices_plus_edit) == 1:
                                            #     specific_plus = True
                                            # else:
                                            #     specific_plus = False
                                            if aa == aa_edited and aa != "codon_incomplete":
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
                                            consequenc = 'splice_site'

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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'
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
                                                        consequenc = 'utr'

                                            codons.append(codon)
                                            codons_edited.append(codon_edited)
                                            aas.append(aa)
                                            aas_edited.append(aa_edited)

                                            consequence.append(consequenc)

                                            aa_positions.append(aa)

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
                                        # codons.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # codons_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # aas.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # aas_edited.append('ANNOTATIONWARNING') # why does this happen? maybe different exon lenghts? annotations issues?
                                        # this happens due to off-target effects in introns or the UTRs
                                        codon = "intron"
                                        codon_edited = "intron"
                                        aa = "intron"
                                        aa_edited = "intron"
                                        consequenc = "intron"

                                        if pos < (row["Start"] - 2):
                                            if row["Strand"] == '+' and exon_number == first_transcript_exon:
                                                codon = "5prime_UTR"
                                                codon_edited = "5prime_UTR"
                                                aa = "5prime_UTR"
                                                aa_edited = "5prime_UTR"
                                                consequenc = "utr"
                                            elif row["Strand"] == '-' and exon_number == last_transcript_exon:
                                                codon = "3prime_UTR"
                                                codon_edited = "3prime_UTR"
                                                aa = "3prime_UTR"
                                                aa_edited = "3prime_UTR"
                                                consequenc = "utr"

                                        elif (row["End"] + 2) <= pos:
                                            if row["Strand"] == '-' and exon_number == first_transcript_exon:
                                                codon = "5prime_UTR"
                                                codon_edited = "5prime_UTR"
                                                aa = "5prime_UTR"
                                                aa_edited = "5prime_UTR"
                                                consequenc = "utr"
                                            elif row["Strand"] == '+' and exon_number == last_transcript_exon:
                                                codon = "3prime_UTR"
                                                codon_edited = "3prime_UTR"
                                                aa = "3prime_UTR"
                                                aa_edited = "3prime_UTR"
                                                consequenc = "utr"

                                        codons.append(codon)
                                        codons_edited.append(codon_edited)
                                        aas.append(aa)
                                        aas_edited.append(aa_edited)

                                        consequence.append(consequenc)

                                        aa_positions.append(aa)
                                        synonymous = False

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
                            distance_to_center_variant_unused, \
                            distance_to_center = shared.analyze_guide(shared.revcom(possible_guide),
                                                                      edit_window_start,
                                                                      edit_window_end,
                                                                      edit_window_start_plus,
                                                                      edit_window_end_plus,
                                                                      bes[be]['fwd']['REF'],
                                                                      None,
                                                                      distance_median_dict,
                                                                      fiveprimepam)
                            # the following lines have been removed from the statement above:
                            # quality_scores_variant_unused, \
                            # distance_to_center, \
                            # quality_scores_all = shared.analyze_guide(shared.revcom(possible_guide),
                                                                    #   distance_median_dict,
                                                                    #   quality_scores_dict) # _unused variables just for testing purposes; can be removed later

                            # edit_window_bases_unused = shared.revcom(edit_window_bases_unused)
                            # edit_window_plus_bases_unused = shared.revcom(edit_window_plus_bases_unused)

                            # if  (edit_window_bases_unused != edit_window_bases or
                            #     edit_window_plus_bases_unused != edit_window_plus_bases): # just for testing purposes; can be removed later
                            #     raise Exception("output of analyze_guide() wrong in reverse search")

                            # guide_is_missense = False
                            # guide_is_synonymous = False
                            # guide_is_stopgain = False
                            # guide_is_startlost = False
                            # guide_is_stoplost = False

                            # for i, aa in enumerate(aas):
                            #     aa_edited = aas_edited[i]
                            #     if aa in non_stop_aas and aa_edited in non_stop_aas and aa != aa_edited:
                            #         guide_is_missense = True
                            #         # consequence.append('missense')
                            #     elif aa in non_stop_aas and aa_edited in non_stop_aas and aa == aa_edited:
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     elif aa != 'Stop' and aa_edited == 'Stop':
                            #         guide_is_stopgain = True
                            #         # consequence.append('stopgain')
                            #     elif aa == 'Stop' and aa_edited != 'Stop':
                            #         guide_is_stoplost = True
                            #         # consequence.append('stoplost')
                            #     elif aa == 'Stop' and aa_edited == 'Stop': # can this happen?
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     elif aa == 'StartM' and aa_edited !='M': # StartM not used in aas_edited yet; Any edit in StartM would be a startlost
                            #         guide_is_startlost = True
                            #         # consequence.append('startlost')
                            #     elif aa == 'StartM' and aa_edited =='M': # can't happen
                            #         guide_is_synonymous = True
                            #         # consequence.append('synonymous')

                            #     # else: # should not happen; only necessary, if using consequence list
                            #     #     consequence.append('complex')

                            # consequence = []
                            # if synonymous:
                            #     consequence.append('synonymous')
                            # if includes_splice_site:
                            #     consequence.append('splice_site')
                            # if specific:
                            #     consequence.append('specific')
                            # if guide_is_missense:
                            #     consequence.append('missense')
                            # if guide_is_stopgain:
                            #     consequence.append('stopgain')
                            # if guide_is_stoplost:
                            #     consequence.append('stoplost')
                            # if guide_is_startlost:
                            #     consequence.append('startlost')
                            # if not consequence:
                            #     consequence.append('complex')

                            # if not ((filter_synonymous) and (not synonymous) or
                            #         (filter_splice_site) and (not includes_splice_site) or
                            #         (filter_specific) and (not specific) or
                            #         (filter_missense) and (not guide_is_missense) or
                            #         (filter_nonsense) and (not guide_is_stopgain) or
                            #         (filter_stoplost) and (not guide_is_stoplost) or
                            #         (filter_startlost) and (not guide_is_startlost)):

                            if fiveprimepam:
                                possible_guide = shared.revcom(possible_guide)
                                possible_guide_with_pam = shared.revcom(possible_guide_with_pam)
                                edit_window_bases = shared.revcom(edit_window_bases)
                                edit_window_plus_bases = shared.revcom(edit_window_plus_bases)
                                safety_region = shared.revcom(safety_region)
                                edit_string = edit_string[::-1]
                                edit_pos_string = edit_pos_string[::-1]
                                distance_to_center = distance_to_center[::-1]

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
                            consequence_output.append(consequence)
                            direction_output.append('-' if not fiveprimepam else '+') # reverse
                            codon_ref.append(codons) # list
                            aa_ref.append(aas) # list
                            aa_pos.append(aa_positions)
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
                            distance_to_center_output.append(distance_to_center)
                            # quality_scores_all_output.append(quality_scores_all)
                            guide_starts_output.append(str(startpos_rev + len(pamsite) + 1))
                            guide_ends_output.append(str(startpos_rev + guidelength + len(pamsite) + 1))
                            chroms_output.append(str(chrom))

                    else:
                        if bes[be]['rev']['REF'] in edit_window_plus_bases:
                            ne_plus = False
                        else:
                            ne_plus = True

                        if fiveprimepam:
                            possible_guide = shared.revcom(possible_guide)
                            possible_guide_with_pam = shared.revcom(possible_guide_with_pam)
                            edit_window_bases = shared.revcom(edit_window_bases)
                            edit_window_plus_bases = shared.revcom(edit_window_plus_bases)
                            safety_region = shared.revcom(safety_region)

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
                        direction_output_ne.append('-' if not fiveprimepam else '+') # reverse
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

    aa_pos_to_modify_very_first = ['aa_pos']

    transcript_cols_to_modify_first = ['exon_number',
                                    'transcript',
                                    'first_transcript_exon',
                                    'last_transcript_exon']

    variant_cols_to_modify_first = ['codon_ref',
                                    'aa_ref',
                                    'codon_edit',
                                    'aa_edit']

    variant_cols_to_modify_second = ['variant']

    transcript_cols_to_modify_second = ['synonymous_specific',
                                        'codon_ref',
                                        'aa_ref',
                                        'codon_edit',
                                        'aa_edit',
                                        'splice_site_included',
                                        'exon_number',
                                        'transcript',
                                        'first_transcript_exon',
                                        'last_transcript_exon']

    consequences_to_modify = ['consequence']

    # editing guides found
    if guide_output: # if no guide is found at all (very unlikely)
        sgrnas = pl.DataFrame({"variant": variant_output,
                            "base_change": be_output,
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
                            "synonymous_specific": synonymous_output,
                            "consequence": consequence_output,
                            "strand": direction_output,
                            "codon_ref": codon_ref,
                            "aa_ref": aa_ref,
                            "aa_pos": aa_pos,
                            "codon_edit": codon_edit,
                            "aa_edit": aa_edit,
                            "splice_site_included": splice_site_included,
                            "originally_intended_ALT": "NA_for_genes",
                            "ref_match": "NA_for_genes",
                            "off_target_bases": edit_string_output,
                            "edited_positions": edit_pos_string_output,
                            "specificity":  "NA_for_genes",
                            "distance_to_center_variant":  "NA_for_genes",
                            # "efficiency_scores_variant":  "NA_for_genes",
                            "distance_to_center": distance_to_center_output,
                            # "efficiency_scores_all": quality_scores_all_output,
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
        sgrnas = sgrnas.with_columns(synonymous_specific = sgrnas['synonymous_specific'].cast(str))
        sgrnas = sgrnas.with_columns(splice_site_included = sgrnas['splice_site_included'].cast(str))
        sgrnas = sgrnas.with_columns(variant = sgrnas['variant'])
        sgrnas = sgrnas.with_columns(codon_ref = sgrnas['codon_ref'])
        sgrnas = sgrnas.with_columns(aa_ref = sgrnas['aa_ref'])
        sgrnas = sgrnas.with_columns(aa_pos = sgrnas['aa_pos'])
        sgrnas = sgrnas.with_columns(codon_edit = sgrnas['codon_edit'])
        sgrnas = sgrnas.with_columns(aa_edit = sgrnas['aa_edit'])
        # sgrnas = sgrnas.with_columns(variant = sgrnas['variant'].list.join('&'))
        # sgrnas = sgrnas.with_columns(codon_ref = sgrnas['codon_ref'].list.join('&'))
        # sgrnas = sgrnas.with_columns(aa_ref = sgrnas['aa_ref'].list.join('&'))
        # sgrnas = sgrnas.with_columns(codon_edit = sgrnas['codon_edit'].list.join('&'))
        # sgrnas = sgrnas.with_columns(aa_edit = sgrnas['aa_edit'].list.join('&'))
        sgrnas = sgrnas.group_by([col for col in sgrnas.columns if col not in ['aa_pos',
                                                                               'exon_number',
                                                                            'transcript',
                                                                            'first_transcript_exon',
                                                                            'last_transcript_exon']],
                                maintain_order=True).agg(pl.all())
        sgrnas = sgrnas.group_by([col for col in sgrnas.columns if col not in ['synonymous_specific',
                                                                            'consequence',
                                                                            'codon_ref',
                                                                            'aa_ref',
                                                                            'aa_pos',
                                                                            'codon_edit',
                                                                            'aa_edit',
                                                                            'splice_site_included',
                                                                            'exon_number',
                                                                            'transcript',
                                                                            'first_transcript_exon',
                                                                            'last_transcript_exon']],
                                maintain_order=True).agg(pl.all())

        sgrnas = sgrnas[sgrnas_cols]

        if blast:

            blast_guides.check_blastdb(ref_genome, False)

            sgrnas = sgrnas.with_row_index('index')
            guides_blast = sgrnas.select('index', 'guide')
            blast_results = blast_guides.guide_blast(guides_blast,
                                                    guidelength,
                                                    ref_genome,
                                                    'genes',
                                                    no_contigs)
            sgrnas = sgrnas.join(blast_results, on='index', how='left')
            sgrnas = sgrnas.drop('index')

        # add vep annotations, if wanted
        if vep:

            variants_vep = sgrnas.select('variant').with_row_index('original_index').explode('variant')
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

            variants_vep_resorted = variants_vep_resorted.group_by('original_index', maintain_order=True).agg(pl.all())

            variants_vep_resorted = variants_vep_resorted.select(pl.exclude(["original_index", "variant"]))
            sgrnas = sgrnas.with_columns(variants_vep_resorted)

            variant_cols_to_modify_second += variants_vep_resorted.columns # if vep for sgrnas_ne put into use, needs to be integrated into sgrnas_ne block

        if filter_synonymous:
            sgrnas = sgrnas.filter(~(pl.col('synonymous_specific').list.join("^").str.contains("false") | pl.col('synonymous_specific').list.join("^").str.contains("False"))) # capitalization inconsistent
        if filter_splice_site:
            sgrnas = sgrnas.filter(pl.col('consequence').list.eval(pl.element().list.join(";")).list.join("^").str.contains("splice_site"))
            # sgrnas = sgrnas.filter(pl.col('splice_site_included').list.join("^").str.contains("true"))
        if filter_specific:
            sgrnas = sgrnas.filter(pl.col('specific'))
        if filter_missense:
            sgrnas = sgrnas.filter(pl.col('consequence').list.eval(pl.element().list.join(";")).list.join("^").str.contains("missense"))
        if filter_nonsense:
            sgrnas = sgrnas.filter(pl.col('consequence').list.eval(pl.element().list.join(";")).list.join("^").str.contains("stopgain"))
        if filter_stoplost:
            sgrnas = sgrnas.filter(pl.col('consequence').list.eval(pl.element().list.join(";")).list.join("^").str.contains("stoplost"))
        if filter_startlost:
            sgrnas = sgrnas.filter(pl.col('consequence').list.eval(pl.element().list.join(";")).list.join("^").str.contains("startlost"))

        if aspect == 'exploded': # maybe also show off target edits in codons in small letters as for bedesigner

            sgrnas = sgrnas.with_columns(
                pl.col(aa_pos_to_modify_very_first).map_batches(lambda s:
                    s.to_frame()
                    .with_row_index("outer_row")
                    .explode(pl.last())
                    .with_row_index("row")
                    .explode(pl.last())
                    .with_row_index()
                    .explode(pl.last())
                    .with_columns(pl.int_ranges(pl.col("index").rle().struct.field("len")).flatten().alias("index"))
                    .group_by("index", "row", "outer_row", maintain_order=True)
                    .agg(pl.last())
                    .group_by("row", "outer_row", maintain_order=True)
                    .agg(pl.last())
                    .group_by("outer_row", maintain_order=True)
                    .agg(pl.last())
                    .select(pl.last())
                    .to_series()))

            sgrnas = sgrnas.with_columns(
                pl.col(aa_pos_to_modify_very_first).list.eval(pl.element().list.eval(pl.element().list.join("~"))))

            sgrnas = sgrnas.with_columns(
                pl.col(transcript_cols_to_modify_first).list.eval(pl.element().list.join("~")).list.join("^"),
                pl.col(variant_cols_to_modify_first + aa_pos_to_modify_very_first + consequences_to_modify).map_batches(lambda s:
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
                pl.col([col for col in transcript_cols_to_modify_second if col not in transcript_cols_to_modify_first + variant_cols_to_modify_first]).list.join("^"),
            ).explode(variant_cols_to_modify_second + variant_cols_to_modify_first + aa_pos_to_modify_very_first + consequences_to_modify)

        elif aspect == 'collapsed':
            sgrnas = sgrnas.with_columns(
                pl.col(aa_pos_to_modify_very_first).list.eval(pl.element().list.eval(pl.element().list.join(";")).list.join("~")).list.join("^"),
                pl.col(transcript_cols_to_modify_first).list.eval(pl.element().list.join("~")).list.join("^"),
                pl.col(variant_cols_to_modify_first).list.eval(pl.element().list.join(";")).list.join("^"),
                pl.col(variant_cols_to_modify_second).list.join(";"),
                pl.col([col for col in transcript_cols_to_modify_second if col not in transcript_cols_to_modify_first + variant_cols_to_modify_first]).list.join("^"),
                pl.col(consequences_to_modify).list.eval(pl.element().list.join(";")).list.join("^")
            )

        sgrnas = sgrnas.unique().sort(by=sgrnas.columns)

        sam_df = sgrnas.select(
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

        sgrnas = sgrnas.drop(['ne_plus',
                            'originally_intended_ALT',
                            'ref_match',
                            'off_target_bases',
                            'specificity',
                            'distance_to_center_variant'])
                            # 'distance_to_center_variant',
                            # 'efficiency_scores_variant'])

        if edit_window_start_plus == 0 and edit_window_end_plus == 0:
            sgrnas = sgrnas.drop(['edit_window_plus',
                                'num_edits_plus',
                                'specific_plus',
                                'safety_region',
                                'num_edits_safety',
                                'additional_in_safety'])

    else:
        sgrnas = pl.DataFrame()
        sam_df = pl.DataFrame()

    # non-editing guides found
    if guide_output_ne:
        sgrnas_ne = pl.DataFrame({"variant":  "NA_for_non_editing_guides",
                                "base_change": be_output_ne,
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
                                "synonymous_specific": "NA_for_non_editing_guides",
                                "consequence": "NA_for_non_editing_guides",
                                "strand": direction_output_ne,
                                "codon_ref": "NA_for_non_editing_guides",
                                "aa_ref": "NA_for_non_editing_guides",
                                "aa_pos": "NA_for_non_editing_guides",
                                "codon_edit": "NA_for_non_editing_guides",
                                "aa_edit": "NA_for_non_editing_guides",
                                "splice_site_included": "NA_for_non_editing_guides",
                                "originally_intended_ALT": "NA_for_genes",
                                "ref_match": "NA_for_genes",
                                "off_target_bases": "NA_for_non_editing_guides",
                                "edited_positions": "NA_for_non_editing_guides",
                                "specificity": "NA_for_genes",
                                "distance_to_center_variant": "NA_for_genes",
                                # "efficiency_scores_variant": "NA_for_genes",
                                "distance_to_center": "NA_for_non_editing_guides",
                                # "efficiency_scores_all": "NA_for_non_editing_guides",
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
        sgrnas_ne = sgrnas_ne.with_columns(synonymous_specific = sgrnas_ne['synonymous_specific'].cast(str))
        sgrnas_ne = sgrnas_ne.with_columns(splice_site_included = sgrnas_ne['splice_site_included'].cast(str))
        sgrnas_ne = sgrnas_ne.group_by([col for col in sgrnas_ne.columns if col not in ['exon_number', # 'aa_pos' shouldn't be necessary here, since it's not set
                                                                                    'transcript',
                                                                                    'first_transcript_exon',
                                                                                    'last_transcript_exon']],
                                        maintain_order=True).agg(pl.all().str.join("~"))
        sgrnas_ne = sgrnas_ne.group_by([col for col in sgrnas_ne.columns if col not in ['synonymous_specific', # 'aa_pos' shouldn't be necessary here, since it's not set
                                                                                    'consequence',
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

            sgrnas_ne = sgrnas_ne.with_row_index('index')
            guides_ne_blast = sgrnas_ne.select('index', 'guide')
            blast_results_ne = blast_guides.guide_blast(guides_ne_blast,
                                                        guidelength,
                                                        ref_genome,
                                                        'genes',
                                                        no_contigs)
            sgrnas_ne = sgrnas_ne.join(blast_results_ne, on='index', how='left')
            sgrnas_ne =sgrnas_ne.drop('index')

        sgrnas_ne = sgrnas_ne.unique().sort(by=sgrnas_ne.columns)

        sam_ne_df = sgrnas_ne.select(
            ['variant', 'guide', 'guide_chrom', 'guide_start', 'strand']
            ) # will always be exploded, but maybe unify somehow different?

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

        sgrnas_ne = sgrnas_ne.drop(['variant',
                                    'num_edits',
                                    'specific',
                                    'num_edits_plus',
                                    'specific_plus',
                                    'num_edits_safety',
                                    'additional_in_safety',
                                    'synonymous_specific',
                                    'consequence',
                                    'codon_ref',
                                    'aa_ref',
                                    'aa_pos',
                                    'codon_edit',
                                    'aa_edit',
                                    'splice_site_included',
                                    'originally_intended_ALT',
                                    'ref_match',
                                    'off_target_bases',
                                    'edited_positions',
                                    'specificity',
                                    'distance_to_center_variant',
                                    'distance_to_center'])
                                    # 'efficiency_scores_variant',
                                    # 'distance_to_center',
                                    # 'efficiency_scores_all'])

        if edit_window_start_plus == 0 and edit_window_end_plus == 0:
            sgrnas_ne = sgrnas_ne.drop(['edit_window_plus',
                                        'safety_region',
                                        'ne_plus'])

    else:
        sgrnas_ne = pl.DataFrame()
        sam_ne_df = pl.DataFrame()

    return (sgrnas, sgrnas_ne, sam_df, sam_ne_df, genes_not_found)


def output_sgrnas(sgrnas, sgrnas_ne, output_file):
    if not sgrnas.is_empty():
        sgrnas.write_csv(output_file + ".tsv", separator='\t')
    else:
        print("No editing guides could be identified: No such TSV file will be written.")

    if not sgrnas_ne.is_empty():
        sgrnas_ne.write_csv(output_file + "_ne.tsv", separator='\t')
    else:
        print("No non-editing guides could be identified: No such TSV file will be written.")


def output_guides_sam(sam_df, sam_ne_df, output_file, ref_genome):

    if not sam_df.is_empty():
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
    else:
        print("No editing guides could be identified: No such BAM file will be written.")

    if not sam_ne_df.is_empty():
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
    else:
        print("No non-editing guides could be identified: No such BAM file will be written.")


if __name__ == "__main__":
    pass
