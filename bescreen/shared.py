import statistics
# from scipy.stats import norm
import os
import polars as pl
import pandas as pd
import numpy as np
import pyranges as pr


def qc_precalc(edit_window_start,
               edit_window_end):

    edit_window_poss = range(edit_window_start, edit_window_end+1)
    median_edit_window = statistics.median(edit_window_poss)
    # max_distance_edit_window = max(median_edit_window-edit_window_start, edit_window_end-median_edit_window) # should always be symetrical

    # twosd_quality = max_distance_edit_window
    # onesd_quality = twosd_quality / 2

    # quality_dict for edit window
    distance_median_dict = dict(zip(edit_window_poss, [abs(possible_edit_pos-median_edit_window) for possible_edit_pos in edit_window_poss]))
    # quality_scores_dict = dict(zip(edit_window_poss, norm.pdf(edit_window_poss, median_edit_window, onesd_quality) * (1 / norm.pdf(median_edit_window, median_edit_window, onesd_quality))))

    # return (distance_median_dict, quality_scores_dict)
    return distance_median_dict


def analyze_guide(guide,
                  edit_window_start,
                  edit_window_end,
                  edit_window_plus_start,
                  edit_window_plus_end,
                  editable_base,
                  position_variant,
                  distance_median_dict,
                  fiveprimepam):
                #   distance_median_dict,
                #   quality_scores_dict):

    edit_window = guide[edit_window_start-1:edit_window_end]
    num_edits = edit_window.count(editable_base)
    specific = True if num_edits == 1 else False # find a way to add and edit_window.index(editable_base) == position_variant-edit_window_start, but only if position_variant is set (bedesigner) and not None (besaturate)

    edit_window_plus = guide[edit_window_start-edit_window_plus_start-1:edit_window_end+edit_window_plus_end]
    num_edits_plus = edit_window_plus.count(editable_base)
    specific_plus = True if num_edits_plus == 1 else False # find a way to add and edit_window_plus.index(editable_base) == position_variant-edit_window_start-edit_window_plus_start, but only if position_variant is set (bedesigner) and not None (besaturate)

    safety_region = guide[edit_window_start-edit_window_plus_start-1:edit_window_start-1] + guide[edit_window_start-1:edit_window_end].lower() + guide[edit_window_end:edit_window_end+edit_window_plus_end] # lowercase letters will never be an editable base
    num_edits_safety = safety_region.count(editable_base)
    additional_in_safety = True if num_edits_safety == 1 else False # check the two comments above

    string_edits = "." * len(guide)
    string_edits = list(string_edits) # is this necessary?
    string_edit_poss = "." * len(guide)
    string_edit_poss = list(string_edit_poss) # is this necessary?

    variant_edits = 0
    all_edits = 0
    variant_edits_pos = []
    all_edits_pos = []

    for j in range(edit_window_start-edit_window_plus_start-1, edit_window_start-1):
        string_edits[j] = ","
        string_edit_poss[j] = ","
        if guide[j] == editable_base:
            string_edits[j] = "+"
            string_edit_poss[j] = "+"

    for j in range(edit_window_start-1, edit_window_end):
        string_edits[j] = ":"
        string_edit_poss[j] = ":"
        if guide[j] == editable_base and j == position_variant:
            string_edits[j] = "V"
            string_edit_poss[j] = str(j+1) if not fiveprimepam else str(len(guide) - j)[::-1]
            variant_edits += 1
            all_edits += 1
            variant_edits_pos.append(j+1)
            all_edits_pos.append(j+1)
        elif guide[j] == editable_base:
            string_edits[j] = "*"
            string_edit_poss[j] = str(j+1) if not fiveprimepam else str(len(guide) - j)[::-1]
            all_edits += 1
            all_edits_pos.append(j+1)

    for j in range(edit_window_end, edit_window_end+edit_window_plus_end):
        string_edits[j] = ","
        string_edit_poss[j] = ","
        if guide[j] == editable_base:
            string_edits[j] = "+"
            string_edit_poss[j] = "+"

    string_edits = "".join(string_edits) # can be removed, if not set to list
    string_edit_poss = "".join(string_edit_poss) # can be removed, if not set to list

    # some QC
    specificty = round(variant_edits / all_edits, 2)

    # median_edit_window = statistics.median(range(edit_window_start, edit_window_end+1))
    # max_distance_edit_window = max(median_edit_window-edit_window_start, edit_window_end-median_edit_window) # should always be symetrical

    # twosd_quality = max_distance_edit_window
    # onesd_quality = twosd_quality / 2

    # # quality for variant
    # distance_median_variant = [abs(edit_pos-median_edit_window) for edit_pos in variant_edits_pos]
    distance_median_variant = [distance_median_dict[edit_pos] for edit_pos in variant_edits_pos]
    # quality_scores_variant = norm.pdf(variant_edits_pos, median_edit_window, onesd_quality) * (1 / norm.pdf(median_edit_window, median_edit_window, onesd_quality))
    # quality_scores_variant = [quality_scores_dict[edit_pos] for edit_pos in variant_edits_pos]

    # # quality for all edits
    # distance_median_all = [abs(edit_pos-median_edit_window) for edit_pos in all_edits_pos]
    distance_median_all = [distance_median_dict[edit_pos] for edit_pos in all_edits_pos]
    # quality_scores_all = norm.pdf(all_edits_pos, median_edit_window, onesd_quality) * (1 / norm.pdf(median_edit_window, median_edit_window, onesd_quality))
    # quality_scores_all = [quality_scores_dict[edit_pos] for edit_pos in all_edits_pos]

    distance_median_variant = [str(number) if not fiveprimepam else str(number)[::-1] for number in distance_median_variant]
    # quality_scores_variant = [str(round(float(number), 2)) for number in quality_scores_variant]
    distance_median_all = [str(number) if not fiveprimepam else str(number)[::-1] for number in distance_median_all]
    # quality_scores_all = [str(round(float(number), 2)) for number in quality_scores_all]

    distance_median_variant = ",".join(distance_median_variant)
    # quality_scores_variant = ",".join(quality_scores_variant)
    distance_median_all = ",".join(distance_median_all)
    # quality_scores_all = ",".join(quality_scores_all)

    return (edit_window,
            num_edits,
            specific,
            edit_window_plus,
            num_edits_plus,
            specific_plus,
            safety_region,
            num_edits_safety,
            additional_in_safety,
            string_edits,
            string_edit_poss,
            specificty,
            distance_median_variant,
            # quality_scores_variant,
            distance_median_all)
            # distance_median_all,
            # quality_scores_all)


def revcom(string):
    return string.replace('A', '*').replace('T', 'A').replace('*', 'T').replace('C', '&').replace('G', 'C').replace('&', 'G').replace('a', '*').replace('t', 'a').replace('*', 't').replace('c', '&').replace('g', 'c').replace('&', 'g')[::-1]


def gtf_to_parquet(gtf_path, parquet_path):
    gtf_file = pr.read_gtf(gtf_path)

    cdss_starts_stops_tsv = gtf_file[
        # (gtf_file.Feature == 'start_codon') | # starts not necessary since part and always start of CDS
        (gtf_file.Feature == 'CDS') |
        ((gtf_file.Feature == 'stop_codon') & (gtf_file.Frame == '0')) # only 'real' stop codons; improve this check!
        ]
    cdss_starts_stops_tsv = cdss_starts_stops_tsv[["gene_name", "transcript_name", "Chromosome", "Strand", "Frame", "Start", "End", "exon_number", "Feature", "tag"]].df
    cdss_starts_stops_tsv = cdss_starts_stops_tsv[cdss_starts_stops_tsv["gene_name"].notna()]
    cdss_starts_stops_tsv["exon_number"] = pd.to_numeric(cdss_starts_stops_tsv["exon_number"])
    cdss_starts_stops_tsv["first_transcript_exon"] = cdss_starts_stops_tsv.groupby(['gene_name', 'transcript_name'])['exon_number'].transform('min')
    cdss_starts_stops_tsv["last_transcript_exon"] = cdss_starts_stops_tsv.groupby(['gene_name', 'transcript_name'])['exon_number'].transform('max')
    cdss_starts_stops_tsv = cdss_starts_stops_tsv.drop_duplicates()
    cdss_tsv = cdss_starts_stops_tsv[cdss_starts_stops_tsv.Feature == 'CDS']
    # starts_tsv = cdss_starts_stops_tsv[cdss_starts_stops_tsv.Feature == 'start_codon'] # starts not necessary since part and always start of CDS
    stops_tsv = cdss_starts_stops_tsv[cdss_starts_stops_tsv.Feature == 'stop_codon']
    cols = cdss_tsv.columns

    # join CDSs with stop codons
    cdss_with_stops_tsv = pd.merge(cdss_tsv, stops_tsv, how='outer', on=['Chromosome',
                                                                         'gene_name',
                                                                         'Strand',
                                                                         'transcript_name',
                                                                         'exon_number',
                                                                         'first_transcript_exon',
                                                                         'last_transcript_exon',
                                                                         'tag'], suffixes=('_from_cdss', '_from_stops'))

    # if a stop codon exon has no corresponding CDS exon, use stop codon exon information
    cdss_with_stops_tsv['isolated_stop'] = np.where(cdss_with_stops_tsv.Start_from_cdss.isna() &
                                                    cdss_with_stops_tsv.End_from_cdss.isna() &
                                                    cdss_with_stops_tsv.Frame_from_cdss.isna() &
                                                    cdss_with_stops_tsv.Feature_from_cdss.isna(),
                                                    True,
                                                    False)

    cdss_with_stops_tsv['Start'] = np.where(cdss_with_stops_tsv.isolated_stop,
                                            cdss_with_stops_tsv.Start_from_stops,
                                            cdss_with_stops_tsv.Start_from_cdss)

    cdss_with_stops_tsv['End'] = np.where(cdss_with_stops_tsv.isolated_stop,
                                          cdss_with_stops_tsv.End_from_stops,
                                          cdss_with_stops_tsv.End_from_cdss)

    cdss_with_stops_tsv['Frame'] = np.where(cdss_with_stops_tsv.isolated_stop,
                                            cdss_with_stops_tsv.Frame_from_stops,
                                            cdss_with_stops_tsv.Frame_from_cdss)

    cdss_with_stops_tsv['Feature'] = np.where(cdss_with_stops_tsv.isolated_stop,
                                              cdss_with_stops_tsv.Feature_from_stops,
                                              cdss_with_stops_tsv.Feature_from_cdss)

    # if a stop codon exon has a corresponding CDS exon, integrate this stop codon exon into the CDS exon
    cdss_with_stops_tsv['Start'] = cdss_with_stops_tsv[['Start_from_cdss', 'Start_from_stops']].apply(np.nanmin, axis=1)
    cdss_with_stops_tsv['End'] = cdss_with_stops_tsv[['End_from_cdss', 'End_from_stops']].apply(np.nanmax, axis=1)
    # the frame depends on the first codon from the read direction, not genomic coordinate; otherwise this would need to be changed to the stop codon frame here as well.

    cdss_with_stops_tsv = cdss_with_stops_tsv[cols].drop('Feature', axis=1).astype({'Start': 'int64', 'End': 'int64'}).drop_duplicates() # .drop_duplicates() might be redundant

    # transform tag to MANE_Select yes or no.
    cdss_with_stops_tsv['tag'] = np.where(cdss_with_stops_tsv.tag == 'MANE_Select', True, False)
    cdss_with_stops_tsv = cdss_with_stops_tsv.rename(columns={'tag':'MANE_Select'})

    cdss_with_stops_tsv = pl.read_parquet(cdss_with_stops_tsv.to_parquet())

    cdss_with_stops_tsv = cdss_with_stops_tsv.with_columns(
        transcript_exon_length = pl.col("End") - pl.col("Start")
    ).with_columns(
        transcript_length = pl.col("transcript_exon_length").sum().over(pl.col("transcript_name"), order_by="exon_number")
        ).with_columns(
        transcript_length_cum = pl.col("transcript_exon_length").cum_sum().over(pl.col("transcript_name"), order_by="exon_number")
        ).with_columns(
        transcript_length_before = pl.col("transcript_length_cum") - pl.col("transcript_exon_length")
        )

    cdss_with_stops_tsv.write_parquet(parquet_path)


def gtf_to_parquet_legacy(gtf_path, parquet_path):
    gtf_file = pr.read_gtf(gtf_path)
    cdss_tsv = gtf_file[gtf_file.Feature == "CDS"]
    cdss_tsv = cdss_tsv[["gene_name", "transcript_name", "Chromosome", "Strand", "Frame", "Start", "End", "exon_number"]].df
    cdss_tsv = cdss_tsv[cdss_tsv["gene_name"].notna()]
    cdss_tsv["exon_number"] = pd.to_numeric(cdss_tsv["exon_number"])
    cdss_tsv["first_transcript_exon"] = cdss_tsv.groupby(['gene_name', 'transcript_name'])['exon_number'].transform('min')
    cdss_tsv["last_transcript_exon"] = cdss_tsv.groupby(['gene_name', 'transcript_name'])['exon_number'].transform('max')
    cdss_tsv = cdss_tsv.drop_duplicates()

    cdss_tsv.to_parquet(parquet_path)


def check_parquet(annotation_file, write_parquet):
    # if annotation_file.endswith(".gtf.gz"):
    #     parquet_file = os.path.splitext(os.path.splitext(annotation_file)[0])[0] + ".cdss.bescreen"
    # elif annotation_file.endswith(".gtf"):
    #     parquet_file = os.path.splitext(annotation_file)[0] + ".cdss.bescreen"
    parquet_file = annotation_file + ".cdss.bescreen.parquet"

    if not os.path.isfile(parquet_file) or write_parquet:
        gtf_to_parquet(annotation_file, parquet_file)

    return parquet_file


def get_offset(strand, frame, indexeditcds):

    """
    idea of offset

    original sequence:  ...TTATGACAT...         pos is python indexed!!! pos = position - 1; end = position - 1 + 3
    strand +:                                   base 1st in codon   base 2nd in codon   base 3rd in codon
        exon sequence:  ...TTATGACAT...         index % 3 =
        frame 0:      ...|TTA|TGA|CAT|...       0: pos:end          1: pos-1:end-1      2: pos-2:end-2
        modulo            012 012 012              pos:pos+3           pos-1:pos+2         pos-2:pos+1
        frame 1:      ...T|TAT|GAC|AT.|..       1: pos:end          2: pos-1:end-1      0: pos-2:end-2
        modulo           0 120 120 12              pos:pos+3           pos-1:pos+2         pos-2:pos+1
        frame 2:      ...TT|ATG|ACA|T..|.       2: pos:end          0: pos-1:end-1      1: pos-2:end-2
        modulo           01 201 201 2              pos:pos+3           pos-1:pos+2         pos-2:pos+1
                                            offset:   0                    -1 (1)              -2 (2)

    strand -:
        exon sequence:  ...ATGTCATAA...
        frame 0:      ...|ATG|TCA|TAA|...
        frame 1:      ...A|TGT|CAT|AA.|..
        frame 2:      ...AT|GTC|ATA|A..|.       pos is from the end!!! posr = position - 3 (pos + 1 - 3); endr = position (pos + 1)
    original sequence:  ...TTATGACAT...         (lenght(exon) - index) % 3 = (read from right to left)
        frame 0:      ...|TTA|TGA|ACT|...       0: posr:endr          1: posr+1:endr+1      2: posr+2:endr+2
        modulo            210 210 210              pos-2:pos+1            pos-1:pos+2             pos:pos+3
        frame 1:      ..|.TT|ATG|ACA|T...       1: posr:endr          2: posr+1:endr+1      0: posr+2:endr+2
        modulo            21 021 021 0             pos-2:pos+1            pos-1:pos+2             pos:pos+3
        frame 2:      .|..T|TAT|GAC|AT...       2: posr:endr          0: posr+1:endr+1      1: posr+2:endr+2
        modulo            2 102 102 10             pos-2:pos+1            pos-1:pos+2             pos:pos+3
                                            offsetr:    0                     +1                    +2
                                            offset:    -2 (2)                 -1 (1)                 0
    """

    if strand == '+':
        if frame == 0:
            if indexeditcds % 3 == 0:
                offset = 0
            elif indexeditcds % 3 == 1:
                offset = 1
            elif indexeditcds % 3 == 2:
                offset = 2
        elif frame == 1:
            if indexeditcds % 3 == 1:
                offset = 0
            elif indexeditcds % 3 == 2:
                offset = 1
            elif indexeditcds % 3 == 0:
                offset = 2
        elif frame == 2:
            if indexeditcds % 3 == 2:
                offset = 0
            elif indexeditcds % 3 == 0:
                offset = 1
            elif indexeditcds % 3 == 1:
                offset = 2
    elif strand == '-':
        if frame == 0:
            if indexeditcds % 3 == 0:
                offset = 2
            elif indexeditcds % 3 == 1:
                offset = 1
            elif indexeditcds % 3 == 2:
                offset = 0
        elif frame == 1:
            if indexeditcds % 3 == 1:
                offset = 2
            elif indexeditcds % 3 == 2:
                offset = 1
            elif indexeditcds % 3 == 0:
                offset = 0
        elif frame == 2:
            if indexeditcds % 3 == 2:
                offset = 2
            elif indexeditcds % 3 == 0:
                offset = 1
            elif indexeditcds % 3 == 1:
                offset = 0
    return offset # check, if this is correct


def replace_str_index(text, index=0, replacement=''):
    return f'{text[:index]}{replacement}{text[index+1:]}'


def sort_variantsdf(variantdf):

    variantdf = variantdf.with_row_index(
        'variantindex'
    ).with_columns(
        pl.col('variant')
        .str.split_exact('_', 3)
        .struct.rename_fields(['chrom', 'pos', 'ref', 'alt_tmp'])
        .alias('fields')
    ).unnest(
        'fields'
    ).with_columns(
        alt = pl.when(pl.col('alt_tmp').is_null()).then(pl.col('ref')).otherwise(pl.col('alt_tmp')),
        ref = pl.when(pl.col('alt_tmp').is_null()).then('alt_tmp').otherwise(pl.col('ref'))
    ).drop(
        'alt_tmp'
    ).with_columns(
        pl.col('pos').cast(pl.Int64)
    ).sort(
        'chrom', 'pos', 'ref', 'alt'
    ).drop(
        'chrom', 'pos', 'ref', 'alt'
    )

    return variantdf


def resort_variantsdf(variantdf_sorted):
    return variantdf_sorted.sort('variantindex').drop('variantindex')


def get_be_presets_dict(be_preset_tsv):
    preset_dict = {}
    preset_df = pl.read_csv(be_preset_tsv, separator="\t")
    preset_df = preset_df.filter(
            ~pl.col("is_duplicate") &
            pl.col("is_usable")
        ).select(
            pl.col("name"),
            pl.col("class"),
            pl.col("window_start"),
            pl.col("window_end"),
            pl.col("plus_start"),
            pl.col("plus_end"),
            pl.col("guide_length"),
            pl.col("pam")
        ).filter(
            ~pl.col("name").is_null(),
            ~pl.col("class").is_null(),
            ~pl.col("window_start").is_null(),
            ~pl.col("window_end").is_null(),
            ~pl.col("guide_length").is_null(),
            ~pl.col("pam").is_null()
        ).with_columns(
            pl.col("plus_start").fill_null(strategy="zero"),
            pl.col("plus_end").fill_null(strategy="zero")
        )

    preset_dict['default'] = {"class": 'both',
                              "window_start": 4,
                              "window_end": 8,
                              "plus_start": 0,
                              "plus_end": 0,
                              "guide_length": 20,
                              "pam": 'NG'}

    for row in preset_df.iter_rows(named=True):
        preset_dict[row['name']] = {"class": row['class'],
                                    "window_start": row['window_start'],
                                    "window_end": row['window_end'],
                                    "plus_start": row['plus_start'],
                                    "plus_end": row['plus_end'],
                                    "guide_length": row['guide_length'],
                                    "pam": row['pam']}

    return preset_dict


codon_sun_three_letters = {
    'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys',
    'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',
    'TTA': 'Leu', 'TCA': 'Ser', 'TAA': 'Stop', 'TGA': 'Stop', # replace 'Stop' with '*'
    'TTG': 'Leu', 'TCG': 'Ser', 'TAG': 'Stop', 'TGG': 'Trp',
    'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',
    'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
    'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
    'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
    'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',
    'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
    'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',
    'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
    'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
    'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
} # check, if this is correct; taken from https://stackoverflow.com/questions/61425474/matching-dictionary-values-and-generating-an-output-whether-they-do-or-do-not-ma

codon_sun_one_letter = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': 'Stop', 'TGA': 'Stop', # replace 'Stop' with '*'
    'TTG': 'L', 'TCG': 'S', 'TAG': 'Stop', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
}

iupac_nt_code = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['T'],
    'W': ['A', 'T'],
    'S': ['C', 'G'],
    'M': ['A', 'C'],
    'K': ['G', 'T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}


if __name__ == "__main__":
    pass
