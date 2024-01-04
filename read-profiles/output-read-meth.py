#!/usr/bin/env python3
__description__ =\
"""
Purpose: From a SAM or BAM file, output methylation.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "It works"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from collections import Counter
import subprocess
import re
import pandas as pd
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'reference',
        type=Path,
        help="path to reference genome (.fna/.fasta)")
    parser.add_argument(
        'input',
        type=Path,
        help="path of bam file (.bam/.sam)")
    parser.add_argument(
        'genomic_region',
        type=str,
        help="genomic coordinates of input region (CHR:START-END)")
    parser.add_argument(
        '-s', '--strand',
        action='store_true',
        help="align to other strand (relative to GENOMIC REFERENCE) -- usually for primers with N-notation, as opposed to P")
    parser.add_argument(
        '-m', '--only-meth',
        action='store_true',
        help="only output methylation percent per read")
    parser.add_argument(
        '-n', '--no-pos-filter',
        action='store_true',
        help="don't filter positions by minimal representation")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not args.input.is_file(): parser.error("Input must be a .bam file!")

    return args
# --------------------------------------------------
def _find_CpG(input_ref_seq: str) -> list:
    """
    Find CpG positions in a given string (should be ref_seq) and returns a list of CpG site indices.
    """

    positions = []
    for i, _ in enumerate(input_ref_seq[:-1]):
        if input_ref_seq[i] in 'Cc' and input_ref_seq[i+1] in 'Gg':
            positions.append(i)
    return positions
def _find_C(input_ref_seq: str) -> list:
    """
    Find all cytostine positions in a given string (should be ref_seq) and returns a list of indices.
    Note: CpG included.
    """
    positions = []
    for i, _ in enumerate(input_ref_seq[:-1]):
        if input_ref_seq[i] in 'Cc':
            positions.append(i)
    return positions
def _find_non_CpG(input_ref_seq: str) -> list:
    """
    Find non-CpG cytosine positions in a given string (should be ref_seq) and returns a list of indices.
    """

    positions = []
    for i, _ in enumerate(input_ref_seq[:-1]):
        if input_ref_seq[i] in 'Cc' and input_ref_seq[i+1] not in 'Gg':
            positions.append(i)
    return positions
def _reverse_complement(input_seq: str) -> str:
    """
    Shorthand reverse complementation, apparently
    """
    trans_table = str.maketrans("ATCG", "TAGC")
    return input_seq.upper()[::-1].translate(trans_table)
def _parse_cigar(cigar: str) -> list:
    """
    Parse CIGAR string into number and operations.
    """
    return [(int(count), op) for count, op in re.findall(r'(\d+)([MID])', cigar)]
def _adapt_sequence_to_reference(input_seq: str, cigar_string: str) -> str:
    """
    Function takes input sequence and parses CIGAR string to rectify insertions and deletions, mostly.

    Parameters:
        input_seq (str): input read
        cigar_str (str): CIGAR string

    Returns:
        (str): rectified sequence, dealing with insertions and deletions
    """

    cigar_ops = _parse_cigar(cigar_string)
    adapted_sequence = []
    seq_index = 0
    for count, op in cigar_ops:
        if op == 'M':
            adapted_sequence.append(input_seq[seq_index:seq_index + count])
            seq_index += count
        elif op == 'I':
            # Insertions - skip characters in sequence
            seq_index += count
        elif op == 'D':
            # Deletions - insert gaps in sequence
            adapted_sequence.append('-' * count)
    return ''.join(adapted_sequence)
def _pandas_row_to_genotype(row):
    genotype = row.replace({1: '1', 0: '0'}).fillna('-')
    return ''.join(genotype)
def _get_alignments(sample_path: Path, region_str: str) -> list:
    """
    Function uses samtools to view SAM/BAM file while filtering to specified region.

    Parameters:
        sample_path (Path): input SAM/BAM file path
        region_str (str): genomic coordinate region for filtering
    
    Returns:
        (list): reads in SAM/BAM file, tab-delimited
    """
    
    direct_output = subprocess.run(['samtools', 'view', sample_path, region_str], capture_output=True)
    return direct_output.stdout.decode().split('\n')
def _get_ref_seq(reference_path: Path, ref_region_str: str) -> str:
    """
    Function uses samtools to get reference sequence from alignment info for a given SAM/BAM read.

    Parameters:
        reference_path (Path): path to reference genome
        ref_region_str (str): region string
    
    Returns:
        (str): reference sequence related to SAM/BAM read
    """
    ref_output = subprocess.run(['samtools', 'faidx', reference_path, ref_region_str], capture_output=True)
    ref_seq = ''.join([line for line in ref_output.stdout.decode().split('\n') if not line.startswith('>')])
    return ref_seq
def iter_print_meth(args: Namespace) -> None:
    """
    Function iterates through reads in a SAM/BAM file.

    Parameters:
        (args)
    
    Returns: None
    """

    reads_list: list = _get_alignments(args.input, args.genomic_region)

    all_positions: list = []
    read_profiles: list = []
    stored_refs: dict = {}

    for line in reads_list:
        if len(line.split('\t')) < 9: continue

        read_id, bit_flag, chrom, start_pos, mapq, cigar_str, _, mate_pos, _, seq, *_ = line.split('\t')
        
        adapted_seq = _adapt_sequence_to_reference(seq, cigar_str)
        seq = adapted_seq

        ref_region_str: str = f"{chrom}:{start_pos}-{int(start_pos)+len(seq)-1}"
        if ref_region_str not in stored_refs:
            ref_seq = _get_ref_seq(args.reference, ref_region_str)
            stored_refs[ref_region_str] = ref_seq
        else:
            ref_seq = stored_refs[ref_region_str]

        seq: str = seq if not args.strand else _reverse_complement(seq)
        ref_seq: str = ref_seq if not args.strand else _reverse_complement(ref_seq)

        C_positions = _find_non_CpG(ref_seq)
        CpG_positions = _find_CpG(ref_seq)

        C_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in C_positions]
        CpG_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in CpG_positions if seq[pos] in 'CT']
        if not CpG_pos_methylated: continue

        if args.only_meth:
            print(sum(CpG_pos_methylated)/len(CpG_pos_methylated))
            continue
        adjusted_CpG_pos = [int(start_pos) + pos for pos in CpG_positions] if not args.strand else [int(start_pos)+len(seq) - pos for pos in CpG_positions]
        all_positions += adjusted_CpG_pos
        read_profile = {key: value for key, value in zip(adjusted_CpG_pos, CpG_pos_methylated)}
        read_profiles.append(read_profile)
    
    if args.only_meth: return None
    count_all_positions = Counter(all_positions)
    max_count = max(count_all_positions.values())
    positions_to_exclude = [k for k, v in count_all_positions.items() if v < 0.1 * max_count]
    profiles_df = pd.DataFrame(read_profiles)
    try:
        filtered_df = profiles_df.drop(columns=positions_to_exclude) if not args.no_pos_filter else profiles_df
    except KeyError:
        filtered_df = profiles_df
    profiles_translated = filtered_df.apply(_pandas_row_to_genotype, axis=1)
    for profile in profiles_translated: print(profile)
# --------------------------------------------------
def main() -> None:
    """ Do the thing. """

    args = get_args()

    iter_print_meth(args)
# --------------------------------------------------
if __name__=="__main__":
    main()