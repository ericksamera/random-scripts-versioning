#!/usr/bin/env python3
__description__ =\
"""
Purpose: Process SAM files to assess methylation.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "lmao it works"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import subprocess
import re
import pandas as pd
from collections import Counter

primers_dict: dict = {
    "DMAP1": ("NC_037330.1:101832510-101832748", True),
    "DNMT3A": ("NC_037338.1:74014960-74015346", True),
    "DNMT3B": ("NC_037340.1:62166550-62166873", True),
    "GNAS": ("NC_037340.1:57520605-57520891", False),
    "H19": ("NC_037356.1:49504607-49504946", True),
    "IGF2R": ("NC_037336.1:96223067-96223367", False),
    "KCNQ1": ("NC_037356.1:48907634-48907882", True),
    "LIF": ("NC_037344.1:69264000-69264354", False),
    "LIFR": ("NC_037347.1:35949053-35949420", True),
    "MEST": ("NC_037331.1:94249893-94250283", False),
    "NNAT": ("NC_037340.1:66465753-66466057", True),
    "PEG10": ("NC_037331.1:12063279-12063673", True),
    "PEG3": ("NC_037345.1:64120688-64120941", False),
    "PLAGL1": ("NC_037336.1:81418732-81419041", False),
    "RTL1": ("NC_037348.1:65778329-65778707", False),
    "SLC2A8": ("NC_037338.1:98154376-98154714", True),
    "SNRPN": ("NC_037348.1:1937284-1937513", True),
    "SUV39H1": ("NC_037357.1:86781615-86781929", False),
    "TXNIP": ("NC_037330.1:21383423-21383784", True),
    "XIST": ("NC_037357.1:77162729-77163121", True),
}
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path of bam file")
    parser.add_argument(
        'input_gene',
        type=str,
        help="path of bam file")

    parser.add_argument(
        '-s', '--sample',
        type=str,
        metavar='SAMPLE',
        default=None,
        help="override the filename with a sample name in outputs")

    parser.add_argument(
        '--out-profiles',
        type=Path,
        metavar='PATH',
        default=None,
        help="output methylation profiles per read and counts (.csv)")
    parser.add_argument(
        '--out-perc-meth',
        type=Path,
        metavar='PATH',
        default=None,
        help="output methylation percent across reads (.csv)")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not args.input_path.is_file(): parser.error("Input must be a .bam file!")
    if (not args.out_profiles) and (not args.out_perc_meth): parser.error("Choose at least one of --out-profiles/--out-perc-meth !")

    return args
# --------------------------------------------------
def _reverse_complement(input_seq: str) -> str:
    """
    Shorthand reverse complementation, apparently
    """
    trans_table = str.maketrans("ATCG", "TAGC")
    return input_seq.upper()[::-1].translate(trans_table)
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

def _parse_cigar(cigar: str) -> list:
    """
    Parse CIGAR string into number and operations.
    """
    return [(int(count), op) for count, op in re.findall(r'(\d+)([MID])', cigar)]

def adapt_sequence_to_reference(input_seq: str, cigar_string: str) -> str:
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

def _get_alignments(sample_path: Path, region_str: str) -> list:
    """
    Function uses samtools to view sam/bam file while filtering to specified region.

    Parameters:
        sample_path (Path): input bam file path
        region_str (str): genomic coordinate region for filtering
    """
    
    direct_output = subprocess.run(['samtools', 'view', sample_path, region_str], capture_output=True)
    return direct_output.stdout.decode().split('\n')

def _get_ref_seq(ref_region_str: str) -> str:
    """
    Function uses samtools to get 
    """
    REFERENCE_PATH: Path = Path("GCF_002263795.2_ARS-UCD1.3_genomic.fna")
    ref_output = subprocess.run(['samtools', 'faidx', REFERENCE_PATH, ref_region_str], capture_output=True)
    ref_seq = ''.join([line for line in ref_output.stdout.decode().split('\n') if not line.startswith('>')])
    return ref_seq

def _pandas_row_to_genotype(row):
    genotype = row.replace({1: '1', 0: '0'}).fillna('-')
    return ''.join(genotype)

def iter_print_meth(sample_path: Path, region_str: str, strand: bool, print_profiles: bool = False) -> (list[float], list[tuple]):
    """
    Function opens a given bam file, performs filtering to target region, and outputs methylation status of reads.

    Parameters:
        sample_path (Path): input bam file path
        region_str (str): genomic coordinate region for filtering
        strand (bool): handles alignment of target strand, i.e., BP1X is the (+) strand = True; BN2X is the (-) strand = False
        print_profiles (bool): print profiles
    
    Returns:
        list[float]: list of methylation percent of reads
    """

    profiles = []
    all_positions = []
    methylation_per_read = []

    reads_list: list = _get_alignments(sample_path, region_str)

    return_result = [None, None]
    for line in reads_list:
        if len(line.split('\t')) < 9: continue

        read_id, bit_flag, chrom, start_pos, mapq, cigar_str, _, mate_pos, _, seq, *_ = line.split('\t')
        
        # deal with CIGAR strings
        adapted_seq = adapt_sequence_to_reference(seq, cigar_str)
        seq = adapted_seq

        ref_region_str: str = f"{chrom}:{start_pos}-{int(start_pos)+len(seq)-1}"
        ref_seq = _get_ref_seq(ref_region_str)

        seq: str = seq if strand else _reverse_complement(seq)
        ref_seq: str = ref_seq if strand else _reverse_complement(ref_seq)
        
        C_positions = _find_non_CpG(ref_seq)
        CpG_positions = _find_CpG(ref_seq)

        C_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in C_positions]
        CpG_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in CpG_positions if seq[pos] in 'CT']
        if not CpG_pos_methylated: continue

        read_methylation_perc: float = ''.join([str(i) for i in CpG_pos_methylated]).count('1')/len(CpG_pos_methylated)

        if print_profiles:
            adjusted_CpG_pos = [int(start_pos) + pos for pos in CpG_positions] if strand else [int(start_pos)+len(seq) - pos for pos in CpG_positions]
            all_positions += adjusted_CpG_pos
            read_profile = {key: value for key, value in zip(adjusted_CpG_pos, CpG_pos_methylated)}        
            profiles.append(read_profile)

        methylation_per_read.append(read_methylation_perc)
    return_result[0] = methylation_per_read

    profiles_with_counts = []
    if print_profiles:
        count_all_positions = Counter(all_positions)
        max_count = max(count_all_positions.values())
        positions_to_exclude = [k for k, v in count_all_positions.items() if v < 0.1 * max_count]

        profiles_df = pd.DataFrame(profiles)
        try:
            filtered_df = profiles_df.drop(columns=positions_to_exclude)
            profiles_translated = filtered_df.apply(_pandas_row_to_genotype, axis=1)
            profiles_with_counts = sorted(Counter(profiles_translated).items(), key=lambda x: x[1])
        except KeyError:
            print(f"ERROR with {sample_path}")
    return_result[1] = profiles_with_counts

    return return_result

def _write_profiles(args, input_list: list) -> None:
    """
    Function writes methylation profiles and counts into long CSV, especially for parsing in R. 
    """
    with open(args.out_profiles, 'w') as output_file:
        for profile, counts in input_list:
            row_info = [profile, counts]
            output_file.write(','.join([str(i) for i in row_info]) + "\n")

def _write_perc_meth(args, input_list: list) -> None:
    """
    Function writes percent methylation across all reads into long CSV, especially for parsing in R.
    """
    with open(args.out_perc_meth, 'w') as output_file:
        for methylation_percent in input_list:
            row_info = [args.sample_name, args.input_gene, methylation_percent]
            output_file.write(','.join([str(i) for i in row_info]) + "\n")

def main():
    """ Do the thing. """

    args = get_args()

    input_path: Path = args.input_path
    input_gene: str = args.input_gene

    args.sample_name = args.sample if args.sample else input_path

    methylation_per_read = iter_print_meth(input_path, *primers_dict[input_gene], True if args.out_profiles else False)

    if args.out_profiles: _write_profiles(args, methylation_per_read[1])
    if args.out_perc_meth: _write_perc_meth(args, methylation_per_read[0])

if __name__=="__main__":
    main()
