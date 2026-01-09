#!/usr/bin/env python3
"""
Extract chromosome 19 and 22 sequences from HPRC assemblies into chunked PanSN FASTA files
using pre-generated contig list files.
"""

import argparse
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

def load_contig_list(contig_list_file):
    """Load contig identifiers (one per line)."""
    contigs = set()
    with open(contig_list_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            contigs.add(line)
    return contigs





def open_output(path, use_bgzip=False):
    """Return a text handle and a closer for gzip/bgzip output."""
    if use_bgzip:
        bgzip_bin = shutil.which('bgzip')
        if bgzip_bin:
            raw_out = open(path, 'wb')
            proc = subprocess.Popen([bgzip_bin, '-c'], stdin=subprocess.PIPE, stdout=raw_out)
            text_handle = io.TextIOWrapper(proc.stdin, encoding='utf-8', newline='\n')

            def closer():
                text_handle.close()
                proc.wait()
                raw_out.close()

            return text_handle, closer
        else:
            print("bgzip not found; falling back to gzip", file=sys.stderr)

    gz_handle = gzip.open(path, 'wt')

    def closer():
        gz_handle.close()

    return gz_handle, closer


def process_single_file(args):
    """Process a single FASTA file and return extracted sequences by chromosome and sample."""
    fasta_file, sample_name, haplotype_id, chr19_contigs, chr22_contigs, file_idx, total_files = args
    
    sequences = {'19': [], '22': []}  # List of (header, seq) tuples
    stats = {'19': 0, '22': 0}
    
    print(f"[{file_idx}/{total_files}] Processing {sample_name} haplotype {haplotype_id}...", file=sys.stderr)
    
    try:
        with gzip.open(fasta_file, 'rt') as in_f:
            current_header = None
            current_seq = []
            current_chrom = None

            for line in in_f:
                line = line.strip()

                if line.startswith('>'):
                    # Flush previous
                    if current_chrom and current_header and current_seq:
                        header_no_gt = current_header[1:]
                        pansn_header = f">{sample_name}#{haplotype_id}#{header_no_gt}"
                        seq_data = ''.join(current_seq)
                        sequences[current_chrom].append((pansn_header, seq_data))
                        stats[current_chrom] += 1

                    # Start new contig
                    current_header = line
                    current_seq = []
                    current_chrom = None

                    # Check if this contig is in our contig lists
                    header_no_gt = current_header[1:]
                    if header_no_gt in chr19_contigs:
                        current_chrom = '19'
                    elif header_no_gt in chr22_contigs:
                        current_chrom = '22'
                    
                else:
                    if current_chrom:
                        current_seq.append(line)

            # Flush last contig
            if current_chrom and current_header and current_seq:
                header_no_gt = current_header[1:]
                pansn_header = f">{sample_name}#{haplotype_id}#{header_no_gt}"
                seq_data = ''.join(current_seq)
                sequences[current_chrom].append((pansn_header, seq_data))
                stats[current_chrom] += 1
    
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}", file=sys.stderr)
    
    return sample_name, haplotype_id, sequences, stats


def extract_chr_sequences(
    input_dir,
    output_dir,
    chr19_list_file=None,
    chr22_list_file=None,
    chunk_size=20,
    use_bgzip=True,
    num_cores=32,
):
    """
    Extract chromosome 19 and 22 sequences from all haplotype assemblies
    and write to chunked PanSN-compliant FASTA files (20 individuals per chunk).
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_path.glob("*_hap[12]_hprc*.fa.gz"))
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fasta_files)} haplotype assemblies")

    # Load contig lists for each chromosome
    chr19_contigs = set()
    chr22_contigs = set()
    if chr19_list_file:
        chr19_contigs = load_contig_list(chr19_list_file)
        print(f"Loaded {len(chr19_contigs)} chr19 contigs from {chr19_list_file}")
    if chr22_list_file:
        chr22_contigs = load_contig_list(chr22_list_file)
        print(f"Loaded {len(chr22_contigs)} chr22 contigs from {chr22_list_file}")

    # Build a mapping of (sample, haplotype) -> list of contigs for each chromosome
    sample_chr19_contigs = defaultdict(list)
    sample_chr22_contigs = defaultdict(list)
    
    for header in chr19_contigs:
        parts = header.split('#')
        if len(parts) >= 2:
            sample = parts[0]
            hap = parts[1]
            sample_chr19_contigs[(sample, hap)].append(header)
    
    for header in chr22_contigs:
        parts = header.split('#')
        if len(parts) >= 2:
            sample = parts[0]
            hap = parts[1]
            sample_chr22_contigs[(sample, hap)].append(header)
    
    # Get unique samples in order
    samples_19 = sorted(set(s for s, h in sample_chr19_contigs.keys()))
    samples_22 = sorted(set(s for s, h in sample_chr22_contigs.keys()))
    
    print(f"\nChr19: {len(samples_19)} samples, {len(chr19_contigs)} total contigs")
    print(f"Chr22: {len(samples_22)} samples, {len(chr22_contigs)} total contigs")

    print(f"\n=== STEP 1: Creating chunk files ===", file=sys.stderr)

    # Prepare chunk files for each chromosome
    chunk_files_19 = {}
    chunk_files_22 = {}
    sequences_by_chunk = {'19': defaultdict(int), '22': defaultdict(int)}
    
    def create_chunk_files(samples, chrom, output_dir, chunk_size, use_bgzip):
        """Create chunk files for a chromosome's samples."""
        chunks = {}
        for chunk_idx, i in enumerate(range(0, len(samples), chunk_size), start=1):
            chunk_samples = samples[i:i+chunk_size]
            chunk_path = output_dir / f"chrom{chrom}_chunk{chunk_idx}.fa.gz"
            out_handle, closer = open_output(str(chunk_path), use_bgzip=use_bgzip)
            chunks[chunk_idx] = {
                'path': chunk_path,
                'handle': out_handle,
                'closer': closer,
                'samples': set(chunk_samples)
            }
            print(f"Created chunk {chunk_idx} for {len(chunk_samples)} samples: {chunk_path}", file=sys.stderr)
        return chunks
    
    if chr19_contigs:
        chunk_files_19 = create_chunk_files(samples_19, '19', output_path, chunk_size, use_bgzip)
    if chr22_contigs:
        chunk_files_22 = create_chunk_files(samples_22, '22', output_path, chunk_size, use_bgzip)
    
    sequences_written = {'19': 0, '22': 0}
    samples_processed = set()

    print(f"\n=== STEP 2: Extracting sequences in parallel ({num_cores} cores) ===", file=sys.stderr)

    # Prepare arguments for parallel processing
    process_args = []
    total_files = len(fasta_files)
    for file_idx, fasta_file in enumerate(fasta_files, start=1):
        match = re.match(r'([^_]+)_hap([12])_', fasta_file.name)
        if not match:
            print(f"Warning: Could not parse filename {fasta_file.name}, skipping")
            continue
        
        sample_name = match.group(1)
        haplotype_id = match.group(2)
        samples_processed.add(sample_name)
        
        process_args.append((
            fasta_file,
            sample_name,
            haplotype_id,
            chr19_contigs,
            chr22_contigs,
            file_idx,
            total_files
        ))

    # Process files in parallel
    print(f"Processing {len(process_args)} files with {num_cores} cores...", file=sys.stderr)
    with Pool(num_cores) as pool:
        results = pool.map(process_single_file, process_args)

    # Write results to chunk files
    print(f"\n=== STEP 3: Writing sequences to chunk files ===", file=sys.stderr)
    
    # Organize sequences by sample and chromosome for efficient writing
    for sample_name, haplotype_id, sequences, stats in results:
        # Write chr19 sequences
        for header, seq in sequences['19']:
            for chunk_idx, chunk_info in chunk_files_19.items():
                if sample_name in chunk_info['samples']:
                    chunk_info['handle'].write(header + '\n')
                    chunk_info['handle'].write(seq + '\n')
                    sequences_by_chunk['19'][chunk_idx] += 1
            sequences_written['19'] += 1
        
        # Write chr22 sequences
        for header, seq in sequences['22']:
            for chunk_idx, chunk_info in chunk_files_22.items():
                if sample_name in chunk_info['samples']:
                    chunk_info['handle'].write(header + '\n')
                    chunk_info['handle'].write(seq + '\n')
                    sequences_by_chunk['22'][chunk_idx] += 1
            sequences_written['22'] += 1

    # Close all chunk handles
    for chunk_idx, chunk_info in chunk_files_19.items():
        chunk_info['closer']()
    for chunk_idx, chunk_info in chunk_files_22.items():
        chunk_info['closer']()

    print("\n=== COMPLETE ===")
    print(f"Samples processed: {len(samples_processed)}")
    
    if chr19_contigs:
        print(f"\nChromosome 19:")
        print(f"  Total sequences written: {sequences_written.get('19', 0)}")
        for chunk_idx in sorted(chunk_files_19.keys()):
            print(f"  Chunk {chunk_idx}: {sequences_by_chunk['19'][chunk_idx]} sequences -> {chunk_files_19[chunk_idx]['path']}")
    
    if chr22_contigs:
        print(f"\nChromosome 22:")
        print(f"  Total sequences written: {sequences_written.get('22', 0)}")
        for chunk_idx in sorted(chunk_files_22.keys()):
            print(f"  Chunk {chunk_idx}: {sequences_by_chunk['22'][chunk_idx]} sequences -> {chunk_files_22[chunk_idx]['path']}")

    return sequences_written.get('22', 0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract chromosome 19 and 22 sequences from HPRC assemblies into chunked PanSN FASTA files")
    parser.add_argument("input_dir", help="Directory with *_hap[12]_hprc*.fa.gz")
    parser.add_argument("output_dir", help="Output directory for chunked FASTA files")
    parser.add_argument("--chr19-list", dest="chr19_list_file", required=True,
                        help="File with chr19 contig headers (one per line)")
    parser.add_argument("--chr22-list", dest="chr22_list_file", required=True,
                        help="File with chr22 contig headers (one per line)")
    parser.add_argument("--chunk-size", type=int, default=20,
                        help="Number of individuals per chunk file (default 20)")
    parser.add_argument("--no-bgzip", dest="use_bgzip", action="store_false",
                        help="Use gzip instead of bgzip for output")

    args = parser.parse_args()

    extract_chr_sequences(
        args.input_dir,
        args.output_dir,
        chr19_list_file=args.chr19_list_file,
        chr22_list_file=args.chr22_list_file,
        chunk_size=args.chunk_size,
        use_bgzip=args.use_bgzip,
        num_cores=24,
    )