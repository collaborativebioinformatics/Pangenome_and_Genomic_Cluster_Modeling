#!/usr/bin/env python3
"""
Federated Pangenome Graph Construction Pipeline
================================================
Dockerized version - runs PGGB and vg via Docker containers
Logs saved to /mnt/shared_vol/graphs/

Steps:
- Step 0: Create subchunks (20 individuals each)
- Step 1: Build local graphs with PGGB
- Step 2: Aggregate graphs with vg combine → MEGAGRAPH
- Step 3: Feedback with minigraph → improved local graphs
"""

import subprocess
import os
import sys
import logging
from datetime import datetime
from pathlib import Path
import gzip

# Configuration
INPUT_DIR = "/mnt/shared_vol/hprc_mini_fasta"
SUBCHUNK_DIR = "/mnt/shared_vol/hprc_mini_fasta/subchunks"
OUTPUT_DIR = "/mnt/shared_vol/graphs"
NUM_INDIVIDUALS = 20  # Number of individuals per subchunk
NUM_THREADS = 8

# Log file in shared volume
LOG_FILE = f"{OUTPUT_DIR}/pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def run_command(cmd, description, timeout=None):
    """Run a shell command and log output"""
    logger.info(f"Starting: {description}")
    logger.info(f"Command: {cmd[:500]}..." if len(cmd) > 500 else f"Command: {cmd}")
    
    try:
        result = subprocess.run(
            cmd, shell=True, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=timeout
        )
        logger.info(f"Completed: {description}")
        if result.stdout:
            for line in result.stdout.split('\n')[-20:]:  # Last 20 lines
                logger.debug(line)
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout: {description}")
        return False
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed: {description}")
        if e.stdout:
            logger.error(f"Error output: {e.stdout[-1000:]}")  # Last 1000 chars
        return False

def step0_create_subchunks():
    """Step 0: Extract first N individuals from each chunk"""
    logger.info("=" * 60)
    logger.info("STEP 0: Creating subchunks")
    logger.info("=" * 60)
    
    os.makedirs(SUBCHUNK_DIR, exist_ok=True)
    
    # Find chr19 chunk files (exclude diploid and sub files)
    chunks = []
    for f in sorted(Path(INPUT_DIR).glob("chrom19_chunk*.fa.gz")):
        if "sub" not in f.name and "diploid" not in f.name:
            chunks.append(f)
    
    logger.info(f"Found {len(chunks)} chr19 chunks to process")
    
    for chunk in chunks:
        chunk_num = chunk.name.replace("chrom19_chunk", "").replace(".fa.gz", "")
        output_file = f"{SUBCHUNK_DIR}/chr19_chunk{chunk_num}_sub{NUM_INDIVIDUALS}.fa.gz"
        
        logger.info(f"Processing {chunk.name} -> {os.path.basename(output_file)}")
        
        # Read and extract first N sequences
        try:
            seq_count = 0
            with gzip.open(chunk, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
                current_header = None
                current_seq = []
                
                for line in fin:
                    if line.startswith('>'):
                        # Write previous sequence if exists
                        if current_header and seq_count < NUM_INDIVIDUALS:
                            fout.write(current_header)
                            fout.write(''.join(current_seq))
                            seq_count += 1
                        
                        if seq_count >= NUM_INDIVIDUALS:
                            break
                        
                        current_header = line
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Write last sequence
                if current_header and seq_count < NUM_INDIVIDUALS:
                    fout.write(current_header)
                    fout.write(''.join(current_seq))
                    seq_count += 1
            
            logger.info(f"  Created {output_file}: {seq_count} sequences")
            
        except Exception as e:
            logger.error(f"Failed to process {chunk.name}: {e}")
            return False
    
    # Create samtools index for each subchunk
    logger.info("Creating samtools indexes...")
    for subchunk in Path(SUBCHUNK_DIR).glob("*.fa.gz"):
        cmd = f"samtools faidx {subchunk}"
        run_command(cmd, f"Indexing {subchunk.name}")
    
    return True

def step1_build_local_graphs():
    """Step 1: Run PGGB on each subchunk"""
    logger.info("=" * 60)
    logger.info("STEP 1: Building local graphs with PGGB")
    logger.info("=" * 60)
    
    subchunks = sorted(Path(SUBCHUNK_DIR).glob("chr19_chunk*_sub*.fa.gz"))
    logger.info(f"Found {len(subchunks)} subchunks to process")
    
    for subchunk in subchunks:
        chunk_name = subchunk.stem.replace(".fa", "")
        output_subdir = f"{OUTPUT_DIR}/{chunk_name}"
        os.makedirs(output_subdir, exist_ok=True)
        
        # Get sequence count for -n parameter
        count_cmd = f"zcat {subchunk} | grep -c '^>'"
        result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
        n_seqs = int(result.stdout.strip()) if result.stdout.strip() else NUM_INDIVIDUALS
        
        logger.info(f"Processing {subchunk.name} with {n_seqs} sequences")
        
        # Run PGGB via Docker
        pggb_cmd = f"""
docker run --rm \
    -v {SUBCHUNK_DIR}:/data:ro \
    -v {output_subdir}:/output \
    ghcr.io/pangenome/pggb:latest \
    pggb -i /data/{subchunk.name} \
         -o /output \
         -n {n_seqs} \
         -t {NUM_THREADS} \
         -p 90 \
         -s 10000
"""
        
        success = run_command(pggb_cmd, f"PGGB on {subchunk.name}", timeout=7200)  # 2 hour timeout per chunk
        
        if not success:
            logger.warning(f"PGGB failed for {subchunk.name}, continuing to next chunk...")
            continue
        
        # Verify output
        gfa_files = list(Path(output_subdir).glob("*.smooth.final.gfa"))
        if gfa_files:
            logger.info(f"  Output GFA: {gfa_files[0].name}")
        else:
            logger.warning(f"  No GFA output found for {subchunk.name}")
    
    return True

def step2_aggregate_graphs():
    """Step 2: Combine all local GFAs into MEGAGRAPH using vg combine"""
    logger.info("=" * 60)
    logger.info("STEP 2: Aggregating graphs into MEGAGRAPH")
    logger.info("=" * 60)
    
    # Find all GFA files from PGGB output
    gfa_files = []
    for chunk_dir in sorted(Path(OUTPUT_DIR).glob("chr19_chunk*_sub*")):
        if chunk_dir.is_dir():
            gfas = list(chunk_dir.glob("*.smooth.final.gfa"))
            if gfas:
                gfa_files.append(gfas[0])
                logger.info(f"Found GFA: {gfas[0]}")
    
    if len(gfa_files) < 2:
        logger.error(f"Need at least 2 GFA files to combine, found {len(gfa_files)}")
        return False
    
    logger.info(f"Combining {len(gfa_files)} GFA files")
    
    # Convert each GFA to VG
    vg_files = []
    for i, gfa in enumerate(gfa_files):
        vg_file = f"temp_{i}.vg"
        vg_path = f"{OUTPUT_DIR}/{vg_file}"
        
        cmd = f"""
docker run --rm \
    -v {gfa.parent}:/input:ro \
    -v {OUTPUT_DIR}:/output \
    quay.io/vgteam/vg:v1.71.0 \
    vg convert /input/{gfa.name} -v > {vg_path}
"""
        
        if run_command(cmd, f"Converting {gfa.name} to VG"):
            if os.path.exists(vg_path) and os.path.getsize(vg_path) > 0:
                vg_files.append(vg_file)
                logger.info(f"  Created {vg_file}")
            else:
                logger.warning(f"  VG file empty or not created: {vg_file}")
        else:
            logger.warning(f"  Failed to convert {gfa.name}")
    
    if len(vg_files) < 2:
        logger.error("Not enough VG files to combine")
        return False
    
    # Combine VG files iteratively
    logger.info(f"Combining {len(vg_files)} VG files into MEGAGRAPH")
    
    # Start with first file
    current_combined = f"{OUTPUT_DIR}/{vg_files[0]}"
    
    for i in range(1, len(vg_files)):
        next_vg = f"{OUTPUT_DIR}/{vg_files[i]}"
        temp_output = f"{OUTPUT_DIR}/temp_combined_{i}.vg"
        
        cmd = f"""
docker run --rm \
    -v {OUTPUT_DIR}:/graphs \
    quay.io/vgteam/vg:v1.71.0 \
    vg combine -p /graphs/{os.path.basename(current_combined)} /graphs/{vg_files[i]} > {temp_output}
"""
        
        if not run_command(cmd, f"Combining graph {i}/{len(vg_files)-1}"):
            logger.error("Failed to combine graphs")
            return False
        
        # Clean up previous combined file if it's a temp file
        if "temp_combined" in current_combined:
            os.remove(current_combined)
        
        current_combined = temp_output
    
    # Rename final combined file to MEGAGRAPH.vg
    os.rename(current_combined, f"{OUTPUT_DIR}/MEGAGRAPH.vg")
    logger.info(f"Created MEGAGRAPH.vg")
    
    # Convert MEGAGRAPH to GFA for minigraph
    convert_cmd = f"""
docker run --rm \
    -v {OUTPUT_DIR}:/graphs \
    quay.io/vgteam/vg:v1.71.0 \
    vg convert -f /graphs/MEGAGRAPH.vg > {OUTPUT_DIR}/MEGAGRAPH.gfa
"""
    
    if not run_command(convert_cmd, "Converting MEGAGRAPH to GFA"):
        return False
    
    logger.info(f"Created MEGAGRAPH.gfa")
    
    # Cleanup temp VG files
    for vg in vg_files:
        vg_path = f"{OUTPUT_DIR}/{vg}"
        if os.path.exists(vg_path):
            os.remove(vg_path)
    
    return True

def step3_feedback_loop():
    """Step 3: Use minigraph to improve local graphs with MEGAGRAPH"""
    logger.info("=" * 60)
    logger.info("STEP 3: Feedback loop with minigraph")
    logger.info("=" * 60)
    
    federated_dir = f"{OUTPUT_DIR}/federated"
    os.makedirs(federated_dir, exist_ok=True)
    
    megagraph = f"{OUTPUT_DIR}/MEGAGRAPH.gfa"
    if not os.path.exists(megagraph):
        logger.error("MEGAGRAPH.gfa not found")
        return False
    
    subchunks = sorted(Path(SUBCHUNK_DIR).glob("chr19_chunk*_sub*.fa.gz"))
    
    for subchunk in subchunks:
        chunk_name = subchunk.stem.replace(".fa", "")
        output_gfa = f"{federated_dir}/{chunk_name}_federated.gfa"
        
        logger.info(f"Feedback for {chunk_name}")
        
        # Decompress subchunk temporarily
        temp_fa = f"/tmp/{chunk_name}.fa"
        decompress_cmd = f"zcat {subchunk} > {temp_fa}"
        run_command(decompress_cmd, f"Decompressing {subchunk.name}")
        
        # minigraph: MEGAGRAPH.gfa + chunk.fa → improved GFA
        cmd = f"minigraph -cxggs {megagraph} {temp_fa} > {output_gfa}"
        
        success = run_command(cmd, f"Minigraph feedback for {chunk_name}")
        
        # Cleanup temp file
        if os.path.exists(temp_fa):
            os.remove(temp_fa)
        
        if not success:
            logger.warning(f"Minigraph failed for {chunk_name}, continuing...")
            continue
        
        # Log stats
        if os.path.exists(output_gfa) and os.path.getsize(output_gfa) > 0:
            count_cmd = f"grep -c '^S' {output_gfa}"
            result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
            logger.info(f"  {os.path.basename(output_gfa)}: {result.stdout.strip()} nodes")
        else:
            logger.warning(f"  Output GFA empty or not created: {output_gfa}")
    
    return True

def print_summary():
    """Print final summary of outputs"""
    logger.info("=" * 60)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 60)
    
    # Subchunks
    subchunks = list(Path(SUBCHUNK_DIR).glob("*.fa.gz"))
    logger.info(f"\nSubchunks created: {len(subchunks)}")
    for s in subchunks:
        logger.info(f"  - {s.name}")
    
    # PGGB outputs
    pggb_dirs = list(Path(OUTPUT_DIR).glob("chr19_chunk*_sub*"))
    logger.info(f"\nPGGB output directories: {len(pggb_dirs)}")
    for d in pggb_dirs:
        gfas = list(d.glob("*.gfa"))
        logger.info(f"  - {d.name}: {len(gfas)} GFA files")
    
    # MEGAGRAPH
    mega_vg = Path(f"{OUTPUT_DIR}/MEGAGRAPH.vg")
    mega_gfa = Path(f"{OUTPUT_DIR}/MEGAGRAPH.gfa")
    logger.info(f"\nMEGAGRAPH:")
    logger.info(f"  - MEGAGRAPH.vg: {'EXISTS' if mega_vg.exists() else 'MISSING'}")
    logger.info(f"  - MEGAGRAPH.gfa: {'EXISTS' if mega_gfa.exists() else 'MISSING'}")
    
    # Federated outputs
    fed_dir = Path(f"{OUTPUT_DIR}/federated")
    if fed_dir.exists():
        fed_gfas = list(fed_dir.glob("*.gfa"))
        logger.info(f"\nFederated outputs: {len(fed_gfas)}")
        for f in fed_gfas:
            logger.info(f"  - {f.name}")

def main():
    """Run the complete federated pangenome pipeline"""
    start_time = datetime.now()
    
    logger.info("=" * 60)
    logger.info("FEDERATED PANGENOME CONSTRUCTION PIPELINE")
    logger.info(f"Started at: {start_time}")
    logger.info(f"Log file: {LOG_FILE}")
    logger.info("=" * 60)
    logger.info(f"Configuration:")
    logger.info(f"  Input directory: {INPUT_DIR}")
    logger.info(f"  Subchunk directory: {SUBCHUNK_DIR}")
    logger.info(f"  Output directory: {OUTPUT_DIR}")
    logger.info(f"  Individuals per subchunk: {NUM_INDIVIDUALS}")
    logger.info(f"  Threads: {NUM_THREADS}")
    
    # Run all steps
    steps = [
        ("Step 0: Create subchunks", step0_create_subchunks),
        ("Step 1: Build local graphs (PGGB)", step1_build_local_graphs),
        ("Step 2: Aggregate graphs (vg combine)", step2_aggregate_graphs),
        ("Step 3: Feedback loop (minigraph)", step3_feedback_loop),
    ]
    
    results = {}
    for step_name, step_func in steps:
        logger.info(f"\n{'='*60}")
        logger.info(f"STARTING: {step_name}")
        logger.info(f"{'='*60}\n")
        
        try:
            success = step_func()
            results[step_name] = "SUCCESS" if success else "FAILED"
        except Exception as e:
            logger.error(f"Exception in {step_name}: {e}")
            results[step_name] = f"ERROR: {e}"
        
        logger.info(f"\n{step_name}: {results[step_name]}")
    
    # Print summary
    print_summary()
    
    end_time = datetime.now()
    duration = end_time - start_time
    
    logger.info("\n" + "=" * 60)
    logger.info("PIPELINE COMPLETED")
    logger.info("=" * 60)
    logger.info(f"Duration: {duration}")
    logger.info(f"Results:")
    for step, result in results.items():
        logger.info(f"  {step}: {result}")
    logger.info(f"\nLogs saved to: {LOG_FILE}")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()
