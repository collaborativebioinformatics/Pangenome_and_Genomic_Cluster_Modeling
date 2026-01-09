#!/usr/bin/env python3
"""
Download HPRC raw sequencing data from S3 URLs listed in a TSV file.
"""

import argparse
import os
import pandas as pd


def s3_to_wget_url(s3_url):
    """
    Convert an S3 URL to an HTTPS URL that works with wget.
    
    Args:
        s3_url (str): S3 URL in format s3://bucket-name/path/to/file
        
    Returns:
        str: HTTPS URL in format https://bucket-name.s3.amazonaws.com/path/to/file
        
    Examples:
        >>> s3_to_wget_url("s3://my-bucket/data/file.txt")
        'https://my-bucket.s3.amazonaws.com/data/file.txt'
        
        >>> s3_to_wget_url("s3://my-bucket/folder/subfolder/data.csv")
        'https://my-bucket.s3.amazonaws.com/folder/subfolder/data.csv'
    """
    if not s3_url.startswith("s3://"):
        raise ValueError("URL must start with 's3://'")
    
    # Remove 's3://' prefix
    path = s3_url[5:]
    
    # Split into bucket and key
    parts = path.split('/', 1)
    bucket = parts[0]
    key = parts[1] if len(parts) > 1 else ""
    
    # Construct HTTPS URL
    https_url = f"https://{bucket}.s3.amazonaws.com/{key}"
    
    return https_url


def download_files(tsv_file, output_dir, url_column='AWS FASTA', filename_column='Filename'):
    """
    Download files from S3 URLs listed in a TSV file.
    
    Args:
        tsv_file: Path to TSV file with download information
        output_dir: Directory where files will be downloaded
        url_column: Column name containing S3 URLs
        filename_column: Column name containing output filenames
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    
    print(f"Found {len(df)} files to download")
    print(f"Output directory: {output_dir}")
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        for idx, row in df.iterrows():
            filename = row[filename_column]
            s3_url = row[url_column]
            https_url = s3_to_wget_url(s3_url)
            
            print(f"Downloading {idx + 1}/{len(df)}: {filename}")
            ret = os.system(f"wget -O {filename} {https_url}")
            
            if ret != 0:
                print(f"Warning: Download failed for {filename}")
    finally:
        os.chdir(original_dir)
    
    print("Download complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download HPRC raw sequencing data from S3 URLs"
    )
    parser.add_argument(
        "tsv_file",
        help="Path to TSV file containing download URLs"
    )
    parser.add_argument(
        "output_dir",
        help="Directory where files will be downloaded"
    )
    parser.add_argument(
        "--url-column",
        default="AWS FASTA",
        help="Column name containing S3 URLs (default: 'AWS FASTA')"
    )
    parser.add_argument(
        "--filename-column",
        default="Filename",
        help="Column name containing output filenames (default: 'Filename')"
    )
    
    args = parser.parse_args()
    
    download_files(
        args.tsv_file,
        args.output_dir,
        url_column=args.url_column,
        filename_column=args.filename_column
    )
