import sys

def load_bed(bed_file):
    dups = {}
    with open(bed_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom not in dups:
                dups[chrom] = []
            dups[chrom].append((start, end))
    # Sort and merge overlapping intervals per chromosome
    for chrom in dups:
        dups[chrom].sort()
        merged = []
        for interval in dups[chrom]:
            if not merged or merged[-1][1] < interval[0]:
                merged.append(interval)
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
        dups[chrom] = merged
    return dups

def purge_fasta(fasta_in, bed_file, fasta_out):
    print("Loading duplication intervals...")
    dups = load_bed(bed_file)
    
    print("Processing FASTA sequences...")
    current_chrom = None
    current_seq = []
    
    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if current_chrom:
                    write_purged_seq(current_chrom, "".join(current_seq), dups, fout)
                current_chrom = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_chrom:
            write_purged_seq(current_chrom, "".join(current_seq), dups, fout)
    print("Done! Cleaned assembly saved to", fasta_out)

def write_purged_seq(chrom, seq, dups, fout):
    if chrom not in dups:
        fout.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fout.write(seq[i:i+60] + "\n")
        return
    
    # If the whole contig is marked as duplicate/haplotig, skip it
    intervals = dups[chrom]
    # Simple strategy: if interval covers almost the whole sequence, drop it, 
    # otherwise mask or slice. For standard purge_dups, fully redundant contigs/haplotigs are dropped.
    # Let's check if the first interval covers the whole thing or large chunk:
    total_dup_len = sum(end - start for start, end in intervals)
    if total_dup_len >= len(seq) * 0.9:
        return # Drop entire redundant contig
        
    fout.write(f">{chrom}\n")
    # For internal overlaps, you can either mask with Ns or slice. 
    # Slicing keeps things clean. Let's write the non-duplicate chunks:
    last_idx = 0
    new_seq_chunks = []
    for start, end in intervals:
        if start > last_idx:
            new_seq_chunks.append(seq[last_idx:start])
        last_idx = max(last_idx, end)
    if last_idx < len(seq):
        new_seq_chunks.append(seq[last_idx:])
    
    final_seq = "".join(new_seq_chunks)
    if len(final_seq) > 200: # Keep only if it retains a valid contig length
        for i in range(0, len(final_seq), 60):
            fout.write(final_seq[i:i+60] + "\n")

if __name__ == "__main__":
    purge_fasta("GCA_035222895.1_eeStrDroe1_p1.0_genomic.fna", "dups.bed", "GCA_035222895.1_eeStrDroe1_p1.0_genomic.purged.fna")
