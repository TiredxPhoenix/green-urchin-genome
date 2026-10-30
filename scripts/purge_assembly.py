import sys

def purge_by_self_alignment(fasta_in, paf_in, fasta_out, overlap_fraction=0.92):
    print("Step 1: Reading contig lengths from FASTA...")
    contig_lens = {}
    current_chrom = None
    current_seq_len = 0
    
    with open(fasta_in) as f:
        for line in f:
            if line.startswith(">"):
                if current_chrom:
                    contig_lens[current_chrom] = current_seq_len
                current_chrom = line.strip().split()[0][1:]
                current_seq_len = 0
            else:
                current_seq_len += len(line.strip())
        if current_chrom:
            contig_lens[current_chrom] = current_seq_len

    print(f"Loaded {len(contig_lens)} contigs.")

    print("Step 2: Analyzing self-alignments in PAF to find redundant contigs...")
    q_intervals = {chrom: [] for chrom in contig_lens}
    
    with open(paf_in) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 11: continue
            qname, qlen, qstart, qend = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            tname, tlen, tstart, tend = parts[5], int(parts[6]), int(parts[7]), int(parts[8])
            
            # CRITICAL FIX: Skip self-hits entirely. We only want cross-contig redundancies.
            if qname == tname:
                continue
                
            # If query is shorter than or equal to target, consider query as the potential haplotig
            if qname in contig_lens and tname in contig_lens:
                if contig_lens[qname] <= contig_lens[tname]:
                    q_intervals[qname].append((qstart, qend))

    # Merge overlapping alignment intervals per query contig to calculate true coverage fraction
    dropped_contigs = set()
    for chrom, intervals in q_intervals.items():
        if not intervals:
            continue
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or merged[-1][1] < start:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        
        total_covered = sum(end - start for start, end in merged)
        length = contig_lens.get(chrom, 0)
        if length > 0 and (total_covered / length) >= overlap_fraction:
            dropped_contigs.add(chrom)

    print(f"Identified {len(dropped_contigs)} redundant haplotig contigs to drop.")

    print("Step 3: Writing clean purged FASTA...")
    kept_count = 0
    kept_len = 0
    
    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        write_flag = False
        for line in fin:
            if line.startswith(">"):
                chrom = line.strip().split()[0][1:]
                if chrom not in dropped_contigs:
                    write_flag = True
                    kept_count += 1
                    fout.write(line)
                else:
                    write_flag = False
            else:
                if write_flag:
                    fout.write(line)
                    kept_len += len(line.strip())

    print(f"Done! Kept {kept_count} contigs. New assembly size: {kept_len / 1e6:.2f} Mb")

if __name__ == "__main__":
    purge_by_self_alignment(
        "GCA_035222895.1_eeStrDroe1_p1.0_genomic.fna", 
        "self.paf", 
        "GCA_035222895.1_eeStrDroe1_p1.0_genomic.purged.fna",
        overlap_fraction=0.92
    )
