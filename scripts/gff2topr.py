#!/usr/bin/env python
## modified script from https://github.com/totajuliusd/topr
import sys
import re

if len(sys.argv) != 6:
    print("Usage: {} <file> <feature_name> <field_name> <biotype_name> <output_filename>".format(sys.argv[0]))
    sys.exit(1)

file_name = sys.argv[1]
feature_name = sys.argv[2]
field_name = sys.argv[3]
biotype_name = sys.argv[4]
output_filename = sys.argv[5]

print('\t'.join([file_name, feature_name, field_name, biotype_name, output_filename]))

try:
    with open(file_name, 'r') as fh:
        genes = {}
        for line in fh:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                if len(columns) < 9:  # Skip malformed lines
                    continue
                    
                chrom, start, end = columns[0], columns[3], columns[4]
                
                if columns[2] == 'exon':
                    # Extract Parent ID for exon
                    parent_match = re.search(r'Parent=([^;]+)', columns[8])
                    if parent_match:
                        parent_id = parent_match.group(1)
                        if parent_id in genes:
                            if "exon_chromstart" in genes[parent_id]:
                                genes[parent_id]["exon_chromstart"].append(start)
                            else:
                                genes[parent_id]["exon_chromstart"] = [start]
                            if "exon_chromend" in genes[parent_id]:
                                genes[parent_id]["exon_chromend"].append(end)
                            else:
                                genes[parent_id]["exon_chromend"] = [end]
                
                elif columns[2] == feature_name:
                    # Extract description and biotype
                    desc_match = re.search(fr'{field_name}=([^;]+)', columns[8])
                    # REMOVE for spontaneum biotype_match = re.search(fr'{biotype_name}=([^;]+)', columns[8])
                    id_match = re.search(r'ID=([^;]+)', columns[8])
                    
                    if id_match and desc_match: # and biotype_match:
                        gene_id = id_match.group(1)
                        gene_type = 'protein_coding' #biotype_match.group(1)
                        if gene_id not in genes:
                            genes[gene_id] = {
                                "gene_start": start,
                                "gene_end": end,
                                "chr": chrom,
                                "gene_name":desc_match.group(1),
                                "biotype": gene_type,
                                "exon_chromstart": [],
                                "exon_chromend": []
                            }

except FileNotFoundError:
    print("File not found: {}".format(file_name))
    sys.exit(1)

# Print the output
with open(output_filename, 'w') as w:
    w.write("\t".join(["chrom", "gene_start", "gene_end", "gene_symbol", "biotype", "exon_chromstart", "exon_chromend"]) + '\n')
    for gene, data in genes.items():
        if not gene:
            continue
        exon_starts = ",".join(data.get("exon_chromstart", []))
        exon_ends = ",".join(data.get("exon_chromend", []))
        w.write("\t".join([
            data["chr"], 
            data["gene_start"], 
            data["gene_end"], 
            data["gene_name"], 
            data["biotype"], 
            exon_starts, 
            exon_ends
        ]) + '\n')
