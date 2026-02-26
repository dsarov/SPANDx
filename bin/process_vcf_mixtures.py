import sys
import os

def get_type(ref, alts):
    """Determine if variant is SNP or INDEL based on length."""
    # Check the first real alt allele (ignore '*')
    for alt in alts:
        if alt == '*': 
            continue
        if len(ref) != len(alt):
            return "INDEL"
    return "SNP"

def parse_snpeff(info_str):
    """
    Parses the SnpEff 'EFF' string.
    Format usually: Effect(Effect_Impact|Functional_Class|Codon_Change|AAChange|Length|Gene_Name|Transcript_BioType...)
    """
    annotations = {
        'Effect': '.', 'Impact': '.', 'Functional_Class': '.', 
        'Codon_change': '.', 'Protein_and_nucleotide_change': '.', 
        'Amino_Acid_Length': '.', 'Gene_name': '.', 'Biotype': '.'
    }
    
    if 'EFF=' not in info_str:
        return annotations

    try:
        # Extract the value after EFF= until the next semicolon
        eff_full = info_str.split('EFF=')[1].split(';')[0]
        # Handle multiple effects (comma separated); take the first one
        first_effect = eff_full.split(',')[0]
        
        # SnpEff Format: Effect(Part|Part|...)
        if '(' in first_effect:
            effect_type = first_effect.split('(')[0]
            details = first_effect.split('(')[1].rstrip(')')
            parts = details.split('|')
            
            # 0: Impact, 1: FuncClass, 2: Codon, 3: AA, 4: Length, 5: Gene, 6: BioType...
            annotations['Effect'] = effect_type
            if len(parts) > 0: annotations['Impact'] = parts[0]
            if len(parts) > 1: annotations['Functional_Class'] = parts[1]
            if len(parts) > 2: annotations['Codon_change'] = parts[2]
            if len(parts) > 3: annotations['Protein_and_nucleotide_change'] = parts[3]
            if len(parts) > 4: annotations['Amino_Acid_Length'] = parts[4]
            if len(parts) > 5: annotations['Gene_name'] = parts[5]
            if len(parts) > 6: annotations['Biotype'] = parts[6]
            
    except Exception:
        pass 
        
    return annotations

def main(vcf_filename):
    if not os.path.exists(vcf_filename):
        print(f"Error: {vcf_filename} not found.")
        return

    with open(vcf_filename, 'r') as f:
        samples = []
        
        # Removed "Noise" from fixed headers
        fixed_headers_start = ["CHROM", "POS", "REF", "ALT", "TYPE"]
        fixed_headers_end = ["Effect", "Impact", "Functional_Class", "Codon_change", 
                             "Protein_and_nucleotide_change", "Amino_Acid_Length", 
                             "Gene_name", "Biotype"]
        
        for line in f:
            if line.startswith('##'):
                continue
            
            # Header Line
            if line.startswith('#'):
                header_parts = line.strip().split('\t')
                samples = header_parts[9:]
                
                # Print Full Header
                full_header = fixed_headers_start + samples + fixed_headers_end
                print("\t".join(full_header))
                continue

            # Data Lines
            parts = line.strip().split('\t')
            
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt_str = parts[4]
            alts = alt_str.split(',')
            info = parts[7]
            format_col = parts[8]
            
            # 1. Determine Type
            var_type = get_type(ref, alts)
            
            # 2. Parse SnpEff
            snpeff = parse_snpeff(info)
            
            # 3. Process Samples
            sample_percentages = []
            
            try:
                format_keys = format_col.split(':')
                ad_index = format_keys.index('AD')
            except ValueError:
                ad_index = -1

            for i, sample_name in enumerate(samples):
                if ad_index == -1:
                    sample_percentages.append(".")
                    continue
                
                sample_data = parts[9 + i]
                if sample_data.startswith('./.') or sample_data == '.':
                    sample_percentages.append(".")
                    continue

                vals = sample_data.split(':')
                if len(vals) <= ad_index:
                    sample_percentages.append(".")
                    continue
                
                ad_val = vals[ad_index]
                if ad_val == '.' or not ad_val:
                    sample_percentages.append(".")
                    continue
                
                # Calculate Percentage
                try:
                    reads = [int(x) for x in ad_val.split(',')]
                    total = sum(reads)
                    
                    if total == 0:
                        sample_percentages.append("0.00%")
                    else:
                        # reads[0] is Ref. reads[1:] are Alts.
                        alt_reads = reads[1:]
                        alt_percents = []
                        for count in alt_reads:
                            p = (count / total) * 100
                            alt_percents.append(f"{p:.2f}%")
                        
                        sample_percentages.append(",".join(alt_percents))

                except ValueError:
                    sample_percentages.append(".")

            # 4. Construct Row (Removed "FALSE" / Noise column)
            row = [chrom, pos, ref, alt_str, var_type]
            row.extend(sample_percentages)
            row.extend([
                snpeff['Effect'],
                snpeff['Impact'],
                snpeff['Functional_Class'],
                snpeff['Codon_change'],
                snpeff['Protein_and_nucleotide_change'],
                snpeff['Amino_Acid_Length'],
                snpeff['Gene_name'],
                snpeff['Biotype']
            ])
            
            print("\t".join(row))

if __name__ == "__main__":
    main("out.annotated.vcf")