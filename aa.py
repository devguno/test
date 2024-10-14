import math

def read_vcf(file_path):
    snps = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, id, ref, alt, qual, filter, info, format = parts[:9]
            genotypes = parts[9:]
            snps.append({
                'pos': int(pos),
                'id': id,
                'genotypes': genotypes
            })
    return snps

def count_genotypes(genotypes):
    counts = {'0/0': 0, '0/1': 0, '1/1': 0}
    for gt in genotypes:
        gt = gt.split(':')[0]
        if gt in counts:
            counts[gt] += 1
    return counts

def calculate_allele_frequencies(counts):
    total = sum(counts.values()) * 2
    p = (counts['0/0'] * 2 + counts['0/1']) / total
    q = 1 - p
    return p, q

def chi_square_test(observed, expected):
    chi_square = sum((o - e) ** 2 / e for o, e in zip(observed, expected) if e != 0)
    return chi_square

def hwe_test(counts):
    n = sum(counts.values())
    p, q = calculate_allele_frequencies(counts)
    
    expected = {
        '0/0': n * (p ** 2),
        '0/1': n * 2 * p * q,
        '1/1': n * (q ** 2)
    }
    
    chi_square = chi_square_test(list(counts.values()), list(expected.values()))
    
    # Degrees of freedom for HWE test is 1
    p_value = 1 - chi_square_cdf(chi_square, 1)
    
    return p_value

def chi_square_cdf(x, k):
    # This is a simplified approximation of the chi-square CDF
    # It's not very accurate for small k or extreme x values
    if x < 0 or k <= 0:
        return 0
    return 1 - math.exp(-x/2) * sum((x/2)**i / math.factorial(i) for i in range(k//2))

def analyze_vcf(file_path):
    snps = read_vcf(file_path)
    p_values = []
    significant_count = 0
    
    for snp in snps:
        counts = count_genotypes(snp['genotypes'])
        p_value = hwe_test(counts)
        p_values.append((snp['pos'], p_value))
        
        if p_value < 1e-5:
            significant_count += 1
    
    return p_values, significant_count

def create_histogram(values, bins=20):
    min_val, max_val = min(values), max(values)
    if min_val == max_val:
        print(f"All values are equal to {min_val}")
        return
    
    bin_size = (max_val - min_val) / bins
    hist = [0] * bins
    
    for value in values:
        if bin_size > 0:
            bin_index = min(int((value - min_val) / bin_size), bins - 1)
        else:
            bin_index = 0
        hist[bin_index] += 1
    
    max_count = max(hist)
    for i, count in enumerate(hist):
        bar = '#' * int(50 * count / max_count) if max_count > 0 else ''
        print(f"{min_val + i*bin_size:.2e} - {min_val + (i+1)*bin_size:.2e}: {bar}")

def create_manhattan_plot(positions, p_values, window_size=50):
    sorted_data = sorted(zip(positions, p_values))
    max_pos = max(positions)
    
    for i in range(0, max_pos, window_size):
        window_data = [(-math.log10(p) if p > 0 else 0) for pos, p in sorted_data if i <= pos < i + window_size]
        if window_data:
            avg_log_p = sum(window_data) / len(window_data)
            bar = '#' * int(50 * avg_log_p / 10) if avg_log_p > 0 else ''
            print(f"{i}: {bar}")

# Main execution
file_path = 'EAS_chr22.vcf'
p_values, significant_count = analyze_vcf(file_path)

print(f"Number of SNPs with HWE p-value < 1e-5: {significant_count}")

# Create histogram of p-values
print("\nHistogram of p-values:")
create_histogram([p for _, p in p_values])

# Create Manhattan plot
print("\nManhattan plot of -log10 p-values:")
create_manhattan_plot([pos for pos, _ in p_values], [p for _, p in p_values])

# Additional analysis
total_snps = len(p_values)
print(f"\nTotal number of SNPs analyzed: {total_snps}")
print(f"Proportion of SNPs with p-value < 1e-5: {significant_count / total_snps:.2%}")

min_p = min(p for _, p in p_values)
max_p = max(p for _, p in p_values)
print(f"\nMinimum p-value: {min_p:.2e}")
print(f"Maximum p-value: {max_p:.2e}")

# Check for any patterns or clusters
print("\nChecking for patterns or clusters:")
threshold = 1e-5
cluster_size = 5
cluster_count = 0
for i in range(len(p_values) - cluster_size + 1):
    if all(p < threshold for _, p in p_values[i:i+cluster_size]):
        cluster_count += 1
print(f"Number of clusters with {cluster_size} consecutive SNPs below p-value threshold: {cluster_count}")