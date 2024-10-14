import math
import matplotlib.pyplot as plt

# Read VCF file
def read_vcf(file_path):
    snps = []
    positions = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                data = line.strip().split('\t')
                chrom = data[0]
                pos = int(data[1])
                genotypes = [x.split(':')[0] for x in data[9:]]
                snps.append(genotypes)
                positions.append(pos)
    return snps, positions

# Calculate allele frequencies and genotype counts
def allele_and_genotype_counts(genotypes):
    n_00, n_01, n_11 = 0, 0, 0
    for g in genotypes:
        if g == '0/0':
            n_00 += 1
        elif g == '0/1':
            n_01 += 1
        elif g == '1/1':
            n_11 += 1
    n_total = n_00 + n_01 + n_11
    p = (2 * n_00 + n_01) / (2 * n_total)  # frequency of reference allele
    q = 1 - p  # frequency of alternative allele
    return n_00, n_01, n_11, p, q, n_total

# Perform chi-square test for HWE
def hwe_chi_square_test(n_00, n_01, n_11, p, q, n_total):
    # Expected counts under HWE
    e_00 = n_total * p * p
    e_01 = 2 * n_total * p * q
    e_11 = n_total * q * q

    # Chi-square test statistic
    chi_square = ((n_00 - e_00) ** 2 / e_00) + ((n_01 - e_01) ** 2 / e_01) + ((n_11 - e_11) ** 2 / e_11)
    return chi_square

# Convert chi-square to p-value
def chi_square_to_pvalue(chi_square):
    # Degrees of freedom = 1 for HWE test
    p_value = math.exp(-0.5 * chi_square) * (1 - chi_square / 2)
    return p_value

# Main function to process SNPs
def process_snps(snps):
    p_values = []
    significant_snp_count = 0
    chi_square_threshold = 10**(-5)

    for genotypes in snps:
        n_00, n_01, n_11, p, q, n_total = allele_and_genotype_counts(genotypes)
        if n_total == 0:
            continue
        chi_square = hwe_chi_square_test(n_00, n_01, n_11, p, q, n_total)
        p_value = chi_square_to_pvalue(chi_square)
        p_values.append(p_value)
        if p_value < chi_square_threshold:
            significant_snp_count += 1

    return p_values, significant_snp_count

# Plot histogram of p-values
def plot_histogram(p_values):
    plt.hist(p_values, bins=50, edgecolor='black')
    plt.xlabel('p-value')
    plt.ylabel('Frequency')
    plt.title('Histogram of p-values')
    plt.show()

# Plot Manhattan plot of -log10 p-values
def plot_manhattan(p_values, positions):
    log_p_values = [-math.log10(p) if p > 0 else 0 for p in p_values]
    plt.scatter(positions, log_p_values, c='blue', s=10)
    plt.xlabel('Genomic Position')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot')
    plt.show()

# Path to your VCF file
file_path = 'EAS_chr22.vcf'

# Run the analysis
snps, positions = read_vcf(file_path)
p_values, significant_snp_count = process_snps(snps)

# Print the result
print(f'Number of SNPs with HWE p-value < 10^-5: {significant_snp_count}')

# Plot results
plot_histogram(p_values)
plot_manhattan(p_values, positions)
