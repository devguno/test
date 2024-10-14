import math
import matplotlib.pyplot as plt

def read_vcf(file_path):
    snps = []
    locations = []  # to store SNP locations
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  
            fields = line.strip().split('\t')
            locations.append(fields[1])  # SNP position is typically in the second column
            genotypes = [gt.split(':')[0] for gt in fields[9:]]  # extract genotype
            snps.append(genotypes)
    return snps, locations


def count_genotypes(genotypes):
    # 각 genotype(0/0, 0/1, 1/1)의 개수 세기
    counts = {'0/0': 0, '0/1': 0, '1/1': 0}
    for gt in genotypes:
        if gt in counts:
            counts[gt] += 1
    return counts['0/0'], counts['0/1'], counts['1/1']

def calculate_allele_frequency(aa, ab, bb):
    # allele frequency 계산
    n = aa + ab + bb
    p = (2 * aa + ab) / (2 * n)  
    q = 1 - p  
    return p, q

def expected_genotypes(n, p, q):
    # HWE에서 allele frequency 계산
    return n * p * p, 2 * n * p * q, n * q * q

def chi_square_test(observed, expected):
    # chi square test
    chi_square = sum((o - e) ** 2 / e for o, e in zip(observed, expected) if e != 0)
    df = 1 
    p_value = 1 - chi_square_cdf(chi_square, df)
    return p_value

def chi_square_cdf(x, df):
    return 1 - math.exp(-x/2) * sum((x/2)**i / math.factorial(i) for i in range(df//2))

def plot_histogram(p_values):
    plt.figure(figsize=(10, 6))
    plt.hist(p_values, bins=50, edgecolor='black')
    plt.title('Histogram of HWE p-values')
    plt.xlabel('p-value')
    plt.ylabel('Frequency')
    plt.savefig('hwe_pvalue_histogram.png')
    plt.close()

def plot_manhattan(locations, p_values):
    plt.figure(figsize=(12, 6))
    plt.scatter(locations, [-math.log10(p) for p in p_values], alpha=0.6)
    plt.title('Manhattan Plot of HWE Test Results')
    plt.xlabel('Location')
    plt.ylabel('-log10(p-value)')
    plt.savefig('hwe_manhattan_plot.png')
    plt.close()

def main():
    file_path = 'EAS_chr22.vcf'
    snps = read_vcf(file_path)  # Only one value to unpack
    
    p_values = []
    significant_snps = 0
    for genotypes in snps:
        aa, ab, bb = count_genotypes(genotypes)
        n = aa + ab + bb
        p, q = calculate_allele_frequency(aa, ab, bb)
        expected = expected_genotypes(n, p, q)
        observed = (aa, ab, bb)
        
        p_value = chi_square_test(observed, expected)
        p_values.append(p_value)
        
        if p_value < 1e-5:
            significant_snps += 1
    
    print(f"Number of SNPs with HWE p-value < 10^-5: {significant_snps}")

    plot_histogram(p_values)


if __name__ == "__main__":
    main()