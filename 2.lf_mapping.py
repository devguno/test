def read_bwt_sequence(file_path):
    with open(file_path, 'r') as file:
        return file.read().strip()

def create_sorted_bwt(bwt):
    return ''.join(sorted(bwt))

def get_char_positions(sequence):
    char_positions = {}
    for i, char in enumerate(sequence):
        if char not in char_positions:
            char_positions[char] = []
        char_positions[char].append(i)
    return char_positions

def reverse_bwt(bwt):
    n = len(bwt)
    sorted_bwt = create_sorted_bwt(bwt)
    
    bwt_char_positions = get_char_positions(bwt)
    sorted_char_positions = get_char_positions(sorted_bwt)
    
    lf_mapping = []
    for char in bwt:
        pos = bwt_char_positions[char].pop(0)
        lf_mapping.append(sorted_char_positions[char].pop(0))
    
    result = [''] * n
    current_index = lf_mapping[0]  
    
    for i in range(n - 1, -1, -1):
        result[i] = sorted_bwt[current_index]
        current_index = lf_mapping[current_index]
    
    return ''.join(result[1:]) + '$'  

def main():
    bwt_sequence = read_bwt_sequence('BWT.txt')
    original_sequence = reverse_bwt(bwt_sequence)
    
    print(f"Reconstructed sequence: {original_sequence}")
    print(f"Last 10 bps: {original_sequence[-11:-1]}")  

if __name__ == "__main__":
    main()