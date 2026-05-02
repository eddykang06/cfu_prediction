# Functions for splitting train/test set according to specified metrics

def drug_reader(x):

    # Single drug
    if len(x) < 11:
        
        # If NDC
        if x[:3] == "NDC":
            return "NDC"
        
        # If another drug
        else:
            
            # Construct of array of 1s and 0s depending on if character
            binary = np.array([int(char.isalpha()) for char in x])

            # Find first 1, then grap drug name
            first_char = np.where(binary == 1)[0][0]
            return x[first_char:first_char+3]
    
    # Multiple drug
    else:

        # Extract drug
        drug_name =  x[2:5] + "+" + x[7:10]
        return drug_name
    
# Parameters : cond_list (list of condition labels, ex: "12CIP1hr-a", "12CIP1hr-b", etc)
# Output     : groups (list of labels according to condition, ex: 0, 0, 0, 1, 1, 1, etc.)
def groups_by_cond(cond_list):

    # Initialize list to store groups
    groups = [0]

    # Loop through conditions, beginning at second position
    for i in range(1, len(cond_list)):
        
        # Previous condition
        prev = cond_list[i-1]   
        curr = cond_list[i]

        # Matching label (apart from last letter)
        if prev[:-1] == curr[:-1]:
            groups.append(groups[i-1])
        
        # New label
        else:
            groups.append(groups[i-1] + 1)
    
    return groups



def train_test_split()