from Binning_Data import get_prefiltered_ids
import json
import numpy as np

def calc_Q_states_percentage(Q_state):
    """calculates the percentage appearance of a secondary structure elements in a protein"""
    
    sec_struc_elements_percent = []
    
    for bmrb_id, sec_struc_info in secondary_structures.items():
        sec_struc_elements_count = sec_struc_info[f"counts_in_{Q_state}"]
        sec_struc_elements_total = len(sec_struc_info[f"sec_structure_{Q_state}"])

        
        elements_percentage = []
        
        for element_count in sec_struc_elements_count:
            element_percentage = round(element_count / sec_struc_elements_total, 2) # maybe 3 is better?
            elements_percentage.append(element_percentage)
            
            
        sec_struc_elements_percent.append(np.array(elements_percentage))
            
    return sec_struc_elements_percent


if __name__=='__main__':

    with open("secondary_structure.json", "r") as infile:
        unfiltered_secondary_structures = {int(bmrb_id): sec_struc_info for bmrb_id, sec_struc_info in json.load(infile).items()}
        secondary_structures = get_prefiltered_ids(unfiltered_secondary_structures)

    Q3_percentages = calc_Q_states_percentage("3_state")
    #with open("Q3_percentage.txt", "w") as outfile:
    #    outfile.writelines("\n".join(map(str, Q3_percentages)))

    with open("Q3_percentage.npy", "wb") as outfile:
        np.save(outfile, Q3_percentages)
        
    Q8_percentages = calc_Q_states_percentage("8_state")
    #with open("Q8_percentage.txt", "w") as outfile:
    #    outfile.writelines("\n".join(map(str, Q8_percentages)))
      
    with open("Q8_percentage.npy", "wb") as outfile:
        np.save(outfile, Q8_percentages)
