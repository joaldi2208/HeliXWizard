import json
import re

def pressure_exceptions(entry):
    """cares about pressure annotation exeptions"""
    written_pressure_terms = ["AMBIENT", "AMBINET", "ATMOSPHERE", "ATMOSPHERIC"]
    for term in written_pressure_terms:
        if term in entry:
            entry = "1"
    return entry


def pressure_units(entry):
    """removes units and converts values to bar"""
    units_to_remove = ["ATM", "BAR"]
    for unit in units_to_remove:
        if unit in entry:
            entry = entry.replace(unit, "")
    if "VARIABLE" in entry:
        entry = "-1"
    if "PA" in entry:
        entry_in_pa = re.search(r"([0-9]*)", entry).group(1)
        entry_in_bar = float(entry_in_pa) / 101300
        return str(entry_in_bar)

    if entry == "":
        entry = "-1"

    return entry
    

def check_multiple_entries(entry):
    """hnt"""
    if ";" in entry:
        multiple_entries = [e.strip() for e in entry.split(";")]

        # add here what is checked
        if multiple_entries[-1] == "":
            multiple_entries = multiple_entries[:-1]
            
        if multiple_entries[0] == "":
            multiple_entries = multiple_entries[1:]

        if len(set(multiple_entries)) == 1:

            return multiple_entries[0]
        else:
            return entry
    else:
        return entry
            
  

def NULL2minus1(measurement_conditions):
    """replace 'NULL' with '-1'"""
    for bmrb_id, condition_type in measurement_conditions.items():

        for condition in condition_type:
        
            if condition == "pressure":
                condition_type[condition] = pressure_exceptions(condition_type[condition])
                condition_type[condition] = pressure_units(condition_type[condition])
                
            condition_type[condition] = check_multiple_entries(condition_type[condition])
            condition_type[condition] = condition_type[condition].replace("NULL", "-1")

    return measurement_conditions




def filter_by_condition(measurement_conditions, condition_range, selected_condition):
    """filters measurement conditions by the given condition and condition range"""
    
    min_val = condition_range[0]
    max_val = condition_range[1]

    out_range = []
    in_range = []

    for bmrb_id, condition_types in measurement_conditions.items():
        
        condition = condition_types[selected_condition]
        # add here why this is needed
        if ";" in condition:
            pass

        elif float(condition) >= min_val and float(condition) <= max_val:
            in_range.append(bmrb_id)

        else:
            out_range.append(bmrb_id)

    return in_range, out_range
    
def toLessConditions(measurement_conditions):
    bmrb_ids = []
    for bmrb_id, condition_types in measurement_conditions.items():
        if "pressure" not in condition_types or "temperature" not in condition_types or "ph" not in condition_types:
            bmrb_ids.append(bmrb_id)

    return bmrb_ids



if __name__ == '__main__':
    # only updated measurement conditions
    # but you have to check if the results have an dssp result as well???????????
    with open("updated_measurement_conditions.json", "r") as infile:
        measurement_conditions = json.load(infile)



    print(measurement_conditions["15645"])
    # define parameter ranges
    pressure_range = (0.9, 1.2)
    ph_range = (5, 8)
    temperature_range = (273, 313) 

    bmrb_no_cond = toLessConditions(measurement_conditions)
    new_measurement_conditions = measurement_conditions.copy()

    for bmrb_id in bmrb_no_cond:
        del new_measurement_conditions[bmrb_id]
        
    new_measurement_conditions = NULL2minus1(new_measurement_conditions)

    in_temperature_range, out_temperature_range = filter_by_condition(new_measurement_conditions, temperature_range, "temperature")
    print("number in temperature range: ", len(in_temperature_range))
    print("number not in temperature range: ",len(out_temperature_range))
    print("--------------------")
    
    in_pressure_range, out_pressure_range = filter_by_condition(new_measurement_conditions, pressure_range, "pressure")
    print("number in pressure range: ",len(in_pressure_range))
    print("number not in pressure range: ",len(out_pressure_range))
    print("--------------------")
    
    in_ph_range, out_ph_range = filter_by_condition(new_measurement_conditions, ph_range, "ph")
    print("number in ph range: ",len(in_ph_range))
    print("number not in ph range: ",len(out_ph_range))
    print("--------------------")
   
    overlapp = set(in_temperature_range) & set(in_pressure_range) & set(in_ph_range)
    print("number in overlapp range: ",len(overlapp))
    print("###################")

    not_in_list = out_ph_range + out_pressure_range + out_temperature_range + bmrb_no_cond
    not_in = set(out_ph_range) | set(out_pressure_range) | set(out_temperature_range) | set(bmrb_no_cond)
    print("number not all conditions: ",len(bmrb_no_cond))
    print("number not in in general: ",len(not_in_list))
    print("number not in in general set: ",len(not_in))
    print("--------------------")
    print()
    print("###################")
    print()

    not_in = set(out_ph_range) | set(out_pressure_range) | set(out_temperature_range) | set(bmrb_no_cond)

    #print(len(out_ph_range + out_pressure_range + out_temperature_range))

    with open("secondary_structure.json", "r") as infile:
        secondary_structure = json.load(infile)

    dssp = list(secondary_structure)
    print("dssp: ",len(dssp))
    with_dssp = set(dssp) & overlapp
    print("overlap with dssp: ", len(with_dssp))
    # write to a file



    
    #with open("not_inclueded_bmrbs.txt", "w") as outfile:
    #    groups = [out_ph_range, out_pressure_range, out_temperature_range, bmrb_no_cond]
    #    groups_name = ["out_ph_range", "out_pressure_range", "out_temperature_range", "bmrb_no_cond"]
    #    for group, group_name in zip(groups, groups_name):
    #        outfile.write(f"{group_name}\n")
    #        outfile.write("---------------------------\n")
    #        keys = list(group)
    #        for entry in keys:
    #             outfile.write(f"{entry}...{measurement_conditions[entry]}\n")
                
