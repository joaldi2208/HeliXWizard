
import pandas as pd
import re
from functools import partial
import numpy as np




###########################################
#        GENERAL FUNCTIONS
###########################################

def remove_spaces(entry):
    return entry.replace(" ","")

def eliminate_redundancy(entry):
    values = entry.split(";")
    unique_values = ",".join(set(values))
    return unique_values



###########################################
#       PRESSURE TRANSFORMATION
###########################################

def uniform_pressure(sep_val_unit_rex, pressure):
    """hnt"""
    value, unit = determine_unit(sep_val_unit_rex, pressure)
    if value != None and unit != None:
        only_unit = remove_words(unit)
        p_in_pa = pressure2pascal(value, only_unit)
        #print(p_in_pa)
        #if p_in_pa < 100000 or p_in_pa > 101325:
            #print(p_in_pa, value, unit)
    return value, unit

def determine_unit(sep_val_unit_rex, pressure):
    units = ["TORR", "ATM", "PA", "MMHG", "BAR", "MBAR"]
    match = sep_val_unit_rex.search(pressure)
    return match.groups()

def remove_words(unit):
    words_to_remove = ["AMBIENT", "ATMOSPHERIC", "VARIABLE"] # check variable value in literature
    for word in words_to_remove:
        unit = unit.replace(word, "101.325")
    if unit == "MBARMBAR":
        unit = unit.replace("MBARMBAR", "MBAR")
    elif unit == "PAATM":
        unit = unit.replace("PAATM", "PA")
    return unit

def pressure2pascal(value, unit):
    if unit == "MBAR":
       return float(value) * 100
    elif unit == "BAR":
       return float(value) * 100000
    elif unit == "ATM":
       return float(value) * 101325
    elif unit == "MMHG":
       return float(value) * 133.322
    elif unit == "TORR":
       return float(value) * 133.322
    elif unit == "PA":
       return float(value)
    else:
      return -1          

def replace_words():
    pass

    

        
###########################################
#       PH, TEMPERATURE TRANSFORMATION
###########################################
    
def calc_mean_values(entry):
    splitted_entries = entry.split(",")
    while "NULL" in splitted_entries:
        splitted_entries.remove("NULL")
    if "-" in entry:
        entries = []
        for value in splitted_entries:
            if "-" in value:
                mean_range = calc_mean_range(value)
                entries.append(mean_range)            
            else:  	        
                entries.append(value)
        return np.mean(list(map(float, entries)))
    elif len(splitted_entries) > 1:
        return np.mean(list(map(float, splitted_entries)))
    else:
        return entry
	
	    	 
def calc_mean_range(entry):
    value1, value2 = entry.split("-")
    mean_range = ( float(value1) + float(value2) ) / 2
    return mean_range

def convert2float(entry):
    try:
        float_entry = float(entry)
        return float_entry
    except ValueError:
        mean_values = calc_mean_values(entry)
        return mean_values




if __name__== "__main__":
    sep_val_unit_rex = re.compile("([\d.]+)?([A-Z]+)?")
    uniform_pressure_rex = partial(uniform_pressure, sep_val_unit_rex)
    #chemical_shifts = pd.read_pickle("chemical-shifts.pkl.gz",compression="gzip")
    struc_and_cond = pd.read_pickle("struc-and-cond.pkl.gz",compression="gzip")
    struc_and_cond = struc_and_cond.fillna(value="NONE")
    #parameters = ["Temp", "Pressure", "PH", "Ionic_strength"]
    parameters = ["Pressure"]
    for para in parameters:
        para_no_spaces = struc_and_cond[para].apply(lambda x: remove_spaces(x))
        para_no_double = para_no_spaces.apply(lambda x: eliminate_redundancy(x))
        print(para_no_double)
    pressure_in_pa = para_no_double.apply(lambda x: uniform_pressure_rex(x))
    #print(pressure_in_pa)
    for i in pressure_in_pa:
        print(i)
    #entries_converted = para_no_double.apply(lambda x: convert2float(x))
    #for i in entries_converted:
    #    print(i)
