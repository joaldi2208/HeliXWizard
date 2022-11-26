
import pandas as pd
import re
from functools import partial
import numpy as np

# write functions for one entry and use the lambda apply function
def remove_spaces(entry):
    return entry.replace(" ","")

def eliminate_redundancy(entry):
    values = entry.split(";")
    unique_values = ",".join(set(values))
    return unique_values

def uniform_pressure(sep_val_unit_rex, pressure):
    """hnt"""
    def determine_unit(sep_val_unit_rex, pressure):
        units = ["TORR", "ATM", "PA", "MMHG", "BAR", "MBAR"]
        match = sep_val_unit_rex.search(pressure)
        return match.groups()
    def remove_words(unit):
        words_to_remove = ["AMBIENT", "ATMOSPHERIC", "VARIABLE"]
        for word in words_to_remove:
             unit = unit.replace(word, "")
        #if unit == "MBARMBAR":
        #    unit = unit.replace("MBARMBAR", "MBAR")
        #elif unit == "PAATM":
        #    unit = unit.replace("PAATM", "PA")
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
    value, unit = determine_unit(sep_val_unit_rex, pressure)
    if value != None and unit != None:
        only_unit = remove_words(unit)
        p_in_pa = pressure2pascal(value, only_unit)
        #print(p_in_pa)
        if p_in_pa < 100000 or p_in_pa > 101325:
            print(p_in_pa, value, unit)
    return value, unit

def replace_words():
    pass

def pressure2pascal():
    pass
    

def calc_mean_value(entry):
    try:
        value1, value2 = entry.split("-")
        mean_value = ( float(value1) + float(value2) ) / 2
        return mean_value
    except:
        return entry
        
    

if __name__== "__main__":
    sep_val_unit_rex = re.compile("([\d.]+)?([A-Z]+)?")
    uniform_pressure_rex = partial(uniform_pressure, sep_val_unit_rex)
    #chemical_shifts = pd.read_pickle("chemical-shifts.pkl")
    struc_and_cond = pd.read_pickle("struc-and-cond.pkl")
    struc_and_cond = struc_and_cond.fillna(value="NONE")
    parameters = ["Temp", "Pressure", "PH", "Ionic_strength"]
    for para in parameters:
        para_no_spaces = struc_and_cond[para].apply(lambda x: remove_spaces(x))
        para_no_double = para_no_spaces.apply(lambda x: eliminate_redundancy(x))
        print(para_no_double)
    #pressure_in_pa = para_no_double.apply(lambda x: uniform_pressure_rex(x))
    
