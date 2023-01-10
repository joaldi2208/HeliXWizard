import pandas as pd

import pycurl
import io
import certifi

import requests

from collections import Counter

import re

import json
import pickle

from tqdm import tqdm

import argparse


def read_in_db_ids()-> pd.DataFrame:
    """Reads the 'BMRB_PDB_ids.csv' file from the same directory. It shows BMRB ids and there corresponding PDB ids. The file can be found here: https://bmrb.io/search/ under 'Matched submitted BMRB-PDB entries'. """
    
    bmrb_pdb_ids = pd.read_csv("BMRB_PDB_ids.csv",
                               names=["bmrb_ids", "pdb_ids"],
                               encoding='latin-1')

    # There a possibly more than one entry for one bmrb_id, but it seems that the first one
    # corresponds to the right pdb_id mentioned on the top of the bmrb webpage
    bmrb_pdb_ids.drop_duplicates(subset="bmrb_ids", keep="first", inplace=True, ignore_index=True)

    return bmrb_pdb_ids


def download_webpages(url: str) -> str: 
    """Downloads the N-HSQC spectra from from the BMRB database."""
    
    buffer = io.BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.perform()
    c.close()
    
    body = buffer.getvalue() 
    decoded_body = body.decode("iso-8859-1")
    
    return decoded_body


def give_corrsponding_ids(bmrb_pdb_ids: pd.DataFrame, query_ids: list, ref_ids="bmrb_ids", target_ids="pdb_ids")-> list:
    """Returns the corresponding ids to either a list of PDB or BMRB ids"""

    corresponding_ids = []

    for index, ref_id in enumerate(bmrb_pdb_ids[ref_ids]):
        if ref_id in query_ids:
            corresponding_id = bmrb_pdb_ids[target_ids][index]
            corresponding_ids.append(corresponding_id)
            
    return corresponding_ids




###############################################
#       CHEMICAL SHIFT INFORMATION
##############################################


def get_N_HSQC_spectra(bmrb_pdb_ids: pd.DataFrame)-> dict:
    """build an url based on the given BMRB ids from the dataframe and calls a function to download this data. The downloaded data is then stored as dataframe in a dictionary with the BMRB ids as keys"""

    url_start = "https://api.bmrb.io/current/entry/"
    url_end = "/simulate_hsqc?format=csv&filter=backbone"

    chemical_shifts = {}
    bmrb_ids = bmrb_pdb_ids.bmrb_ids
    for index, bmrb_id in tqdm(enumerate(bmrb_ids), total=len(bmrb_ids), desc="Downloading N-HSQC spectra..."):
        
        url = url_start + str(bmrb_id) + url_end
        N_HSQC_spectra = download_webpages(url)
        
        if N_HSQC_spectra[0:3] == "seq":
            N_HSQC_spectra_IOobject = io.StringIO(N_HSQC_spectra)
            N_HSQC_spectra_table = pd.read_csv(N_HSQC_spectra_IOobject)
            
            chemical_shifts[bmrb_id] = N_HSQC_spectra_table
            
    return chemical_shifts




###############################################
#       SECONDARY STRUCTURE INFORMATION
##############################################


def get_sec_struc_counts(sec_structure_8_state: str, sec_structure_3_state:str)-> list:
    """counts the appearance of symbols in secondary structure string"""
    
    counter_8_state = Counter(sec_structure_8_state)
    avail_8_sec_structure = ["H","B","E","G","I","T","S"]
    counts_in_8_state = [ counter_8_state[sym] for sym in avail_8_sec_structure ]

    counter_3_state = Counter(sec_structure_3_state)
    avail_3_sec_structure = ["C","E","H"]
    counts_in_3_state = [ counter_3_state[sym] for sym in avail_3_sec_structure ]
    
    return tuple(counts_in_8_state), tuple(counts_in_3_state)



def download_DSSP(pdb_ids: list)-> list:
    """downloads dssp files from 'https://dssp.bellstedt-lab.ch' and extract amino acid sequence nad secondary structure information"""

    dssp_files = []
    
    for pdb_id in tqdm(pdb_ids, total=len(pdb_ids), desc="Downloading DSSP..."):
        
        pdb_id = pdb_id.lower()
        dssp = requests.get(f"https://dssp.bellstedt-lab.ch/get/{pdb_id}")
        dssp_file = dssp.text

        dssp_files.append(dssp_file)

    return dssp_files


def get_secondary_structure(bmrb_ids: list, dssp_files: list)-> dict:
    """reads secondary structure information from the downloaded DSSP files."""

    secondary_structure = {}
    
    for bmrb_id, dssp_file in zip(bmrb_ids, dssp_files):

        lines = dssp_file.split("\n")

        if len(lines) >= 9:
            
            AA_sequence = lines[4]
            sec_structure_8_state = lines[6]
            sec_structure_3_state = lines[8]

            counts_in_8_state, counts_in_3_state = get_sec_struc_counts(sec_structure_8_state, sec_structure_3_state)

            secondary_structure_summary = [("AA_sequence", AA_sequence),
                                           ("sec_structure_8_state", sec_structure_8_state),
                                           ("counts_in_8_state", counts_in_8_state),
                                           ("sec_structure_3_state", sec_structure_3_state),
                                           ("counts_in_3_state", counts_in_3_state)]
        
            secondary_structure[bmrb_id] = dict(secondary_structure_summary)
    
    return secondary_structure




###############################################
#       MEASUREMENT CONDITION INFORMATION
##############################################

def find_conditions(webpage, database, bmrb_id=None):
    """finds the measurement conditions on the BMRB webpage based on a given condition id or on the PDB header file. The given database name decides which regex to use"""

    if database == "BMRB":
        if "Page not found" in webpage:
            webpage = download_webpages(f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str")
            
        condition_id = re.search(r"condition.*?_label.*? \$([A-Za-z0-9_.-]*)", webpage).group(1)
        condition_lines = re.findall(rf"save_{condition_id}([0-9a-zA-Z, \n\r\t_./';:\"=\]\[<~>*?()%+#-]*)(?=save)", webpage)[0]
        conditions = re.findall(r"\n *'?(temperature|pH|pressure|ionic strength)'? *([A-Z0-9.]*)", condition_lines)
        
        return conditions
       

    elif database == "PDB":
        conditions = re.findall(r"(TEMPERATURE|PH|PRESSURE|IONIC STRENGTH).*?: (.*)", webpage)
        return conditions

    
def get_bmrb_measurement_conditions(bmrb_ids):
    """copies the measurement conditions from the BMRB webpage"""

    measurement_conditions = {}
    
    for bmrb_id in tqdm(bmrb_ids, total=len(bmrb_ids), desc="Downloading BMRB conditions..."):

        bmrb_webpage = download_webpages(f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_21.str")
        
        conditions = find_conditions(bmrb_webpage, "BMRB", bmrb_id)
        conditions_lowercase = [(parameter.lower(), value) for parameter, value in conditions]

        measurement_conditions[bmrb_id] = dict(conditions_lowercase)

    return measurement_conditions
        

def get_pdb_measurement_conditions(pdb_ids, bmrb_ids, measurement_conditions):
    """copies the measurement conditions from the PDB file header and adds it to the existing measurement conditions"""

    for bmrb_id, pdb_id in tqdm(zip(bmrb_ids, pdb_ids), total=len(bmrb_ids), desc="Downloading PDB conditions..."):

        pdb_webpage = download_webpages(f"https://files.rcsb.org/header/{pdb_id}.pdb")

        conditions = find_conditions(pdb_webpage, "PDB")

        for condition in conditions:
            parameter = condition[0].lower()
            value = condition[1].strip()

            measurement_conditions[str(bmrb_id)].setdefault(parameter, value)

            if measurement_conditions[str(bmrb_id)][parameter] == "" or measurement_conditions[str(bmrb_id)][parameter] == ".":
                measurement_conditions[str(bmrb_id)].update({parameter: value})

    return measurement_conditions
        
    
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--do", type=str, nargs="+")
    
    args = parser.parse_args().do
    
    bmrb_pdb_ids = read_in_db_ids()

    
    if "shift" in args:
        
        chemical_shifts = get_N_HSQC_spectra(bmrb_pdb_ids)
        with open("chemical_shifts.pkl", "wb") as outfile:
            pickle.dump(chemical_shifts, outfile)
    
        found_bmrb_ids = chemical_shifts.keys()
        with open("found_bmrb_ids.txt", "w") as outfile:
            for bmrb_id in found_bmrb_ids:
                outfile.write(f"{bmrb_id},")
                # delete last comma or
                #replace all of it with the expressions used in the filter scripts


    elif "shift" not in args:
        
        with open("found_bmrb_ids.txt", "r") as infile:
            found_bmrb_ids_strings = infile.read().replace("\n","").split(",")
            found_bmrb_ids = [int(bmrb_id) for bmrb_id in found_bmrb_ids_strings]

            
    found_pdb_ids = give_corrsponding_ids(bmrb_pdb_ids,found_bmrb_ids)

    
    if "dssp" in args:
        
        dssp_files = download_DSSP(found_pdb_ids)
        secondary_structure = get_secondary_structure(found_bmrb_ids, dssp_files)
        with open("secondary_structure.json", "w") as outfile:
            json.dump(secondary_structure, outfile, indent=4, separators=(", ", ": "), sort_keys=True)
            

    if "bmrb" in args:

        measurement_conditions = get_bmrb_measurement_conditions(found_bmrb_ids)
        with open("measurement_conditions.json", "w") as outfile:
            json.dump(measurement_conditions, outfile, indent=4, separators=(", ", ": "), sort_keys=True)

            
    if "pdb" in args:

        with open("measurement_conditions.json", "r") as infile:
            measurement_conditions = json.load(infile)

        updated_measurement_conditions = get_pdb_measurement_conditions(found_pdb_ids, found_bmrb_ids, measurement_conditions)
        with open("updated_measurement_conditions.json", "w") as outfile:
            json.dump(updated_measurement_conditions, outfile, indent=4, separators=(", ", ": "), sort_keys=True)
