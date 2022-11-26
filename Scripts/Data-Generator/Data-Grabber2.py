# before using the script, download the data with the following command 
# rsync -rlpt -v -z --delete --port=33444 \
# rsync.rcsb.org::ftp_data/structures/divided/nmr_chemical_shifts/
# into a folder with a relativ path to this file: ../../Data


#H = alpha-helix
#B = beta-bridge residue
#E = extended strand (in beta ladder)
#G = 3/10-helix
#I = 5-helix
#T = H-bonded turn
#S = bend
 
import os
import gzip
import re
import pandas as pd
import io
import pycurl
import certifi
import requests
import json
import time
from functools import partial
import numpy as np
from collections import Counter
from ordered_set import OrderedSet


def get_paths():
    """paths"""
    path2files = []
    for root, _, filenames in os.walk("../../Data/"):
        path2files.extend([os.path.join(root, name) for name in filenames])
    return path2files


def get_data(path2files, get_pdb_id_rex, get_chemical_shifts_rex, get_sec_struc_rex):
    """get data"""
    files_not_found = []
    df_per_id = []
    cond_struc_depot = []
    for i, files in enumerate(path2files):
        try:
            with gzip.open(files, "rt") as f:
                if i in list(range(50,6000,50)):
                    print(i)
                nmrstar_file = f.read()
                chemical_shifts_short = get_chemical_shifts_rex(nmrstar_file)
                if type(chemical_shifts_short) != str:
                    pdb_id = get_pdb_id_rex(nmrstar_file)
                    if "Error" not in pdb_id:
                        # mach die pdb_id doch uppercase
                        # du solltest auch schauen anhand eines bsp ob das mit duplicate drop klappt
                        # steht im pdb file nach header
                        cond = get_cond(pdb_id)
                        if "Error" not in cond:
                            counts_sec_struc, seq_length = get_sec_struc_rex(pdb_id)
                            if "Error" not in counts_sec_struc: 
                                cond_and_struc = [pdb_id] + counts_sec_struc + list(cond.values()) + [seq_length]
                                cond_struc_depot.append(cond_and_struc)
                                df_per_id.append(chemical_shifts_short)
                            else:
                                files_not_found.append(counts_sec_struc[0])
                        else:
                            files_not_found.append(cond)
                    else:
                        files_not_found.append(pdb_id)
                else:
                    files_not_found.append(chemical_shifts_short)
        except FileNotFoundError:
            files_not_found.append("Path Error")
    cond_struc_table = from_list2table(cond_struc_depot)
    return pd.concat(df_per_id), cond_struc_table, files_not_found



def get_chemical_shifts(headers_rex, nmrstar_file):
    """ids"""
    #def is_protein(chemical_shifts):
    #    """RNA has a one letter id while proteins have three"""
    #    amino_acid_code = [
    #        "ALA",
    #        "ASP",
    #        "CYS",
    #        "GLU",
    #        "PHE",
    #        "GLY",
    #        "HIS",
    #        "ILE",
    #        "LYS",
    #        "LEU",
    #        "MET",
    #        "ASN",
    #        "PRO",
    #        "GLN",
    #        "ARG",
    #        "SER",
    #        "THR",
    #        "VAL",
    #        "TRP",
    #        "TYR"
    #        ]
    #    if chemical_shifts["Comp_ID"][0] in amino_acid_code:
    #        return True
        
    def create_dataframe(last_header, nmrstar_file, n):
        chemical_shifts_range = nmrstar_file[(nmrstar_file.find(last_header, nmrstar_file.find(last_header)+n)+len(last_header)) : ]
        chemical_shifts_range = chemical_shifts_range[: (chemical_shifts_range.find("stop_"))]
        stream = io.StringIO()
        stream.write(chemical_shifts_range)
        stream.seek(0)
        chemical_shifts = pd.read_csv(stream, delim_whitespace=True, names=unique_headers, quotechar="'", dtype={"Entry_ID": object})
        stream.close()
        length_comp_id = chemical_shifts.Comp_ID.str.len()
        if length_comp_id.min() == 3:
            return chemical_shifts
        else:
            return pd.DataFrame(columns=unique_headers) # Protein Error
        

    headers = headers_rex.findall(nmrstar_file)
    if headers:
        df_in_file = []
        unique_headers = list(OrderedSet(headers))
        runtimes = len(headers) / len(unique_headers)
        last_header = headers[-1]
        for n in range(int(runtimes)):
            chemical_shifts = create_dataframe(last_header, nmrstar_file, n)
            df_in_file.append(chemical_shifts)
        if len(df_in_file) > 1: ###
            print(nmrstar_file[:9]) ###
            n_rows = [df.shape[0] for df in df_in_file]
            max_n_row = max(n_rows)
            index_max_rows = n_rows.index(max_n_row)
            df_in_file = df_in_file.append(df_in_file[index_max_rows])
        ##if df_in_file:
            ##combined_chemical_shifts = pd.concat(df_in_file)
        headers_of_interest = ["Comp_index_ID","Comp_ID","Atom_ID","Atom_type","Val","Entry_ID"]
        chemical_shifts_short = chemical_shifts[chemical_shifts.columns.intersection(headers_of_interest)]
        if chemical_shifts_short.shape[0] != 0:
            if "Entry_ID" in chemical_shifts_short.columns:
                return chemical_shifts_short
            else:
                return "Entry_ID Error"
        else:
            return "Protein Error"
        #else:
        #    return "Protein Error"
    else:
        return "Column Header Error"
    
    
def get_pdb_id(pdb_id_rex, nmrstar_file):
    pdb_id = pdb_id_rex.search(nmrstar_file)
    try:
        return pdb_id.group(1)
    except AttributeError:
        return "PDB ID Error"


def get_sec_struc(sec_struc_rex, pdb_id):
    """dssp"""
    def load_dssp(pdb_id):
        """loads the dssp file from the api with the pdb id"""
        formdata = {"data": pdb_id}
        try:
            url_create = ("https://www3.cmbi.umcn.nl/xssp/api/create/pdb_id/dssp/")
            r = requests.post(url_create, data=formdata)
            r.raise_for_status()
            job_id = json.loads(r.text)['id']
            ready = False
            while not ready:
                url_status = ("https://www3.cmbi.umcn.nl/xssp/api/status/pdb_id/dssp/{}/".format(job_id))
                r = requests.get(url_status)
                r.raise_for_status()
                status = json.loads(r.text)['status']
                if status == "SUCCESS":
                    ready = True
                elif status in ["FAILURE", "REVOKED"]:
                    #raise Exception(json.loads(r.text)["message"])
                    return "DSSP Error"
                else:
                    time.sleep(2)
            else:  
                url_result = ("https://www3.cmbi.umcn.nl/xssp/api/result/pdb_id/dssp/{}/".format(job_id))
                r = requests.get(url_result)
                r.raise_for_status()
                result = json.loads(r.text)['result']
                return result

        except requests.exceptions.HTTPError:
            return "DSSP Error"

        
    def get_sec_struc_counts(sec_struc):
        """counts the appearance of symbols in sec struc list"""
        c = Counter(sec_struc)
        avail_sec_struc = ["H","B","E","G","I","T","S"]
        counts_in_order = [ c[sym] for sym in avail_sec_struc ]
        return counts_in_order

    dssp_file = load_dssp(pdb_id)
    if dssp_file != "DSSP Error":
        sec_struc = sec_struc_rex.findall(dssp_file)
        counts_in_order = get_sec_struc_counts(sec_struc)
        return counts_in_order, len(sec_struc)
    else:
        return ["DSSP Error"], "_"


def get_cond(pdb_id):
    """conditions"""
    
    def load_pdb_header(pdb_id):
        """loads the pbd header in a buffer"""
        buffer = io.BytesIO()
        c = pycurl.Curl()
        c.setopt(c.URL, f"https://files.rcsb.org/header/{pdb_id}.pdb")
        c.setopt(c.WRITEDATA, buffer)
        c.setopt(c.CAINFO, certifi.where())
        c.perform()
        c.close()
        body = buffer.getvalue() # byte string
        pdb_header = body.decode("iso-8859-1")
        return pdb_header
    
    def find_cond(pdb_header_lines):
        """find the values for the condition parameters"""
        cond_ensemble = {}
        cond_types = ["Temp", "PH", "Ionic_Strengh", "Pressure"]
        catch_count = 0
        line_count = 0
        line = pdb_header_lines[line_count] #here is a mistake, maybe empty file after 2900
        catch = False
        while "SAMPLE CONTENTS" not in line:
            if catch == True:
                values = line[delimiter_location+1 :]
                if ":" in line:
                    cond = {cond_types[catch_count] : values.strip()}
                    cond_ensemble.update(cond)
                    catch_count += 1
                else: 
                    cond_ensemble[cond_types[catch_count-1]] += values.strip()
            if "EXPERIMENT TYPE" in line and "NMR" in line:
                catch = True
                delimiter_location = line.find(":")
            line_count += 1
            line = pdb_header_lines[line_count]
        return cond_ensemble

    def make_iterable(pdb_header):
        """returns file like object to iterate through"""
        stream = io.StringIO()
        stream.write(pdb_header)
        stream.seek(0)
        pdb_header_lines = stream.readlines()
        stream.close()
        return pdb_header_lines
        
    pdb_header = load_pdb_header(pdb_id)
    if "The requested URL was not found on this server" in pdb_header:
        return "PDB Request Error"
    else:
        pdb_header_lines = make_iterable(pdb_header)
        cond = find_cond(pdb_header_lines)
        return cond

def from_list2table(cond_struc_depot):
    headers = ["Entry_ID","H","B","E","G","I","T","S","Temp","PH","Ionic_strength","Pressure","Seq_length"]
    cond_struc_table = pd.DataFrame(cond_struc_depot, columns = headers)
    return cond_struc_table

                
if __name__ == "__main__":
    headers_rex = re.compile("_Atom_chem_shift\.([^\s]*)")
    pdb_id_rex = re.compile("data_([^\s]*)")
    sec_struc_rex = re.compile("[A-Z]\ [A-Za-z]\ \ (.)")
    get_pdb_id_rex = partial(get_pdb_id, pdb_id_rex)
    get_sec_struc_rex = partial(get_sec_struc, sec_struc_rex)
    get_chemical_shifts_rex = partial(get_chemical_shifts, headers_rex)
    path2files = get_paths()
    all_chemical_shifts, all_struc_and_cond, files_not_found = get_data(path2files, get_pdb_id_rex, get_chemical_shifts_rex, get_sec_struc_rex)
    all_chemical_shifts.to_pickle("chemical-shifts.pkl")
    all_struc_and_cond.to_pickle("struc-and-cond.pkl")
    with open("files-not-found.txt", "w") as outfile:
        for error in files_not_found:
            outfile.write(f"{error}\n")
