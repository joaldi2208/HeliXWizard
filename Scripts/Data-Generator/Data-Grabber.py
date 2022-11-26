# before using the script, download the data with the following command 
# rsync -rlpt -v -z --delete --port=33444 \
# rsync.rcsb.org::ftp_data/structures/divided/nmr_chemical_shifts/ ../Data

import pycurl
import certifi
from io import BytesIO

import gzip
import os
from collections import namedtuple
from itertools import chain
from functools import partial

from tempfile import gettempdir
import biotite.sequence.io.genbank as gb
import biotite.database.entrez as entrez
import biotite

import re
# part 1
def get_paths():
    path2files = []
    for root, _, filenames in os.walk("../../Data/"):
        path2files.extend([os.path.join(root, name) for name in filenames])
    return path2files

# part 2
def get_ids(path2files):
    entry_ids = []
    files_not_found = 0
    for files in path2files:
        entry_id = namedtuple("entry_id", ["pdb", "bmrb"])
        try:
            with gzip.open(files, "rt") as f:
                for line in f:
                    if "_Entry.ID" in line:
                        entry_id.pdb = line.split()[1]
                    elif "_Entry.Assigned_BMRB_ID" in line: 
                        entry_id.bmrb = line.split()[1]
                        if entry_id.bmrb != "1" and entry_id.bmrb != "?":
                            entry_ids.append(entry_id) # bmrb-id is last interesting line in file
                            break # no need to read entire file
        except FileNotFoundError:
            files_not_found += 1
    return entry_ids, files_not_found


# part 4
def get_cond(cond_rex, bmrb_id):
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, f"https://bmrb.io/data_library/summary/index.php?bmrbId={bmrb_id}")
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.perform()
    c.close()
    body = buffer.getvalue() # byte string
    bmrb_file = body.decode("iso-8859-1")
    #print(bmrb_file)
    
    conditions = cond_rex.findall(str(bmrb_file))
   
    return conditions


    
def get_csv_HSQC(bmrb_id):
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, f"https://api.bmrb.io/current/entry/{bmrb_id}/simulate_hsqc?format=csv&filter=backbone")
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.perform()
    c.close()
    body = buffer.getvalue() # byte string
    csv_HSQC = body.decode("iso-8859-1")
    return csv_HSQC


def get_chain_type(chain_rex, pdb_id):
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, f"https://files.rcsb.org/header/{pdb_id}.pdb")
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.perform()
    c.close()
    body = buffer.getvalue() # byte string
    pdb_header = body.decode("iso-8859-1")
    
    chain_type = chain_rex.search(str(pdb_header)).group(1)
    return chain_type


def get_2nd_str_data(loc_rex, len_rex, type_rex, conditions_rex, gb_file):
    sec_str_loc = loc_rex.findall(str(gb_file))
    sec_str_type = type_rex.findall(str(gb_file))
    sec_str_len = len_rex.findall(str(gb_file))[1]
    if sec_str_len != None:
        return sec_str_loc, sec_str_type, sec_str_len
     

def get_data(entry_ids, get_chain_type_rex, get_2nd_str_data_rex, get_cond_rex):
    # you should save the data here in a list and then return the list
    data = []
    files_not_found = 0
    for entry_id in entry_ids:
        csv_HSQC = get_csv_HSQC(entry_id.bmrb)
        chain_type = get_chain_type_rex(entry_id.pdb)
        if "error" not in csv_HSQC:
            
            conditions = get_cond_rex(entry_id.bmrb)
            
            try:
                filename = entrez.fetch(f"{entry_id.pdb}_{chain_type}", gettempdir(), "gb", "protein", "gb")
                gb_file = gb.GenBankFile.read(filename)
                sec_str_loc, sec_str_type, sec_str_len = get_2nd_str_data_rex(gb_file)

            except biotite.database.RequestError:
                files_not_found += 1
            finally:
                data.append((entry_id, chain_type, sec_str_loc, sec_str_type, sec_str_len))
                print("n: ", len(data))
                print("pdb_id: ", entry_id.pdb)
                print("bmrb_id: ", entry_id.bmrb)
                print("chain type: ", chain_type)
                print("sec_str_loc: ", sec_str_loc)
                print("sec_str_type: ", sec_str_type)
                print("sec_str_len: ", sec_str_len)
                for cond, value in conditions:
                    print(cond, value) # not as it should be right now
                # https://bmrb.io/data_library/summary/index.php?bmrbId=30300 interesting example
                # example 3 shows more then one condition
                print("files_not_found: ", files_not_found)
                print("------------------------------------------------")
               
    return data
        

if __name__ == "__main__":
    chain_type_rex = re.compile("CHAIN:\ (\w)")
    loc_rex = re.compile("SecStr\ *(\d*)\.\.(\d*)")
    len_rex = re.compile("source\ *\d*..(\d*)")
    type_rex = re.compile("sec_str_type=\"([^\"]*)")
    cond_rex = re.compile("(pH:|pressure:|temperature:|ionic strength:)[^\d]*([\d.]*)")

    get_chain_type_rex = partial(get_chain_type, chain_type_rex)
    get_2nd_str_data_rex = partial(get_2nd_str_data, loc_rex, len_rex, type_rex, cond_rex)
    get_cond_rex = partial(get_cond, cond_rex)
    
    path2files = get_paths() # from data dir, should be changed later
    entry_ids, path_files_not_found = get_ids(path2files)
    print("number entry_ids: ", len(entry_ids))
    print("path_files_not_found: ", path_files_not_found)
    print("++++++++++++++++++++++++++++++++++++++++++")
    
    data = get_data(entry_ids, get_chain_type_rex, get_2nd_str_data_rex, get_cond_rex)
    # what about pdb_id:  6V1N, https://www.ncbi.nlm.nih.gov/Structure/pdb/6V1N
    print("++++++++++++++++++++++++++++++++++++++++++")
    print(len(data))
