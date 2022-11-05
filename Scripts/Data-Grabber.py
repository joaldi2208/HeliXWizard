
# rsync -rlpt -v -z --delete --port=33444 \
# rsync.rcsb.org::ftp_data/structures/divided/nmr_chemical_shifts/ ../Data
# before using the script, download the data with the command above
import pycurl
import certifi
from io import BytesIO

import gzip
import os
from collections import namedtuple
from itertools import chain

from tempfile import gettempdir
import biotite.sequence.io.genbank as gb
import biotite.database.entrez as entrez

# part 1
path2files = []
for root, _, filenames in os.walk("../Data/"):
    path2files.extend([os.path.join(root, name) for name in filenames])


# part 2
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
        #print(files)

print(files_not_found)
print(len(entry_ids))

#https://api.bmrb.io/current/entry/34549/simulate_hsqc?format=csv&filter=backbone

# part 3
for entry in entry_ids:
    pdb_id = entry.pdb
    bmrb_id = entry.bmrb
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, f"https://api.bmrb.io/current/entry/{bmrb_id}/simulate_hsqc?format=csv&filter=backbone")
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.perform()
    c.close()
    body = buffer.getvalue() # byte string
    csv_HSQC = body.decode("iso-8859-1")
    z = 0
    if "error" not in csv_HSQC:
        #print(csv_HSQC)
        z += 1
        filename = entrez.fetch(f"{pdb_id}_A", gettempdir(), "gb", "protein", "gb")
        gb_file = gb.GenBankFile.read(filename)
        annotation = gb.get_annotation(gb_file, include_only=["SecStr"])
        print(annotation)
    
    # TODO: How to save the data???

