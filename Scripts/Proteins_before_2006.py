from Data_Downloads import download_webpages
from datetime import datetime
from tqdm import tqdm
import re


def get_deposition_date(found_bmrb_ids):
    """greps deposition date from the bmrb webpage"""

    deposition_dates = []
    for bmrb_id in tqdm(found_bmrb_ids, total=len(found_bmrb_ids), desc="Dowloading BMRB Deposition Date"):
        
        bmrb_webpage = download_webpages(f"https://bmrb.io/data_library/summary/index.php?bmrbId={bmrb_id}")
        if "Page not found" in bmrb_webpage or "only exists in NMR-STAR 2.0" in bmrb_webpage:
            bmrb_webpage = download_webpages(f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_21.str")
            deposition_date = re.search(r"_Submission_date *(\d{4}-\d{2}-\d{2})", bmrb_webpage).group(1)
        else:
            deposition_date = re.search(r"Deposition date:(.|\n)*?(\d{4}-\d{2}-\d{2})", bmrb_webpage).group(2)
            
        deposition_dateformat = datetime.strptime(deposition_date, "%Y-%m-%d")
        deposition_dates.append(deposition_dateformat)

    return deposition_dates

def filter_only_before_2006(deposition_dates, found_brmb_ids):
    """return bmrb ids for those entries with a deposition date newer then 2005-09-20."""

    bmrb_ids_before_2006 = []
    for bmrb_id, deposition_dateformat in zip(found_bmrb_ids, deposition_dates):
        if deposition_dateformat < datetime.strptime("2005-09-20", "%Y-%m-%d"):
            bmrb_ids_before_2006.append(bmrb_id)
            
    return bmrb_ids_before_2006
            
            
if __name__=="__main__":
    with open("found_bmrb_ids.txt", "r") as infile:
        found_bmrb_ids_strings = infile.read().replace("\n","").split(",")
        found_bmrb_ids = [int(bmrb_id) for bmrb_id in found_bmrb_ids_strings]

    deposition_dates = get_deposition_date(found_bmrb_ids)
    print(len(deposition_dates))
    with open("deposition_dates.txt", "w") as outfile:
        bmrb_id_with_date = list(zip(deposition_dates, found_bmrb_ids))
        outfile.writelines(",".join(map(str, bmrb_id_with_date)))
            
    only_before_2006 = filter_only_before_2006(deposition_dates, found_bmrb_ids)
    print(len(only_before_2006))
    with open("proteins_before_2006.txt", "w") as outfile:
        outfile.writelines(",".join(map(str, only_before_2006)))

    
    
