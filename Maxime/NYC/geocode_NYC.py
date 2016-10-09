import sys
sys.path.append("..")
from geocode import *
import time
import pandas
import re
import numpy as np

doWrite = False

raw_data_dir = "NYC_data/raw_data/"
raw_csv_files = os.listdir(raw_data_dir)
raw_csv_files = [rf for rf in raw_csv_files if rf.endswith(".csv")]
output_file_dir = "NYC_data/geocoded/"
out_header = ["ADDRESS", "ZIP CODE", "lat", "lng"]
cache = {}
for raw_filename in raw_csv_files:
    basename = raw_filename.split(".csv")[0]
    #borough = re.match("rollingsales_(\w+)", basename).groups()[0]

    out_filename = output_file_dir+basename+"_geocoded.csv"
    already_geocoded = set([])
    try:
        with open(out_filename, "r") as out_file:
            existing=True
            csv_existing = csv.reader(out_file)
            header=next(csv_existing)
            assert header==out_header
            for row in csv_existing:
                row_dict = dict(zip(header, row))
                address = row_dict["ADDRESS"]
                zipcode = row_dict["ZIP CODE"]
                lat = row_dict["lat"]
                lng = row_dict["lng"]
                cache[address, zipcode] = {"lat": lat, "lng": lng}
    except FileNotFoundError:
            existing=False

    with open(out_filename, "a") as out_file:
        csv_writer = csv.writer(out_file)
        if existing:
            pass
        else:
            csv_writer.writerow(out_header)
        with open(raw_data_dir+raw_filename, "r") as raw_file:
            csv_file = pandas.read_csv(raw_file, header=4, dtype={"ZIP CODE": str})
            for (i,row) in csv_file.iterrows():
                address = row.ADDRESS
                zipcode = row["ZIP CODE"]
                if not isinstance(address,str):
                    continue
                if not zipcode:
                    continue
                address_full = ",".join([address, zipcode, "New York", "NY", "USA"])
                if (address, zipcode) in cache:
                    print("Address in cache: %s" % address)
                    continue
                try:
                    latlng = geocode_address(address_full)
                except:
                    print("failed getting latlng for %s" % address_full)
                    time.sleep(1)
                    continue
                    #raise
                print("%s %s: %f,%f" % (address, zipcode, latlng["lat"], latlng["lng"]))
                csv_writer.writerow([address, zipcode, latlng["lat"], latlng["lng"]])
                cache[address, zipcode] = latlng
