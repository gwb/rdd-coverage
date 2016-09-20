from geocode import *
import time
import pandas

raw_data_dir = "NYC_data/raw_data/"
raw_csv_files = os.listdir(raw_data_dir)
raw_csv_files = [rf for rf in raw_csv_files if rf.endswith(".csv")]
output_file_dir = "NYC_data/geocoded/"
out_header = ["Address", "lat", "lng"]
cache = {}
for raw_filename in raw_csv_files:
    basename = raw_filename.split(".csv")[0]
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
                address = row_dict["Address"]
                lat = row_dict["lat"]
                lng = row_dict["lng"]
                cache[address] = {"lat": lat, "lng": lng}
    except FileNotFoundError:
            existing=False

    with open(out_filename, "a") as out_file:
        csv_writer = csv.writer(out_file)
        if existing:
            pass
        else:
            csv_writer.writerow(out_header)
        with open(raw_data_dir+raw_filename, "r") as raw_file:
            csv_file = pandas.read_csv(raw_file, header=4)
            for (i,row) in csv_file.iterrows():
                address = row.ADDRESS
                if not isinstance(address,str):
                    continue
                if address in cache:
                    print("Address in cache: %s" % address)
                    continue
                address_full = ",".join([address, "New York", "NY", "USA"])
                try:
                    latlng = geocode_address(address_full)
                except:
                    print("failed getting latlng for %s" % address_full)
                    time.sleep(1)
                    continue
                    #raise
                print("%s: %f,%f" % (address, latlng["lat"], latlng["lng"]))
                csv_writer.writerow([address, latlng["lat"], latlng["lng"]])
                cache[address] = latlng
