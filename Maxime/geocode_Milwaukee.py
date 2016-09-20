from geocode import *

raw_data_dir = "Milwaukee_data/raw_data/"
raw_csv_files = os.listdir(raw_data_dir)
raw_csv_files = [rf for rf in raw_csv_files if rf.endswith(".csv")]
output_file_dir = "Milwaukee_data/geocoded/"
out_header = ["Taxkey", "Address", "lat", "lng"]
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
                taxkey = row_dict["Taxkey"]
                address = row_dict["Address"]
                lat = row_dict["lat"]
                lng = row_dict["lng"]
                already_geocoded.add(taxkey)
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
            csv_file = csv.reader(raw_file)
            header = next(csv_file)
            for row in csv_file:
                row_dict = dict(zip(header, row))
                taxkey = row_dict["Taxkey"]
                if taxkey in already_geocoded:
                    print("skipping %s" % taxkey)
                    continue
                address = row_dict["Address"]
                if address in cache:
                    print("Address in cache: %s" % address)
                    csv_writer.writerow([row_dict["Taxkey"], address, cache[address]["lat"], cache[address]["lng"]])
                    continue
                address_full = ",".join([address, "Milwaukee", "Wisconsin", "USA"])
                latlng = geocode_address(address_full)
                print("%s: %f,%f" % (address, latlng["lat"], latlng["lng"]))
                csv_writer.writerow([row_dict["Taxkey"], address, latlng["lat"], latlng["lng"]])
                cache[address] = latlng
