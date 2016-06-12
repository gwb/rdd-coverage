import sys
import urllib
import base64
import hashlib
import hmac
#import json
import unicodedata
import time
import csv
import requests
import os

#This script uses Harvard's Google Maps geocoding service, allowing for geocoding up to 250,000 addresses per day, per machine.
#Setup:  1)Save this script into the same folder your input addresses are in.
#        2)Format your input addresses into a tab-delimited text file without headers.
#        3)Change "geocoding_input.txt" below to the name of your input file name.
inputfile = r".\geocoding_input.txt"
#        4)Open a command prompt, and change directories into the folder where this file is.  At the command
#        prompt, type in:  google_geocoding_harvard.py  The script will run, geocoding each input address, and
#        outputting the results into a file named "geocoding_output.txt".
outputfile = r".\geocoding_output.txt"

google_url = b"http://maps.googleapis.com"
geocoding_endpoint = b"/maps/api/geocode/json?"
client = b"gme-harvarduniversity1"
channel = b""
from geocoding_private import privateKey

field1 = "ID"
field2 = "In_Address"
field3 = "In_City"
field4 = "In_State"
field5 = "In_Country"
field6 = "Address_Matched"
field7 = "City_Matched"
field8 = "State_Matched"
field9 = "Country_Matched"
field10 = "Location_Type"
field11 = "Latitude"
field12 = "Longitude"

def request_geojson(address):
   address_full = ",".join([address, "Milwaukee", "Wisconsin", "USA"])
   #decode the private key into its binary format
   decodeKey = base64.urlsafe_b64decode(privateKey)
   urltosign = geocoding_endpoint + encodedParams
   #create a signature using the private key and the url encoded, string using HMAC SHA1. This signature will be binary.
   signature = hmac.new(decodeKey, urltosign, hashlib.sha1)
   #encode the binary signature into base64 for use within a URL
   encodedsignature = base64.urlsafe_b64encode(signature.digest())
   signedurl = google_url + geocoding_endpoint + encodedParams + "&signature=" + encodedsignature
   payload = {"address": address, "client": client, "signature": encodedsignature}
   r = requests.get(signedurl, params=payload)
   data_json = r.json()
   return data_json

def request_geojson(address):
   address_full = ",".join([address, "Milwaukee", "Wisconsin", "USA"])
   #Generate valid signature
   encodedParams = urllib.parse.urlencode({"address":address_full, "client": client}
           ).encode("ascii","ignore")
   #decode the private key into its binary format
   decodeKey = base64.urlsafe_b64decode(privateKey)
   urltosign = geocoding_endpoint + encodedParams
   #create a signature using the private key and the url encoded, string using HMAC SHA1. This signature will be binary.
   signature = hmac.new(decodeKey, urltosign, hashlib.sha1)
   #encode the binary signature into base64 for use within a URL
   encodedsignature = base64.urlsafe_b64encode(signature.digest())
   signedurl = google_url + geocoding_endpoint + encodedParams + b"&signature=" + encodedsignature
   data = requests.get(signedurl)
   data_json = data.json()
   r = requests.get(signedurl)
   #print("get ", r.url)
   data_json = r.json()
   return data_json

def parse_geojson(geojson):
   latlng = geojson["results"][0]["geometry"]["location"]
   return latlng

def geocode_address(address):
   geojson = request_geojson(address)
   latlng = parse_geojson(geojson)
   return latlng

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
                latlng = geocode_address(address)
                print("%s: %f,%f" % (address, latlng["lat"], latlng["lng"]))
                csv_writer.writerow([row_dict["Taxkey"], address, latlng["lat"], latlng["lng"]])
                cache[address] = latlng

def _old():
    for line in f_in:
       fields = line.strip().replace("\"", "").split('\t')
       address = "%s+%s,%s,%s" % (fields[2], fields[1], fields[3], fields[4])
       address = unicodedata.normalize('NFKD', address.decode("utf-8", "replace")).encode('ascii', 'ignore')
       address = address.replace("n/a", "").replace(" ", "+")
       city_name = "N/A"
       admin1 = "N/A"
       country = "N/A"
       streetNum = "N/A"
       street = "N/A"
       for result in data_json["results"]:
          for component in result["address_components"]:
             if "locality" in component["types"]:
                city_name = component["long_name"]
          if city_name != "N/A":
             break
       for result in data_json["results"]:
          for component in result["address_components"]:
             if "administrative_area_level_1" in component["types"]:
                admin1 = component["long_name"]
          if admin1 != "N/A":
             break
       for result in data_json["results"]:
          for component in result["address_components"]:
             if "country" in component["types"]:
                country = component["long_name"]
          if country != "N/A":
             break
       for result in data_json["results"]:
          for component in result["address_components"]:
             if "street_number" in component["types"]:
                streetNum = component["long_name"]
          if streetNum != "N/A":
             break
       for result in data_json["results"]:
          for component in result["address_components"]:
             if "route" in component["types"]:
                street = component["long_name"]
          if street != "N/A":
             break
       p_city = unicodedata.normalize('NFKD', unicode(city_name)).encode('ascii', 'ignore')
       p_admin1 = unicodedata.normalize('NFKD', unicode(admin1)).encode('ascii', 'ignore')
       p_country = unicodedata.normalize('NFKD', unicode(country)).encode('ascii', 'ignore')
       p_address = unicodedata.normalize('NFKD', unicode(streetNum)).encode('ascii', 'ignore') + " " + unicodedata.normalize('NFKD', unicode(street)).encode('ascii', 'ignore')
       print("Processed" + "  " + fields[0] + "  "+  p_address)
       try:
          f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (fields[0], fields[1], fields[2], fields[3], fields[4], p_address, p_city, p_admin1 , p_country, data_json["results"][0]["geometry"]["location_type"], data_json["results"][0]["geometry"]["location"]["lat"], data_json["results"][0]["geometry"]["location"]["lng"]))
       #except: continue
       except: f_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (fields[0], fields[1], fields[2], fields[3], fields[4], "not found"))
       time.sleep(.5)
       #print data.read()

    print("finished")
    f_out.flush()
    f_in.close()
    f_out.close()
