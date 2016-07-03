import sys
import urllib
import base64
import hashlib
import hmac
import unicodedata
import time
import csv
import requests
import os

#This script uses Harvard's Google Maps geocoding service, allowing for geocoding up to 250,000 addresses per day, per machine.

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
