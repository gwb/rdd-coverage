import sys, urllib, base64, hashlib, hmac, json, unicodedata, time
from geocoding_private import privateKey

#This script uses Harvard's Google Maps geocoding service, allowing for geocoding up to 250,000 addresses per day, per machine.
#Setup:  1)Save this script into the same folder your input addresses are in.
#        2)Format your input addresses into a tab-delimited text file without headers.
#        3)Change "geocoding_input.txt" below to the name of your input file name.
inputfile = r".\geocoding_input.txt"
#        4)Open a command prompt, and change directories into the folder where this file is.  At the command
#        prompt, type in:  google_geocoding_harvard.py  The script will run, geocoding each input address, and
#        outputting the results into a file named "geocoding_output.txt".
outputfile = r".\geocoding_output.txt"

google_url = "http://maps.googleapis.com"
geocoding_endpoint = "/maps/api/geocode/json?"
client = "gme-harvarduniversity1"
channel = ""

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

f_in = open(inputfile, 'r')
f_out = open(outputfile, 'w')
f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (field1, field2, field3, field4, field5, field6, field7, field8, field9, field10, field11, field12))
for line in f_in:
   fields = line.strip().replace("\"", "").split('\t')
   address = "%s+%s,%s,%s" % (fields[2], fields[1], fields[3], fields[4])
   address = unicodedata.normalize('NFKD', address.decode("utf-8", "replace")).encode('ascii', 'ignore')
   address = address.replace("n/a", "").replace(" ", "+")
   #Generate valid signature
   encodedParams = urllib.urlencode({"address":address, "client": client})
   #decode the private key into its binary format
   decodeKey = base64.urlsafe_b64decode(privateKey)
   urltosign = geocoding_endpoint + encodedParams
   #create a signature using the private key and the url encoded, string using HMAC SHA1. This signature will be binary.
   signature = hmac.new(decodeKey, urltosign, hashlib.sha1)
   #encode the binary signature into base64 for use within a URL
   encodedsignature = base64.urlsafe_b64encode(signature.digest())
   signedurl = google_url + geocoding_endpoint + encodedParams + "&signature=" + encodedsignature
   data = urllib.urlopen(signedurl)
   data_json = json.loads(data.read())
   city_name = "N/A"
   admin1 = "N/A"
   country = "N/A"
   streetNum = "N/A"
   street = "N/A"
   for i in range(len(data_json["results"])):
      for component in data_json["results"][i]["address_components"]:
         if "locality" in component["types"]:
            city_name = component["long_name"]
      if city_name != "N/A":
         break
   for i in range(len(data_json["results"])):
      for component in data_json["results"][i]["address_components"]:
         if "administrative_area_level_1" in component["types"]:
            admin1 = component["long_name"]
      if admin1 != "N/A":
         break
   for i in range(len(data_json["results"])):
      for component in data_json["results"][i]["address_components"]:
         if "country" in component["types"]:
            country = component["long_name"]
      if country != "N/A":
         break
   for i in range(len(data_json["results"])):
      for component in data_json["results"][i]["address_components"]:
         if "street_number" in component["types"]:
            streetNum = component["long_name"]
      if streetNum != "N/A":
         break
   for i in range(len(data_json["results"])):
      for component in data_json["results"][i]["address_components"]:
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

print "finished"
f_out.flush()
f_in.close()
f_out.close()
