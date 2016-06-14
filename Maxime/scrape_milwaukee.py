# coding: utf-8

import lxml
import requests
import lxml.html
import csv


resid_template = "%d_RVS_Dist%d"
condo_template = "%d_Condominium"
aptmt_template = "%d_Apartments"


url_head = "http://assessments.milwaukee.gov/SalesData/"
csv_directory = "Milwaukee_data/"


def scrape_page(page):
    r=requests.get(url_head + page + ".htm")
    tree=lxml.html.fromstring(r.text)
    table=tree.cssselect("table")[0]
    data_rows = table.cssselect("tr.the_tdm")
    colspan_row = table.cssselect("tr.the_colspan")[0]
    header_row = colspan_row.getnext()
    # The only identifying feature of the header row is its background color. 
    # Let's use it to check we've got the correct row:
    assert header_row.attrib["bgcolor"] == "#cc9966"
    header = [cell.text_content().strip() for cell in header_row.cssselect("td")]
    table_data = [[cell.text_content().strip() 
        for cell in row.cssselect("td")] 
        for row in data_rows]

    with open(csv_directory + page + ".csv", "w") as writefile:
        w=csv.writer(writefile)
        w.writerow(header)
        w.writerows(table_data)

for year in (2014,2015,2016):
    condo_page = condo_template % year
    scrape_page(condo_page)
    aptmt_page = aptmt_template % year
    scrape_page(aptmt_page)
    for district in range(1,16):
        resid_page = resid_template % (year, district)
        scrape_page(resid_page)

