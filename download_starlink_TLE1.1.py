#for downloading TLEs on a regular timescale

import pandas as pd
import requests
from datetime import date
import os.path

today = date.today()
# dd/mm/YY
d1 = today.strftime("%d/%m/%Y")
# print("d1 =", d1)

url = 'https://www.celestrak.com/NORAD/elements/starlink.txt'
req = requests.get(url)
url_content = req.content
txt_file = open('starlink_TLEdata.txt', 'wb')
txt_file.write(url_content)
# txt_file.close()

save_path = 'C:\\Users\\rasta\\Documents\\TLEs'
file_name = 'starlink_TLEdata.txt'
completeName = os.path.join(save_path, file_name)

file1 = open(completeName, "wb")
file1.write(url_content)

file1.close()
txt_file.close()

# file2 = open('starlink_TLEdata%s.txt'% d1, "w")
# file2.close()