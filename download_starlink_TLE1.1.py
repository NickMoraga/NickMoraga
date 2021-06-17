#for downloading TLEs on a regular timescale

import requests
from datetime import date
import os.path

# creates date
today = date.today()
# dd/mm/YY
d1 = today.strftime("%d-%m-%Y")


# takes Url and downloads into bytes(?)
url = 'https://www.celestrak.com/NORAD/elements/starlink.txt'
req = requests.get(url)
url_content = req.content

#adds timstamp to name, and directs to the right file
save_path = 'C:\\Users\\rasta\\Documents\\TLEs'
save_file = f'starlink_TLEdata-{d1}.txt'
completeName = os.path.join(save_path, file_name)
with open(completeName, 'wb') as txt_file:

    txt_file.write(url_content)

    # txt_file = open('starlink_TLEdata.txt', 'wb')
    # txt_file.write(url_content)
    # txt_file.close()

    save_path = 'C:\\Users\\rasta\\Documents\\TLEs'
    file_name = 'starlink_TLEdata.txt'
    completeName = os.path.join(save_path, file_name)

    file1 = open(completeName, "wb")
    file1.write(url_content)

    file1.close()

# txt_file.close()
