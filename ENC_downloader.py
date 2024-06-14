from tkinter import *
from tkinter import filedialog

from lxml import etree

import requests
from bs4 import BeautifulSoup

from pathlib import Path

from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO


enc_url = 'https://www.charts.noaa.gov/ENCs/ENCsIndv.shtml?msclkid=bcd1d4a1d09011ec84ec21ed4ec896d3'


# File browser
def openFile():
    global filepath
    filepath = filedialog.askopenfilename(title='Open XML',
                                          filetypes=(('XML Files','*.xml'),
                                                     ('All Files','*.*')))
    file = open(filepath, 'r')
    entry1.insert(END, filepath)
    file.close()

# Scans XML for ENC Name, then matches ENC Name with ZIP file found in enc_url
# The zipfile.extractall() function below is currently set to my own local directory for ease of testing--feel free to
#   set it to your own directory
def scanFile():
    tree = etree.parse(filepath)
    root = tree.getroot()

    enc_list = []

    global scan_result
    for child in root.iter("{*}" + "ENC"):
        for grandchild in child.iter("{*}" + "name"):
            scan_result = grandchild.text
            enc_list.insert(1, scan_result)

    reqs = requests.get(enc_url)
    soup = BeautifulSoup(reqs.text, 'html.parser')
    soup.find_all('a')

    for link in soup.find_all('a'):
        for enc in enc_list:
            if enc + '.zip' in link['href']:
                url = 'https://www.charts.noaa.gov/ENCs/' + link['href']
                zip_open = urlopen(url)
                zipfile = ZipFile(BytesIO(zip_open.read()))
                zipfile.extractall(str(Path.home() / "Downloads") + '/' + link['href'].replace('.zip', ''))
                print('Download Complete')

# Primitive GUI
window = Tk()
window.title('ENC Downloader')
window.geometry('500x500')

entry1 = Entry(window, width=70, bg='white')
entry1.pack()

button = Button(text='Browse', command=openFile)
button.pack()

button = Button(text='Scan', command=scanFile)
button.pack()

window.mainloop()