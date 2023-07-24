import sys
import requests
from bs4 import BeautifulSoup

pathway = sys.argv[1]
url = f'https://www.kegg.jp/entry/{pathway}'

# Send a GET request to the URL and retrieve the page content
response = requests.get(url)
html_code = response.text

soup = BeautifulSoup(html_code, 'html.parser')

# Find the 'Gene' section using its CSS class
gene_section = soup.find(class_='th31 deft tal vtop', string='Gene')

gene_info = gene_section.find_next_sibling()
gene_text = gene_info.get_text()

html_parts = gene_text.split(';')

for element in html_parts[:-1]:
    print(element.split()[-1])
