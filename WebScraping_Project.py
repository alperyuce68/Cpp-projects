# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 20:03:24 2024

@author: Asus
"""

import requests  # Unused in this script, typically used for HTTP requests
from bs4 import BeautifulSoup  # For parsing HTML and XML documents
import csv  # For reading from and writing to CSV files
from selenium import webdriver  # To automate web browser interaction
from selenium.webdriver.common.by import By  # For locating elements within a webpage
from selenium.webdriver.support.ui import WebDriverWait  # To wait for certain conditions in the browser
from selenium.webdriver.support import expected_conditions as EC  # To specify what you are waiting for in a browser
import re  # For regular expression operations
import os  # For interacting with the operating system

querry=input("Enter company name, keyword or job title:")  # Asks the user for a job title or keyword
location=input("Enter desired location:")  # Asks the user for a location
base_url = "https://de.indeed.com/jobs?q={query}&l={location1}".format(query=querry, location1=location) # This line formats the base URL with the user-provided query and location.

driver = webdriver.Chrome()  # Initializes the Chrome WebDriver
driver.get(base_url)  # Navigates to the specified URL in the web browser

drivers = driver.page_source  # Retrieves the HTML source of the page

soup = BeautifulSoup(drivers, 'html.parser')  # Parses the HTML content using BeautifulSoup

job_postings = soup.find_all('span', id=re.compile("jobTitle"))  # Finds all job postings
job_names = soup.find_all('span', attrs={'data-testid':"company-name"})  # Finds all company names
job_locations = soup.find_all('div', attrs={'data-testid': "text-location"})  # Finds all job locations

jobs = []  # Initializes a list to store job data
for job_posting, job_name, job_location in zip(job_postings, job_names, job_locations):
    jobs.append([job_posting.text, job_name.text, job_location.text])  # Appends job data to the list

driver.quit()  # Closes the web browser

file_path = 'C:\\Users\Asus\Desktop\jobs_finder.csv'  # Sets the file path
file_exists = os.path.isfile(file_path)  # Checks if the file already exists

with open(file_path, mode='a', newline='', encoding='utf-8') as file:  # Opens the CSV file
    writer = csv.writer(file)
    writer.writerow(['Job Title', 'Company Name', 'Location'])  # Writes the header
    for start in range(20, 100, 10):  # Iterates through page numbers
        url = f'{base_url}&start={start}'  # Updates the URL with the new page number
        for job in jobs:  # Iterates through the jobs list
            writer.writerow(job)  # Writes each job to the CSV file

print(f'Job data have been saved to {file_path}')  # Prints a confirmation message