#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:38:11 2017

Developed and tested in : 
- Python version 3.6.3
- Linux Ubuntu version 16.04 LTS (64-bit)

@author: Nikola Knezevic
"""

import os
import requests
import json
from collections import OrderedDict

############################# PARAMETERS #############################
# API key for Bot                                                    #
api_key="Here copy your Bot's API key."                              #
# list that represents json file for search obj                      #
search_obj=[("ra",""), ("dec",""), ("radius",""), ("units",""),      #
            ("objname",""), ("internal_name","")]                    #
# list that represents json file for get obj                         #
get_obj=[("objname",""), ("photometry","0"), ("spectra","1")]        #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
url_tns_api="https://wis-tns.weizmann.ac.il/api/get"                 #
url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"     #
######################################################################

############################# DIRECTORIES ############################
# current working directory                                          #
cwd=os.getcwd()                                                      #
# directory for downloaded files                                     #
download_dir=os.path.join(cwd,'downloaded_files')                    #
######################################################################

########################## API FUNCTIONS #############################
# function for changing data to json format                          #
def format_to_json(source):                                          #
    # change data to json format and return                          #
    parsed=json.loads(source,object_pairs_hook=OrderedDict)          #
    result=json.dumps(parsed,indent=4)                               #
    return result                                                    #
#--------------------------------------------------------------------#
# function for search obj                                            #
def search(url,json_list):                                           #
  try:                                                               #
    # url for search obj                                             #
    search_url=url+'/search'                                         #
    # change json_list to json format                                #
    json_file=OrderedDict(json_list)                                 #
    # construct the list of (key,value) pairs                        #
    search_data=[('api_key',(None, api_key)),                        #
                 ('data',(None,json.dumps(json_file)))]              #
    # search obj using request module                                #
    response=requests.post(search_url, files=search_data)            #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for get obj                                               #
def get(url,json_list):                                              #
  try:                                                               #
    # url for get obj                                                #
    get_url=url+'/object'                                            #
    # change json_list to json format                                #
    json_file=OrderedDict(json_list)                                 #
    # construct the list of (key,value) pairs                        #
    get_data=[('api_key',(None, api_key)),                           #
                 ('data',(None,json.dumps(json_file)))]              #
    # get obj using request module                                   #
    response=requests.post(get_url, files=get_data)                  #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for downloading file                                      #
def get_file(url):                                                   #
  try:                                                               #
    # take filename                                                  #
    filename=os.path.basename(url)                                   #
    # downloading file using request module                          #
    response=requests.post(url, files=[('api_key',(None, api_key))], #
                           stream=True)                              #
    # saving file                                                    #
    path=os.path.join(download_dir,filename)                         #
    if response.status_code == 200:                                  #
        with open(path, 'wb') as f:                                  #
            for chunk in response:                                   #
                f.write(chunk)                                       #
        print ('File : '+filename+' is successfully downloaded.')    #
    else:                                                            #
        print ('File : '+filename+' was not downloaded.')            #
        print ('Please check what went wrong.')                      #
  except Exception as e:                                             #
    print ('Error message : \n'+str(e))                              #
######################################################################

# EXAMPLE

# api key for your Bot
api_key="4aea2e3d2f25e2d93be7f2dc45b31340170a9779"

# Comment/Uncomment sections for testing the various examples:

"""

# search obj (here an example of cone search)
search_obj=[("ra","15:57:28"), ("dec","+30:03:39"), ("radius","5"), ("units","arcsec"),
            ("objname",""), ("internal_name","")]                    
response=search(url_tns_sandbox_api,search_obj)
if None not in response:
    # Here we just display the full json data as the response
    json_data=format_to_json(response.text)
    print (json_data)
else:
    print (response[1])


# get obj
get_obj=[("objname","2017A"), ("photometry","1"), ("spectra","1")]
response=get(url_tns_sandbox_api,get_obj)
if None not in response:
    # Here we just display the full json data as the response
    json_data=format_to_json(response.text)
    print (json_data)
else:
    print (response[1])


# get ascii file (according to the "asciifile" name obtained from the get obj reply)
ascii_file_url="https://sandbox-tns.weizmann.ac.il/system/files/uploaded/"\
               "Padova-Asiago/tns_2017A_2457777.69_Ekar_AFOSC_Padova-Asiago.txt"
get_file(ascii_file_url)

# get fits file (according to the "fitsfile" name obtained from the get obj reply)
fits_file_url="https://sandbox-tns.weizmann.ac.il/system/files/uploaded/"\
              "Padova-Asiago/tns_2017A_2457777.69_Ekar_AFOSC_Padova-Asiago.fits"
get_file(fits_file_url)

"""




