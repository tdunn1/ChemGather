# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 12:01:20 2023

@author: Tim
"""

import requests
import time

class pubchem_access(object):
    def __init__(self, chem_id, data_type, id_type, max_rec = 101, sim_search = None):
        self.chem_id = chem_id
        self.data_type = data_type
        self.id_type = id_type
        self.max_rec = max_rec
        self.sim_search = sim_search
        
    def __repr__(self):
        return f'PubChem search results for {self.id_type} {self.chem_id}'
    
    def pubchem_req(self):
        if self.sim_search:
            base = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{self.sim_search}'
        else:
            base = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure'
        form = 'JSON'
        full_search = f'{base}/{self.id_type}/{self.chem_id}/{self.data_type}/{form}?MaxRecords={self.max_rec}'
        req = requests.get(full_search)
        while req.status_code == 202:
            time.sleep(5)
        
        return req
    
    def get_smis(self, group_up:bool, k:int):
        org_res = self.pubchem_req()
        res_dict = org_res.json()
        cids = res_dict['IdentifierList']['CID']
        smis = []
        if group_up:
            cids = self.frac_list(cids, k)
        for cid in cids:
            base = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
            req_input = 'compound/cid'
            req_output = 'property/CanonicalSMILES/txt'
            full_req = f'{base}/{req_input}/{cid}/{req_output}'
            re = requests.get(full_req)
            out_put = re.text
            out_put = out_put.split('\n')
            smis.extend(out_put[:-1]) #Remove empty space
            time.sleep(5)
            
        return smis
    
    def frac_list(self, cids_list, k):
        n = len(cids_list)//k
        r = len(cids_list)%k
        set_list = []
        for i in range(n):
            cid_nums = ','.join(map(str,cids_list[k*i:k*(i + 1)]))
            set_list.append(cid_nums)
            
        if r:
            cid_nums = ','.join(map(str,cids_list[-r:]))
            set_list.append(cid_nums)
        
        return set_list
            
            
if __name__ == '__main__':
    cid = 101
    test = pubchem_access(cid,'cids','cid', sim_search = 'fastsimilarity_2d')
    print(test.get_smis(group_up = True, k = 5))
