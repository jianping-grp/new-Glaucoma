import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", 'glaucoma.settings')
import django
django.setup()
import pandas as pd

from compounds import models

def upload_structure(file):
    df = pd.read_csv(file)
    for idx, line in df.iterrows():
        print idx
        structure = line['smiles']
        jchem_iupac = line['jchem_iupac']
        generic_name = line['generic_name']
        synonyms = line['synonyms']
        products = line['products']
        id = line['id']
        drugbank_object = models.DrugBankID.objects.get(id=id)
        # drugbank_object.smiles = structure
        drugbank_object.synonyms = synonyms
        drugbank_object.generic_name = generic_name
        # drugbank_object.jchem_iupac = jchem_iupac
        # drugbank_object.products = products
        drugbank_object.save()

if __name__ == '__main__':
    file_path = '/home/jianping/django_test/new-Glaucoma-master/compounds/data/drugbank_structure_1.csv'
    upload_structure(file_path)
