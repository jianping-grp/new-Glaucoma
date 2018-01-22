import os
import pandas as pd

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()

from compounds.models import *
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned

def extract_pathway_info(cell):
    pathways = cell.split('\n') if isinstance(cell, str) else []
    for pathway in pathways:
        yield pathway.split(':')

def upload_drugs(line):
    drug_name = line['drug_name']
    drug_state = line['drug_state']
    iupac = line['IUPAC']
    drugbank_id = line['drugbank_id']
    cas = line['cas']
    cid = line['cid']
    smiles = line['smiles']
    pathways = extract_pathway_info(line['pathway'])


    drug, created = Drug.objects.get_or_create(
        drug_name=drug_name,
        drug_state=drug_state,
        IUPAC=iupac,
        drugbank_id=drugbank_id,
        cas=cas,
        cid=cid,
        smiles=smiles
    )

    for pathway in pathways:
        try:
            p, created = Pathway.objects.get_or_create(
                pathway_name=pathway[0],
                descripor=pathway[1]
            )
            p.drugs.add(drug)
            p.save()
        except pathway.DoesNotExist:
            print("Can't find %s:%s in database" %(pathway[0],pathway[1]))
        except pathway.MultipleObjectsReturned:
            print("return multiple pathways")

if __name__ == '__main__':
    drug_file = '/home/jianping/Desktop/weiyu/process/drugs/drugs_to_glaucoma.xlsx'
    df = pd.read_excel(drug_file, sheet_name=0)
    for idx, line in df.iterrows():
        print(idx)
        upload_drugs(line)

