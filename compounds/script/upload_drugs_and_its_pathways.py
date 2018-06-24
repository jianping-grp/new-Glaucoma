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
    pathways = line['pathway'].split('\n') if not isinstance(line['pathway'], float) else []

    drug = Drug.objects.get(
    # drug, created = Drug.objects.get_or_create(
    #     drug_name=drug_name,
    #     drug_state=drug_state,
    #     IUPAC=iupac,
    #     drugbank_id=drugbank_id,
        cas=cas,
        cid=cid,
        # smiles=smiles
    )

    for pathway in pathways:
        el = pathway.split(':')
        pathway_name = el[0].strip()
        descriptor = el[1].strip()
        pathway, created = Pathway.objects.get_or_create(
            pathway_name=pathway_name,
            descripor=descriptor
        )
        pathway.drugs.add(drug)
        pathway.save()
        # try:
        #     # p = Pathway.objects.get(
        #     p, created = Pathway.objects.get_or_create(
        #         pathway_name=pathway_name,
        #         descripor=descriptor
        #     )
        #     p.drugs.add(drug)
        #     p.save()
        # except pathway.DoesNotExist:
        #     print("Can't find %s:%s in database" %(pathway[0],pathway[1]))
        # except pathway.MultipleObjectsReturned:
        #     print("return multiple pathways")

if __name__ == '__main__':
    drug_file = '/home/jianping/Desktop/weiyu/process/drugs/drugs_to_glaucoma.xlsx'
    df = pd.read_excel(drug_file, sheet_name=0)
    for idx, line in df.iterrows():
        print(idx)
        upload_drugs(line)