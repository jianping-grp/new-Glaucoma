import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()
import pandas as pd
from compounds.models import *
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned

def upload_target(line):
    drug_name = line['drug_name']
    drug_state = line['drug_state']
    entry_name_list = line['entry_name'].split(', ') if isinstance(line['entry_name'], unicode) else []
    uniprot_accession_list = line['uniprot_accession'].split(', ') if isinstance(line['uniprot_accession'], unicode) else []
    uniprot_accession_list = [el.strip() for el in uniprot_accession_list]
    gene_list = line['Gene'].split(', ') if isinstance(line['Gene'], unicode) else []
    protein_description_list = line['protein_description'].split(', ') if isinstance(line['protein_description'], unicode) else []
    target_info_list = zip(entry_name_list, uniprot_accession_list, gene_list, protein_description_list)

    drug, created = Drug.objects.get_or_create(drug_name=drug_name, drug_state=drug_state)
    for target_info in target_info_list:
        # print target_info
        t, created = Target.objects.get_or_create(
            entry_name=target_info[0],
            # uniprot_accession=target_info[1],
            # gene=target_info[2],
            # protein_description=target_info[3]
        )
        t.uniprot_accession = target_info[1]
        t.gene = target_info[2]
        t.protein_description = target_info[3]
        t.save()
        drug.targets.add(t)
        drug.save()

if __name__ == '__main__':
    target_file = '/home/jianping/Desktop/weiyu/process/drugs/drug_related_protein.xlsx'
    df = pd.read_excel(target_file, sheet_name=0)
    for idx, line in df.iterrows():
        print(idx)
        # print type(line['entry_name']) is not float
        upload_target(line)
