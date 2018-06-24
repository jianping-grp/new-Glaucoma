import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()
import pandas as pd
from compounds.models import *
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned

def extract_pathway(cell):
    pathways = cell.split('\n') if isinstance(cell, unicode) else []
    for pathway in pathways:
        yield pathway.split(':')

def upload_target_other_info(line):
    entry_name = line['entry_name']
    uniprot_accession = line['uniprot_accession']
    gene = line['gene']
    description = line['description']
    chembl_id = line['chembl_id']
    type = line['type']
    keggid = line['keggid']
    pdbid = line['pdbid']
    pathway_list = line['pathways'].split('\n') if not isinstance(line['pathways'], float) else []

    target, created = Target.objects.get_or_create(
        entry_name=entry_name,
        # uniprot_accession=uniprot_accession,
        # protein_description=description,
        # gene=gene
    )
    target.uniprot_accession = uniprot_accession
    target.protein_description = description
    target.gene = gene
    target.type = type,
    target.chemblid = chembl_id
    target.keggid = keggid
    target.pdbid = pdbid
    target.save()
    for pathway_info in pathway_list:
        el = pathway_info.split(':')
        pathway_name = el[0].strip()
        descriptor = el[1].strip()
        try:
            p, create = Pathway.objects.get_or_create(pathway_name=pathway_name, descripor=descriptor)
            p.targets.add(target)
            p.save()
        except:
            print(pathway_info)
if __name__ == '__main__':
    target_file = '/home/jianping/Desktop/weiyu/process/target/target_other_info.xlsx'
    df = pd.read_excel(target_file)
    for idx, line in df.iterrows():
        print(idx)
        upload_target_other_info(line)

