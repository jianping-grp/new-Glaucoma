import os
import pandas as pd

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()
from compounds.models import *

def upload_small_molecule_all_info(line):
    target_chembl_id = line['target_chembl_id']
    molecule_chembl_id = line['molecule_chembl_id']
    molecule_smile = line['canonical_smiles']

    activity_pchembl_value = line['pchembl_value']
    activity_standard_value = line['standard_value']
    standard_type = line['standard_type']
    standard_units = line['standard_units']
    activity_doc_id = line['doc_id']
    activity_assay_id = line['assay_assay_id']

    assay_chembl_id = line['assay_chembl_id']
    assay_description = line['description']
    assay_assay_id = line['assay_assay_id']
    doc_chembl_id = line['doc_chembl_id']
    doc_doi = line['doi']
    docs_title = line['title']
    doc_pubmed_id = line['pubmed_id']

    target = Target.objects.get(chemblid=target_chembl_id)

    compound, created = ChEMBL_small_molecule_all_info.objects.get_or_create(
        molecule_chembl_id=molecule_chembl_id,
        molecule_smile=molecule_smile,
        activity_pchembl_value=activity_pchembl_value,
        activity_standard_value=activity_standard_value,
        activity_standard_type=standard_type,
        activity_standard_units=standard_units,
        activity_doc_id=activity_doc_id,
        activity_assay_id=activity_assay_id,
        assay_chembl_id=assay_chembl_id,
        assay_description=assay_description,
        assay_assay_id=assay_assay_id,
        doc_chembl_id=doc_chembl_id,
        doc_doi=doc_doi,
        docs_title=docs_title,
        doc_pubmed_id=doc_pubmed_id,
    )
    target.chembl_small_molecules_all_infos.add(compound)
    target.save()

if __name__ == '__main__':
    file = '/home/jianping/Desktop/weiyu/process/chembl/all_info.csv'
    df = pd.read_csv(file)
    for idx, row in df.iterrows():
        print(idx)
        upload_small_molecule_all_info(row)





