import os
import pandas as pd
os.environ.setdefault("DJANGO_SETTINGS_MODULE", 'glaucoma.settings')
import django
django.setup()

from compounds.models import *

def upload_assay(line):
    molecule_chembl_id = line['molecule_chembl_id']
    assay_chembl_id = line['assay_chembl_id'].split('\n') if line['assay_chembl_id'] else []
    descripion = line['description'].split('\n') if isinstance(line['description'], unicode) else []
    assay_assay_id = line['assay_assay_id'].split('\n') if line['assay_assay_id'] else []

    assay_info_set = set(zip(assay_chembl_id, descripion, assay_assay_id))
    small_molecule = ChEMBLSmallMolecule.objects.get(chembl_id=molecule_chembl_id)

    for assay_info in assay_info_set:
        assay, created = Assay.objects.get_or_create(
            chembl_id=assay_info[0],
            description=assay_info[1],
            assay_id=assay_info[2]
        )
        assay.chembl_small_molecules.add(small_molecule)
        assay.save()

if __name__ == '__main__':
    file = '/home/jianping/Desktop/weiyu/process/chembl/assay.csv'
    df = pd.read_csv(file)
    for idx, line in df.iterrows():
        print(idx)
        upload_assay(line)

