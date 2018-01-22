import os
import pandas as pd
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()

from compounds.models import *

def upload_small_molecule_chembl(line):
    target_chembl_id = line['target_chembl_id']
    molecule_chembl_ids = line['molecule_chembl_id'].split('\n') if line['molecule_chembl_id'] else []
    smiles = line['canonical_smiles'].split('\n') if line['canonical_smiles'] else []

    molecule_info_set = set(zip(molecule_chembl_ids, smiles))

    target = Target.objects.get(chemblid=target_chembl_id)
    for molecule in molecule_info_set:
        molecule, create = ChEMBLSmallMolecule.objects.get_or_create(
            chembl_id=molecule[0],
            canonical_smiles=molecule[1]
        )
        target.chembl_small_molecules.add(molecule)
        target.save()

if __name__ == '__main__':
    file = '/home/jianping/Desktop/weiyu/process/chembl/small_molecule.csv'
    df = pd.read_csv(file)
    for idx, row in df.iterrows():
        print(idx)
        upload_small_molecule_chembl(row)




