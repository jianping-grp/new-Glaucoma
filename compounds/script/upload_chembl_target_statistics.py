import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", 'glaucoma.settings')
import django
django.setup()
import pandas as pd

from compounds import models

def upload_target_statistics(line):
    target_id = line['target_id']
    reference_count = line['doi_count']
    compound_count = line['molecule_chembl_id_count']
    bioactivity_count = line['bio_count']

    target_object = models.Target.objects.get(id=target_id)
    target_object.reference_count = reference_count
    target_object.compound_count = compound_count
    target_object.bioactivity_count = bioactivity_count
    target_object.save()

if __name__ == '__main__':
    file_path = '/home/jianping/django_test/new-Glaucoma-master/compounds/data/chembl_target_statistics.csv'
    df = pd.read_csv(file_path)
    for idx, line in df.iterrows():
        print idx
        upload_target_statistics(line)

