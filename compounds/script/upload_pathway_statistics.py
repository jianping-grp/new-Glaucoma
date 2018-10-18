import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", 'glaucoma.settings')
import django
django.setup()
import pandas as pd
from compounds import models


def upload_pathway_statistics():
    for pathway in models.Pathway.objects.all():
        print pathway
        # print pathway
        target_count = pathway.targets.count()
        drug_count = pathway.drugs.count()
        # print targets_count
        pathway.drug_count = drug_count
        pathway.target_count = target_count
        pathway.save()
if __name__ == '__main__':
    upload_pathway_statistics()
