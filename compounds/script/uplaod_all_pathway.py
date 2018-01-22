import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
import pandas as pd
django.setup()

from compounds.models import *

def upload_pathways(line):
    pathway_name = line['ID']
    descriptor = line['Name']
    p, created = Pathway.objects.get_or_create(
        pathway_name=pathway_name,
        descripor=descriptor
    )

if __name__ == '__main__':
    pathway_file = '/home/jianping/Desktop/weiyu/pathway/path.xlsx'
    df = pd.read_excel(pathway_file, sheet_name=0)
    for idx, line in df.iterrows():
        print(idx)
        upload_pathways(line)
