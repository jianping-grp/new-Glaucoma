import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
import pandas as pd
django.setup()

from compounds.models import Target

for target in Target.objects.all():
    print target.id
    description = target.protein_description
    # target.update(target_description=description.strip())
    Target.objects.filter(id=target.id).update(protein_description=description.strip())