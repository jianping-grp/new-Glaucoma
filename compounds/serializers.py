from . import models
from dynamic_rest import serializers
from rest_framework.permissions import AllowAny

class DrugSerializer(serializers.DynamicModelSerializer):
    targets = serializers.DynamicRelationField('TargetSerializer', many=True, deferred=True, embed=True)
    pathway_set = serializers.DynamicRelationField('PathwaySerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.Drug
        exclude = ['bfp', 'mol_block', 'mol']

class TargetSerializer(serializers.DynamicModelSerializer):
    drug_set = serializers.DynamicRelationField('DrugSerializer', many=True, deferred=True, embed=True)
    drugbankids = serializers.DynamicRelationField('DrugBankIDSerializer', many=True, deferred=True, embed=True)
    chembl_small_molecules_all_infos = serializers.DynamicRelationField('ChEMBL_small_molecule_all_info_Serializer', many=True, deferred=True, embed=True)
    pathway_set = serializers.DynamicRelationField('PathwaySerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.Target
        exclude = []

class PathwaySerializer(serializers.DynamicModelSerializer):
    drugs = serializers.DynamicRelationField('DrugSerializer', many=True, deferred=True, embed=True)
    targets = serializers.DynamicRelationField('TargetSerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.Pathway
        exclude = []

class DrugBankIDSerializer(serializers.DynamicModelSerializer):
    target_set = serializers.DynamicRelationField('TargetSerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.DrugBankID
        exclude = []

class ChEMBL_small_molecule_all_info_Serializer(serializers.DynamicModelSerializer):
    target_set = serializers.DynamicRelationField('TargetSerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.ChEMBL_small_molecule_all_info
        exclude = []

class FeedbackSerializer(serializers.DynamicModelSerializer):
    class Meta:
        model = models.Feedback
        exclude = []
